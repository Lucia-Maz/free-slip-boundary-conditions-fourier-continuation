#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Puente numérico para [P-02]: decaimiento de una condición inicial de banda ancha.

La forma del espectro horizontal sale de los cuatro PIV en t ~= 11,1 s guardados en
``salidas/tablas/decaimiento.json``. Se construye un campo horizontal solenoidal con
fases deterministas y el perfil vertical fundamental del par rígido--free-slip,

    u_h(x,y,z) = sin(pi z / 2h) grad_perp psi(x,y),   w = 0.

Se corre la MISMA forma espectral a amplitud experimental (delta ~= 9) y a amplitud
1e-4 veces menor. La segunda es el control lineal medido con el propio código. Los dos
casos usan fondo no-deslizante, tope free-slip plano y periodicidad horizontal.

La malla 16^2 sólo es una prueba de la cadena y trunca el espectro. La 32^2 contiene los
anillos medidos 1--10, pero el extremo alto queda cerca del corte de dealiasing: es un
diagnóstico local, no el resultado convergido. La convergencia 64^2/128^2 va al clúster.
Dentro de una asignación de SLURM se lanza con ``srun --mpi=pmix``; fuera de ella se usa
el ``mpirun`` del prefijo ``SPECTER_TOOLCHAIN``.

Uso local, siempre secuencial:

    ulimit -v 2500000
    /home/lucia/miniforge3/envs/piv-dt/bin/python \
        verificacion/decaimiento_banda_ancha.py --resolucion 16 --prueba
    /home/lucia/miniforge3/envs/piv-dt/bin/python \
        verificacion/decaimiento_banda_ancha.py --resolucion 16
    /home/lucia/miniforge3/envs/piv-dt/bin/python \
        verificacion/decaimiento_banda_ancha.py --resolucion 32

No modifica las puertas de aceptación. Compila una copia temporal de SPECTER y escribe
``salidas/tablas/decaimiento_banda_ancha.json`` después de cada corrida.
"""

import argparse
import importlib.util
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(AQUI)
CODIGO = os.path.join(RAIZ, "codigo")
sys.path.insert(0, CODIGO)
from celda import CAMPANA_02_06_25 as CELDA  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "fase2", os.path.join(AQUI, "test_aceptacion_fase2.py"))
fase2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fase2)

ENTRADA = os.path.join(RAIZ, "salidas", "tablas", "decaimiento.json")
SALIDA = os.path.join(RAIZ, "salidas", "tablas", "decaimiento_banda_ancha.json")

NZ, CZ = 64, 25
NZ_FIS = NZ - CZ
DT = 5.0e-4
CSTEP = 100
LZ_REF, NU_REF = 0.5, 1.0e-2
T_OBJETIVO_S = 11.133333333333333
ANILLOS_OBJETIVO = tuple(range(1, 11))
ANILLOS_AJUSTE = tuple(range(3, 11))  # mismos ocho bins de [V-18]
SEMILLA_FASES = 20260921
FACTOR_LINEAL = 1.0e-4
INTERVALO_CAMPOS_MEMORIAS = 3.2
N_INTERVALOS = 5


def espectro_piv_limpio(e):
    """Mismo piso por modo que 07/13, y corrección de la separación n entre cuadros."""
    k = np.asarray(e["k_1_m"], float)
    energia = np.asarray(e["E"], float)
    cuenta = np.asarray(e.get("modos_por_anillo", k), float)
    k_nyquist = math.pi / float(e["paso_grilla_m"])
    banda = (k > 0.35 * k_nyquist) & (k < 0.95 * k_nyquist)
    piso = float(np.median((energia / np.maximum(cuenta, 1.0))[banda]))
    return k, np.maximum(energia - piso * cuenta, 0.0) / float(e["n"]) ** 2


def objetivo_desde_piv():
    """Pesos por anillo y delta en t~=11 s, con procedencia en el propio JSON."""
    with open(ENTRADA) as fh:
        datos = json.load(fh)
    registros = sorted(datos["registros"])
    espectros, velocidades = [], []
    dk = None
    for registro in registros:
        r = datos["registros"][registro]
        e = min(r["espectros"], key=lambda x: abs(x["t_s"] - T_OBJETIVO_S))
        k, energia = espectro_piv_limpio(e)
        if dk is None:
            dk = float(k[1] - k[0])
        pesos = []
        for n in ANILLOS_OBJETIVO:
            mascara = (k / dk >= n) & (k / dk < n + 1)
            pesos.append(float(energia[mascara].sum()))
        pesos = np.asarray(pesos)
        pesos /= pesos.sum()
        espectros.append(pesos)
        p = min(r["calibracion"], key=lambda x: abs(x["t_s"] - T_OBJETIVO_S))
        velocidades.append(float(p["u_rms_superficie_m_s"]))

    pesos = np.mean(espectros, axis=0)
    pesos /= pesos.sum()
    # Los centros (n+1/2)dk son exactamente los guardados por el binneado del PIV.
    k_rms_fis = math.sqrt(sum(w * ((n + 0.5) * dk) ** 2
                               for n, w in zip(ANILLOS_OBJETIVO, pesos)))
    u_sup = float(np.median(velocidades))
    u_prom = u_sup / CELDA.sesgo_piv_superficie()
    delta = u_prom * k_rms_fis * CELDA.h ** 2 / CELDA.fluido.nu
    return {
        "registros": registros,
        "t_objetivo_s": T_OBJETIVO_S,
        "dk_fisico_1_m": dk,
        "pesos_anillos": {str(n): float(w)
                           for n, w in zip(ANILLOS_OBJETIVO, pesos)},
        "k_rms_fisico_1_m": k_rms_fis,
        "u_rms_superficie_mediana_m_s": u_sup,
        "delta_objetivo": delta,
    }


def modos_discretos(nxy, objetivo):
    """Una mitad independiente del plano de Fourier real, con energía normalizada."""
    k_componente_max = nxy // 3
    anillo_max = min(max(ANILLOS_OBJETIVO), k_componente_max)
    pesos_obj = {int(k): v for k, v in objetivo["pesos_anillos"].items()
                 if int(k) <= anillo_max}
    norma = sum(pesos_obj.values())
    pesos_obj = {k: v / norma for k, v in pesos_obj.items()}
    rng = np.random.default_rng(SEMILLA_FASES)
    modos = []
    for n in sorted(pesos_obj):
        candidatos = []
        # kx>0 representa una pareja conjugada completa. ky puede tener ambos signos.
        for kx in range(1, k_componente_max + 1):
            for ky in range(-k_componente_max, k_componente_max + 1):
                if int(math.floor(math.hypot(kx, ky))) == n:
                    candidatos.append((kx, ky))
        if not candidatos:
            raise RuntimeError("el anillo %d no tiene modos en la malla %d" % (n, nxy))
        peso_modo = pesos_obj[n] / len(candidatos)
        for kx, ky in candidatos:
            k2 = float(kx * kx + ky * ky)
            # psi=A cos(k.r+phi): <|grad psi|^2>=A^2 k^2/2=peso_modo.
            amplitud_psi = math.sqrt(2.0 * peso_modo / k2)
            modos.append({"kx": kx, "ky": ky, "anillo": n, "k2": k2,
                          "peso": peso_modo, "amplitud_psi": amplitud_psi,
                          "fase": float(rng.uniform(0.0, 2.0 * math.pi))})
    k_rms = math.sqrt(sum(m["peso"] * m["k2"] for m in modos))
    return modos, pesos_obj, k_rms


def condicion_inicial(modos, nxy):
    """Fuente Fortran de la condición inicial solenoidal de banda ancha."""
    lineas = [
        "! Banda ancha derivada del PIV en t ~= 11.1 s ([V-18]).",
        "! u_h=sin(pi*z/2Lz)*grad_perp(psi), w=0; fases deterministas.",
        "      vx = 0; vy = 0; vz = 0",
        "",
        "      IF (ista .eq. 1) THEN",
        "         DO k = 1,nz-Cz",
        "            rmt = SIN(0.5_GP*pi*z(k)/Lz)",
    ]
    for m in modos:
        ix = m["kx"] + 1
        iy = m["ky"] + 1 if m["ky"] >= 0 else nxy + m["ky"] + 1
        # Coeficiente complejo de A cos(k.r+phi): (A/2) exp(i phi).
        coef = 0.5 * m["amplitud_psi"]
        c, s = math.cos(m["fase"]), math.sin(m["fase"])
        lineas += [
            "            rmq = %.16e_GP*u0*nx*ny" % coef,
            ("            vx(k,%d,%d) = vx(k,%d,%d) + im*ky(%d)*rmq*"
             "(%.16e_GP+im*%.16e_GP)*rmt")
            % (iy, ix, iy, ix, iy, c, s),
            ("            vy(k,%d,%d) = vy(k,%d,%d) - im*kx(%d)*rmq*"
             "(%.16e_GP+im*%.16e_GP)*rmt")
            % (iy, ix, iy, ix, ix, c, s),
        ]
    lineas += [
        "         ENDDO",
        "      ENDIF",
        "",
        "      CALL fftp1d_real_to_complex_z(planfc,vx,MPI_COMM_WORLD)",
        "      CALL fftp1d_real_to_complex_z(planfc,vy,MPI_COMM_WORLD)",
        "      CALL fftp1d_real_to_complex_z(planfc,vz,MPI_COMM_WORLD)",
        "",
    ]
    return "\n".join(lineas)


def compilar_liviano():
    """Como [V-16]: -O1 evita el pico de memoria de gfortran en esta laptop."""
    original = fase2._parchear_makefile_in

    def envoltura(ruta, solver="HD"):
        original(ruta, solver)
        with open(ruta) as fh:
            texto = fh.read()
        with open(ruta, "w") as fh:
            fh.write(texto.replace("-O2", "-O1"))

    fase2._parchear_makefile_in = envoltura


def correr_largo(bindir, etiqueta, *, lz, nu, u0, pasos, timeout=86400):
    """Equivalente a fase2.correr con timeout apto para 16 memorias modales."""
    dirrun = os.path.join(os.path.dirname(bindir), "corridas", etiqueta)
    os.makedirs(os.path.join(dirrun, "out"), exist_ok=True)
    os.makedirs(os.path.join(dirrun, "in"), exist_ok=True)
    with open(os.path.join(dirrun, "parameter.inp"), "w") as fh:
        fh.write(fase2.PLANTILLA_INP.format(
            lz=lz, dt=DT, step=pasos, cstep=CSTEP, u0=u0, nu=nu,
            vparam0=0.0, vparam1=1.0, bcsta="noslip", bcend=fase2.CADENA_LIBRE))
    enlace = os.path.join(os.path.dirname(dirrun), "tables")
    if not os.path.exists(enlace):
        os.symlink(os.path.join(os.path.dirname(bindir), "tables"), enlace)
    binario = os.path.join(bindir, "HD")
    if os.environ.get("SLURM_JOB_ID"):
        comando = ["srun", "--mpi=pmix", "--ntasks=1", binario]
    else:
        comando = [os.path.join(fase2.TOOLCHAIN, "bin", "mpirun"), "-np", "1",
                   "--oversubscribe", binario]
    timeout = int(os.environ.get("SPECTER_RUN_TIMEOUT_S", timeout))
    inicio = time.monotonic()
    proc = subprocess.run(
        comando,
        cwd=dirrun, env=fase2._entorno(), capture_output=True, text=True,
        timeout=timeout)
    if proc.returncode != 0:
        raise RuntimeError("la corrida %s falló (rc=%d):\n%s" %
                           (etiqueta, proc.returncode,
                            (proc.stdout + proc.stderr)[-3000:]))
    return dirrun, time.monotonic() - inicio


def leer_campo(dirrun, indice, nxy):
    ext = "%04d" % indice
    carpeta = os.path.join(dirrun, "out")

    def leer(nombre):
        return np.fromfile(os.path.join(carpeta, "%s.%s.out" % (nombre, ext)),
                           dtype="<f8").reshape((nxy, nxy, NZ_FIS), order="F")

    return leer("vx"), leer("vy"), leer("vz")


def diagnostico_campo(dirrun, indice, nxy, lz, nu, dk_fisico, factor_tasa):
    u, v, w = leer_campo(dirrun, indice, nxy)
    z = np.loadtxt(os.path.join(dirrun, "z.txt"))
    us, vs = u[:, :, -1], v[:, :, -1]
    fu, fv = np.fft.fft2(us), np.fft.fft2(vs)
    potencia = (np.abs(fu) ** 2 + np.abs(fv) ** 2) / float(nxy * nxy) ** 2
    kx = np.fft.fftfreq(nxy, d=1.0 / nxy)
    ky = np.fft.fftfreq(nxy, d=1.0 / nxy)
    K = np.hypot(*np.meshgrid(kx, ky, indexing="ij"))
    anillo = np.floor(K).astype(int)
    energia = np.bincount(anillo.ravel(), weights=potencia.ravel(),
                          minlength=int(anillo.max()) + 1)
    energia[0] = 0.0
    total = float(energia.sum())
    k_rms = math.sqrt(float((K ** 2 * potencia).sum()) / float(potencia.sum()))
    ub = np.trapezoid(u, z, axis=2) / lz
    vb = np.trapezoid(v, z, axis=2) / lz
    rms_sup = math.sqrt(float(np.mean(us ** 2 + vs ** 2)))
    rms_prom = math.sqrt(float(np.mean(ub ** 2 + vb ** 2)))
    rms_3d = math.sqrt(float(np.mean(u ** 2 + v ** 2 + w ** 2)))
    dzu0 = np.gradient(u, z, axis=2, edge_order=2)[:, :, 0]
    dzv0 = np.gradient(v, z, axis=2, edge_order=2)[:, :, 0]
    alfa_eff = (nu * np.mean(dzu0 * ub + dzv0 * vb)
                / (lz * np.mean(ub ** 2 + vb ** 2)))
    alfa = nu * math.pi ** 2 / (4.0 * lz ** 2)
    return {
        "indice": indice,
        "rms_superficie_codigo": rms_sup,
        "rms_promedio_codigo": rms_prom,
        "rms_3d_codigo": rms_3d,
        "sup_sobre_prom": rms_sup / rms_prom,
        "k_rms_codigo": k_rms,
        "k_rms_fisico_1_m": k_rms * dk_fisico,
        "delta": rms_prom * k_rms * lz ** 2 / nu,
        "alfa_eff_sobre_alfa": float(alfa_eff / alfa),
        "w2_sobre_v2": float(np.mean(w ** 2) / np.mean(u ** 2 + v ** 2 + w ** 2)),
        "energia_superficie_codigo": total,
        "energia_por_anillo": {str(i): float(e) for i, e in enumerate(energia) if e > 0},
        "factor_tasa_codigo_a_fisico": factor_tasa,
    }


def ajuste_exponencial(t, amplitud):
    t, amplitud = np.asarray(t), np.asarray(amplitud)
    coef = np.polyfit(t, np.log(amplitud), 1)
    residuo = np.log(amplitud) - np.polyval(coef, t)
    return -float(coef[0]), float(np.std(residuo))


def analizar_corrida(dirrun, nxy, config, tiempo_pared):
    campos = []
    for i in range(N_INTERVALOS + 1):
        indice = i + 1
        d = diagnostico_campo(
            dirrun, indice, nxy, config["Lz"], config["nu"],
            config["objetivo"]["dk_fisico_1_m"], config["factor_tasa_codigo_a_fisico"])
        d["t_codigo"] = i * config["tstep"] * DT
        d["t_fisico_s"] = d["t_codigo"] / config["factor_tasa_codigo_a_fisico"]
        campos.append(d)

    # Se descarta t=0: el primer campo usado está a 3,2 memorias, como [D-42].
    usados = campos[1:]
    tt = [x["t_codigo"] for x in usados]
    lsup, dsup = ajuste_exponencial(tt, [x["rms_superficie_codigo"] for x in usados])
    lprom, _ = ajuste_exponencial(tt, [x["rms_promedio_codigo"] for x in usados])
    tasas = []
    for n in sorted(config["pesos_usados"], key=int):
        ee = [x["energia_por_anillo"].get(str(n), 0.0) for x in usados]
        if min(ee) <= 0:
            continue
        lam, disp = ajuste_exponencial(tt, np.sqrt(ee))
        tasas.append({"anillo": int(n),
                      "k_centro_fisico_1_m": (int(n) + 0.5)
                                               * config["objetivo"]["dk_fisico_1_m"],
                      "lambda_codigo_1_t": lam,
                      "lambda_fisico_1_s": lam * config["factor_tasa_codigo_a_fisico"],
                      "dispersion_log_amplitud": disp})

    puntos = [p for p in tasas if p["anillo"] in ANILLOS_AJUSTE]
    ajuste_k2 = None
    if len(puntos) >= 3:
        kval = np.asarray([p["k_centro_fisico_1_m"] for p in puntos])
        lam = np.asarray([p["lambda_fisico_1_s"] for p in puntos])
        X = np.column_stack([np.ones(len(kval)), kval ** 2])
        coef, *_ = np.linalg.lstsq(X, lam, rcond=None)
        residuo = lam - X @ coef
        den = float((lam - lam.mean()) @ (lam - lam.mean()))
        ajuste_k2 = {"alfa_intercepto_1_s": float(coef[0]),
                     "nu_pendiente_m2_s": float(coef[1]),
                     "R2": 1.0 - float(residuo @ residuo) / den if den > 0 else None,
                     "anillos": [p["anillo"] for p in puntos]}

    balance_t, balance_e, _ = fase2.leer_balance(dirrun)
    mascara = balance_t >= config["tstep"] * DT - 1e-12
    lbal, dbal = ajuste_exponencial(balance_t[mascara], np.sqrt(balance_e[mascara]))
    diag = np.loadtxt(os.path.join(dirrun, "freeslip_diagnostic.txt"))
    return {
        "tiempo_pared_s": tiempo_pared,
        "lambda_superficie_codigo": lsup,
        "lambda_superficie_fisico_1_s": lsup * config["factor_tasa_codigo_a_fisico"],
        "dispersion_superficie": dsup,
        "lambda_promedio_codigo": lprom,
        "lambda_promedio_fisico_1_s": lprom * config["factor_tasa_codigo_a_fisico"],
        "lambda_balance_codigo": lbal,
        "lambda_balance_fisico_1_s": lbal * config["factor_tasa_codigo_a_fisico"],
        "dispersion_balance": dbal,
        "tasas_por_anillo": tasas,
        "ajuste_lambda_k2": ajuste_k2,
        "campos": campos,
        "residuo_freeslip_max": float(np.max(diag[:, -1])),
    }


def guardar(resultado):
    os.makedirs(os.path.dirname(SALIDA), exist_ok=True)
    with open(SALIDA, "w") as fh:
        json.dump(resultado, fh, indent=1, ensure_ascii=False)


def comparar(corridas):
    if not {"lineal", "experimental"}.issubset(corridas):
        return None
    li = {x["anillo"]: x for x in corridas["lineal"]["tasas_por_anillo"]}
    nl = {x["anillo"]: x for x in corridas["experimental"]["tasas_por_anillo"]}
    comun = sorted(set(li) & set(nl))
    return {
        "delta_lambda_no_lineal_menos_lineal_por_anillo": [
            {"anillo": n, "k_centro_fisico_1_m": nl[n]["k_centro_fisico_1_m"],
             "delta_lambda_1_s": (nl[n]["lambda_fisico_1_s"]
                                    - li[n]["lambda_fisico_1_s"])}
            for n in comun
        ],
        "cociente_lambda_superficie": (corridas["experimental"]["lambda_superficie_fisico_1_s"]
                                        / corridas["lineal"]["lambda_superficie_fisico_1_s"]),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--resolucion", type=int, choices=(16, 32, 64, 128), required=True)
    ap.add_argument("--prueba", action="store_true",
                    help="cadena corta: seis campos en 2 memorias; no es física")
    ap.add_argument("--solo", choices=("lineal", "experimental", "ambas"),
                    default="ambas")
    ap.add_argument("--repetir", action="store_true")
    args = ap.parse_args()

    objetivo = objetivo_desde_piv()
    modos, pesos_usados, k_rms = modos_discretos(args.resolucion, objetivo)
    dk = objetivo["dk_fisico_1_m"]
    lz = CELDA.h * dk
    nu = NU_REF * (lz / LZ_REF) ** 2
    factor_tasa = CELDA.fluido.nu * dk ** 2 / nu
    t_mem = lz ** 2 / (2.0 * math.pi ** 2 * nu)
    if args.prueba:
        tstep = int(round(0.4 * t_mem / DT))
        n_intervalos = 3
    else:
        tstep = int(round(INTERVALO_CAMPOS_MEMORIAS * t_mem / DT))
        n_intervalos = N_INTERVALOS
    # analizar_corrida usa la constante global: la prueba también necesita 6 archivos.
    # Para no duplicar caminos, la prueba escribe seis a intervalos cortos.
    if args.prueba:
        n_intervalos = N_INTERVALOS
    # SPECTER escribe al comenzar un paso que coincide con tstep, no después del paso
    # terminal. El +1 garantiza el sexto campo cuando step es múltiplo de tstep.
    pasos = n_intervalos * tstep + 1
    u0_exp = (objetivo["delta_objetivo"] * nu
              / ((2.0 / math.pi) * k_rms * lz ** 2))
    u0s = {"lineal": u0_exp * FACTOR_LINEAL, "experimental": u0_exp}
    config = {
        "resolucion": args.resolucion, "nz": NZ, "cz": CZ,
        "Lz": lz, "nu": nu, "dt": DT, "pasos": pasos, "tstep": tstep,
        "intervalo_campos_memorias": (tstep * DT / t_mem),
        "n_intervalos": n_intervalos, "t_mem_codigo": t_mem,
        "factor_tasa_codigo_a_fisico": factor_tasa,
        "bordes": {"xy": "periodicos", "z0": "noslip", "zL": "freeslip plano",
                    "presion_z": "Neumann-Neumann"},
        "objetivo": objetivo,
        "anillos_disponibles": sorted(pesos_usados),
        "pesos_usados": {str(k): float(v) for k, v in pesos_usados.items()},
        "k_rms_codigo_inicial": k_rms,
        "u0s": u0s,
        "semilla_fases": SEMILLA_FASES,
        "n_modos_independientes": len(modos),
        "prueba": bool(args.prueba),
    }
    clave = "n%d%s" % (args.resolucion, "_prueba" if args.prueba else "")
    if os.path.exists(SALIDA):
        with open(SALIDA) as fh:
            resultado = json.load(fh)
    else:
        resultado = {"descripcion": __doc__.split("\n\n")[0], "resoluciones": {}}
    previo = resultado["resoluciones"].get(clave)
    if previo and previo.get("config") != config and not args.repetir:
        raise RuntimeError("hay una salida %s con otra configuración; usar --repetir" % clave)
    if not previo or args.repetir:
        resultado["resoluciones"][clave] = {"config": config, "corridas": {}}
    bloque = resultado["resoluciones"][clave]
    etiquetas = ("lineal", "experimental") if args.solo == "ambas" else (args.solo,)
    faltan = [e for e in etiquetas if args.repetir or e not in bloque["corridas"]]
    if not faltan:
        print("%s ya está calculado: %s" % (clave, ", ".join(etiquetas)))
        return 0

    fase2.NX = fase2.NY = args.resolucion
    fase2.NZ, fase2.CZ, fase2.NZFIS = NZ, CZ, NZ_FIS
    fase2.IC_SHEAR = condicion_inicial(modos, args.resolucion)
    assert "tstep = 1000000" in fase2.PLANTILLA_INP
    fase2.PLANTILLA_INP = fase2.PLANTILLA_INP.replace(
        "tstep = 1000000", "tstep = %d" % tstep)
    compilar_liviano()

    scratch_base = os.environ.get("SPECTER_SCRATCH")
    temporal = scratch_base is None
    scratch = (tempfile.mkdtemp(prefix="banda_ancha_") if temporal
               else os.path.join(scratch_base, "banda_ancha_%s" % clave))
    os.makedirs(scratch, exist_ok=True)
    print("config %s: delta=%.3f, Lz=%.5f, nu=%.3e, k_rms=%.3f, "
          "%d modos, %d pasos" %
          (clave, objetivo["delta_objetivo"], lz, nu, k_rms, len(modos), pasos),
          flush=True)
    print("scratch:", scratch, flush=True)
    bindir = fase2.construir(os.path.join(scratch, "arbol"))
    print("compilado", flush=True)
    try:
        for etiqueta in faltan:
            print("=== %s  u0=%.6e" % (etiqueta, u0s[etiqueta]), flush=True)
            dirrun, pared = correr_largo(
                bindir, "%s_%s" % (clave, etiqueta), lz=lz, nu=nu,
                u0=u0s[etiqueta], pasos=pasos)
            bloque["corridas"][etiqueta] = analizar_corrida(
                dirrun, args.resolucion, config, pared)
            bloque["comparacion"] = comparar(bloque["corridas"])
            guardar(resultado)
            r = bloque["corridas"][etiqueta]
            print("    pared %.1f min; lambda_sup=%.5f 1/s; ajuste k2=%s" %
                  (pared / 60.0, r["lambda_superficie_fisico_1_s"],
                   r["ajuste_lambda_k2"]), flush=True)
    finally:
        if temporal and os.environ.get("SPECTER_KEEP") != "1":
            shutil.rmtree(scratch, ignore_errors=True)

    print("escrito", os.path.relpath(SALIDA, RAIZ))
    return 0


if __name__ == "__main__":
    sys.exit(main())
