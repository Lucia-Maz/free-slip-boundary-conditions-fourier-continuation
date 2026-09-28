#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Campos de superficie con la superficie deformable de orden N = 2 ([D-51], [D-52]) a la
escala física de la celda.

La corrida es la de [V-19] con otro tope:
- condición inicial de banda ancha, con el espectro del PIV en t ≈ 11 s y amplitud
  experimental;
- fondo no deslizante;
- tope `freesurface` con fsorder = 2, en vez de free-slip plano;
- g y σ/ρ físicos pasados a unidades de código con la misma escala que [V-19]. Los valores
  salen de salidas/tablas/superficie_deformable_escalas.json, con su procedencia.

La superficie arranca plana (fsamp = 0), así que al principio se ajusta al campo de presión
del flujo e irradia ondas de gravedad-capilaridad débilmente amortiguadas. Guarda, en pocos
instantes, la elevación η(x,y) y la velocidad en la superficie, en unidades físicas:

  salidas/tablas/fase6_campos_superficie_n{N}_m{M}.json   configuración, escalas, series
  salidas/campos/fase6_campos_superficie_n{N}_m{M}.npz    η, u, v y ω_z en z = h
(N la malla horizontal, M la duración en memorias modales).

Sólo numpy y json: corre en Sakura con el módulo python/3.13.13. La figura la hace
numerico/fase6/figura_campos_superficie.py, en la laptop.

Uso (laptop, secuencial, con ulimit -v 2500000; en Sakura, vía SLURM con
numerico/fase6/campos_superficie.sbatch):

  $PY numerico/fase6/campos_superficie.py --nxy 16 --prueba          # cadena, un minuto
  $PY numerico/fase6/campos_superficie.py --nxy 16 --comparar-np     # np=1 contra np=2
  $PY numerico/fase6/campos_superficie.py --nxy 32 --memorias 1.5    # la de respaldo
  $PY numerico/fase6/campos_superficie.py --nxy 64 --memorias 1.5    # la de Sakura
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
RAIZ = os.path.dirname(os.path.dirname(AQUI))
VERIF = os.path.join(RAIZ, "verificacion")
sys.path.insert(0, VERIF)
sys.path.insert(0, os.path.join(RAIZ, "codigo"))


def _modulo(nombre, ruta):
    spec = importlib.util.spec_from_file_location(nombre, ruta)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


BA = _modulo("banda_ancha", os.path.join(VERIF, "decaimiento_banda_ancha.py"))
F2 = BA.fase2
H6 = _modulo("puerta6", os.path.join(VERIF, "test_aceptacion_fase6.py"))  # sólo la plantilla
CELDA = BA.CELDA
ESCALAS = os.path.join(RAIZ, "salidas", "tablas", "superficie_deformable_escalas.json")
SAL_J = os.path.join(RAIZ, "salidas", "tablas", "fase6_campos_superficie_n%d_m%g%s.json")
SAL_N = os.path.join(RAIZ, "salidas", "campos", "fase6_campos_superficie_n%d_m%g%s.npz")

NZ, CZ = BA.NZ, BA.CZ
NZ_FIS = NZ - CZ
DT = BA.DT
CSTEP = 20


def escalas(nxy):
    """Unidades de código de [V-19] y los parámetros físicos de la superficie en ellas."""
    objetivo = BA.objetivo_desde_piv()
    modos, pesos, k_rms = BA.modos_discretos(nxy, objetivo)
    dk = objetivo["dk_fisico_1_m"]            # 1/L_c: k de código entero <-> k_fis = k dk
    lz = CELDA.h * dk
    nu = BA.NU_REF * (lz / BA.LZ_REF) ** 2
    t_c = nu / (CELDA.fluido.nu * dk ** 2)    # segundos por unidad de tiempo de código
    e = json.load(open(ESCALAS))
    g, s_rho = float(e["g_m_s2"]), float(e["sigma_sobre_rho_m3_s2"])
    g_c = g * t_c ** 2 * dk
    s_c = s_rho * t_c ** 2 * dk ** 3
    u0 = (objetivo["delta_objetivo"] * nu / ((2.0 / math.pi) * k_rms * lz ** 2))
    t_mem = lz ** 2 / (2.0 * math.pi ** 2 * nu)
    return {"objetivo": objetivo, "modos": modos, "k_rms_codigo": k_rms, "dk_1_m": dk,
            "Lz": lz, "nu": nu, "t_c_s": t_c, "L_c_m": 1.0 / dk, "g_codigo": g_c,
            "sigma_rho_codigo": s_c, "g_m_s2": g, "sigma_rho_m3_s2": s_rho,
            "procedencia_g_sigma": e.get("procedencia_g_rho_sigma", ESCALAS),
            "u0": u0, "t_mem_codigo": t_mem}


def estabilidad(esc, nxy, dt):
    """Chequeo previo del paso: difusión explícita y ondas de superficie con RK2."""
    dz = esc["Lz"] / (NZ_FIS - 1)
    D = esc["nu"] * dt / dz ** 2
    peor = None
    for k in range(1, nxy // 2):
        w2 = (esc["g_codigo"] * k + esc["sigma_rho_codigo"] * k ** 3) * math.tanh(k * esc["Lz"])
        x = math.sqrt(w2) * dt
        crece = x ** 4 / 8.0                  # |G_RK2|-1 para una oscilación pura
        amort = 2.0 * esc["nu"] * k ** 2 * dt   # amortiguamiento viscoso de la onda, 2 nu k^2
        c = crece / amort
        if peor is None or c > peor["cociente"]:
            peor = {"k": k, "omega_dt": x, "crecimiento_rk2": crece, "amortiguamiento": amort,
                    "cociente": c}
    return {"nu_dt_dz2": D, "ondas_peor_caso": peor}


def condicion_inicial_mpi(modos, nxy):
    """La de [V-19] (BA.condicion_inicial) pero correcta con varios procesos MPI: aquélla
    escribe todos los modos dentro de IF (ista .eq. 1), con índices globales de kx, así que
    con más de un proceso los modos de los otros se pierden o pisan memoria. Acá cada modo
    lo escribe el proceso que tiene su columna kx."""
    lineas = [
        "! Banda ancha derivada del PIV en t ~= 11.1 s ([V-18]); version MPI.",
        "! u_h=sin(pi*z/2Lz)*grad_perp(psi), w=0; fases deterministas.",
        "      vx = 0; vy = 0; vz = 0",
        "",
        "      DO k = 1,nz-Cz",
        "         rmt = SIN(0.5_GP*pi*z(k)/Lz)",
    ]
    for m in modos:
        ix = m["kx"] + 1
        iy = m["ky"] + 1 if m["ky"] >= 0 else nxy + m["ky"] + 1
        coef = 0.5 * m["amplitud_psi"]
        c, s = math.cos(m["fase"]), math.sin(m["fase"])
        lineas += [
            "         IF ( %d .ge. ista .AND. %d .le. iend ) THEN" % (ix, ix),
            "            rmq = %.16e_GP*u0*nx*ny" % coef,
            ("            vx(k,%d,%d) = vx(k,%d,%d) + im*ky(%d)*rmq*"
             "(%.16e_GP+im*%.16e_GP)*rmt") % (iy, ix, iy, ix, iy, c, s),
            ("            vy(k,%d,%d) = vy(k,%d,%d) - im*kx(%d)*rmq*"
             "(%.16e_GP+im*%.16e_GP)*rmt") % (iy, ix, iy, ix, ix, c, s),
            "         ENDIF",
        ]
    lineas += [
        "      ENDDO",
        "",
        "      CALL fftp1d_real_to_complex_z(planfc,vx,MPI_COMM_WORLD)",
        "      CALL fftp1d_real_to_complex_z(planfc,vy,MPI_COMM_WORLD)",
        "      CALL fftp1d_real_to_complex_z(planfc,vz,MPI_COMM_WORLD)",
        "",
    ]
    return "\n".join(lineas)


def construir(destino, nxy, esc):
    if os.path.exists(destino):
        shutil.rmtree(destino)
    shutil.copytree(os.path.join(RAIZ, "numerico", "SPECTER-trabajo"), destino,
                    ignore=shutil.ignore_patterns(".git", "*.o", "*.mod"))
    F2.NX = F2.NY = nxy
    F2.NZ, F2.CZ, F2.NZFIS = NZ, CZ, NZ_FIS
    F2._parchear_makefile_in(os.path.join(destino, "src", "Makefile.in"))
    if os.environ.get("SPECTER_O1"):          # la laptop: gfortran -O2 pica en memoria
        ruta = os.path.join(destino, "src", "Makefile.in")
        texto = open(ruta).read().replace("-O2", "-O1")
        open(ruta, "w").write(texto)
    with open(os.path.join(destino, "src", "initialv.f90"), "w") as fh:
        fh.write(condicion_inicial_mpi(esc["modos"], nxy))
    src = os.path.join(destino, "src")
    env = F2._entorno()
    env["PWD"] = src                          # GHOME = $(PWD), ver CLAUDE.md
    log = subprocess.run(["make", "main"], cwd=src, env=env, capture_output=True, text=True)
    binario = os.path.join(destino, "bin", "HD")
    if log.returncode != 0 or not os.path.exists(binario):
        raise RuntimeError("la compilación falló:\n%s" % (log.stdout + log.stderr)[-3000:])
    return os.path.join(destino, "bin")


def correr(bindir, etiqueta, esc, *, pasos, tstep, nproc, fsorder=2):
    dirrun = os.path.join(os.path.dirname(bindir), "corridas", etiqueta)
    if os.path.exists(dirrun):
        shutil.rmtree(dirrun)
    os.makedirs(os.path.join(dirrun, "out"))
    os.makedirs(os.path.join(dirrun, "in"))
    plantilla = H6.PLANTILLA.replace("tstep = 1000000", "tstep = %d" % tstep)
    assert "tstep = %d" % tstep in plantilla
    with open(os.path.join(dirrun, "parameter.inp"), "w") as fh:
        fh.write(plantilla.format(
            lz=esc["Lz"], dt=DT, step=pasos, cstep=CSTEP, u0=esc["u0"], nu=esc["nu"],
            vparam0=0.0, vparam1=0.0, bcsta="noslip", bcend="freesurface",
            fsgrav=esc["g_codigo"], fstens=esc["sigma_rho_codigo"], fsamp=0.0, fskx=1,
            fsky=0, fscoup=2, fsorder=fsorder, fsmean=0.0, fstol=1e-10, fsmaxit=30, fsgen=0))
    enlace = os.path.join(os.path.dirname(dirrun), "tables")
    if not os.path.exists(enlace):
        os.symlink(os.path.join(os.path.dirname(bindir), "tables"), enlace)
    binario = os.path.join(bindir, "HD")
    if os.environ.get("SLURM_JOB_ID"):
        comando = ["srun", "--mpi=pmix", "--ntasks=%d" % nproc, binario]
    else:
        comando = [os.path.join(F2.TOOLCHAIN, "bin", "mpirun"), "-np", str(nproc),
                   "--oversubscribe", binario]
    t0 = time.monotonic()
    proc = subprocess.run(comando, cwd=dirrun, env=F2._entorno(), capture_output=True,
                          text=True, timeout=int(os.environ.get("SPECTER_RUN_TIMEOUT_S", 86400)))
    avisos = [l.strip() for l in (proc.stdout + proc.stderr).splitlines()
              if "WARNING" in l or "ERROR" in l]
    if proc.returncode != 0:
        raise RuntimeError("la corrida %s falló (rc=%d):\n%s" % (
            etiqueta, proc.returncode, (proc.stdout + proc.stderr)[-3000:]))
    return dirrun, time.monotonic() - t0, avisos


def leer_eta(dirrun, indice, nxy):
    """fs_eta.NNNN.dat: los coeficientes (ky,kx) con la normalización nx*ny (v_fs_output)."""
    ruta = os.path.join(dirrun, "out", "fs_eta.%04d.dat" % indice)
    c = np.fromfile(ruta, dtype="<c16").reshape((nxy // 2 + 1, nxy)).T   # (ky, kx)
    return np.fft.irfft2(c, s=(nxy, nxy)).T                                # (x, y)


def extraer(dirrun, esc, nxy, n_campos, tstep):
    """Planos de superficie en unidades físicas: η [m], u, v [m/s], ω_z [1/s]."""
    L, T = esc["L_c_m"], esc["t_c_s"]
    kk = np.fft.fftfreq(nxy, d=1.0 / nxy)
    KX, KY = np.meshgrid(kk, kk, indexing="ij")
    out = {"t_codigo": [], "t_s": [], "eta_m": [], "u_m_s": [], "v_m_s": [], "omega_z_1_s": []}
    for i in range(1, n_campos + 2):
        if not os.path.exists(os.path.join(dirrun, "out", "vx.%04d.out" % i)):
            break
        u, v, _ = BA.leer_campo(dirrun, i, nxy)
        us, vs = u[:, :, -1], v[:, :, -1]
        wz = np.real(np.fft.ifft2(1j * KX * np.fft.fft2(vs) - 1j * KY * np.fft.fft2(us)))
        t = (i - 1) * tstep * DT
        out["t_codigo"].append(t)
        out["t_s"].append(t * T)
        out["eta_m"].append(leer_eta(dirrun, i, nxy) * L)
        out["u_m_s"].append(us * L / T)
        out["v_m_s"].append(vs * L / T)
        out["omega_z_1_s"].append(wz / T)
    return {k: np.asarray(v) for k, v in out.items()}


def resumen(campos, dirrun, esc):
    d = np.atleast_2d(np.loadtxt(os.path.join(dirrun, "freesurface_order_diagnostic.txt")))
    dz = esc["Lz"] / (NZ_FIS - 1)
    instantes = []
    for j in range(len(campos["t_s"])):
        eta = campos["eta_m"][j]
        instantes.append({
            "t_s": float(campos["t_s"][j]),
            "eta_rms_m": float(np.sqrt(np.mean(eta ** 2))),
            "eta_max_abs_m": float(np.max(np.abs(eta))),
            "eta_medio_m": float(np.mean(eta)),
            "u_rms_superficie_m_s": float(np.sqrt(np.mean(campos["u_m_s"][j] ** 2
                                                          + campos["v_m_s"][j] ** 2))),
            "omega_z_rms_1_s": float(np.sqrt(np.mean(campos["omega_z_1_s"][j] ** 2))),
            "max_eta_sobre_dz": float(np.max(np.abs(eta)) / (dz * esc["L_c_m"])),
        })
    serie = {"t_s": (d[:, 0] * esc["t_c_s"]).tolist(),
             "eta_rms_m": (d[:, 8] * esc["L_c_m"]).tolist(),
             "pasadas": d[:, 10].tolist(), "residuo": d[:, 11].tolist(),
             "max_eta_sobre_dz": d[:, 12].tolist()}
    return instantes, serie


def comparar_np(args, esc):
    """El camino general con 1 y 2 procesos MPI: η y el diagnóstico tienen que coincidir."""
    scratch = tempfile.mkdtemp(prefix="fase6_np_")
    try:
        b = construir(os.path.join(scratch, "m"), args.nxy, esc)
        res = {}
        for n in (1, 2):
            d, s, _ = correr(b, "np%d" % n, esc, pasos=101, tstep=50, nproc=n)
            res[n] = (np.loadtxt(os.path.join(d, "freesurface_order_diagnostic.txt")),
                      leer_eta(d, 3, args.nxy))
            print("np=%d: %.0f s" % (n, s), flush=True)
        dd = float(np.max(np.abs(res[1][0] - res[2][0])[:, [1, 2, 7, 8]]))
        de = float(np.max(np.abs(res[1][1] - res[2][1])) / np.max(np.abs(res[1][1])))
        print("np=1 contra np=2: diagnóstico %.2e (absoluto), eta %.2e (relativo)" % (dd, de))
        return de < 1e-10
    finally:
        shutil.rmtree(scratch)


def main():
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("--nxy", type=int, choices=(16, 32, 64, 128), required=True)
    ap.add_argument("--memorias", type=float, default=1.5,
                    help="duración en memorias modales h^2/(2 pi^2 nu)")
    ap.add_argument("--campos", type=int, default=6, help="instantes guardados, además del 0")
    ap.add_argument("--np", type=int, default=int(os.environ.get("SLURM_NTASKS", 1)))
    ap.add_argument("--prueba", action="store_true", help="cadena corta; no es física")
    ap.add_argument("--comparar-np", action="store_true")
    args = ap.parse_args()

    esc = escalas(args.nxy)
    est = estabilidad(esc, args.nxy, DT)
    print("escalas: L_c = %.3f mm, t_c = %.3f s, g = %.1f, sigma/rho = %.3f (código); "
          "nu dt/dz^2 = %.3f; ondas: crecimiento RK2 / amortiguamiento = %.2f en k = %d" % (
              esc["L_c_m"] * 1e3, esc["t_c_s"], esc["g_codigo"], esc["sigma_rho_codigo"],
              est["nu_dt_dz2"], est["ondas_peor_caso"]["cociente"], est["ondas_peor_caso"]["k"]),
          flush=True)
    if est["ondas_peor_caso"]["cociente"] > 0.5 or est["nu_dt_dz2"] > 0.19:
        raise RuntimeError("el paso no es estable con margen: %s" % est)
    if args.comparar_np:
        return 0 if comparar_np(args, esc) else 1

    memorias = 0.2 if args.prueba else args.memorias
    campos = 3 if args.prueba else args.campos
    pasos_tot = int(round(memorias * esc["t_mem_codigo"] / DT))
    tstep = max(1, pasos_tot // campos)
    pasos = campos * tstep + 1                # SPECTER escribe al empezar el paso
    scratch_base = os.environ.get("SPECTER_SCRATCH")
    scratch = scratch_base or tempfile.mkdtemp(prefix="fase6_campos_")
    os.makedirs(scratch, exist_ok=True)
    t0 = time.monotonic()
    b = construir(os.path.join(scratch, "m%d" % args.nxy), args.nxy, esc)
    t_comp = time.monotonic() - t0
    dirrun, t_run, avisos = correr(b, "campos", esc, pasos=pasos, tstep=tstep, nproc=args.np)
    campos_sup = extraer(dirrun, esc, args.nxy, campos, tstep)
    instantes, serie = resumen(campos_sup, dirrun, esc)
    config = {k: v for k, v in esc.items() if k not in ("modos", "objetivo")}
    config.update({"nxy": args.nxy, "nz": NZ, "cz": CZ, "dt": DT, "pasos": pasos,
                   "tstep": tstep, "memorias": memorias, "fsorder": 2, "np": args.np,
                   "prueba": bool(args.prueba), "estabilidad": est,
                   "delta_objetivo": esc["objetivo"]["delta_objetivo"],
                   "bordes": {"xy": "periódicos", "z0": "noslip", "zL": "freesurface, N = 2"}})
    resultado = {"descripcion": __doc__.strip().splitlines()[0], "config": config,
                 "tiempo_compilacion_s": t_comp, "tiempo_corrida_s": t_run,
                 "avisos": avisos, "instantes": instantes, "serie_diagnostico": serie}
    sufijo = "_prueba" if args.prueba else ""
    nombre = (args.nxy, memorias, sufijo)
    os.makedirs(os.path.dirname(SAL_N % nombre), exist_ok=True)
    np.savez_compressed(SAL_N % nombre, **campos_sup)
    json.dump(resultado, open(SAL_J % nombre, "w"), indent=1)
    print("->", SAL_J % nombre)
    print("corrida: %.0f s; %d campos; eta_rms final %.3e m; pasadas <= %d; residuo <= %.1e; "
          "avisos: %d" % (t_run, len(instantes), instantes[-1]["eta_rms_m"],
                          int(max(serie["pasadas"])), max(serie["residuo"]), len(avisos)))
    if not scratch_base:
        shutil.rmtree(scratch)
    return 0


if __name__ == "__main__":
    sys.exit(main())
