#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Etapa 1 del barrido en delta: ¿se aparta la tasa de decaimiento de la predicción
lineal cuando el término no lineal deja de ser chico?

QUÉ CONTESTA
============
[P-01] dejó dos flecos: el factor geométrico O(1) de la estimación de la distorsión del
perfil vertical, que puse en 1 sin calcularlo, y si el sesgo sobre alfa es de orden delta
o delta^2. Los dos se cierran midiendo cuánto se aparta la tasa medida de la lineal a
medida que sube la amplitud.

CÓMO, Y POR QUÉ ASÍ
===================
La referencia lineal NO se escribe a mano: se MIDE, corriendo el mismo binario, la misma
malla y la misma condición inicial a amplitud minúscula, donde el término no lineal es
despreciable por construcción. Así la referencia incluye automáticamente el k discreto de
la malla, el sesgo del integrador y el filtro de dealiasing, y el cociente
lambda_eff/lambda_lin es un número puro: cualquier apartamiento es no linealidad y nada
más. Escribir la fórmula analítica habría mezclado el efecto buscado con los errores de
discretización.

Dos resoluciones horizontales, y no es opcional: [H-09] mostró que con superficie libre el
término no lineal se desestabiliza sobre mallas donde el canal aguanta. Un apartamiento
que NO converge con la resolución es numérico, no físico, y ese es el confundido que hay
que descartar.

Se usa la condición inicial vparam1=1 de la puerta de la Fase 2 (un modo solenoidal con
dependencia en x), sin escribir Fortran nuevo. No cumple el no-deslizamiento del fondo, así
que hay un transitorio inicial; se descarta y se ajusta sobre una ventana temprana, con
el flujo todavía fuerte, que es donde delta es grande y donde habría algo que ver.

Uso:
    /home/lucia/miniforge3/envs/piv-dt/bin/python verificacion/barrido_delta_etapa1.py
"""

import importlib.util
import json
import math
import os
import sys

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(AQUI)

_spec = importlib.util.spec_from_file_location(
    "fase2", os.path.join(AQUI, "test_aceptacion_fase2.py"))
fase2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fase2)

SALIDA = os.path.join(RAIZ, "salidas", "tablas", "barrido_delta_etapa1.json")

# --- geometría y física de las corridas ------------------------------------------
# Los números de onda de SPECTER son ENTEROS (specter.fpp:777-780, ky(j)=j-1), o sea
# que el Lx=1 del archivo de parámetros son 2*pi de largo físico y el modo fundamental
# es k=1. Verificado comparando <v^2> del código contra la fórmula analítica de la
# condición inicial: 4.750e-3 contra 4.750e-3 con u0=0.1. Por eso Lz NO lleva 2*pi:
# es la longitud literal.
LZ = 0.5               # eps = k*h = 0.5, prácticamente el 0.53 de la celda medida
NU = 1.0e-2
NZ, CZ = 64, 25        # 39 puntos físicos en z: de sobra para sin(pi z / 2h)
DT = 5.0e-4            # ~1/6 del límite de estabilidad viscosa explícita
RESOLUCIONES = (16, 32)
# delta ~ 15.5*u0 para esta condición inicial: el barrido cubre 0.5 a 15, que es el
# rango que tiene el experimento (de ~20 al arranque a ~1 al final del registro)
U0S = (0.03, 0.1, 0.2, 0.4, 0.7, 1.0)
U0_REFERENCIA = 1.0e-4  # referencia lineal: delta ~ 1.5e-3
VPARAM1 = 1.0

# --- condición inicial propia -----------------------------------------------------
# La de la puerta NO sirve acá y conviene decir por qué: es psi ~ cos(kx x) sin(pi z/Lz)
# en el plano (x,z), que es autofunción del laplaciano 2D, de modo que J(psi, lap psi)
# se anula idénticamente y el término no lineal es CERO. Además, en un flujo 2D en
# (x,z) la continuidad integrada en z da d/dx <u> = -[w] = 0, o sea que el promedio
# vertical no puede tener estructura horizontal, que es justo la variable de la
# clausura. La prueba tiene que ser tridimensional.
#
# Se usa entonces  u_h = F(z) grad_perp(psi(x,y)),  w = 0,  con
#     F(z)  = sin(pi z / 2Lz)          el modo fundamental del par rígido-libre
#     psi   = cos(k x)cos(k y) + b cos(2k x)cos(k y)
# El campo es horizontalmente solenoidal por construcción (u = d psi/dy, v = -d psi/dx)
# y tiene w = 0 exacto, así que cumple impermeabilidad en las dos caras. La suma de dos
# modos con |k|^2 distinto (2k^2 y 5k^2) rompe la degeneración: psi ya no es autofunción
# del laplaciano horizontal y J(psi, lap psi) no se anula.
B_SEGUNDO_MODO = 0.6

IC_CELULAR = r"""
! Condicion inicial celular 3D para el barrido en delta.
!   u_h = F(z) grad_perp(psi),  w = 0
!   F(z) = sin(pi z / 2 Lz);  psi = cos(k x)cos(k y) + b cos(2k x)cos(k y)
      vx = 0; vy = 0; vz = 0

      IF (ista .eq. 1) THEN
         DO k = 1,nz-Cz
            rmt = SIN(0.5_GP*pi*z(k)/Lz)
            ! modo 1: kx(2), +-ky(2)
            rmq = 0.25_GP*u0*nx*ny
            vx(k,2,2)  = vx(k,2,2)  + im*ky(2)*rmq*rmt
            vy(k,2,2)  = vy(k,2,2)  - im*kx(2)*rmq*rmt
            vx(k,ny,2) = vx(k,ny,2) - im*ky(2)*rmq*rmt
            vy(k,ny,2) = vy(k,ny,2) - im*kx(2)*rmq*rmt
            ! modo 2: kx(3), +-ky(2), con amplitud relativa b
            rmq = 0.25_GP*BSEG*u0*nx*ny
            vx(k,2,3)  = vx(k,2,3)  + im*ky(2)*rmq*rmt
            vy(k,2,3)  = vy(k,2,3)  - im*kx(3)*rmq*rmt
            vx(k,ny,3) = vx(k,ny,3) - im*ky(2)*rmq*rmt
            vy(k,ny,3) = vy(k,ny,3) - im*kx(3)*rmq*rmt
         ENDDO
      ENDIF

      CALL fftp1d_real_to_complex_z(planfc,vx,MPI_COMM_WORLD)
      CALL fftp1d_real_to_complex_z(planfc,vy,MPI_COMM_WORLD)
      CALL fftp1d_real_to_complex_z(planfc,vz,MPI_COMM_WORLD)
"""

# La ventana de ajuste va TEMPRANO y es corta, y eso es una decisión de diseño, no
# una economía. Con esta condición inicial la energía cae unas dos décadas por unidad
# de tiempo; si se ajustara sobre la cola, el ajuste quedaría a amplitud minúscula, o
# sea en el régimen LINEAL, que es justo donde no hay nada que medir. Se ajusta
# entonces sobre [0,25 T, 0,70 T]: después del transitorio de arranque —la condición
# inicial no cumple el no-deslizamiento del fondo y el código se lo impone por
# inyección— y con el flujo todavía fuerte. Dentro de esa ventana delta cae un factor
# ~2,5, así que se informa su valor en el medio y también el rango.
# La condición inicial ES el modo fundamental del par rígido-libre, así que en el caso
# lineal no hay transitorio: decae como una exponencial pura desde t=0. La ventana no
# necesita ser tardía, y conviene que sea corta para que delta no cambie mucho adentro:
# en [0.25 T, 0.70 T] la energía cae un 25 % y delta un 16 %.
T_FINAL = 3.0
PASOS = int(round(T_FINAL / DT))
CSTEP = 30
VENT = (0.25, 0.70)


def tasa_en_ventana(t, e, t0, t1):
    """lambda del ajuste de log(<v^2>) ~ -2*lambda*t en [t0, t1], y su dispersión."""
    m = (t >= t0) & (t <= t1) & (e > 0)
    if m.sum() < 5:
        return float("nan"), float("nan")
    coef = np.polyfit(t[m], np.log(e[m]), 1)
    disp = float(np.std(np.log(e[m]) - np.polyval(coef, t[m])))
    return -float(coef[0]) / 2.0, disp


def _compilar_liviano():
    """Baja la optimización de gfortran de -O2 a -O1 sólo para estas corridas.

    No es una decisión física: compilar SPECTER con -O2 es lo que hace picar la memoria
    en esta máquina, y ya mató dos corridas. Con -O1 el binario es algo más lento pero
    idéntico en resultados, y estas corridas son chicas. La puerta de la Fase 2 sigue
    compilando con -O2, porque no se la toca.
    """
    orig = fase2._parchear_makefile_in

    def envoltura(ruta, solver="HD"):
        orig(ruta, solver)
        with open(ruta) as f:
            txt = f.read()
        with open(ruta, "w") as f:
            f.write(txt.replace("-O2", "-O1"))

    fase2._parchear_makefile_in = envoltura


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--resolucion", type=int, default=None,
                    help="correr una sola malla, para no tener dos compilaciones "
                         "vivas en la misma sesión de memoria")
    args = ap.parse_args()
    resoluciones = (args.resolucion,) if args.resolucion else RESOLUCIONES
    _compilar_liviano()

    scratch = os.environ.get("SPECTER_SCRATCH")
    temporal = scratch is None
    if temporal:
        import tempfile
        scratch = tempfile.mkdtemp(prefix="delta1_")
    print("scratch: %s" % scratch, flush=True)
    print("Lz=%.2f  nu=%.1e  dt=%.1e  pasos=%d  t_final=%.3f  ventana=[%.3f, %.3f]"
          % (LZ, NU, DT, PASOS, T_FINAL, VENT[0]*T_FINAL, VENT[1]*T_FINAL), flush=True)

    # la malla vertical es la misma en todas las corridas
    fase2.NZ, fase2.CZ = NZ, CZ
    fase2.NZFIS = NZ - CZ
    # se reemplaza la condición inicial de la puerta por la celular 3D
    # rmt y rmq ya están declaradas REAL(KIND=GP) en el ámbito donde se
    # incluye initialv.f90 (specter.fpp:104-105), así que se usan tal cual.
    fase2.IC_SHEAR = IC_CELULAR.replace("BSEG", "%.4f_GP" % B_SEGUNDO_MODO)

    res = {"Lz": LZ, "nu": NU, "dt": DT, "pasos": PASOS, "t_final": T_FINAL,
           "nz": NZ, "cz": CZ, "ord": fase2.ORD, "u0s": list(U0S),
           "u0_referencia": U0_REFERENCIA, "resoluciones": list(RESOLUCIONES),
           "corridas": {}}

    # Se ANEXA sobre lo ya calculado y se guarda DESPUÉS DE CADA MALLA. Cada malla
    # cuesta varios minutos y esta máquina ya mató dos corridas por falta de memoria;
    # perder una malla entera porque el guardado estaba al final es un defecto.
    if os.path.exists(SALIDA):
        try:
            with open(SALIDA) as fh:
                previo = json.load(fh)
            if all(previo.get(c) == res[c] for c in
                   ("Lz", "nu", "dt", "pasos", "nz", "cz", "u0s")):
                res["corridas"].update(previo.get("corridas", {}))
                res.setdefault("lambda_referencia", {}).update(
                    previo.get("lambda_referencia", {}))
                print("anexando sobre: %s" % ", ".join(sorted(res["corridas"])),
                      flush=True)
            else:
                print("[aviso] el JSON previo es de otros parámetros; se ignora",
                      flush=True)
        except (ValueError, OSError) as exc:
            print("[aviso] no se pudo leer %s (%s)" % (SALIDA, exc), flush=True)

    for nxy in resoluciones:
        if "n%d" % nxy in res["corridas"]:
            print("\n=== malla %dx%d ya calculada, se saltea" % (nxy, nxy), flush=True)
            continue
        fase2.NX = fase2.NY = nxy
        bindir = fase2.construir(os.path.join(scratch, "arbol_%d" % nxy))
        print("\n=== malla %dx%dx%d compilada" % (nxy, nxy, NZ), flush=True)

        filas = []
        for u0 in (U0_REFERENCIA,) + U0S:
            etiqueta = "n%d_u%.0e" % (nxy, u0)
            _, _, dirrun = fase2.correr(
                bindir, etiqueta, bcsta="noslip", bcend=fase2.CADENA_LIBRE,
                vparam1=VPARAM1, dt=DT, step=PASOS, cstep=CSTEP, lz=LZ, nu=NU, u0=u0)
            t, e, w = fase2.leer_balance(dirrun)

            t0, t1 = VENT[0] * T_FINAL, VENT[1] * T_FINAL
            lam, disp = tasa_en_ventana(t, e, t0, t1)
            # k efectivo y amplitud en el MEDIO de la ventana. <w^2>/<v^2> es un
            # cociente, así que no depende de la normalización del dominio extendido.
            tm = 0.5 * (t0 + t1)
            im = int(np.argmin(np.abs(t - tm)))
            i0 = int(np.argmin(np.abs(t - t0)))
            i1 = int(np.argmin(np.abs(t - t1)))
            k_ef = float(np.sqrt(w[im] / e[im]))
            filas.append({
                "u0": u0, "lambda": lam, "dispersion": disp,
                "k_efectivo": k_ef,
                "v2_inicial": float(e[0]), "v2_medio": float(e[im]),
                "v2_ventana": [float(e[i0]), float(e[i1])],
                "v2_final": float(e[-1]),
                "crece": bool(e[-1] > e[0]),
            })
            print("   u0=%.1e  lambda=%.6e  disp=%.2e  k_ef=%.2f  "
                  "<v2>: %.3e -> %.3e%s"
                  % (u0, lam, disp, k_ef, e[0], e[-1],
                     "   [CRECE: inestable]" if e[-1] > e[0] else ""), flush=True)
            del t, e, w

        ref = filas[0]["lambda"]
        for f in filas:
            f["lambda_sobre_referencia"] = f["lambda"] / ref
            # delta a menos de una constante de forma O(1): U ~ sqrt(<v^2>) en la
            # ventana de ajuste, k del propio campo, h = Lz
            # delta CALIBRADO. Para esta condición inicial <v^2> = <F^2>|grad_perp
            # psi|^2_rms = (1/2)|grad_perp psi|^2_rms, y la velocidad que entra en la
            # clausura es el PROMEDIO VERTICAL, U = <F>|grad_perp psi|_rms con
            # <F> = 2/pi. De ahí U = (2/pi) sqrt(2 <v^2>), y delta = U k h^2 / nu.
            def _delta(v2):
                U = (2.0 / math.pi) * math.sqrt(2.0 * max(v2, 0.0))
                return U * f["k_efectivo"] * LZ ** 2 / NU
            f["delta"] = _delta(f["v2_medio"])
            f["delta_rango"] = [_delta(v) for v in f["v2_ventana"]]
        res["corridas"]["n%d" % nxy] = filas
        res.setdefault("lambda_referencia", {})["n%d" % nxy] = ref
        os.makedirs(os.path.dirname(SALIDA), exist_ok=True)
        with open(SALIDA, "w") as f:
            json.dump(res, f, indent=1)
        print("   guardado parcial (%s)" % ", ".join(sorted(res["corridas"])),
              flush=True)

    os.makedirs(os.path.dirname(SALIDA), exist_ok=True)
    with open(SALIDA, "w") as f:
        json.dump(res, f, indent=1)

    if len(res["corridas"]) < len(RESOLUCIONES):
        print("\nfalta(n) la(s) malla(s): %s"
              % ", ".join(str(n) for n in RESOLUCIONES
                          if "n%d" % n not in res["corridas"]))
        print("escrito %s" % os.path.relpath(SALIDA, RAIZ))
        return 0

    print("\nRESUMEN  lambda/lambda_lineal")
    print("  u0        " + "".join("  %dx%d      " % (n, n) for n in RESOLUCIONES))
    for i, u0 in enumerate(U0S):
        fila = "  %.1e " % u0
        for n in RESOLUCIONES:
            f = res["corridas"]["n%d" % n][i + 1]
            fila += "  %9.6f%s" % (f["lambda_sobre_referencia"],
                                   "*" if f["crece"] else " ")
        print(fila)
    print("  (* = la energía creció: corrida inestable, no usable)")
    print("\nescrito %s" % os.path.relpath(SALIDA, RAIZ))

    if temporal and os.environ.get("SPECTER_KEEP") != "1":
        import shutil
        shutil.rmtree(scratch, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
