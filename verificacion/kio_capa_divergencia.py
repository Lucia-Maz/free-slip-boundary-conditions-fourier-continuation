#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
El análisis de capa límite de divergencia de Karniadakis, Israeli y Orszag (1991),
sección 2.3, aplicado al esquema de Fontana et al. (2020) que usa SPECTER.

QUÉ SE ESTÁ PONIENDO A PRUEBA
=============================
KIO §2.3 muestra que en su esquema —viscosidad IMPLÍCITA— la divergencia obedece

    Q - c·nu·dt·grad^2(Q) = 0 ,     Q = div(v) ,

o sea un operador de Helmholtz cuyas soluciones homogéneas decaen como exp(-s/l) con
l = sqrt(c·nu·dt): esa es la capa límite numérica. El valor de Q en la pared lo fija el
error de la condición de contorno de la presión, y el error de velocidad resulta un orden
en dt MENOR que el de la divergencia de borde.

En el esquema de Fontana el término viscoso es EXPLÍCITO —verificado en
`src/include/hd/hd_rkstep2.f90`: vx = C1 + dt*(nu*lap(vx) - NL + f)/o— así que al tomar
la divergencia del paso de proyección

    v^{n+1} = v* - dt·grad(p) ,   grad^2(p) = div(v*)/dt

queda  div(v^{n+1}) = div(v*) - dt·grad^2(p) = 0  IDÉNTICAMENTE, sin ningún operador
diferencial actuando sobre Q. No hay ecuación de Helmholtz, no hay solución homogénea, y
por lo tanto no hay capa límite de divergencia. El error de partición tiene que aparecer
en otro lado: en la condición de contorno tangencial.

PREDICCIONES FALSABLES, y son tres
==================================
  P1. <|div v|^2> es INDEPENDIENTE de dt, en los dos tipos de borde. Si hubiera capa de
      divergencia a la KIO, tendría que crecer con dt.
  P2. Con no-deslizamiento, la velocidad de deslizamiento residual escala como dt^(2·ORD)
      en <|v_t|^2>, o sea dt^ORD en rms. Ahí va a parar el error de partición.
  P3. Con superficie libre, el residuo de tensión es INDEPENDIENTE de dt, porque la
      condición se transmite por el rotor y la proyección resta un gradiente ([D-24]).
      O sea que con tope libre el error de partición no aparece NI en la divergencia NI
      en la condición de contorno.

Si P1 fallara, la lectura de que "es el mismo error contabilizado en otra variable" sería
incorrecta y habría que rehacerla.

Uso:
    /home/lucia/miniforge3/envs/piv-dt/bin/python verificacion/kio_capa_divergencia.py
"""

import importlib.util
import json
import math
import os
import sys

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(AQUI)

# Se reutiliza la maquinaria de la puerta de la Fase 2: compilar, correr, leer
# diagnósticos. La puerta está fijada por hash y NO se modifica; sólo se importa.
_spec = importlib.util.spec_from_file_location(
    "fase2", os.path.join(AQUI, "test_aceptacion_fase2.py"))
fase2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fase2)

SALIDA = os.path.join(RAIZ, "salidas", "tablas", "kio_capa_divergencia.json")

# Campo 3D con dependencia en x, para que la presión y el término no lineal sean no
# triviales y la condición de contorno se ejercite con kx distinto de cero.
VPARAM1 = 1.0
T_FINAL = 0.08
DTS = (8.0e-4, 4.0e-4, 2.0e-4, 1.0e-4)


def escalada(x, y):
    """Exponente p tal que y ~ dt^p, por ajuste log-log."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = np.isfinite(y) & (y > 0)
    if ok.sum() < 2:
        return float("nan")
    return float(np.polyfit(np.log(x[ok]), np.log(y[ok]), 1)[0])


def main():
    scratch = os.environ.get("SPECTER_SCRATCH")
    temporal = scratch is None
    if temporal:
        import tempfile
        scratch = tempfile.mkdtemp(prefix="kio_")
    print("scratch: %s" % scratch, flush=True)

    bindir = fase2.construir(os.path.join(scratch, "arbol"))
    print("compilado", flush=True)

    res = {"dts": list(DTS), "t_final": T_FINAL, "ord": fase2.ORD,
           "malla": [fase2.NX, fase2.NY, fase2.NZ], "nu": fase2.NU,
           "corridas": {}}

    for etiqueta, bcsta, bcend, archivo in (
            ("noslip", "noslip", "noslip", "noslip_diagnostic.txt"),
            ("freeslip", "noslip", fase2.CADENA_LIBRE, "freeslip_diagnostic.txt")):
        filas = []
        for dt in DTS:
            pasos = int(round(T_FINAL / dt))
            _, _, dirrun = fase2.correr(
                bindir, "%s_dt%.0e" % (etiqueta, dt), bcsta=bcsta, bcend=bcend,
                vparam1=VPARAM1, dt=dt, step=pasos, cstep=max(pasos // 4, 1))
            d = fase2.leer_diagnostico(dirrun, archivo)
            if d.ndim == 1:
                d = d[None, :]
            ult = d[-1]
            # columnas: t, <|div|^2>, col2|z=0, col3|z=Lz, <|v_n|^2>|z=0, |z=Lz
            filas.append({
                "dt": dt, "pasos": pasos,
                "div2": float(ult[1]),
                "borde_z0": float(ult[2]), "borde_zLz": float(ult[3]),
                "vn2_z0": float(ult[4]), "vn2_zLz": float(ult[5]),
            })
            print("  %-9s dt=%.1e  <div^2>=%.3e  borde(z=0)=%.3e  borde(z=Lz)=%.3e"
                  % (etiqueta, dt, ult[1], ult[2], ult[3]), flush=True)
        res["corridas"][etiqueta] = filas

    # --- exponentes -----------------------------------------------------------
    dts = [f["dt"] for f in res["corridas"]["noslip"]]
    res["exponentes"] = {
        "noslip_div2": escalada(dts, [f["div2"] for f in res["corridas"]["noslip"]]),
        "noslip_slip2_z0": escalada(dts, [f["borde_z0"] for f in res["corridas"]["noslip"]]),
        "noslip_slip2_zLz": escalada(dts, [f["borde_zLz"] for f in res["corridas"]["noslip"]]),
        "freeslip_div2": escalada(dts, [f["div2"] for f in res["corridas"]["freeslip"]]),
        "freeslip_tension_zLz": escalada(
            dts, [f["borde_zLz"] for f in res["corridas"]["freeslip"]]),
    }

    # --- el espesor que tendría la capa, si existiera --------------------------
    dz = fase2.LZ / (fase2.NZFIS - 1)
    res["capa"] = {
        "dz": dz,
        "l_sobre_dz": {("%.1e" % dt): math.sqrt(fase2.NU * dt) / dz for dt in DTS},
        "dt_estabilidad_viscosa": 2.0 * dz ** 2 / (math.pi ** 2 * fase2.NU),
    }

    os.makedirs(os.path.dirname(SALIDA), exist_ok=True)
    with open(SALIDA, "w") as f:
        json.dump(res, f, indent=1)

    e = res["exponentes"]
    print("\nEXPONENTES  (y ~ dt^p)")
    print("  P1  <div^2>, noslip      p = %+.3f   (se espera ~0)" % e["noslip_div2"])
    print("  P1  <div^2>, freeslip    p = %+.3f   (se espera ~0)" % e["freeslip_div2"])
    print("  P2  <|v_t|^2> en z=0     p = %+.3f   (se espera %d)"
          % (e["noslip_slip2_z0"], 2 * fase2.ORD))
    print("  P2  <|v_t|^2> en z=Lz    p = %+.3f   (se espera %d)"
          % (e["noslip_slip2_zLz"], 2 * fase2.ORD))
    print("  P3  residuo de tensión   p = %+.3f   (se espera ~0)"
          % e["freeslip_tension_zLz"])
    print("\nESPESOR DE LA CAPA, si existiera:  l = sqrt(nu*dt),  dz = %.4f" % dz)
    for dt in DTS:
        print("   dt=%.1e ->  l/dz = %.3f" % (dt, math.sqrt(fase2.NU * dt) / dz))
    print("  (la estabilidad viscosa explícita exige dt <~ %.1e, o sea l/dz <~ %.2f)"
          % (res["capa"]["dt_estabilidad_viscosa"],
             math.sqrt(fase2.NU * res["capa"]["dt_estabilidad_viscosa"]) / dz))
    print("\nescrito %s" % os.path.relpath(SALIDA, RAIZ))

    if temporal and os.environ.get("SPECTER_KEEP") != "1":
        import shutil
        shutil.rmtree(scratch, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
