#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Corridas de SPECTER con la superficie deformable, para figuras y chequeos adicionales.

La puerta verificacion/test_aceptacion_fase5.py decide si la implementación pasa; este
script guarda las series completas para poder MIRARLAS (salidas/tablas/
superficie_deformable_corridas.json) y agrega tres chequeos que la puerta no tiene,
porque se pensaron después de congelarla:

  E1  desacople toroidal: un modo con velocidad horizontal perpendicular a k (w = 0,
      p = 0) decae igual con 'freesurface' que con 'freeslip', y eta queda en cero.
  E2  reinicio: N pasos + reinicio desde el binario + N pasos da lo mismo que 2N pasos
      seguidos (ejercita fs_eta.NNNN.dat).
  E3  amplitud: el desvío respecto de la referencia LINEAL crece con la amplitud, como
      corresponde al término no lineal del volumen, y a amplitud chica es el piso
      numérico de F4.

Usa el arnés y la referencia de la puerta (importados, no copiados). Corre todo en
secuencia, un proceso MPI por vez salvo el caso np=2. Del orden de 15 minutos en la
laptop, ~150 MB.

Uso:
    ulimit -v 2500000
    /home/lucia/miniforge3/envs/piv-dt/bin/python numerico/fase5/corridas_superficie_deformable.py
"""

import json
import os
import shutil
import sys
import tempfile

import numpy as np

RAIZ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(RAIZ, "verificacion"))
import test_aceptacion_fase5 as G                              # noqa: E402
import test_aceptacion_fase2 as F2                             # noqa: E402

SALIDA = os.path.join(RAIZ, "salidas", "tablas", "superficie_deformable_corridas.json")

# Condición inicial toroidal agregada a la de la puerta (vparam1 = 3). El factor 1/2
# hace que vy = u0 sin(pi z / 2Lz) cos(x): el modo kx = 1 guarda la mitad del coseno.
IC_TOROIDAL = """         ELSEIF ( int(vparam1) .eq. 3 ) THEN
            DO k = 1,nz-Cz
               vy(k,1,2) = u0*nx*ny*SIN(0.5_GP*pi*z(k)/Lz)/2
            ENDDO
         ENDIF
      ENDIF
"""


def construir(scratch):
    ic = G.IC_FS.replace("         ENDIF\n      ENDIF\n", IC_TOROIDAL, 1)
    if "vparam1) .eq. 3" not in ic:
        raise RuntimeError("no se pudo agregar la condición inicial toroidal")
    G.IC_FS = ic
    return G.construir(os.path.join(scratch, "main"), "main")


def serie(dirrun, dt, cstep):
    t, r, d = G.serie_eta(dirrun, dt, cstep)
    return t, r, d


def main():
    scratch = os.environ.get("SPECTER_SCRATCH") or tempfile.mkdtemp(prefix="fase5_corr_")
    b = construir(scratch)
    res = {"malla": [G.NX, G.NY, G.NZ, G.CZ], "ORD": G.ORD, "nu": G.NU, "g": G.GRAV,
           "A": G.AMP}

    # Trayectorias de F4 (series completas) y el caso con 2 procesos
    res["trayectorias"] = {}
    for nombre, kx, ky, tens in G.CASOS:
        K = float(np.hypot(kx, ky))
        pasos = int(round(G.T_F4 / G.DT_F4))
        _, _, d = G.correr(b, "tr_" + nombre, step=pasos, fskx=kx, fsky=ky, fstens=tens)
        t, r, dd = serie(d, G.DT_F4, G.CSTEP_F4)
        ref, err = G.referencia_eta(K, G.GRAV + tens * K * K, G.CSTEP_F4 * G.DT_F4,
                                    len(t) - 1)
        var = {}
        for etq, kw in (("sin_2nu_dwdz", dict(visc=False)),
                        ("freeslip_plano", dict(tension="plana"))):
            var[etq] = G._trayectoria(200, K, G.GRAV + tens * K * K,
                                      G.CSTEP_F4 * G.DT_F4, len(t) - 1, **kw)[0][-1].real
        res["trayectorias"][nombre] = {
            "K": K, "tension": tens, "t": t.tolist(), "specter": r.real.tolist(),
            "referencia": ref.real.tolist(), "error_referencia": err,
            "max_error": float(np.max(np.abs(r.real - ref.real))),
            "tension_superficie": dd[:, 5].tolist(), "tension_fondo": dd[:, 6].tolist(),
            "variantes": {k: v.tolist() for k, v in var.items()}}
        print(nombre, "max error", res["trayectorias"][nombre]["max_error"], flush=True)

    # Tensión tangencial: exacto contra rezagado, tres dt
    res["tension"] = []
    for coup in (2, 1):
        for dt in (1.0e-3, 5.0e-4, 2.5e-4):
            pasos = int(round(0.2 / dt))
            _, _, d = G.correr(b, "tau_c%d_%g" % (coup, dt), dt=dt, step=pasos,
                               cstep=max(1, pasos // 10), fscoup=coup)
            dd = G.diag_fs(d)
            m = dd[len(dd) // 2:]
            res["tension"].append({"fscoup": coup, "dt": dt,
                                   "residuo_relativo": float(np.max(m[:, 5] / m[:, 6])),
                                   "discordancia_W": float(np.max(dd[:, 8]))})
            print("tension", res["tension"][-1], flush=True)

    # Orden temporal
    res["orden"] = []
    series = []
    for dt, cstep in ((1.0e-3, 40), (5.0e-4, 80), (2.5e-4, 160)):
        _, _, d = G.correr(b, "ord_%g" % dt, dt=dt, step=int(round(4.0 / dt)),
                           cstep=cstep)
        series.append(serie(d, dt, cstep)[1])
    n = min(len(s) for s in series)
    for i, dt in enumerate((1.0e-3, 5.0e-4)):
        res["orden"].append({"dt": dt, "dif_con_dt_mitad":
                             float(np.max(np.abs(series[i][:n] - series[i + 1][:n])))})
    print("orden", res["orden"], flush=True)

    # E1: desacople toroidal
    e1 = {}
    for top in ("freesurface", "freeslip"):
        _, _, d = G.correr(b, "e1_" + top, bcend=top, dt=F2.DT_BASE, step=6000, cstep=20,
                           u0=1.0, vparam1=3.0, fsamp=0.0)
        t, e, w = F2.leer_balance(d)
        e1[top] = (d, e)
    lam, discrep, disp = F2.tasa_de_corrida(e1["freesurface"][0])
    lam_fs, _, _ = F2.tasa_de_corrida(e1["freeslip"][0])
    ref = F2.tasas_referencia("neumann", 1)[0] + G.NU * 1.0
    res["E1_toroidal"] = {
        "lambda_freesurface": lam, "lambda_freeslip": lam_fs, "referencia_FD": ref,
        "error_relativo": abs(lam - ref) / ref,
        "max_dif_relativa_energia": float(np.max(np.abs(e1["freesurface"][1]
                                                        - e1["freeslip"][1])
                                                 / e1["freeslip"][1])),
        "max_rms_eta": float(np.max(np.abs(G.diag_fs(e1["freesurface"][0])[:, 4])))}
    print("E1", res["E1_toroidal"], flush=True)

    # E2: reinicio. SPECTER escribe el binario 0001 al empezar (t=0) y el 0002 al empezar
    # el paso N+1, así que el reinicio desde t = N dt es stat = 2.
    N = 2000
    base = dict(cstep=50, fskx=1, fsky=1, fstens=0.2)
    _, _, d_dir = G.correr(b, "e2_directo", step=2 * N, **base)
    plantilla = G.PLANTILLA
    try:
        G.PLANTILLA = plantilla.replace("tstep = 1000000", "tstep = %d" % N)
        _, _, d1 = G.correr(b, "e2_tramo1", step=N + 1, **base)
        d2 = os.path.join(os.path.dirname(b), "corridas", "e2_tramo2")
        os.makedirs(os.path.join(d2, "in"), exist_ok=True)
        for f in os.listdir(os.path.join(d1, "out")):
            shutil.copy(os.path.join(d1, "out", f), os.path.join(d2, "in", f))
        G.PLANTILLA = G.PLANTILLA.replace("stat = 0", "stat = 2")
        _, _, d2 = G.correr(b, "e2_tramo2", step=2 * N, **base)
    finally:
        G.PLANTILLA = plantilla
    a, c = G.diag_fs(d_dir), G.diag_fs(d2)
    ia = {int(round(x / G.DT_F4)): i for i, x in enumerate(a[:, 0])}
    comunes = [(ia[int(round(x / G.DT_F4))], j) for j, x in enumerate(c[:, 0])
               if int(round(x / G.DT_F4)) in ia]
    res["E2_reinicio"] = {
        "tiempos_comunes": len(comunes),
        "t_inicio_tramo2": float(c[0, 0]),
        "max_dif_eta_sobre_A": float(max(abs(a[i, 1] - c[j, 1]) for i, j in comunes))
        / G.AMP,
        "archivos_tramo1": sorted(os.listdir(os.path.join(d1, "out")))}
    print("E2", res["E2_reinicio"], flush=True)

    # E3: amplitud
    res["E3_amplitud"] = []
    for A in (1e-6, 1e-4, 1e-3, 1e-2):
        _, _, d = G.correr(b, "e3_%g" % A, step=int(round(4.0 / G.DT_F4)), fsamp=A)
        dd = G.diag_fs(d)
        r = dd[:, 1] / A
        ref, _ = G.referencia_eta(1.0, G.GRAV, G.CSTEP_F4 * G.DT_F4, len(r) - 1)
        res["E3_amplitud"].append({"A": A, "max_desvio": float(np.max(np.abs(r - ref.real)))})
        print("E3", res["E3_amplitud"][-1], flush=True)

    os.makedirs(os.path.dirname(SALIDA), exist_ok=True)
    with open(SALIDA, "w") as f:
        json.dump(res, f, indent=1)
    print("escrito", SALIDA)
    if "SPECTER_SCRATCH" not in os.environ and os.environ.get("SPECTER_KEEP") != "1":
        shutil.rmtree(scratch, ignore_errors=True)


if __name__ == "__main__":
    main()
