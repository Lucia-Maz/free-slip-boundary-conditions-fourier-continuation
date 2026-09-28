#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Las dos escaleras de la superficie de orden N, fuera de la puerta congelada, con sus mismas
funciones (se importa verificacion/test_aceptacion_fase6.py; no se edita). [H-25], [H-26].

  profundidad  elevación media uniforme eta0 y onda lineal de amplitud 1e-6, modo (1,1),
               malla 8x8xNZ: error de eta(t) contra la referencia MAC a profundidad 1 + eta0,
               con el piso de malla (fsorder = 1 plano contra MAC a profundidad 1) restado.
               Se compara con el truncamiento NO viscoso esperado para cada N (derivación de
               esta sesión): omega_N^2 = (g k + gam k^3) S_N/C_N, con S_N, C_N las trazas de
               orden N-1 de sinh y cosh.
  no_lineal    onda estacionaria 2D (modo (1,0), 16x1x128) de amplitud eps, contra la
               referencia de geometría exacta (verificacion/referencia_fase6.py), con la
               métrica de H7. N = 1 y N = 3: N = 2 aborta para eta < -0,7 dz ([H-24]).
  figura       f_fs6_escaleras.png desde los dos JSON.

Uso (laptop, en secuencia, con ulimit -v 2500000):
  PY=/home/lucia/miniforge3/envs/piv-dt/bin/python
  $PY numerico/fase6/escaleras.py profundidad --nz 256 --dt 2.5e-4 --T 2    # ~25 min, 300 MB
  $PY numerico/fase6/escaleras.py no_lineal                                  # ~5 min, 320 MB
  $PY numerico/fase6/escaleras.py figura
Compila SPECTER en un directorio temporal (o en --scratch) desde numerico/SPECTER-trabajo.
"""

import argparse
import json
import math
import os
import shutil
import sys
import tempfile
import time

import numpy as np

RAIZ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(RAIZ, "verificacion"))
import test_aceptacion_fase6 as H          # noqa: E402  (congelado: sólo se importa)
import test_aceptacion_fase5 as G5         # noqa: E402

TABLAS = os.path.join(RAIZ, "salidas", "tablas")
FIG = os.path.join(RAIZ, "salidas", "figuras", "f_fs6_escaleras")
J_PROF = os.path.join(TABLAS, "fase6_escalera_profundidad.json")
J_NL = os.path.join(TABLAS, "fase6_escalera_no_lineal.json")
ETA0 = (0.01, 0.02, 0.04)


def truncamiento(N, eta0, T, k=math.sqrt(2.0), h=1.0, g=H.GRAV, gam=H.TENS):
    """max_t |cos(omega_N t) - cos(omega t)| en [0, T], límite lineal no viscoso."""
    s, c = math.sinh(k * h), math.cosh(k * h)
    S = sum((k * eta0)**m / math.factorial(m) * (s if m % 2 == 0 else c) for m in range(N))
    C = sum((k * eta0)**m / math.factorial(m) * (c if m % 2 == 0 else s) for m in range(N))
    w_n = math.sqrt((g * k + gam * k**3) * S / C)
    w_x = math.sqrt((g * k + gam * k**3) * math.tanh(k * (h + eta0)))
    t = np.linspace(0.0, T, 4001)
    return float(np.max(np.abs(np.cos(w_n * t) - np.cos(w_x * t))))


def _correr(b, lab, **kw):
    d = os.path.join(os.path.dirname(b), "corridas", lab)
    if os.path.exists(d):
        shutil.rmtree(d)
    t0 = time.time()
    H.correr(b, lab, **kw)
    return d, time.time() - t0


def profundidad(args):
    scratch = args.scratch or tempfile.mkdtemp(prefix="fase6_prof_")
    b = H.construir(os.path.join(scratch, "m3d_%d" % args.nz), "main", (8, 8, args.nz))
    dt, T = args.dt, args.T
    cs = int(round(0.04 / dt))
    pasos = int(round(T / dt)) + 1
    K = math.sqrt(2.0)
    G = H.GRAV + H.TENS * K * K
    kw = dict(dt=dt, cstep=cs, step=pasos, fskx=1, fsky=1, fstens=H.TENS)
    d0, s0 = _correr(b, "piso", fsorder=1, **kw)
    t, r0, _ = G5.serie_eta(d0, dt, cs)
    piso = r0.real - H.mac_eta(K, G, 1.0, cs * dt, len(t) - 1)
    res = {"nz": args.nz, "dt": dt, "T": T, "eta0": list(ETA0),
           "piso_max": float(np.max(np.abs(piso))), "N": {}}
    print("nz %d dt %g T %g: piso de malla %.2e (%.0f s)" % (args.nz, dt, T, res["piso_max"], s0),
          flush=True)
    for N in (1, 2, 3):
        err, esp, pas, rmax = [], [], [], []
        for e0 in ETA0:
            d, s = _correr(b, "N%d_%g" % (N, e0), fsorder=N, fsmean=e0, fsgen=1, **kw)
            _, r, _, _, dd = H.serie_modos(d, dt, cs)
            ref = H.mac_eta(K, G, 1.0 + e0, cs * dt, len(t) - 1)
            err.append(float(np.max(np.abs((r.real / H.AMP_LIN - ref) - piso))))
            esp.append(truncamiento(N, e0, T))
            pas.append(int(dd[:, 10].max()))
            rmax.append(float(dd[:, 11].max()))
            print("  N=%d eta0=%g: error %.3e (truncamiento esperado %.2e); pasadas <= %d, "
                  "residuo <= %.1e (%.0f s)" % (N, e0, err[-1], esp[-1], pas[-1], rmax[-1], s),
                  flush=True)
        res["N"][str(N)] = {"error": err, "esperado": esp, "ordenes": H.ordenes(err),
                            "pasadas_max": pas, "residuo_max": rmax}
    os.makedirs(TABLAS, exist_ok=True)
    json.dump(res, open(J_PROF, "w"), indent=1)
    print("->", J_PROF)
    if not args.scratch:
        shutil.rmtree(scratch)


def no_lineal(args):
    scratch = args.scratch or tempfile.mkdtemp(prefix="fase6_nl_")
    b = H.construir(os.path.join(scratch, "m2d"), "main", H.MALLA_2D)
    ctx = {}
    pasos = int(round(H.T_FIN / H.DT)) + 1
    kw = dict(step=pasos, fskx=1, fsky=0, fstens=H.TENS)
    d0, _ = _correr(b, "piso", fsorder=1, fsgen=1, fsamp=H.AMP_LIN, **kw)
    S0 = H.serie_modos(d0)
    R0 = H._referencia(ctx, H.AMP_LIN)
    piso1 = (S0[1] - R0["eta1"]) / H.AMP_LIN
    res = {"eps": list(H.EPS_ONDA), "malla": list(H.MALLA_2D), "T": H.T_FIN, "dt": H.DT, "N": {}}
    for N in (1, 3):
        err, pas, rmax, rel = [], [], [], []
        for eps in H.EPS_ONDA:
            d, s = _correr(b, "N%d_%g" % (N, eps), fsorder=N, fsgen=1, fsamp=eps, **kw)
            S = H.serie_modos(d)
            err.append(H._error_estado(S[:4], H._referencia(ctx, eps), piso1, eps))
            dd = S[4]
            pas.append(int(dd[:, 10].max()))
            rmax.append(float(dd[:, 11].max()))
            rel.append(float(dd[:, 12].max()))
            print("  N=%d eps=%g: error %.3e; pasadas <= %d, residuo <= %.1e, max|eta|/dz %.2f "
                  "(%.0f s)" % (N, eps, err[-1], pas[-1], rmax[-1], rel[-1], s), flush=True)
        res["N"][str(N)] = {"error": err, "ordenes": H.ordenes(err), "pasadas_max": pas,
                            "residuo_max": rmax, "max_eta_sobre_dz": rel}
        print("N=%d: órdenes %s" % (N, " ".join("%.2f" % q for q in res["N"][str(N)]["ordenes"])),
              flush=True)
    os.makedirs(TABLAS, exist_ok=True)
    json.dump(res, open(J_NL, "w"), indent=1)
    print("->", J_NL)
    if not args.scratch:
        shutil.rmtree(scratch)


def figura(args):
    import matplotlib
    matplotlib.use("Agg")
    sys.path.insert(0, os.path.expanduser("~/Desktop/PhD/Num/Codigos"))
    try:
        import estetica  # noqa: F401  (rcParams compartidos)
        import matplotlib.pyplot as plt
        plt.rcParams["font.serif"] = ["Times New Roman", "Nimbus Roman", "Liberation Serif"]
        plt.rcParams["font.family"] = "serif"
        ancho = estetica.width_in
    except ImportError:
        import matplotlib.pyplot as plt
        ancho = 7.0
    P, NL = json.load(open(J_PROF)), json.load(open(J_NL))
    fig, ax = plt.subplots(1, 2, figsize=(ancho, 0.40 * ancho))
    col = {"1": "C0", "2": "C1", "3": "C3"}
    x = np.array(P["eta0"])
    a = ax[0]
    for N, d in P["N"].items():
        a.loglog(x, d["error"], "o-", color=col[N], lw=1.2, label="SPECTER, N = %s" % N)
        a.loglog(x, d["esperado"], "--", color=col[N], lw=0.9, alpha=0.8)
    a.axhline(P["piso_max"], color="0.6", lw=0.8, ls=":")
    a.text(x[0], P["piso_max"] * 1.3, "piso de malla (restado)", fontsize=7, color="0.4")
    a.set_xlabel(r"$\eta_0$")
    a.set_ylabel(r"$\max_t|\eta_1(t)/A-\eta_{1,\mathrm{ref}}(t)|$")
    a.set_title(r"profundidad, $n_z=%d$, $T=%g$ (trazos: truncamiento esperado)" % (P["nz"], P["T"]),
                fontsize=8)
    a.legend(fontsize=7)
    a = ax[1]
    e = np.array(NL["eps"])
    for N, d in NL["N"].items():
        a.loglog(e, d["error"], "o-", color=col[N], lw=1.2,
                 label="N = %s: órdenes %s" % (N, ", ".join("%.2f" % q for q in d["ordenes"])))
        pend = int(N)
        a.loglog(e, d["error"][-1] * (e / e[-1])**pend, "--", color=col[N], lw=0.8, alpha=0.7)
    a.set_xlabel(r"$\varepsilon$")
    a.set_ylabel("error de los armónicos 1–3, relativo a ε")
    a.set_title("onda estacionaria no lineal 2D, contra la geometría exacta", fontsize=8)
    a.legend(fontsize=7)
    fig.savefig(FIG + ".png", dpi=200)
    fig.savefig(FIG + ".pdf")
    print("->", FIG + ".png")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("profundidad")
    p.add_argument("--nz", type=int, default=256)
    p.add_argument("--dt", type=float, default=2.5e-4)
    p.add_argument("--T", type=float, default=2.0)
    p.add_argument("--scratch", default=None)
    q = sub.add_parser("no_lineal")
    q.add_argument("--scratch", default=None)
    sub.add_parser("figura")
    args = ap.parse_args()
    {"profundidad": profundidad, "no_lineal": no_lineal, "figura": figura}[args.cmd](args)


if __name__ == "__main__":
    main()
