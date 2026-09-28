#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Figuras de la superficie deformable, desde salidas/tablas/superficie_deformable_corridas.json
(numerico/fase5/corridas_superficie_deformable.py). No corre SPECTER.

  f_sd1  eta(t) de SPECTER contra la referencia MAC, los dos casos de la puerta, con el
         error y dos modelos equivocados como escala
  f_sd2  tensión tangencial en la superficie, exacta contra rezagada, y orden temporal
"""

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                      # noqa: E402
import numpy as np                                   # noqa: E402

RAIZ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATOS = os.path.join(RAIZ, "salidas", "tablas", "superficie_deformable_corridas.json")
FIG = os.path.join(RAIZ, "salidas", "figuras")


def main():
    d = json.load(open(DATOS))
    os.makedirs(FIG, exist_ok=True)

    titulos = {"A_x_gravedad": "A: modo (1,0), gravedad",
               "B_oblicuo_capilar": "B: modo (1,1), gravedad + capilaridad"}
    fig, ax = plt.subplots(2, 2, figsize=(11, 6), sharex=True,
                           gridspec_kw={"height_ratios": [2, 1]})
    for c, (nombre, tr) in enumerate(d["trayectorias"].items()):
        t = np.array(tr["t"])
        s, r = np.array(tr["specter"]), np.array(tr["referencia"])
        a = ax[0, c]
        a.plot(t, r, "-", color="0.6", lw=3, label="referencia MAC (Richardson)")
        a.plot(t, s, "-", color="C0", lw=1.2, label="SPECTER")
        a.plot(t, tr["variantes"]["sin_2nu_dwdz"], ":", color="C3",
               label="modelo sin 2ν ∂w/∂z")
        a.plot(t, tr["variantes"]["freeslip_plano"], "--", color="C1",
               label="freeslip plano reusado")
        a.set_title(titulos.get(nombre, nombre), fontsize=10)
        a.set_ylabel("η(t) / A")
        a.legend(fontsize=7)
        b = ax[1, c]
        b.semilogy(t[1:], np.abs(s - r)[1:], color="C0", label="|SPECTER − ref|")
        b.axhline(1e-4, color="k", ls="--", lw=0.8, label="tolerancia F4")
        b.axhline(tr["error_referencia"], color="0.5", ls=":", lw=0.8,
                  label="error estimado de la ref.")
        b.set_xlabel("t (unidades del código)")
        b.set_ylabel("error")
        b.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "f_sd1_trayectorias.png"), dpi=150)

    fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
    for coup, mk, lab in ((2, "o-", "acoplamiento exacto (fscoup=2)"),
                          (1, "s--", "rezagado, control negativo (fscoup=1)")):
        sel = [x for x in d["tension"] if x["fscoup"] == coup]
        ax[0].loglog([x["dt"] for x in sel], [x["residuo_relativo"] for x in sel], mk,
                     label=lab)
    dts = np.array([1e-3, 2.5e-4])
    ref = [x for x in d["tension"] if x["fscoup"] == 1][0]
    ax[0].loglog(dts, ref["residuo_relativo"] * (dts / ref["dt"]) ** 2, ":", color="0.5",
                 label="∝ dt² (en cuadrados)")
    ax[0].set_xlabel("dt")
    ax[0].set_ylabel("⟨|τ_t|²⟩ sup. / ⟨|∂v_t/∂z|²⟩ fondo")
    ax[0].set_title("Tensión tangencial en la superficie", fontsize=10)
    ax[0].legend(fontsize=7)
    o = d["orden"]
    ax[1].loglog([x["dt"] for x in o], [x["dif_con_dt_mitad"] for x in o], "o-",
                 label="max |η_dt − η_dt/2| / A")
    dt0 = o[0]["dt"]
    x = np.array([o[-1]["dt"], dt0])
    ax[1].loglog(x, o[0]["dif_con_dt_mitad"] * (x / dt0) ** 2, ":", color="0.5", lw=2,
                 label="∝ dt² (pasando por el primer punto)")
    ax[1].set_xlabel("dt")
    ax[1].set_title("Orden temporal (caso A, ORD=2)", fontsize=10)
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "f_sd2_tension_y_orden.png"), dpi=150)
    print("figuras escritas en", FIG)


if __name__ == "__main__":
    main()
