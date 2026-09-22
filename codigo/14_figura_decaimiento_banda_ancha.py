#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Figura de la prueba puente SPECTER--PIV con condición inicial de banda ancha.

Lee solamente las salidas compactas. Si ya existe convergencia en 64² o 128², usa la
resolución más fina que tenga completos el control lineal y el caso experimental.

Uso:
    /home/lucia/miniforge3/envs/piv-dt/bin/python \
        codigo/14_figura_decaimiento_banda_ancha.py
"""

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt              # noqa: E402
import numpy as np                            # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(AQUI)
BANDA = os.path.join(RAIZ, "salidas", "tablas", "decaimiento_banda_ancha.json")
PIV = os.path.join(RAIZ, "salidas", "tablas", "decaimiento_espectral.json")
FIGURA = os.path.join(RAIZ, "salidas", "figuras", "f10_decaimiento_banda_ancha")

AZUL = "#2878a8"
NARANJA = "#d97706"
GRIS = "#606770"


def cargar():
    with open(BANDA) as fh:
        banda = json.load(fh)
    with open(PIV) as fh:
        piv = json.load(fh)
    completos = []
    for clave, bloque in banda["resoluciones"].items():
        if bloque["config"].get("prueba"):
            continue
        if {"lineal", "experimental"}.issubset(bloque["corridas"]):
            completos.append((int(bloque["config"]["resolucion"]), clave, bloque))
    if not completos:
        raise RuntimeError("no hay una resolución completa en %s" % BANDA)
    return max(completos), piv


def serie(corrida):
    p = corrida["tasas_por_anillo"]
    return (np.asarray([x["k_centro_fisico_1_m"] for x in p]),
            np.asarray([x["lambda_fisico_1_s"] for x in p]))


def main():
    (nxy, clave, bloque), piv = cargar()
    lineal = bloque["corridas"]["lineal"]
    no_lineal = bloque["corridas"]["experimental"]
    kp = np.asarray([x["k_1_m"] for x in piv["ajuste_modal_k2"]["puntos"]])
    lp = np.asarray([x["lambda_1_s"] for x in piv["ajuste_modal_k2"]["puntos"]])
    ep = np.asarray([x["error_formal_1_s"] for x in piv["ajuste_modal_k2"]["puntos"]])
    kl, ll = serie(lineal)
    kn, ln = serie(no_lineal)

    plt.rcParams.update({"font.size": 10.2, "axes.spines.top": False,
                         "axes.spines.right": False})
    fig, (ax, ad) = plt.subplots(1, 2, figsize=(11.7, 4.45),
                                 gridspec_kw={"wspace": 0.28})

    kfit = np.linspace(70.0, 260.0, 300)
    nominal = piv["parametros_nominales"]
    ax.plot(kfit, nominal["alfa_1_s"] + nominal["nu_m2_s"] * kfit ** 2,
            color="black", lw=1.35, ls="--", label=r"nominal $\alpha+\nu k^2$")
    ax.errorbar(kp, lp, yerr=ep, fmt="o", ms=4.7, color=GRIS, mfc="white",
                capsize=2.2, lw=1, label="PIV, t ≥ 10 s")
    ax.plot(kl, ll, "s", ms=4.5, color=AZUL, label="SPECTER, control lineal")
    ax.plot(kn, ln, "o", ms=4.5, color=NARANJA,
            label=r"SPECTER, banda ancha ($\delta_0=8{,}89$)")
    for corrida, color in ((lineal, AZUL), (no_lineal, NARANJA)):
        ajuste = corrida["ajuste_lambda_k2"]
        ax.plot(kfit, ajuste["alfa_intercepto_1_s"]
                + ajuste["nu_pendiente_m2_s"] * kfit ** 2,
                color=color, lw=1.2, alpha=0.75)
    ax.set(xlim=(72, 258), ylim=(0.042, 0.151), xlabel=r"$k$ [m$^{-1}$]",
           ylabel=r"tasa modal $\lambda$ [s$^{-1}$]")
    ax.grid(alpha=0.18)
    ax.legend(frameon=False, fontsize=8.5, loc="upper left")
    ax.text(0.03, 0.05,
            ("ajustes SPECTER:\n"
             r"lineal: $\alpha=%.4f$, $\nu=%.3f\times10^{-6}$" "\n"
             r"banda: $\alpha=%.4f$, $\nu=%.3f\times10^{-6}$")
            % (lineal["ajuste_lambda_k2"]["alfa_intercepto_1_s"],
               lineal["ajuste_lambda_k2"]["nu_pendiente_m2_s"] * 1e6,
               no_lineal["ajuste_lambda_k2"]["alfa_intercepto_1_s"],
               no_lineal["ajuste_lambda_k2"]["nu_pendiente_m2_s"] * 1e6),
            transform=ax.transAxes, va="bottom", fontsize=8.25,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82,
                  "pad": 2.2})

    lineal_por_k = {round(k, 6): y for k, y in zip(kl, ll)}
    delta_piv = np.asarray([y - lineal_por_k[round(k, 6)] for k, y in zip(kp, lp)])
    comp = bloque["comparacion"]["delta_lambda_no_lineal_menos_lineal_por_anillo"]
    kc = np.asarray([x["k_centro_fisico_1_m"] for x in comp])
    dc = np.asarray([x["delta_lambda_1_s"] for x in comp])
    piv_fit = piv["ajuste_modal_k2"]
    delta_intercepto_piv = (piv_fit["alfa_espectral_1_s"]
                            - nominal["alfa_1_s"])
    error_intercepto_piv = piv_fit["alfa_espectral_sem_entre_registros_1_s"]
    ad.axhspan(delta_intercepto_piv - error_intercepto_piv,
               delta_intercepto_piv + error_intercepto_piv,
               color=GRIS, alpha=0.11)
    ad.axhline(0.0, color="black", lw=1.0)
    ad.errorbar(kp, delta_piv, yerr=ep, fmt="o", ms=4.7, color=GRIS,
                mfc="white", capsize=2.2, lw=1, label="PIV − control lineal")
    ad.plot(kc, dc, "o-", ms=4.3, lw=1.15, color=NARANJA,
            label="banda ancha − control lineal")
    ajuste_l = lineal["ajuste_lambda_k2"]
    ajuste_n = no_lineal["ajuste_lambda_k2"]
    delta_intercepto_num = (ajuste_n["alfa_intercepto_1_s"]
                             - ajuste_l["alfa_intercepto_1_s"])
    ad.set(xlim=(72, 258), ylim=(-0.032, 0.018), xlabel=r"$k$ [m$^{-1}$]",
           ylabel=r"cambio de tasa $\Delta\lambda$ [s$^{-1}$]")
    ad.grid(alpha=0.18)
    ad.legend(frameon=False, fontsize=8.5, loc="upper left")
    ad.text(0.03, 0.05,
            (r"$\Delta\alpha_{PIV}=%.4f\pm%.4f$ s$^{-1}$" "\n"
             r"$\Delta\alpha_{SPECTER}=%.4f$ s$^{-1}$" "\n"
             r"$\lambda_{sup}^{banda}/\lambda_{sup}^{lin}=%.3f$")
            % (delta_intercepto_piv, error_intercepto_piv,
               delta_intercepto_num,
               bloque["comparacion"]["cociente_lambda_superficie"]),
            transform=ad.transAxes, va="bottom", fontsize=8.5,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.84,
                  "pad": 2.2})

    nzfis = bloque["config"]["nz"] - bloque["config"]["cz"]
    fig.suptitle("Puente numérico PIV–SPECTER: banda ancha a resolución "
                 r"$%d^2\times%d$" % (nxy, nzfis), y=1.01, fontsize=12)
    os.makedirs(os.path.dirname(FIGURA), exist_ok=True)
    fig.savefig(FIGURA + ".png", dpi=180, bbox_inches="tight", pad_inches=0.18)
    fig.savefig(FIGURA + ".pdf", bbox_inches="tight", pad_inches=0.18)
    plt.close(fig)
    print("resolución usada:", clave)
    print("escrito", os.path.relpath(FIGURA + ".png", RAIZ))
    print("escrito", os.path.relpath(FIGURA + ".pdf", RAIZ))


if __name__ == "__main__":
    main()
