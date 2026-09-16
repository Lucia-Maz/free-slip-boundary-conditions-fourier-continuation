#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Figuras de la etapa 1 del barrido en delta, para el PDF.

Van en una sola versión (fondo claro): el destino es un PDF, no una página que siga el
tema de quien mira. Paleta: slots 1 y 2 de la paleta categórica de referencia, en orden;
con dos series la validación por pares adyacentes de esa paleta aplica directamente.

Uso:
    /home/lucia/miniforge3/envs/piv-dt/bin/python codigo/10_figuras_barrido_delta.py --caso forzado
    /home/lucia/miniforge3/envs/piv-dt/bin/python codigo/10_figuras_barrido_delta.py --caso ventana

Un JSON por caso ([D-42]): salidas/tablas/barrido_delta_etapa1_<caso>.json.
"""

import importlib.util
import json
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt        # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(AQUI)
sys.path.insert(0, AQUI)

_spec = importlib.util.spec_from_file_location(
    "figs", os.path.join(AQUI, "07_figuras_decaimiento.py"))
figs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(figs)

TABLAS = os.path.join(RAIZ, "salidas", "tablas")
FIGS = os.path.join(RAIZ, "salidas", "figuras")

# barra de error del alfa medido en el experimento: 0,00488 / 0,06687
BARRA_EXPERIMENTAL = 0.0730


def guardar(fig, nombre):
    os.makedirs(FIGS, exist_ok=True)
    ruta = os.path.join(FIGS, nombre + ".pdf")
    fig.savefig(ruta, bbox_inches="tight", pad_inches=0.2)
    plt.close(fig)
    return ruta


def fig_apartamiento(d):
    t = figs.estilo("claro")
    fig, ax = plt.subplots(figsize=(7.0, 4.3))
    caso = d.get("caso", "")

    ax.axhspan(1 - BARRA_EXPERIMENTAL, 1 + BARRA_EXPERIMENTAL,
               color=t["grilla"], alpha=0.8, lw=0, zorder=0)
    ax.axhline(1.0, color=t["tinta2"], lw=1.3, ls="--", zorder=1)

    for i, clave in enumerate(sorted(d["corridas"])):
        filas = [f for f in d["corridas"][clave] if not f["crece"]]
        malos = [f for f in d["corridas"][clave] if f["crece"]]
        c = t["series"][i]
        n = clave[1:]
        x = [f["delta"] for f in filas]
        y = [f["lambda_sobre_referencia"] for f in filas]
        ax.plot(x, y, "o-", color=c, ms=6.5, mec=t["fondo"], mew=1.5, zorder=3)
        # alfa_eff/alfa medido sobre el campo: la fricción de fondo sin clausura
        if filas and "alfa_eff_sobre_alfa" in filas[0]:
            ax.plot(x, [f["alfa_eff_sobre_alfa"] for f in filas], "s--", color=c,
                    ms=5.5, mfc=t["fondo"], mec=c, mew=1.5, lw=1.1, zorder=3)
        if malos:
            ax.plot([f["delta"] for f in malos],
                    [f["lambda_sobre_referencia"] for f in malos],
                    "x", color=c, ms=9, mew=2, zorder=3)
        if x:
            ax.annotate("%s×%s" % (n, n), (x[-1], y[-1]), xytext=(7, 0),
                        textcoords="offset points", color=c, fontsize=10,
                        va="center", fontweight="semibold")

    ax.set_xscale("log")
    ax.set_xlabel("δ medido en la ventana de ajuste  (k horizontal del propio campo)")
    ax.set_ylabel("cociente contra la referencia lineal")
    ax.set_title("¿Se aparta la fricción al crecer el término no lineal?  caso %s, ε = %.2f"
                 % (caso, d.get("eps", float("nan"))))
    ax.annotate("círculos llenos: λ/λ_lin · cuadrados vacíos: α_eff/α del campo\n"
                "banda gris: ±7 %, la barra del α medido en el experimento",
                (0.02, 0.05), xycoords="axes fraction", color=t["tinta2"], fontsize=9)
    return guardar(fig, "f7_barrido_delta_%s" % caso if caso else "f7_barrido_delta")


def fig_decaimientos(d):
    t = figs.estilo("claro")
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    clave = "n%d" % max(d["resoluciones"])
    filas = d["corridas"][clave]
    tf = d["t_final"]
    ax.axvspan(tf * 2 / 3, tf, color=t["grilla"], alpha=0.8, lw=0, zorder=0)
    muestras = [filas[0]] + [filas[i] for i in (1, len(filas) // 2, len(filas) - 1)]
    for i, f in enumerate(muestras[:4]):
        c = t["series"][i % len(t["series"])]
        ax.plot([0, tf], [1.0, f["v2_final"] / f["v2_inicial"]], color=c, lw=0)
    ax.set_yscale("log")
    ax.set_xlabel("tiempo (unidades de código)")
    ax.set_ylabel("⟨v²⟩ / ⟨v²⟩(0)")
    ax.set_title("Ventana de ajuste")
    return guardar(fig, "f8_decaimientos")


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--caso", default="forzado")
    args = ap.parse_args()
    with open(os.path.join(TABLAS, "barrido_delta_etapa1_%s.json" % args.caso)) as f:
        d = json.load(f)
    r1 = fig_apartamiento(d)
    print("escrito %s" % os.path.relpath(r1, RAIZ))

    # resumen numérico que va al texto del PDF
    print("\nlambda / lambda_lineal")
    for clave in sorted(d["corridas"]):
        filas = d["corridas"][clave]
        print("  %s: " % clave + "  ".join(
            "%.5f%s" % (f["lambda_sobre_referencia"], "*" if f["crece"] else "")
            for f in filas[1:]))
        if "alfa_eff_sobre_alfa" in filas[1]:
            print("  %s  alfa_eff/alfa: " % (" " * len(clave)) + "  ".join(
                "%.5f" % f["alfa_eff_sobre_alfa"] for f in filas[1:]))
        buenos = [f for f in filas[1:] if not f["crece"]]
        if buenos:
            ap = max(abs(f["lambda_sobre_referencia"] - 1) for f in buenos)
            print("     apartamiento máximo entre las estables: %.3f %%" % (100 * ap))
            print("     delta cubierto: %.3g a %.3g"
                  % (min(f["delta"] for f in buenos),
                     max(f["delta"] for f in buenos)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
