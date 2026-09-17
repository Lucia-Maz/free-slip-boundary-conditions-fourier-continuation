#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Figuras del barrido en delta, para el PDF y la página de entrega: f7 (etapa 1,
[V-16]) y f8 (etapa 1b, la velocidad de superficie, [V-17]).

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
    """PDF para el informe y PNG para la página HTML de la entrega."""
    os.makedirs(FIGS, exist_ok=True)
    ruta = os.path.join(FIGS, nombre + ".pdf")
    fig.savefig(ruta, bbox_inches="tight", pad_inches=0.2)
    fig.savefig(os.path.join(FIGS, nombre + ".png"), dpi=150, bbox_inches="tight",
                pad_inches=0.2)
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
    ax.set_xlim(left=0.1)
    ax.set_xlabel("δ medido en la ventana de ajuste  (k horizontal del propio campo)")
    ax.set_ylabel("cociente contra la referencia lineal")
    ax.set_title("¿Se aparta la fricción al crecer el término no lineal?  caso %s, ε = %.2f"
                 % (caso, d.get("eps", float("nan"))))
    ax.annotate("círculos llenos: λ/λ_lin · cuadrados vacíos: α_eff/α del campo\n"
                "banda gris: ±7 %, la barra del α medido en el experimento",
                (0.02, 0.05), xycoords="axes fraction", color=t["tinta2"], fontsize=9)
    return guardar(fig, "f7_barrido_delta_%s" % caso if caso else "f7_barrido_delta")


def fig_superficie(ds, d1):
    """Etapa 1b ([V-17]): la tasa de la velocidad de superficie contra la del promedio
    vertical. `ds` es el JSON con sufijo _superficie; `d1` el de la etapa 1, del que sale
    qué puntos convergieron entre mallas (la superficie sólo se midió en 16²).

    Cada tasa se divide por LA MISMA tasa en la referencia lineal (λ_sup/λ_sup,lin, etc.),
    no por la λ de balance.txt: sobre la propia referencia los dos estimadores difieren
    0,05 % (ventana) y 0,5 % (forzado), y esa normalización lo quita ([V-17]).
    """
    t = figs.estilo("claro")
    caso = ds.get("caso", "")
    filas = ds["corridas"]["n16"]
    ref = filas[0]
    filas = [f for f in filas[1:] if not f["crece"]]

    # convergencia heredada de la etapa 1: |λ16 − λ32| / λ32 sobre λ/λ_lin
    conv = {}
    if d1 is not None and "n32" in d1["corridas"]:
        l16 = {round(f["u0"], 5): f["lambda_sobre_referencia"] for f in d1["corridas"]["n16"]}
        for f in d1["corridas"]["n32"]:
            a = l16.get(round(f["u0"], 5))
            if a is not None:
                conv[round(f["u0"], 5)] = abs(a - f["lambda_sobre_referencia"]) / f["lambda_sobre_referencia"] < 0.015
    ok = [f for f in filas if conv.get(round(f["u0"], 5), True)]
    no = [f for f in filas if not conv.get(round(f["u0"], 5), True)]

    c_sup, c_prom = t["series"][0], t["series"][1]
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(7.0, 6.4), sharex=True,
                                  gridspec_kw={"height_ratios": [3, 2], "hspace": 0.12})

    ax.axhspan(1 - BARRA_EXPERIMENTAL, 1 + BARRA_EXPERIMENTAL,
               color=t["grilla"], alpha=0.8, lw=0, zorder=0)
    ax.axhline(1.0, color=t["tinta2"], lw=1.3, ls="--", zorder=1)
    x = [f["delta"] for f in ok]
    ax.plot(x, [f["lambda_promedio"] / ref["lambda_promedio"] for f in ok], "s--",
            color=c_prom, ms=5.5, mfc=t["fondo"], mec=c_prom, mew=1.5, lw=1.1, zorder=3,
            label="⟨u⟩, promedio vertical (lo que predice la clausura)")
    ax.plot(x, [f["lambda_superficie"] / ref["lambda_superficie"] for f in ok], "o-",
            color=c_sup, ms=6.5, mec=t["fondo"], mew=1.5, zorder=4,
            label="u(h), superficie (lo que ajusta el PIV)")
    if no:
        xn = [f["delta"] for f in no]
        ax.plot(xn, [f["lambda_promedio"] / ref["lambda_promedio"] for f in no], "x",
                color=c_prom, ms=8, mew=1.8, zorder=3)
        ax.plot(xn, [f["lambda_superficie"] / ref["lambda_superficie"] for f in no], "x",
                color=c_sup, ms=8, mew=1.8, zorder=4)
    ax.legend(loc="upper left", fontsize=9.5)
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi + 0.22 * (hi - lo))      # aire para la leyenda
    ax.set_ylabel("tasa / tasa lineal")
    ax.set_title("¿Decae igual la superficie que el promedio vertical?  caso %s, ε = %.2f"
                 % (caso, ds.get("eps", float("nan"))))
    ax.annotate(("cruces: no convergido entre mallas en [V-16], sólo la dirección\n"
                 if no else "")
                + "banda gris: ±7 %, la barra del α medido en el experimento",
                (0.02, 0.05), xycoords="axes fraction", color=t["tinta2"], fontsize=9)

    # u(h)/<u> al principio y al final de la ventana: flecha de inicio a fin
    ax2.axhline(math.pi / 2, color=t["tinta2"], lw=1.3, ls="--", zorder=1)
    for f in ok + no:
        a, b = f["sup_sobre_prom_inicio"], f["sup_sobre_prom_fin"]
        al = 1.0 if f in ok else 0.4        # no convergido: atenuado
        ax2.annotate("", (f["delta"], b), (f["delta"], a),
                     arrowprops=dict(arrowstyle="-|>", color=c_sup, lw=1.6, alpha=al,
                                     shrinkA=0, shrinkB=0), zorder=3)
        ax2.plot([f["delta"]], [a], "o", color=c_sup, ms=5, mec=t["fondo"], mew=1.2,
                 alpha=al, zorder=4)
    ax2.annotate("π/2: el perfil fundamental", (0.02, 0.86), xycoords="axes fraction",
                 color=t["tinta2"], fontsize=9)
    ax2.set_ylabel("u(h) / ⟨u⟩")
    ax2.set_xscale("log")
    ax2.set_xlim(left=0.1)
    ax2.set_xlabel("δ medido en la ventana de ajuste  (k horizontal del propio campo)")
    ax2.annotate("flecha: del inicio al final de la ventana", (0.02, 0.06),
                 xycoords="axes fraction", color=t["tinta2"], fontsize=9)
    return guardar(fig, "f8_superficie_delta_%s" % caso if caso else "f8_superficie_delta")


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--caso", default="forzado")
    args = ap.parse_args()
    with open(os.path.join(TABLAS, "barrido_delta_etapa1_%s.json" % args.caso)) as f:
        d = json.load(f)
    r1 = fig_apartamiento(d)
    print("escrito %s" % os.path.relpath(r1, RAIZ))
    ruta_s = os.path.join(TABLAS, "barrido_delta_etapa1_%s_superficie.json" % args.caso)
    if os.path.exists(ruta_s):
        with open(ruta_s) as f:
            ds = json.load(f)
        r2 = fig_superficie(ds, d)
        print("escrito %s" % os.path.relpath(r2, RAIZ))

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
