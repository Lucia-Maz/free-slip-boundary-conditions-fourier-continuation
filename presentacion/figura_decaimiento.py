#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Figura de los decaimientos para la presentación: la misma medición que f1 de
codigo/07_figuras_decaimiento.py ([V-13]), en la paleta azul de las diapositivas, con
las pendientes que predicen las dos condiciones de contorno encima.

Lee salidas/tablas/decaimiento.json (no hace falta el disco externo) y usa el mismo
ajuste (ajuste_exponencial de 07) y los mismos alfa de codigo/celda.py.

Uso, desde la raíz del repositorio:
    /home/lucia/miniforge3/envs/piv-dt/bin/python presentacion/figura_decaimiento.py
"""

import importlib.util
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt              # noqa: E402
import matplotlib.ticker                      # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(AQUI)
sys.path.insert(0, os.path.join(RAIZ, "codigo"))
from celda import CAMPANA_02_06_25 as CELDA  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "figs07", os.path.join(RAIZ, "codigo", "07_figuras_decaimiento.py"))
figs07 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(figs07)

AZUL = "#1F4E79"
AZUL_CLARO = "#7FA7CF"
AZUL_FONDO = "#EAF1F8"
TINTA = "#1E2A38"
T_CORTE = 10.0
# k del flujo en la ventana del ajuste ([P-02]): los bins donde pica el espectro limpio en
# t > 10 s, de salidas/tablas/decaimiento_resumen.json ("k_por_instante", k_pico_1_m). De
# 20 espectros, 17 pican en 59,0 / 82,6 / 129,8 m^-1; los otros tres son tempranos (t <= 17 s):
# dos todavía en el fundamental de la red, 294,9, y uno en 176,9.
K_VENTANA = (58.98, 129.80)


def main():
    with open(os.path.join(RAIZ, "salidas", "tablas", "decaimiento.json")) as f:
        d = json.load(f)

    plt.rcParams.update({"font.family": "Inter", "font.size": 10,
                         "axes.edgecolor": AZUL, "axes.labelcolor": TINTA,
                         "xtick.color": TINTA, "ytick.color": TINTA,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(4.9, 3.25))
    ax.axvspan(-2, T_CORTE, color=AZUL_FONDO, lw=0, zorder=0)
    ax.text(T_CORTE / 2 - 1, 0.105, "forzado\nencendido", ha="center", va="bottom",
            color=AZUL, fontsize=8.5)

    tt, uu = [], []
    for med, reg, mk in zip(*zip(*sorted(d["registros"].items())), "osD^"):
        cal = reg["calibracion"]
        tc = np.array([p["t_s"] for p in cal])
        uc = 1e3 * np.array([p["u_rms_superficie_m_s"] for p in cal])
        ax.plot(tc, uc, "-", color=AZUL_CLARO, lw=1.0, zorder=2)
        ax.plot(tc, uc, mk, color=AZUL_CLARO, ms=6, mec="white", mew=1.0, zorder=3)
        tt += list(tc[tc >= T_CORTE]); uu += list(uc[tc >= T_CORTE])

    a, ea, u0 = figs07.ajuste_exponencial(tt, uu)
    u10 = u0 * np.exp(-a * T_CORTE)
    xx = np.linspace(T_CORTE, 45, 60)
    coma = lambda x, n=3: ("%.*f" % (n, x)).replace(".", ",")
    # Lo medido es una tasa total, lambda; la predicción comparable es alfa + nu k^2 con el k
    # que el flujo tiene en la ventana (bandas), no alfa pelado (k -> 0, línea de trazos).
    nu = CELDA.fluido.nu
    a_libre, a_rigida = CELDA.alfa(), CELDA.alfa("noslip-noslip")
    lam = lambda a_, k: a_ + nu * k ** 2
    # la figura muestra las pendientes alfa (k -> 0); los intervalos alfa + nu k^2 van en la
    # recta de la slide (presentacion.tex), con los valores que imprime este script
    ax.plot(xx, u0 * np.exp(-a * xx), color=AZUL, lw=2.4, zorder=5,
            label="tasa medida: %s ± %s s⁻¹" % (coma(a), coma(ea)))
    ax.plot(xx, u10 * np.exp(-a_libre * (xx - T_CORTE)), color=AZUL, lw=1.4,
            ls=(0, (4, 3)), zorder=4, label="superficie libre, α: %s s⁻¹" % coma(a_libre))
    xr = np.linspace(T_CORTE, 30, 50)
    ax.plot(xr, u10 * np.exp(-a_rigida * (xr - T_CORTE)), color=AZUL, lw=1.4,
            ls=(0, (1, 2)), zorder=4, label="tapa rígida, α: %s s⁻¹" % coma(a_rigida))
    ax.legend(loc="upper right", frameon=False, fontsize=8.5, handlelength=2.4,
              title="pendientes desde t = 10 s", title_fontsize=8, labelcolor=AZUL)

    ax.set_yscale("log")
    ax.set_ylim(0.1, 6)
    ax.set_xlim(-2, 50)
    ax.set_yticks([0.1, 0.3, 1, 3])
    ax.set_yticklabels(["0,1", "0,3", "1", "3"])
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel("tiempo [s]")
    ax.set_ylabel("velocidad rms de superficie [mm/s]")
    fig.tight_layout()
    salida = os.path.join(AQUI, "figuras", "decaimiento_azul.pdf")
    fig.savefig(salida)
    print("tasa medida, ensemble t>10 s = %.5f ± %.5f (V-13: 0,06687 ± 0,00488)" % (a, ea))
    print("alfa (k->0): libre %.5f, rígida %.5f" % (a_libre, a_rigida))
    print("lambda = alfa + nu k^2, k = %.1f-%.1f 1/m: libre %.4f-%.4f, rígida %.4f-%.4f" % (
        K_VENTANA[0], K_VENTANA[1], lam(a_libre, K_VENTANA[0]), lam(a_libre, K_VENTANA[1]),
        lam(a_rigida, K_VENTANA[0]), lam(a_rigida, K_VENTANA[1])))
    print("escribí", salida)
    return 0


if __name__ == "__main__":
    sys.exit(main())
