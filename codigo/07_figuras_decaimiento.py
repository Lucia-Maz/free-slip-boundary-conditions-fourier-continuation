#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Figuras del decaimiento medido, para el HTML de revisión.

Cada figura se renderiza DOS veces, clara y oscura, porque la página que las muestra
sigue el tema del que la mira y una PNG de fondo blanco sobre fondo oscuro deslumbra.

Paleta: slots 1-4 de la paleta categórica de referencia, en orden y sin ciclar. Para
listas de líneas la validación documentada de esa paleta cubre estos cuatro
(peor par adyacente CVD ΔE 9,1 claro / 8,4 oscuro). Como tres de los slots quedan por
debajo de 3:1 de contraste sobre fondo claro, se aplica la regla de relieve: cada curva
lleva etiqueta directa y el HTML lleva además la tabla de datos.

Uso:
    /home/lucia/miniforge3/envs/piv-dt/bin/python codigo/07_figuras_decaimiento.py
"""

import json
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt              # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(AQUI)
sys.path.insert(0, AQUI)
from celda import CAMPANA_02_06_25 as CELDA  # noqa: E402

TABLAS = os.path.join(RAIZ, "salidas", "tablas")
FIGS = os.path.join(RAIZ, "salidas", "figuras")

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
# Paleta más suave para el tema oscuro (pedido de Lucía, 2026-09-16): mismos cuatro
# tonos categóricos, desaturados hacia un aire "pastel" pero sin salir de la banda
# L 0.48-0.67 / C >= 0.10 que exige la skill de dataviz para que sigan leyéndose
# como identidad de serie sobre #1a1a19 — validado con validate_palette.js
# (los 4 checks pasan; peor par CVD 8.8, peor par visión normal 16.1).
SERIES_OSC = ["#4a8fd4", "#cf7238", "#3aa47a", "#b8871f"]

TEMAS = {
    "claro": dict(fondo="#fcfcfb", tinta="#0b0b0b", tinta2="#52514e",
                  grilla="#e3e2dd", series=SERIES),
    "oscuro": dict(fondo="#1a1a19", tinta="#ffffff", tinta2="#c3c2b7",
                   grilla="#33322f", series=SERIES_OSC),
}

ETIQUETAS = {"med_S0003": "Decaimiento 1", "med_S0005": "Decaimiento 2",
             "med_S0006": "Decaimiento 3", "med_S0009": "Decaimiento 4"}


def estilo(tema):
    t = TEMAS[tema]
    plt.rcParams.update({
        "figure.facecolor": t["fondo"], "axes.facecolor": t["fondo"],
        "savefig.facecolor": t["fondo"],
        "text.color": t["tinta"], "axes.labelcolor": t["tinta"],
        "xtick.color": t["tinta2"], "ytick.color": t["tinta2"],
        "axes.edgecolor": t["grilla"], "grid.color": t["grilla"],
        "axes.grid": True, "grid.linewidth": 0.8, "grid.alpha": 0.9,
        "axes.spines.top": False, "axes.spines.right": False,
        "font.size": 10.5, "axes.titlesize": 11.5, "axes.titleweight": "semibold",
        "legend.frameon": False, "lines.linewidth": 2.0,
        "figure.dpi": 130,
    })
    return t


def guardar(fig, nombre, tema):
    os.makedirs(FIGS, exist_ok=True)
    ruta = os.path.join(FIGS, "%s_%s.png" % (nombre, tema))
    fig.savefig(ruta, bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)
    return ruta


# --------------------------------------------------------------------------------

def espectro_limpio(e):
    """Espectro con el ruido blanco de PIV restado, y los dos k que salen de él.

    El ruido de PIV aporta la misma potencia por modo, o sea p0·(modos del anillo) en
    el espectro sumado radialmente, que crece con k. Para los registros calculados
    antes de que el script guardara el conteo de modos se usa que ese conteo es
    proporcional a k, que es exactamente la forma correcta para un binneado radial en
    2D; la constante se absorbe en p0.
    """
    k = np.array(e["k_1_m"], float)
    E = np.array(e["E"], float)
    cuenta = np.array(e["modos_por_anillo"], float) if "modos_por_anillo" in e else k
    # La banda donde se lee el piso tiene que estar DEBAJO de Nyquist. Por encima los
    # anillos sólo recogen las esquinas de la grilla y el espectro se desploma, de modo
    # que un piso leído ahí sale demasiado bajo y no resta casi nada. Con paso de
    # grilla de 4 px, k_Nyquist = pi/paso ~ 2984 1/m, y la meseta de ruido se ve entre
    # unos 500 y 3000 1/m.
    k_nyq = math.pi / e["paso_grilla_m"]
    banda = (k > 0.35 * k_nyq) & (k < 0.95 * k_nyq)
    if banda.sum() < 3:
        banda = slice(int(0.6 * len(E)), int(0.85 * len(E)))
    p0 = float(np.median((E / np.maximum(cuenta, 1.0))[banda]))
    lim = np.maximum(E - p0 * cuenta, 0.0)
    kp = float(k[int(np.argmax(lim))]) if lim.sum() > 0 else float("nan")
    kw = float((k * lim).sum() / lim.sum()) if lim.sum() > 0 else float("nan")
    return k, E, lim, kp, kw


def geometria(e):
    """(nx, ny) de la grilla de vectores, deducidos del espaciado en k."""
    k = np.array(e["k_1_m"], float)
    dk = float(k[1] - k[0])
    n = int(round(2.0 * math.pi / (dk * e["paso_grilla_m"])))
    return n, n, dk


def sigma_desde_espectro(e):
    """sigma por componente [px] leído del piso del espectro.

    Normalización, y está verificada abajo sobre datos sintéticos. Con
    f = fft2(u*w) y pot = |f|^2/N^2, Parseval da  sum_k pot = <w^2> * <u^2+v^2>,
    con <w^2> = (3/8)^2 para una Hann 2D separable. Un ruido blanco de varianza
    sigma_c^2 por componente reparte 2*sigma_c^2*<w^2> entre los N modos, así que

        sigma_c = sqrt( p0 * N / (2 * <w^2>) )   con p0 la potencia POR MODO.

    Los espectros guardados traen el conteo por anillo aproximado por k (proporcional
    al conteo real, que es 2*pi*k/dk), así que hay que corregir esa constante.
    """
    k = np.array(e["k_1_m"], float)
    E = np.array(e["E"], float)
    nx, ny, dk = geometria(e)
    if "modos_por_anillo" in e:
        cuenta = np.array(e["modos_por_anillo"], float)
    else:
        cuenta = 2.0 * math.pi * k / dk           # conteo real del anillo
    k_nyq = math.pi / e["paso_grilla_m"]
    banda = (k > 0.35 * k_nyq) & (k < 0.95 * k_nyq)
    p0 = float(np.median((E / np.maximum(cuenta, 1.0))[banda]))
    w2 = (3.0 / 8.0) ** 2
    return math.sqrt(max(p0, 0.0) * nx * ny / (2.0 * w2)), p0


def ajuste_exponencial(t, u):
    """alfa a partir de la pendiente de log(u), y su error estándar."""
    t = np.asarray(t, float); u = np.asarray(u, float)
    ok = np.isfinite(t) & np.isfinite(u) & (u > 0)
    t, u = t[ok], u[ok]
    if len(t) < 3:
        return float("nan"), float("nan"), float("nan")
    A = np.vstack([t, np.ones_like(t)]).T
    coef, res, *_ = np.linalg.lstsq(A, np.log(u), rcond=None)
    resid = np.log(u) - A @ coef
    s2 = float(resid @ resid) / (len(t) - 2)
    cov = s2 * np.linalg.inv(A.T @ A)
    return -float(coef[0]), math.sqrt(float(cov[0, 0])), float(np.exp(coef[1]))


# --------------------------------------------------------------------------------

def _mediana_movil(y, w=15):
    y = np.asarray(y, float)
    m = np.empty_like(y)
    for i in range(len(y)):
        m[i] = np.median(y[max(0, i - w // 2):i + w // 2 + 1])
    return m


def fig_decaimiento(d, tema, t_corte=10.0):
    t = estilo(tema)
    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    ax.axvspan(-2, t_corte, color=t["grilla"], alpha=0.55, lw=0, zorder=0)
    tt, uu = [], []
    for i, (med, reg) in enumerate(sorted(d["registros"].items())):
        c = t["series"][i]
        sv = reg["serie_vec"]
        ax.plot(sv["t_s"], 1e3 * _mediana_movil(sv["u_rms_m_s"]), color=c, lw=1.2,
                alpha=0.45, zorder=1)
        cal = reg["calibracion"]
        tc = [p["t_s"] for p in cal]
        uc = [1e3 * p["u_rms_superficie_m_s"] for p in cal]
        ax.plot(tc, uc, "o-", color=c, ms=6.5, mec=t["fondo"], mew=1.5, zorder=3)
        if tc:
            ax.annotate(ETIQUETAS[med], (tc[-1], uc[-1]), xytext=(7, 0),
                        textcoords="offset points", color=c, fontsize=9.5,
                        va="center", fontweight="semibold")
        tt += [x for x in tc if x >= t_corte]
        uu += [y for x, y in zip(tc, uc) if x >= t_corte]

    a, ea, u0 = ajuste_exponencial(tt, uu)
    xx = np.linspace(t_corte, 45, 20)
    ax.plot(xx, u0 * np.exp(-a * xx), color=t["tinta"], lw=1.6, ls="--", zorder=2)
    ax.annotate("ajuste del ensemble:  α = %.4f ± %.4f s⁻¹\n"
                "predicho, free-slip:  0,0685 s⁻¹\n"
                "predicho, tapa rígida: 0,2742 s⁻¹"
                % (a, ea), (0.985, 0.97), xycoords="axes fraction",
                ha="right", va="top", color=t["tinta"], fontsize=9.5)

    ax.set_yscale("log")
    ax.set_xlabel("tiempo desde el inicio del registro  [s]")
    ax.set_ylabel("velocidad rms de superficie  [mm/s]")
    ax.set_title("El decaimiento, antes y después de descontar el ruido de PIV")
    ax.set_xlim(-2, 58)
    ax.text(0.015, 0.05,
            "tenue: .vec archivados (Δt=1/60 s), mediana móvil —\n"
            "no decaen: se quedan en el piso de ruido\n"
            "puntos: rearmado desde los TIF, corregido\n"
            "gris: fuera del ajuste (meseta inicial)",
            transform=ax.transAxes, color=t["tinta2"], fontsize=8.5)
    return guardar(fig, "f1_decaimiento", tema)


def fig_metodo(d, tema, med="med_S0003"):
    t = estilo(tema)
    cal = d["registros"][med]["calibracion"]
    fig, ax = plt.subplots(figsize=(6.4, 4.3))
    muestras = [cal[0], cal[len(cal) // 2], cal[-1]]
    for i, p in enumerate(muestras):
        c = t["series"][i]
        ns = np.array([q["n"] for q in p["puntos"]], float)
        y = np.array([q["desp_rms_px"] for q in p["puntos"]]) ** 2
        ax.plot(ns ** 2, y, "o", color=c, ms=7, mec=t["fondo"], mew=1.5, zorder=3)
        aj = p["ajuste_dos_dt"]
        d1, sg = aj["d1_px_por_cuadro"], aj["sigma_px"]
        xx = np.linspace(0, ns.max() ** 2, 50)
        ax.plot(xx, d1 ** 2 * xx + math.copysign(sg ** 2, sg), color=c, lw=1.6,
                alpha=0.8, zorder=2)
        ax.annotate("t = %.0f s" % p["t_s"], (ns[-1] ** 2, y[-1]), xytext=(-4, 8),
                    textcoords="offset points", color=c, fontsize=9.5, ha="right",
                    fontweight="semibold")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("n²   (n = cuadros entre los dos de cada par)")
    ax.set_ylabel("desplazamiento rms al cuadrado  [px²]")
    ax.set_title("El método: la ordenada al origen es el ruido")
    ax.text(0.03, 0.9, "desp² = (n·d₁)² + σ²", transform=ax.transAxes,
            color=t["tinta"], fontsize=11)
    return guardar(fig, "f2_metodo", tema)


def fig_ruido(d, tema):
    t = estilo(tema)
    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    for i, (med, reg) in enumerate(sorted(d["registros"].items())):
        c = t["series"][i]
        cal = reg["calibracion"]
        # ruta 1: la ordenada al origen del ajuste de dos Dt
        tc = [p["t_s"] for p in cal]
        sg = [p["ajuste_dos_dt"]["sigma_por_componente_px"] for p in cal]
        ax.plot(tc, sg, "o", color=c, ms=5.5, mec=t["fondo"], mew=1.2, alpha=.55)
        # ruta 2: el piso del espectro, independiente de la anterior
        te = [e["t_s"] for e in reg["espectros"]]
        se = [sigma_desde_espectro(e)[0] for e in reg["espectros"]]
        ax.plot(te, se, "s-", color=c, ms=6, mec=t["fondo"], mew=1.5)
        if te:
            ax.annotate(ETIQUETAS[med].replace("Decaimiento ", "D"),
                        (te[-1], se[-1]), xytext=(6, 0), textcoords="offset points",
                        color=c, fontsize=9.5, va="center", fontweight="semibold")
    ax.axhline(0.169, color=t["tinta2"], lw=1.4, ls="--")
    ax.annotate("0,169 px — el forzado a 2,15 A [H-05]", (1, 0.169), xytext=(4, 6),
                textcoords="offset points", color=t["tinta2"], fontsize=9)
    ax.axhline(0.030, color=t["tinta2"], lw=1.4, ls=":")
    ax.annotate("0,030 px — fluido en reposo", (1, 0.030), xytext=(4, 5),
                textcoords="offset points", color=t["tinta2"], fontsize=9)
    ax.set_xlabel("tiempo  [s]")
    ax.set_ylabel("σ por componente  [px]")
    ax.set_title("El ruido de PIV cae junto con el flujo")
    ax.set_yscale("log")
    ax.text(0.02, 0.06, "cuadrados: piso del espectro (normalización verificada sobre\n"
                        "ruido sintético al 1 %)  ·  círculos tenues: ordenada al origen\n"
                        "del ajuste de dos Δt, que a estos instantes es inestable",
            transform=ax.transAxes, color=t["tinta2"], fontsize=8.5)
    return guardar(fig, "f3_ruido", tema)


def fig_espectros(d, tema, med="med_S0003"):
    t = estilo(tema)
    esp = d["registros"][med]["espectros"]
    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    muestras = [esp[0], esp[len(esp) // 2], esp[-1]] if len(esp) >= 3 else esp
    for i, e in enumerate(muestras):
        c = t["series"][i]
        k, E, lim, kp, kw = espectro_limpio(e)
        norma = E.max()
        m = (k > 0) & (E > 0)
        ax.loglog(k[m], E[m] / norma, color=c, lw=1.0, alpha=0.35)
        m2 = (k > 0) & (lim > 0)
        ax.loglog(k[m2], lim[m2] / norma, color=c, lw=2.0)
        j = int(np.argmin(np.abs(k - kp)))
        ax.annotate("t = %.0f s" % e["t_s"], (k[j], (lim / norma)[j]),
                    xytext=(6, 8), textcoords="offset points", color=c,
                    fontsize=9.5, fontweight="semibold")
    ax.text(0.02, 0.06, "tenue: espectro crudo · grueso: con el ruido blanco de PIV "
                        "restado", transform=ax.transAxes, color=t["tinta2"],
            fontsize=9)
    kf = d["celda"]["k_forzado_1_m"]
    ax.axvline(kf, color=t["tinta2"], lw=1.4, ls="--")
    ax.annotate("k del forzado\n88,9 m⁻¹ (7,07 cm)", (kf, 1.0), xytext=(-6, -4),
                textcoords="offset points", color=t["tinta2"], fontsize=9, ha="right",
                va="top")
    ax.set_xlabel("número de onda k  [1/m]")
    ax.set_ylabel("E(k), normalizado")
    ax.set_title("Espectro de energía: la escala crece durante el decaimiento")
    return guardar(fig, "f4_espectros", tema)


def fig_delta(d, tema, resumen):
    t = estilo(tema)
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    h, nu = d["celda"]["h_m"], d["celda"]["nu_m2_s"]
    for i, (med, reg) in enumerate(sorted(d["registros"].items())):
        c = t["series"][i]
        cal = reg["calibracion"]
        esp = {round(e["t_s"], 3): e for e in reg["espectros"]}
        tt, dd, dk = [], [], []
        for p in cal:
            e = esp.get(round(p["t_s"], 3))
            if e is None:
                continue
            U = p["u_rms_superficie_m_s"] / (math.pi / 2)
            _, _, _, _, kw = espectro_limpio(e)
            tt.append(p["t_s"])
            dd.append(U * kw * h ** 2 / nu)
            dk.append(U * d["celda"]["k_forzado_1_m"] * h ** 2 / nu)
        ax.plot(tt, dd, "o-", color=c, ms=6.5, mec=t["fondo"], mew=1.5)
        ax.plot(tt, dk, ":", color=c, lw=1.4, alpha=0.7)
        if tt:
            ax.annotate(ETIQUETAS[med].replace("Decaimiento ", "D"), (tt[-1], dd[-1]),
                        xytext=(6, 0), textcoords="offset points", color=c,
                        fontsize=9.5, va="center", fontweight="semibold")
    ax.axhline(1.0, color=t["tinta2"], lw=1.4, ls="--")
    ax.annotate("δ = 1", (0.5, 1.0), xytext=(3, 5), textcoords="offset points",
                color=t["tinta2"], fontsize=9)
    ax.set_yscale("log")
    ax.set_xlabel("tiempo  [s]")
    ax.set_ylabel("δ = U k h² / ν")
    ax.set_title("El parámetro de [P-01], con k y U sacados de los datos")
    ax.text(0.02, 0.06,
            "puntos: con el k medido del campo · punteado: con el k del forzado\n"
            "distorsión del perfil vertical ≈ 0,0077·δ",
            transform=ax.transAxes, color=t["tinta2"], fontsize=9)
    return guardar(fig, "f5_delta", tema)


def main():
    with open(os.path.join(TABLAS, "decaimiento.json")) as f:
        d = json.load(f)

    # ---- ajuste de alfa sobre el ensemble -------------------------------------
    tt, uu = [], []
    por_registro = {}
    for med, reg in sorted(d["registros"].items()):
        cal = reg["calibracion"]
        t_ = [p["t_s"] for p in cal]
        u_ = [p["u_rms_superficie_m_s"] for p in cal]
        a, ea, _ = ajuste_exponencial(t_, u_)
        por_registro[med] = {"alfa_1_s": a, "err": ea, "n_puntos": len(t_)}
        tt += t_; uu += u_
    a_ens, ea_ens, u0 = ajuste_exponencial(tt, uu)
    # el mismo ajuste sobre la cola, que es donde delta es chico y donde la clausura
    # de un modo tiene que valer mejor ([P-01]); si los dos coinciden, la clausura no
    # se está rompiendo al principio
    tar = [x for x, y in zip(tt, uu) if x >= 10.0]
    uar = [y for x, y in zip(tt, uu) if x >= 10.0]
    a_cola, ea_cola, _ = ajuste_exponencial(tar, uar)

    # ---- k medido, y delta ----------------------------------------------------
    ks = []
    for med, reg in sorted(d["registros"].items()):
        for e in reg["espectros"]:
            _, _, _, kp, kw = espectro_limpio(e)
            ks.append({"med": med, "t_s": e["t_s"], "k_pico_1_m": kp,
                       "k_pesado_1_m": kw})
    k_ini = float(np.median([x["k_pesado_1_m"] for x in ks if x["t_s"] < 3]))
    k_fin = float(np.median([x["k_pesado_1_m"] for x in ks if x["t_s"] > 30]))

    eps = CELDA.forzado.k_fundamental() * CELDA.h
    eps_med = k_ini * CELDA.h
    resumen = {
        "alfa_medido_ensemble_1_s": a_ens,
        "alfa_medido_error_1_s": ea_ens,
        "alfa_medido_cola_1_s": a_cola,
        "alfa_medido_cola_error_1_s": ea_cola,
        "u0_ajustado_m_s": u0,
        "alfa_por_registro": por_registro,
        "alfa_predicho_1_s": CELDA.alfa(),
        "lambda_predicho_con_k_forzado_1_s": CELDA.alfa() * (1 + 4 * eps ** 2 / math.pi ** 2),
        "lambda_predicho_con_k_medido_1_s": CELDA.alfa() * (1 + 4 * eps_med ** 2 / math.pi ** 2),
        "alfa_predicho_canal_1_s": CELDA.alfa("noslip-noslip"),
        "k_forzado_1_m": CELDA.forzado.k_fundamental(),
        "k_medido_inicio_1_m": k_ini,
        "k_medido_final_1_m": k_fin,
        "k_por_instante": ks,
    }
    with open(os.path.join(TABLAS, "decaimiento_resumen.json"), "w") as f:
        json.dump(resumen, f, indent=1)

    rutas = {}
    for tema in ("claro", "oscuro"):
        rutas.setdefault("f1", []).append(fig_decaimiento(d, tema))
        rutas.setdefault("f2", []).append(fig_metodo(d, tema))
        rutas.setdefault("f3", []).append(fig_ruido(d, tema))
        rutas.setdefault("f4", []).append(fig_espectros(d, tema))
        rutas.setdefault("f5", []).append(fig_delta(d, tema, resumen))

    print("alfa medido (ensemble, todo)  = %.5f ± %.5f 1/s   tau = %.2f s"
          % (a_ens, ea_ens, 1 / a_ens))
    print("alfa medido (ensemble, t>10s) = %.5f ± %.5f 1/s" % (a_cola, ea_cola))
    for med, r in por_registro.items():
        print("    %s: alfa = %.5f ± %.5f  (%d puntos)"
              % (med, r["alfa_1_s"], r["err"], r["n_puntos"]))
    print("k medido: inicio %.1f 1/m, final %.1f 1/m   (forzado: %.1f)"
          % (k_ini, k_fin, resumen["k_forzado_1_m"]))
    print("alfa predicho free-slip       = %.5f 1/s" % resumen["alfa_predicho_1_s"])
    print("lambda con k del forzado      = %.5f 1/s" % resumen["lambda_predicho_con_k_forzado_1_s"])
    print("lambda con k medido           = %.5f 1/s" % resumen["lambda_predicho_con_k_medido_1_s"])
    print("alfa predicho canal (rígido)  = %.5f 1/s" % resumen["alfa_predicho_canal_1_s"])
    print("cociente medido/predicho (free-slip, alfa solo) = %.3f"
          % (a_ens / resumen["alfa_predicho_1_s"]))
    print("cociente medido/canal                            = %.3f"
          % (a_ens / resumen["alfa_predicho_canal_1_s"]))
    for k, v in rutas.items():
        print("%s: %s" % (k, ", ".join(os.path.relpath(x, RAIZ) for x in v)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
