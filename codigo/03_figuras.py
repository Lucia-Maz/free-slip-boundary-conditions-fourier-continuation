"""Figuras del barrido de Δt, y el ajuste del que sale σ.

Produce dos figuras a partir de `salidas/tablas/barrido_dt_*.json`:

  fig1  U_med² contra n⁻², que debe ser una recta mientras el modelo valga.
        Ordenada al origen = velocidad verdadera; pendiente = σ².
  fig2  la función de costo: error relativo contra Δt, con las dos ramas que la
        componen (ruido, que cae con Δt, y pérdida de correlación, que crece).

El ajuste se hace sólo sobre los n donde el modelo puede valer, y ese recorte es una
decisión explícita: se excluyen los n con fracción válida por debajo de un umbral,
porque ahí ya no se está midiendo el mismo campo.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

AQUI = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Valor publicado contra el que se contrasta, para ventana de 16x16 px.
# Foucaut, Carlier & Stanislas 2004, DOI 10.1088/0957-0233/15/6/003, tabla 1.
SIGMA_FOUCAUT_PX = 0.087
FRAC_VALIDA_MINIMA = 0.90


def ajustar(n, desp, mascara, escala, fps, n_boot=2000, semilla=0):
    """Ajusta desp(n) = sqrt( (U·s/fps)²n² + σ² ) con pesos RELATIVOS.

    Por qué no un ajuste lineal sin pesos. El modelo es lineal en n² con ordenada al
    origen σ², pero el rango de n² va de 1 a 1024: mínimos cuadrados sin pesos queda
    dominado por el punto más grande y σ² sale como diferencia de números grandes —
    en estos datos da **negativa**, que es física imposible. La incerteza de un
    desplazamiento medido es aproximadamente proporcional al desplazamiento, no
    constante, así que corresponde pesar por el valor. Ver [D-17].

    La barra de error de σ sale por bootstrap sobre los puntos, porque es la única
    forma honesta de mostrar cuán poco constriñe este diseño al ruido: la señal a n
    chico es comparable a σ, y a n grande σ es despreciable y no aporta información.
    """
    x, y = n[mascara], desp[mascara]

    def modelo(nn, c, sig):
        return np.sqrt((c * nn) ** 2 + sig ** 2)

    def una_vez(xx, yy):
        from scipy.optimize import curve_fit
        p0 = [yy[-1] / xx[-1], 0.15]
        try:
            (c, sig), _ = curve_fit(modelo, xx, yy, p0=p0, sigma=yy,
                                    absolute_sigma=False, maxfev=20000)
            return abs(c), abs(sig)
        except Exception:
            return np.nan, np.nan

    c, sig = una_vez(x, y)
    rng = np.random.default_rng(semilla)
    muestras = np.array([una_vez(*(lambda i: (x[i], y[i]))(
        rng.integers(0, len(x), len(x)))) for _ in range(n_boot)])
    muestras = muestras[np.isfinite(muestras).all(axis=1)]

    pred = modelo(x, c, sig)
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    # sigma despejado punto a punto, que muestra la dispersion cruda
    por_punto = np.sqrt(np.clip(y ** 2 - (c * x) ** 2, 0, None))
    return {
        "U_verdadera_m_s": float(c * fps / escala),
        "U_verdadera_ic68_m_s": [float(np.percentile(muestras[:, 0], 16) * fps / escala),
                                 float(np.percentile(muestras[:, 0], 84) * fps / escala)],
        "sigma_px": float(sig),
        "sigma_ic68_px": [float(np.percentile(muestras[:, 1], 16)),
                          float(np.percentile(muestras[:, 1], 84))],
        "sigma_por_punto_px": {int(nn): round(float(ss), 4) for nn, ss in zip(x, por_punto)},
        "R2": 1 - ss_res / ss_tot if ss_tot > 0 else float("nan"),
        "pendiente": float(c ** 2), "ordenada": float(sig ** 2),
        "n_usados": x.tolist(), "n_bootstrap": int(len(muestras)),
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--carpeta", default="med_S0008")
    args = p.parse_args()

    ruta = os.path.join(AQUI, "salidas", "tablas", f"barrido_dt_{args.carpeta}.json")
    with open(ruta, encoding="utf-8") as fh:
        d = json.load(fh)
    fps, escala = d["fps"], d["escala_px_por_m"]

    f = d["filas"]
    n = np.array([r["n"] for r in f], float)
    dt = np.array([r["dt_s"] for r in f])
    desp = np.array([r["desp_rms_px"] for r in f])
    val = np.array([r["frac_valida"] for r in f])
    # el desplazamiento se convierte a velocidad acá, no en el barrido ([D-09])
    U = desp / (dt * escala)

    # El ajuste usa el desplazamiento medido en pixeles, no en m/s: sigma es una
    # propiedad de la correlacion. U_med^2 = U_true^2 + (sigma/(dt*s))^2
    mascara = val >= FRAC_VALIDA_MINIMA
    aj = ajustar(n, desp, mascara, escala, fps)
    aj["sigma_foucaut_px"] = SIGMA_FOUCAUT_PX
    aj["desviacion_relativa_vs_foucaut"] = (
        abs(aj["sigma_px"] - SIGMA_FOUCAUT_PX) / SIGMA_FOUCAUT_PX)
    aj["frac_valida_minima"] = FRAC_VALIDA_MINIMA
    aj["n_excluidos_por_validez"] = n[~mascara].tolist()
    # control: el punto n=1 tiene el brazo de palanca mas largo en 1/n^2 y podria
    # dominar el ajuste por si solo. Se repite sin el para ver cuanto cambia.
    aj["sin_n1"] = ajustar(n, desp, mascara & (n > 1), escala, fps)

    # ---- figura 1: la recta, en desplazamiento -----------------------------
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    ax.plot(n[mascara] ** 2, desp[mascara] ** 2, "o", color="#1f4e79",
            label="usado en el ajuste")
    if (~mascara).any():
        ax.plot(n[~mascara] ** 2, desp[~mascara] ** 2, "o", mfc="none", color="#b03030",
                label=f"validez < {FRAC_VALIDA_MINIMA:.0%}, excluido")
    xr = np.linspace(0, (n ** 2).max() * 1.05, 100)
    ax.plot(xr, aj["pendiente"] * xr + aj["ordenada"], "-", color="#e8a05c",
            label=(f"ajuste: U={aj['U_verdadera_m_s']*1e3:.3f} mm/s, "
                   f"σ={aj['sigma_px']:.3f} px (R²={aj['R2']:.4f})"))
    ax.axhline(aj["ordenada"], ls=":", color="grey", lw=1)
    ax.annotate(f"σ² = {aj['ordenada']:.4f} px²  →  σ = {aj['sigma_px']:.3f} px",
                (0.02, aj["ordenada"]), xycoords=("axes fraction", "data"),
                textcoords="offset points", xytext=(0, 6), fontsize=8, color="#444")
    ax.set_xlabel(r"$n^2$   (n = cuadros entre las dos imágenes)")
    ax.set_ylabel(r"desplazamiento$^2$   [px$^2$]")
    ax.set_title(f"{d['descripcion']} — {d['carpeta']}\n"
                 r"$d^2 = (U\,s/\mathrm{fps})^2 n^2 + \sigma^2$", fontsize=10)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=.3)
    # recuadro con la zona de n chico, donde vive el ruido
    eje = ax.inset_axes([0.58, 0.12, 0.38, 0.34])
    ch = mascara & (n <= 12)
    eje.plot(n[ch] ** 2, desp[ch] ** 2, "o", color="#1f4e79", ms=4)
    xr2 = np.linspace(0, 150, 50)
    eje.plot(xr2, aj["pendiente"] * xr2 + aj["ordenada"], "-", color="#e8a05c")
    eje.axhline(aj["ordenada"], ls=":", color="grey", lw=1)
    eje.set_title("n ≤ 12", fontsize=7)
    eje.tick_params(labelsize=6)
    fig.tight_layout()
    f1 = os.path.join(AQUI, "salidas", "figuras", f"fig1_recta_{args.carpeta}.png")
    fig.savefig(f1, dpi=150)
    plt.close(fig)

    # ---- figura 2: la función de costo --------------------------------------
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    Ut = aj["U_verdadera_m_s"]
    err_ruido = np.abs(U - Ut) / Ut
    ax.plot(dt, 100 * err_ruido, "o-", color="#1f4e79", label="error relativo en U")
    ax2 = ax.twinx()
    ax2.plot(dt, 100 * val, "s--", color="#b03030", ms=4, label="vectores válidos")
    ax2.set_ylabel("vectores válidos [%]", color="#b03030")
    ax2.tick_params(axis="y", colors="#b03030")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\Delta t$ [s]")
    ax.set_ylabel("error relativo en U [%]", color="#1f4e79")
    ax.tick_params(axis="y", colors="#1f4e79")
    # los Δt que Lucía usó realmente
    for nn, etiqueta in ((1, "usado: dt=1/60"), (5, "usado: dt5")):
        ax.axvline(nn / fps, color="grey", ls=":", lw=1)
        ax.annotate(etiqueta, (nn / fps, ax.get_ylim()[1]), rotation=90,
                    fontsize=7, va="top", ha="right", color="grey")
    ax.set_title("Las dos ramas del costo: ruido a Δt chico, decorrelación a Δt grande",
                 fontsize=10)
    ax.grid(alpha=.3)
    fig.tight_layout()
    f2 = os.path.join(AQUI, "salidas", "figuras", f"fig2_costo_{args.carpeta}.png")
    fig.savefig(f2, dpi=150)
    plt.close(fig)

    destino = os.path.join(AQUI, "salidas", "tablas", f"ajuste_{args.carpeta}.json")
    with open(destino, "w", encoding="utf-8") as fh:
        json.dump(aj, fh, indent=1, ensure_ascii=False)

    ic = aj["U_verdadera_ic68_m_s"]; ics = aj["sigma_ic68_px"]
    print(f"U verdadera        = {aj['U_verdadera_m_s']*1e3:.4f} mm/s"
          f"   IC68 [{ic[0]*1e3:.4f}, {ic[1]*1e3:.4f}]")
    print(f"sigma              = {aj['sigma_px']:.4f} px"
          f"   IC68 [{ics[0]:.4f}, {ics[1]:.4f}]  ({len(aj['sigma_por_punto_px'])} puntos)")
    print(f"sigma punto a punto: {aj['sigma_por_punto_px']}")
    print(f"sigma de Foucaut   = {SIGMA_FOUCAUT_PX} px  (ventana 16x16, tabla 1)")
    print(f"desviación         = {aj['desviacion_relativa_vs_foucaut']:.1%}")
    print(f"R² del ajuste      = {aj['R2']:.5f}   sobre n = {aj['n_usados']}")
    s1 = aj["sin_n1"]
    print(f"sin el punto n=1   : U = {s1['U_verdadera_m_s']*1e3:.4f} mm/s, "
          f"sigma = {s1['sigma_px']:.4f} px, R² = {s1['R2']:.5f}")
    if aj["n_excluidos_por_validez"]:
        print(f"excluidos por validez < {FRAC_VALIDA_MINIMA:.0%}: n = {aj['n_excluidos_por_validez']}")
    print(f"\n-> {f1}\n-> {f2}\n-> {destino}")


if __name__ == "__main__":
    main()
