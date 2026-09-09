#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
¿Se está filtrando de más? Verificaciones del tratamiento de ruido.

Pregunta de Lucía: hay manera de saber si estoy suavizando de más y perdiendo
información, y ¿se debería poder identificar el filtro en los espectros?

Se separan los pasos que SON filtros de los que no, y se verifica cada uno:

  A. Elegir n grande (Δt largo) para medir la velocidad. NO es un suavizado espacial,
     pero sí promedia en el tiempo sobre n/fps segundos. El test es directo: si
     d1(n) = desp_rms(n)/n no depende de n en el rango usado, no se está perdiendo
     señal. Si cayera monótonamente, sí.

  B. La ventana de interrogación del PIV. Este SÍ es un filtro pasabajos espacial, no
     lo elegí yo, y es el que Lucía pregunta si se ve en el espectro. Se compara el
     espectro medido contra la respuesta sinc^2(kX/2) del marco de Foucaut para dos
     candidatos de X: la ventana final (16 px) y el paso de la grilla (4 px).

  C. La resta del piso de ruido del espectro. Es el único paso donde yo elijo algo que
     puede borrar señal. Se verifica de dos maneras independientes:
       C1. sobre datos SINTÉTICOS de sigma conocido, que es la única forma de saber si
           el estimador está bien normalizado;
       C2. contra el sigma que sale del ajuste de dos Δt, que es una medición
           completamente independiente sobre los mismos cuadros.
     Y se reporta cuánto se mueve el k pesado si se cambia la banda donde se lee el
     piso: si se moviera mucho, el número no sería del dato sino de la elección.

Lo que NO es un filtro, y conviene decirlo: la mediana móvil de la figura 1 es sólo
cosmética, se aplica a las líneas tenues y no entra en ningún número.

Uso:
    /home/lucia/miniforge3/envs/piv-dt/bin/python codigo/09_verificacion_filtrado.py
"""

import importlib.util
import json
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(AQUI)
sys.path.insert(0, AQUI)

_spec = importlib.util.spec_from_file_location(
    "dec", os.path.join(AQUI, "06_decaimiento.py"))
dec = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dec)

_spec2 = importlib.util.spec_from_file_location(
    "figs", os.path.join(AQUI, "07_figuras_decaimiento.py"))
figs = importlib.util.module_from_spec(_spec2)
_spec2.loader.exec_module(figs)

TABLAS = os.path.join(RAIZ, "salidas", "tablas")
SALIDA = os.path.join(TABLAS, "verificacion_filtrado.json")

VENTANA_FINAL_PX = 16.0
PASO_GRILLA_PX = 4.0


# --------------------------------------------------------------------------------
# C1. El estimador del piso, sobre ruido sintético de sigma conocido
# --------------------------------------------------------------------------------

def verificar_sobre_sintetico(nx=253, sigmas=(0.05, 0.1, 0.2, 0.4), semilla=0):
    """Ruido blanco puro de sigma conocido -> ¿lo recupera el estimador?"""
    rng = np.random.default_rng(semilla)
    filas = []
    for s in sigmas:
        campo = {"u": rng.normal(0.0, s, (nx, nx)),
                 "v": rng.normal(0.0, s, (nx, nx)),
                 "s2n": np.full((nx, nx), 10.0)}
        e = dec.espectro(campo, PASO_GRILLA_PX / 3800.0, 1.05)
        rec, _ = figs.sigma_desde_espectro(e)
        filas.append({"sigma_puesto_px": s, "sigma_recuperado_px": rec,
                      "error_relativo": rec / s - 1.0})
        del campo, e
    return filas


# --------------------------------------------------------------------------------

def main():
    d = json.load(open(os.path.join(TABLAS, "decaimiento.json")))
    res = {}

    # ---- C1 ------------------------------------------------------------------
    res["sintetico"] = verificar_sobre_sintetico()

    # ---- A: d1(n) contra n ---------------------------------------------------
    a_filas = []
    for med in sorted(d["registros"]):
        for p in d["registros"][med]["calibracion"]:
            pts = {q["n"]: q for q in p["puntos"]
                   if q["frac_valida"] >= 0.90 and q["desp_rms_px"] < 4.0}
            grandes = sorted(n for n in pts if n >= 8)
            if len(grandes) < 2:
                continue
            d1 = {n: pts[n]["desp_rms_px"] / n for n in grandes}
            ref = d1[max(grandes)]
            a_filas.append({
                "med": med, "t_s": p["t_s"], "n_usados": grandes,
                "d1_por_n": {str(n): d1[n] for n in grandes},
                "dispersion_relativa": (max(d1.values()) - min(d1.values())) / ref,
                "convergencia_ultimos_dos": (
                    abs(d1[grandes[-1]] - d1[grandes[-2]]) / ref),
            })
    res["d1_vs_n"] = a_filas
    disp = [f["dispersion_relativa"] for f in a_filas]
    conv = [f["convergencia_ultimos_dos"] for f in a_filas]
    res["d1_vs_n_resumen"] = {
        "n_instantes": len(a_filas),
        "dispersion_mediana": float(np.median(disp)) if disp else float("nan"),
        "dispersion_maxima": float(np.max(disp)) if disp else float("nan"),
        "convergencia_mediana": float(np.median(conv)) if conv else float("nan"),
        "convergencia_maxima": float(np.max(conv)) if conv else float("nan"),
    }

    # ---- C2: sigma espectral contra sigma del ajuste -------------------------
    c_filas = []
    for med in sorted(d["registros"]):
        reg = d["registros"][med]
        aj = {round(p["t_s"], 3): p for p in reg["calibracion"]}
        for e in reg["espectros"]:
            s_esp, p0 = figs.sigma_desde_espectro(e)
            p = aj.get(round(e["t_s"], 3))
            s_fit = p["ajuste_dos_dt"]["sigma_por_componente_px"] if p else float("nan")
            c_filas.append({"med": med, "t_s": e["t_s"], "sigma_espectral_px": s_esp,
                            "sigma_ajuste_px": s_fit, "p0": p0})
    res["sigma_cruzado"] = c_filas
    pos = [(f["sigma_espectral_px"], f["sigma_ajuste_px"]) for f in c_filas
           if f["sigma_ajuste_px"] > 0]
    res["sigma_cruzado_resumen"] = {
        "n_con_ajuste_positivo": len(pos),
        "cociente_mediano": float(np.median([a / b for a, b in pos])) if pos else None,
        "sigma_espectral_mediano_px": float(np.median(
            [f["sigma_espectral_px"] for f in c_filas])),
    }

    # ---- C3: sensibilidad del k pesado a la banda del piso -------------------
    def k_pesado_con_banda(e, lo, hi):
        k = np.array(e["k_1_m"], float); E = np.array(e["E"], float)
        _, _, dk = figs.geometria(e)
        cuenta = (np.array(e["modos_por_anillo"], float)
                  if "modos_por_anillo" in e else 2 * math.pi * k / dk)
        kn = math.pi / e["paso_grilla_m"]
        b = (k > lo * kn) & (k < hi * kn)
        p0 = float(np.median((E / np.maximum(cuenta, 1.0))[b]))
        lim = np.maximum(E - p0 * cuenta, 0.0)
        return float((k * lim).sum() / lim.sum()) if lim.sum() > 0 else float("nan")

    bandas = [(0.25, 0.95), (0.35, 0.95), (0.45, 0.95), (0.35, 0.80), (0.55, 0.95)]
    e0 = d["registros"]["med_S0003"]["espectros"][0]
    e1 = d["registros"]["med_S0003"]["espectros"][-1]
    res["sensibilidad_banda"] = [
        {"banda": list(b),
         "k_pesado_t0": k_pesado_con_banda(e0, *b),
         "k_pesado_tfinal": k_pesado_con_banda(e1, *b)} for b in bandas]

    # ---- B: la respuesta de la ventana, contra el espectro medido ------------
    k = np.array(e0["k_1_m"], float)
    _, _, dk = figs.geometria(e0)
    cuenta = 2 * math.pi * k / dk
    por_modo = np.array(e0["E"], float) / np.maximum(cuenta, 1.0)
    kn = math.pi / e0["paso_grilla_m"]
    banda = (k > 0.3 * kn) & (k < 0.98 * kn)

    def sinc2(kk, X_px):
        X = X_px / 3800.0
        a = kk * X / 2.0
        return np.where(a == 0, 1.0, (np.sin(a) / np.where(a == 0, 1, a)) ** 2)

    ref = float(np.median(por_modo[banda]))
    med_norm = por_modo[banda] / ref
    res["respuesta_ventana"] = {
        "k_corte_foucaut_ventana16_1_m": 2.8 / (VENTANA_FINAL_PX / 3800.0),
        "primer_cero_sinc_ventana16_1_m": 2 * math.pi / (VENTANA_FINAL_PX / 3800.0),
        "primer_cero_sinc_paso4_1_m": 2 * math.pi / (PASO_GRILLA_PX / 3800.0),
        "k_nyquist_grilla_1_m": kn,
        "caida_medida_en_la_banda": float(med_norm[-1] / med_norm[0]),
        "caida_predicha_sinc2_ventana16": float(
            sinc2(k[banda][-1], VENTANA_FINAL_PX) / sinc2(k[banda][0], VENTANA_FINAL_PX)),
        "caida_predicha_sinc2_paso4": float(
            sinc2(k[banda][-1], PASO_GRILLA_PX) / sinc2(k[banda][0], PASO_GRILLA_PX)),
    }

    with open(SALIDA, "w") as f:
        json.dump(res, f, indent=1, ensure_ascii=False)

    # ---- colapso de A: exceso medido contra exceso predicho por el ruido ----
    colapso = []
    for med in sorted(d["registros"]):
        reg = d["registros"][med]
        esp = {round(e["t_s"], 3): e for e in reg["espectros"]}
        for p in reg["calibracion"]:
            e = esp.get(round(p["t_s"], 3))
            if e is None:
                continue
            sig = figs.sigma_desde_espectro(e)[0] * math.sqrt(2.0)   # las dos componentes
            pts = {q["n"]: q for q in p["puntos"] if q["frac_valida"] >= 0.90}
            grandes = [n for n in pts if n >= 8 and pts[n]["desp_rms_px"] < 4.0]
            if not grandes:
                continue
            nref = max(grandes)
            d1v = pts[nref]["desp_rms_px"] / nref
            for n, q in sorted(pts.items()):
                if n == nref:
                    continue
                pred = (sig / (n * d1v)) ** 2
                medd = (q["desp_rms_px"] / (n * d1v)) ** 2 - 1.0
                if medd > 0:
                    colapso.append({"med": med, "t_s": p["t_s"], "n": n,
                                    "predicho": pred, "medido": medd})
    res["colapso_ruido"] = colapso
    if colapso:
        coc = [f["medido"] / f["predicho"] for f in colapso]
        res["colapso_resumen"] = {
            "n_puntos": len(colapso),
            "cociente_mediano": float(np.median(coc)),
            "cociente_q25": float(np.quantile(coc, .25)),
            "cociente_q75": float(np.quantile(coc, .75)),
        }

    # ---- figura --------------------------------------------------------------
    for tema in ("claro", "oscuro"):
        t = figs.estilo(tema)
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.3))

        # izquierda: el exceso de d1 a n chico, ¿es sólo ruido?
        # Si desp(n)^2 = (n d1)^2 + sigma^2 y nada más, entonces
        #    (d1(n)/d1_verdadero)^2 - 1  =  (sigma / (n d1_verdadero))^2
        # exactamente, sin parámetros libres: pendiente 1 en log-log. Cualquier
        # pérdida de señal por promediar en el tiempo se vería como puntos POR
        # DEBAJO de la recta a n grande.
        for i, med in enumerate(sorted(d["registros"])):
            c = t["series"][i]
            for f in [x for x in colapso if x["med"] == med]:
                ax1.plot(f["predicho"], f["medido"], "o", color=c, ms=5.5,
                         mec=t["fondo"], mew=1.0, alpha=0.8)
        lim = [3e-4, 30]
        ax1.plot(lim, lim, color=t["tinta2"], lw=1.4, ls="--")
        ax1.set_xscale("log"); ax1.set_yscale("log")
        ax1.set_xlim(*lim); ax1.set_ylim(*lim)
        ax1.set_xlabel("predicho sólo por ruido:  (σ / n·d₁)²")
        ax1.set_ylabel("medido:  (d₁(n)/d₁)² − 1")
        ax1.set_title("A · ¿el Δt corto sube por ruido, o el largo pierde señal?")
        ax1.text(0.035, 0.955,
                 "la recta es y = x, sin ajustar nada\n"
                 "por ENCIMA ⇒ el n grande perdió señal\n"
                 "por debajo ⇒ σ del espectro es cota superior",
                 transform=ax1.transAxes, color=t["tinta2"], fontsize=8.5,
                 va="top")

        # derecha: piso del espectro contra las respuestas candidatas
        ax2.plot(k[banda], med_norm, color=t["series"][0], lw=2.0)
        ax2.plot(k[banda], sinc2(k[banda], VENTANA_FINAL_PX)
                 / sinc2(k[banda][0], VENTANA_FINAL_PX), color=t["series"][1],
                 lw=1.8, ls="--")
        ax2.plot(k[banda], sinc2(k[banda], PASO_GRILLA_PX)
                 / sinc2(k[banda][0], PASO_GRILLA_PX), color=t["series"][2],
                 lw=1.8, ls=":")
        ax2.set_yscale("log")
        ax2.set_ylim(1e-2, 3)
        ax2.set_xlabel("número de onda k  [1/m]")
        ax2.set_ylabel("potencia por modo, normalizada")
        ax2.set_title("B · ¿se ve la ventana en el espectro?")
        ax2.annotate("medido", (k[banda][-1], med_norm[-1]), xytext=(-8, 8),
                     textcoords="offset points", color=t["series"][0], ha="right",
                     fontsize=9.5, fontweight="semibold")
        curva16 = (sinc2(k[banda], VENTANA_FINAL_PX)
                   / sinc2(k[banda][0], VENTANA_FINAL_PX))
        j16 = int(np.argmax(np.where(k[banda] > 1600, curva16, 0)))
        ax2.annotate("sinc², ventana 16 px\n(primer cero en 1492 m⁻¹)",
                     (k[banda][j16], curva16[j16]), xytext=(0, -8),
                     textcoords="offset points", color=t["series"][1],
                     ha="center", va="top", fontsize=9.5)
        ax2.annotate("sinc², paso 4 px", (k[banda][-1],
                     (sinc2(k[banda], PASO_GRILLA_PX)
                      / sinc2(k[banda][0], PASO_GRILLA_PX))[-1]),
                     xytext=(-8, 6), textcoords="offset points",
                     color=t["series"][2], ha="right", fontsize=9.5)
        fig.tight_layout()
        figs.guardar(fig, "f6_filtrado", tema)

    # ---- informe -------------------------------------------------------------
    print("C1 · el estimador del piso, sobre ruido blanco de sigma CONOCIDO")
    for f in res["sintetico"]:
        print("   sigma puesto %.3f px  ->  recuperado %.4f px   (error %+.2f %%)"
              % (f["sigma_puesto_px"], f["sigma_recuperado_px"],
                 100 * f["error_relativo"]))
    print()
    print("A · d1(n) para n >= 8, en los %d instantes con al menos dos n válidos"
          % res["d1_vs_n_resumen"]["n_instantes"])
    print("   dispersión sobre todos los n>=8: mediana %.2f %%, máxima %.2f %%"
          % (100 * res["d1_vs_n_resumen"]["dispersion_mediana"],
             100 * res["d1_vs_n_resumen"]["dispersion_maxima"]))
    print("   entre los DOS n más grandes:     mediana %.2f %%, máxima %.2f %%"
          % (100 * res["d1_vs_n_resumen"]["convergencia_mediana"],
             100 * res["d1_vs_n_resumen"]["convergencia_maxima"]))
    print()
    print("C2 · sigma del espectro contra sigma del ajuste de dos Δt")
    print("   sigma espectral mediano: %.4f px"
          % res["sigma_cruzado_resumen"]["sigma_espectral_mediano_px"])
    print("   cociente mediano espectro/ajuste: %s  (sobre %d instantes con sigma^2>0)"
          % (("%.3f" % res["sigma_cruzado_resumen"]["cociente_mediano"])
             if res["sigma_cruzado_resumen"]["cociente_mediano"] else "—",
             res["sigma_cruzado_resumen"]["n_con_ajuste_positivo"]))
    print()
    if "colapso_resumen" in res:
        cr = res["colapso_resumen"]
        print("A-bis · el exceso a n chico contra lo que predice sólo el ruido")
        print("   %d puntos; cociente medido/predicho: mediana %.2f  (q25 %.2f, q75 %.2f)"
              % (cr["n_puntos"], cr["cociente_mediano"], cr["cociente_q25"],
                 cr["cociente_q75"]))
        print()

    print("C3 · sensibilidad del k pesado a la banda donde se lee el piso")
    for f in res["sensibilidad_banda"]:
        print("   banda %.2f–%.2f k_Nyq :  k(t=0) = %6.1f   k(final) = %6.1f"
              % (f["banda"][0], f["banda"][1], f["k_pesado_t0"], f["k_pesado_tfinal"]))
    print()
    r = res["respuesta_ventana"]
    print("B · la ventana en el espectro")
    print("   k_Nyquist de la grilla de vectores : %.0f 1/m" % r["k_nyquist_grilla_1_m"])
    print("   primer cero de sinc², ventana 16 px: %.0f 1/m" % r["primer_cero_sinc_ventana16_1_m"])
    print("   corte de Foucaut k_c X = 2,8       : %.0f 1/m" % r["k_corte_foucaut_ventana16_1_m"])
    print("   caída del piso medido en la banda  : %.3f" % r["caida_medida_en_la_banda"])
    print("   la que predice sinc² con X = 16 px : %.3f" % r["caida_predicha_sinc2_ventana16"])
    print("   la que predice sinc² con X =  4 px : %.3f" % r["caida_predicha_sinc2_paso4"])
    print("\nescrito %s" % os.path.relpath(SALIDA, RAIZ))
    return 0


if __name__ == "__main__":
    sys.exit(main())
