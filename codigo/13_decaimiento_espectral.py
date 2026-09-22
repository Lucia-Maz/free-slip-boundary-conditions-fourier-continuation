#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Prueba rápida de [P-02]: decaimiento por banda espectral.

Usa solamente ``salidas/tablas/decaimiento.json``; no necesita el disco externo. Los
espectros guardados son de desplazamiento y fueron calculados con separaciones distintas
entre cuadros. Para compararlos en el tiempo se resta el piso de PIV con el mismo criterio
de ``07_figuras_decaimiento.py`` y se divide la energía por n².

La tasa de cada bin se estima con un único coeficiente temporal y una ordenada distinta
para cada uno de los cuatro decaimientos:

    log E_ij(t) = c_ij - 2 lambda_j t.

El factor 2 convierte la pendiente de energía en tasa de amplitud. La prueba contesta dos
preguntas acotadas: si la escala migra hacia k bajos y si las tasas aparentes conservan la
dependencia viscosa k². No identifica por sí sola la fricción de fondo: cada bin también
recibe o entrega energía por transferencia espectral, y el PIV mide la superficie.

Uso:
    /home/lucia/miniforge3/envs/piv-dt/bin/python codigo/13_decaimiento_espectral.py
"""

import json
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt              # noqa: E402
import numpy as np                            # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(AQUI)
sys.path.insert(0, AQUI)
from celda import CAMPANA_02_06_25 as CELDA  # noqa: E402

ENTRADA = os.path.join(RAIZ, "salidas", "tablas", "decaimiento.json")
SALIDA = os.path.join(RAIZ, "salidas", "tablas", "decaimiento_espectral.json")
FIGURA = os.path.join(RAIZ, "salidas", "figuras", "f9_decaimiento_espectral")

T_CORTE = 10.0
T_CORTE_CONTROL = 17.0       # desde acá todos los espectros usan n=32
K_FIT = (80.0, 250.0)       # bins con señal y dispersión chica en los cuatro registros
K_MAX = 600.0
BANDAS = (
    (0.0, 118.0, "59–106 m⁻¹"),
    (118.0, 212.0, "130–200 m⁻¹"),
    (212.0, 354.0, "224–342 m⁻¹"),
    (354.0, 600.0, "366–578 m⁻¹"),
)


def espectro_velocidad_limpio(e):
    """Devuelve (k, E_vel) con piso restado y la dependencia n² eliminada.

    Los JSON actuales son anteriores a que ``06_decaimiento.py`` guardara el número
    exacto de modos por anillo. En ese caso se usa su proporcionalidad con k, igual que
    en ``07_figuras_decaimiento.espectro_limpio``. La constante se absorbe en el piso.
    Las unidades absolutas no hacen falta: sólo se ajustan pendientes y fracciones.
    """
    k = np.asarray(e["k_1_m"], float)
    energia = np.asarray(e["E"], float)
    cuenta = np.asarray(e.get("modos_por_anillo", k), float)
    k_nyquist = math.pi / float(e["paso_grilla_m"])
    banda_ruido = (k > 0.35 * k_nyquist) & (k < 0.95 * k_nyquist)
    if banda_ruido.sum() < 3:
        banda_ruido = slice(int(0.6 * len(k)), int(0.85 * len(k)))
    piso = float(np.median((energia / np.maximum(cuenta, 1.0))[banda_ruido]))
    limpia = np.maximum(energia - piso * cuenta, 0.0)
    return k, limpia / float(e["n"]) ** 2


def ajuste_comun(filas, registros):
    """Ajusta log(E) con intercepto por registro y pendiente temporal común."""
    X, y = [], []
    for registro, tiempo, energia in filas:
        if np.isfinite(energia) and energia > 0:
            X.append([float(registro == r) for r in registros] + [float(tiempo)])
            y.append(math.log(float(energia)))
    X, y = np.asarray(X), np.asarray(y)
    if len(y) <= len(registros) + 1:
        return {"lambda_1_s": float("nan"), "error_formal_1_s": float("nan"),
                "n": int(len(y)), "dispersion_log": float("nan")}
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    residuo = y - X @ coef
    grados = len(y) - X.shape[1]
    s2 = float(residuo @ residuo) / grados
    cov = s2 * np.linalg.inv(X.T @ X)
    return {
        "lambda_1_s": -float(coef[-1]) / 2.0,
        "error_formal_1_s": math.sqrt(float(cov[-1, -1])) / 2.0,
        "n": int(len(y)),
        "dispersion_log": math.sqrt(s2),
    }


def tasa_por_bin(datos, indice, t_corte):
    registros = sorted(datos["registros"])
    filas = []
    for registro in registros:
        for e in datos["registros"][registro]["espectros"]:
            if e["t_s"] < t_corte:
                continue
            _, energia = espectro_velocidad_limpio(e)
            filas.append((registro, e["t_s"], energia[indice]))
    return ajuste_comun(filas, registros)


def ajuste_k2(datos, t_corte):
    primero = next(iter(datos["registros"].values()))["espectros"][0]
    k = np.asarray(primero["k_1_m"], float)
    indices = np.where((k >= K_FIT[0]) & (k <= K_FIT[1]))[0]
    puntos = []
    for j in indices:
        tasa = tasa_por_bin(datos, int(j), t_corte)
        puntos.append({"k_1_m": float(k[j]), **tasa})

    lambdas = np.asarray([p["lambda_1_s"] for p in puntos])
    X = np.column_stack([np.ones(len(indices)), k[indices] ** 2])
    coef, *_ = np.linalg.lstsq(X, lambdas, rcond=None)
    residuo = lambdas - X @ coef
    denominador = float((lambdas - lambdas.mean()) @ (lambdas - lambdas.mean()))
    r2 = 1.0 - float(residuo @ residuo) / denominador

    # La dispersión entre bins subestima la incertidumbre porque todos salen de los
    # mismos campos. La incertidumbre que se informa sale de repetir el ajuste entero
    # por registro y tomar el SEM entre los cuatro resultados.
    por_registro = {}
    coefs = []
    for registro in sorted(datos["registros"]):
        tasas = []
        espectros = [e for e in datos["registros"][registro]["espectros"]
                     if e["t_s"] >= t_corte]
        for j in indices:
            tt, ee = [], []
            for e in espectros:
                _, energia = espectro_velocidad_limpio(e)
                if energia[j] > 0:
                    tt.append(e["t_s"])
                    ee.append(energia[j])
            pendiente = np.polyfit(tt, np.log(ee), 1)[0]
            tasas.append(-float(pendiente) / 2.0)
        c, *_ = np.linalg.lstsq(X, np.asarray(tasas), rcond=None)
        coefs.append(c)
        por_registro[registro] = {"alfa_espectral_1_s": float(c[0]),
                                  "nu_espectral_m2_s": float(c[1])}
    coefs = np.asarray(coefs)
    sem = coefs.std(axis=0, ddof=1) / math.sqrt(len(coefs))

    return {
        "t_corte_s": t_corte,
        "rango_k_ajuste_1_m": list(K_FIT),
        "alfa_espectral_1_s": float(coef[0]),
        "alfa_espectral_sem_entre_registros_1_s": float(sem[0]),
        "nu_espectral_m2_s": float(coef[1]),
        "nu_espectral_sem_entre_registros_m2_s": float(sem[1]),
        "R2_k2": r2,
        "rms_residuo_1_s": float(np.sqrt(np.mean(residuo ** 2))),
        "por_registro": por_registro,
        "puntos": puntos,
    }


def bandas_y_migracion(datos):
    registros = sorted(datos["registros"])
    filas_total = []
    tasas_bandas = []
    for k0, k1, nombre in BANDAS:
        filas = []
        for registro in registros:
            for e in datos["registros"][registro]["espectros"]:
                if e["t_s"] < T_CORTE:
                    continue
                k, energia = espectro_velocidad_limpio(e)
                valor = float(energia[(k >= k0) & (k < k1)].sum())
                filas.append((registro, e["t_s"], valor))
        tasas_bandas.append({"nombre": nombre, "k_1_m": [k0, k1],
                             **ajuste_comun(filas, registros)})

    tiempos = sorted({float(e["t_s"])
                      for r in datos["registros"].values()
                      for e in r["espectros"] if e["t_s"] >= T_CORTE})
    fracciones = []
    for tiempo in tiempos:
        por_registro = []
        krms = []
        for registro in registros:
            e = next(x for x in datos["registros"][registro]["espectros"]
                     if abs(x["t_s"] - tiempo) < 1e-9)
            k, energia = espectro_velocidad_limpio(e)
            usa = k < K_MAX
            total = float(energia[usa].sum())
            por_registro.append([
                float(energia[(k >= k0) & (k < k1)].sum()) / total
                for k0, k1, _ in BANDAS
            ])
            krms.append(math.sqrt(float((k[usa] ** 2 * energia[usa]).sum()) / total))
            filas_total.append((registro, tiempo, total))
        fr = np.asarray(por_registro)
        fracciones.append({
            "t_s": tiempo,
            "fraccion_media": fr.mean(axis=0).tolist(),
            "fraccion_sd_entre_registros": fr.std(axis=0, ddof=1).tolist(),
            "k_rms_medio_1_m": float(np.mean(krms)),
        })

    # Contrafactual: cada bin del primer espectro tardío decae independientemente con
    # alpha + nu k². Compara la amplitud final predicha con la observada.
    pronosticos = {}
    alfa, nu = CELDA.alfa(), CELDA.fluido.nu
    for registro in registros:
        espectros = [e for e in datos["registros"][registro]["espectros"]
                     if e["t_s"] >= T_CORTE]
        k, inicial = espectro_velocidad_limpio(espectros[0])
        _, final = espectro_velocidad_limpio(espectros[-1])
        usa = k < K_MAX
        dt = float(espectros[-1]["t_s"] - espectros[0]["t_s"])
        observado = math.sqrt(float(final[usa].sum() / inicial[usa].sum()))
        propagado = inicial[usa] * np.exp(-2.0 * (alfa + nu * k[usa] ** 2) * dt)
        lineal = math.sqrt(float(propagado.sum() / inicial[usa].sum()))
        pronosticos[registro] = {
            "dt_s": dt,
            "cociente_amplitud_observado": observado,
            "cociente_amplitud_lineal_sin_transferencia": lineal,
            "observado_sobre_lineal": observado / lineal,
        }

    return tasas_bandas, fracciones, ajuste_comun(filas_total, registros), pronosticos


def figura(resultado):
    ajuste = resultado["ajuste_modal_k2"]
    puntos = ajuste["puntos"]
    k = np.asarray([p["k_1_m"] for p in puntos])
    lam = np.asarray([p["lambda_1_s"] for p in puntos])
    err = np.asarray([p["error_formal_1_s"] for p in puntos])

    plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.dpi": 140})
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11.2, 4.4))

    ax0.errorbar(k, lam, yerr=err, fmt="o", color="#2a78d6", ms=6,
                 capsize=2.5, label="bins del PIV")
    xx = np.linspace(0, 270, 300)
    ax0.plot(xx, CELDA.alfa() + CELDA.fluido.nu * xx ** 2, "--", color="#333333",
             label=r"nominal: $\alpha+\nu k^2$")
    ax0.plot(xx, ajuste["alfa_espectral_1_s"] +
             ajuste["nu_espectral_m2_s"] * xx ** 2, color="#eb6834",
             label="ajuste espectral")
    ax0.set(xlabel=r"número de onda $k$  [m$^{-1}$]",
            ylabel=r"tasa de amplitud $\lambda$  [s$^{-1}$]",
            title="La pendiente en k² recupera ν; la ordenada queda baja",
            xlim=(0, 270), ylim=(0.035, 0.16))
    ax0.legend(frameon=True, framealpha=0.92, facecolor="white", fontsize=9)
    ax0.text(0.04, 0.94,
             (r"$\alpha_{spec}=%.4f\pm%.4f$ s$^{-1}$" "\n"
              r"$\nu_{spec}=(%.2f\pm%.2f)\,10^{-6}$ m$^2$/s")
             % (ajuste["alfa_espectral_1_s"],
                ajuste["alfa_espectral_sem_entre_registros_1_s"],
                1e6 * ajuste["nu_espectral_m2_s"],
                1e6 * ajuste["nu_espectral_sem_entre_registros_m2_s"]),
             transform=ax0.transAxes, va="top", fontsize=9)

    colores = ("#2a78d6", "#1baf7a", "#eda100", "#eb6834")
    fracciones = resultado["migracion_espectral"]
    tt = np.asarray([x["t_s"] for x in fracciones])
    for j, ((_, _, nombre), color) in enumerate(zip(BANDAS, colores)):
        media = np.asarray([x["fraccion_media"][j] for x in fracciones])
        sd = np.asarray([x["fraccion_sd_entre_registros"][j] for x in fracciones])
        ax1.plot(tt, media, "o-", color=color, label=nombre)
        ax1.fill_between(tt, np.maximum(media - sd, 0), np.minimum(media + sd, 1),
                         color=color, alpha=0.12, linewidth=0)
    ax1.set(xlabel="tiempo [s]", ylabel="fracción de energía espectral",
            title="La energía migra hacia las escalas grandes", ylim=(0, 0.75))
    ax1.legend(frameon=False, fontsize=8.5, ncol=2)
    ax1.text(0.03, 0.04, "media ± SD entre cuatro decaimientos\n"
             "espectros de superficie; piso de PIV restado",
             transform=ax1.transAxes, fontsize=8.5, color="#555555")

    fig.tight_layout()
    os.makedirs(os.path.dirname(FIGURA), exist_ok=True)
    fig.savefig(FIGURA + ".png", bbox_inches="tight", pad_inches=0.2)
    fig.savefig(FIGURA + ".pdf", bbox_inches="tight", pad_inches=0.2)
    plt.close(fig)


def main():
    with open(ENTRADA) as fh:
        datos = json.load(fh)

    ajuste = ajuste_k2(datos, T_CORTE)
    control = ajuste_k2(datos, T_CORTE_CONTROL)
    tasas_bandas, migracion, total, pronosticos = bandas_y_migracion(datos)
    alfa = ajuste["alfa_espectral_1_s"]
    h_equivalente = math.pi / 2.0 * math.sqrt(CELDA.fluido.nu / alfa)

    resultado = {
        "fecha": "2026-09-21",
        "fuente": os.path.relpath(ENTRADA, RAIZ),
        "pregunta": "si la tasa baja se explica por una transición a un flujo poco turbulento",
        "metodo": {
            "t_corte_s": T_CORTE,
            "normalizacion": "E de desplazamiento dividido por n^2; lambda=-pendiente(log E)/2",
            "piso_PIV": "mediana por modo entre 0.35 y 0.95 de Nyquist, restada por anillo",
            "ajuste_temporal": "pendiente común y una ordenada por cada decaimiento",
            "k_max_migracion_1_m": K_MAX,
            "bandas": [{"nombre": nombre, "k_1_m": [k0, k1]}
                       for k0, k1, nombre in BANDAS],
        },
        "parametros_nominales": {
            "alfa_1_s": CELDA.alfa(),
            "nu_m2_s": CELDA.fluido.nu,
            "h_m": CELDA.h,
        },
        "tasa_total_espectral": total,
        "tasas_por_banda": tasas_bandas,
        "ajuste_modal_k2": ajuste,
        "control_solo_n32_t_mayor_17_s": control,
        "h_equivalente_si_todo_el_intercepto_fuera_friccion_m": h_equivalente,
        "migracion_espectral": migracion,
        "pronostico_lineal_sin_transferencia": pronosticos,
        "conclusion": [
            "La energia migra hacia k bajos en los cuatro decaimientos.",
            "Las tasas aparentes entre 83 y 248 1/m recuperan el coeficiente nu*k^2.",
            "El deficit restante es aproximadamente independiente de k y no es el limite lineal esperado para h=6 mm.",
            "Los bins no son modos aislados: transferencia espectral, perfil vertical de superficie y ruido PIV correlacionado impiden identificar el intercepto directamente con alfa.",
        ],
    }

    os.makedirs(os.path.dirname(SALIDA), exist_ok=True)
    with open(SALIDA, "w") as fh:
        json.dump(resultado, fh, indent=1, ensure_ascii=False)
    figura(resultado)

    print("tasa total espectral = %.5f ± %.5f 1/s"
          % (total["lambda_1_s"], total["error_formal_1_s"]))
    print("lambda(k) = (%.5f ± %.5f) + (%.3f ± %.3f)e-6 k²"
          % (ajuste["alfa_espectral_1_s"],
             ajuste["alfa_espectral_sem_entre_registros_1_s"],
             1e6 * ajuste["nu_espectral_m2_s"],
             1e6 * ajuste["nu_espectral_sem_entre_registros_m2_s"]))
    print("control t>=17 s: alfa=%.5f, nu=%.3fe-6, R²=%.3f"
          % (control["alfa_espectral_1_s"], 1e6 * control["nu_espectral_m2_s"],
             control["R2_k2"]))
    print("h equivalente (interpretación no única) = %.2f mm" % (1e3 * h_equivalente))
    print("escrito", os.path.relpath(SALIDA, RAIZ))
    print("figura", os.path.relpath(FIGURA + ".png", RAIZ))
    return 0


if __name__ == "__main__":
    sys.exit(main())
