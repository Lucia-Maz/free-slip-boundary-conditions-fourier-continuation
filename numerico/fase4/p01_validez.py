#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
[P-01] Validez de la reducción bidimensional en el régimen medido.

Acompaña a teoria/P01_validez_reduccion_2D.md. La derivación está allá; acá está todo
lo que se puede *calcular*, para que ninguna afirmación cuantitativa del texto sea una
afirmación a secas.

Cuatro cosas, en este orden:

  1. El número de onda del forzado, sacado de la transformada del patrón real de imanes
     (discos de diámetro finito, polaridad alternada en red cuadrada), no de una escala
     elegida a ojo. Esto reemplaza el l = 2 cm heredado, que no correspondía a ninguna
     longitud de la celda.
  2. El espectro vertical y los acoplamientos entre modos, por cuadratura, contrastados
     contra las formas cerradas donde existen.
  3. Los parámetros adimensionales de la reducción evaluados sobre la celda, y la
     distorsión del perfil vertical que predicen.
  4. Cómo evolucionan durante un decaimiento, que es lo que decide si hay o no una
     ventana temporal en la que la clausura vale.

Nada de esto pasa de unos pocos MB de memoria.

Uso:
    /home/lucia/miniforge3/envs/piv-dt/bin/python numerico/fase4/p01_validez.py
"""

import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "codigo"))

from celda import CAMPANA_02_06_25, VELOCIDADES_MEDIDAS  # noqa: E402

RAIZ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SALIDA = os.path.join(RAIZ, "salidas", "tablas", "p01_validez.json")


# --------------------------------------------------------------------------------
# 1. El número de onda del forzado, medido sobre el patrón
# --------------------------------------------------------------------------------

def espectro_forzado(forzado, celdas=8, puntos_por_celda=64):
    """Transformada del patrón de polaridad de los imanes.

    El modelo del campo es crudo a propósito: B_z vale +1 o -1 dentro de cada disco y 0
    fuera. Con una corriente horizontal uniforme la fuerza de Lorentz J x B es
    proporcional a ese patrón rotado 90 grados, así que el contenido espectral del
    forzado es el del patrón. La red fija las posiciones de los picos y la forma del
    imán sólo pesa sus amplitudes.

    Lo que esta función NO decide es si el fundamental de la red es donde el forzado
    pone la energía: el argmax modo a modo cae en el fundamental para cualquier red,
    también para una de fuentes puntuales donde el fundamental lleva pocos por ciento.
    Eso lo decide `fraccion_energia_en_fundamental` junto con el espectro medido de
    los campos forzados ([H-16]): con imanes que ocupan dos tercios del paso las dos
    cosas coinciden.

    Devuelve el pico dominante, el número de onda medio pesado por energía, y la
    orientación del pico, que es lo que decide si la contribución es diagonal.
    """
    a = forzado.paso
    n = celdas * puntos_por_celda
    L = celdas * a
    x = (np.arange(n) + 0.5) * (L / n)
    X, Y = np.meshgrid(x, x, indexing="ij")

    campo = np.zeros((n, n))
    radio = forzado.diametro / 2.0
    for i in range(celdas):
        for j in range(celdas):
            cx, cy = (i + 0.5) * a, (j + 0.5) * a
            dentro = (X - cx) ** 2 + (Y - cy) ** 2 <= radio ** 2
            campo[dentro] = 1.0 if (i + j) % 2 == 0 else -1.0

    esp = np.abs(np.fft.fft2(campo)) ** 2
    kx = 2.0 * np.pi * np.fft.fftfreq(n, d=L / n)
    KX, KY = np.meshgrid(kx, kx, indexing="ij")
    K = np.hypot(KX, KY)

    esp[0, 0] = 0.0                     # el modo medio no es forzado con estructura
    i, j = np.unravel_index(np.argmax(esp), esp.shape)
    k_pico = float(K[i, j])
    peso = esp / esp.sum()
    k_medio = float((K * peso).sum())
    # fracción de la energía del patrón que está en el armónico fundamental (los
    # cuatro picos equivalentes por simetría de la red)
    frac_fundamental = float(esp[np.isclose(K, k_pico, rtol=1e-9)].sum() / esp.sum())

    del campo, esp, peso, K, KX, KY, X, Y
    return {
        "k_pico_1_m": k_pico,
        "k_pico_cerrado_1_m": forzado.k_fundamental(),
        "error_relativo_vs_cerrado": abs(k_pico - forzado.k_fundamental())
                                     / forzado.k_fundamental(),
        "longitud_de_onda_m": 2.0 * math.pi / k_pico,
        "orientacion_pico": [int(round(KX_ij)) for KX_ij in
                             (kx[i] * a / math.pi, kx[j] * a / math.pi)],
        "k_medio_pesado_1_m": k_medio,
        "fraccion_energia_en_fundamental": frac_fundamental,
        "celdas": celdas,
        "puntos_por_celda": puntos_por_celda,
    }


# --------------------------------------------------------------------------------
# 2. Espectro vertical y acoplamientos entre modos
# --------------------------------------------------------------------------------

def modos_verticales(m_max=6, n=200001):
    """Autovalores y acoplamientos del problema vertical con fondo rígido y tope libre.

    Base: F_m(z) = sin((2m+1) pi z / 2h), que cumple F(0)=0 y F'(h)=0, con autovalor
    lambda_m = nu ((2m+1) pi / 2h)^2. Todo se calcula por cuadratura sobre z/h en [0,1]
    y se contrasta contra las formas cerradas.

    Los acoplamientos son c_m = <F_0^2 F_m> / <F_m^2>: cuánto del término no lineal
    generado por el modo fundamental cae sobre el modo m. Son los que fijan el tamaño
    de la distorsión del perfil.
    """
    z = np.linspace(0.0, 1.0, n)
    F = [np.sin((2 * m + 1) * np.pi * z / 2.0) for m in range(m_max + 1)]

    def prom(v):
        return float(np.trapezoid(v, z))

    F0 = F[0]
    prom_F0 = prom(F0)                       # 2/pi
    beta = prom(F0 ** 2) / prom_F0 ** 2      # <F^2> con la normalización <F> = 1
    modos = []
    for m in range(m_max + 1):
        Fm = F[m]
        modos.append({
            "m": m,
            "lambda_sobre_lambda0": float((2 * m + 1) ** 2),
            "acoplamiento_c_m": prom(F0 ** 2 * Fm) / prom(Fm ** 2),
            "promedio_vertical": prom(Fm),
            "derivada_en_el_fondo": float((2 * m + 1) * np.pi / 2.0),
        })
    del F, F0
    return {
        "modos": modos,
        "beta_medido": beta,
        "beta_cerrado": math.pi ** 2 / 8.0,
        "sesgo_superficie_medido": 1.0 / prom_F0,
        "sesgo_superficie_cerrado": math.pi / 2.0,
        "c1_cerrado": -4.0 / (15.0 * math.pi) / 0.5,
        "puntos_cuadratura": n,
    }


# --------------------------------------------------------------------------------
# 3. Los adimensionales sobre la celda
# --------------------------------------------------------------------------------

def adimensionales(celda, U, k):
    """U es la velocidad PROMEDIADA EN LA VERTICAL, no la de superficie."""
    h, nu = celda.h, celda.fluido.nu
    eps = k * h
    Re_h = U * h / nu
    return {
        "U_promedio_vertical_m_s": U,
        "epsilon_kh": eps,
        "Re_h": Re_h,
        "Re_l": U / (nu * k),
        "delta": Re_h * eps,                       # = U k h^2 / nu
        "lambda_lin_sobre_alfa": 1.0 + 4.0 * eps ** 2 / math.pi ** 2,
    }


def distorsion(delta, modos, m_max=6):
    """Amplitud relativa de los modos verticales m>=1 excitados por el término no
    lineal del modo fundamental, en balance cuasi-estacionario.

    lambda_m a_m ~ c_m k U a_0  =>  a_m/a_0 ~ c_m (k U / lambda_0) (lambda_0/lambda_m)
    y k U / lambda_0 = (4/pi^2) delta.

    ESTO ES UNA ESTIMACIÓN DE ORDEN, no un cálculo: el factor geométrico O(1) del
    término no lineal, que depende de la estructura horizontal del campo y de la
    proyección de la presión, se toma igual a 1. Ver las salvedades del documento.
    """
    lam = 4.0 / math.pi ** 2 * delta
    out = {}
    total = 0.0
    for mo in modos["modos"]:
        if mo["m"] == 0:
            continue
        if mo["m"] > m_max:
            break
        a = abs(mo["acoplamiento_c_m"]) * lam / mo["lambda_sobre_lambda0"]
        out["a%d_sobre_a0" % mo["m"]] = a
        total += a ** 2
    out["distorsion_rms"] = math.sqrt(total)
    out["k_U_sobre_lambda0"] = lam
    return out


# --------------------------------------------------------------------------------
# 4. Cómo evoluciona durante un decaimiento
# --------------------------------------------------------------------------------

def ventana_temporal(celda, delta_0, umbral):
    """Instante a partir del cual delta(t) = delta_0 exp(-alfa t) cae bajo el umbral.

    Es exacto en la medida en que el decaimiento sea el de la clausura, que es
    justamente lo que se está poniendo a prueba: sirve como estimación autoconsistente,
    no como demostración.
    """
    if delta_0 <= umbral:
        return 0.0
    return math.log(delta_0 / umbral) / celda.alfa()


def pico_medido_forzado(t_max=10.0):
    """El pico del espectro de los campos medidos en la meseta forzada, si existe.

    Lee `salidas/tablas/decaimiento.json` (lo produce `codigo/06_decaimiento.py`, que
    necesita el disco). La meseta inicial de cada registro es el estado forzado
    ([D-30]), así que el pico ahí es la escala de inyección medida, independiente de
    cualquier modelo del imán. Devuelve None si el archivo no está.
    """
    ruta = os.path.join(RAIZ, "salidas", "tablas", "decaimiento.json")
    if not os.path.exists(ruta):
        return None
    with open(ruta) as f:
        d = json.load(f)
    picos, anchos = [], []
    for reg in d["registros"].values():
        for e in reg["espectros"]:
            if e["t_s"] < t_max:
                picos.append(e["k_pico_1_m"])
                anchos.append(e["k_1_m"][1] - e["k_1_m"][0])
    if not picos:
        return None
    return {"k_1_m": float(np.median(picos)), "n": len(picos),
            "k_min_1_m": float(min(picos)), "k_max_1_m": float(max(picos)),
            "ancho_bin_1_m": float(np.median(anchos)),
            "registros": sorted(d["registros"]), "t_max_s": t_max}


def main():
    c = CAMPANA_02_06_25
    esp = espectro_forzado(c.forzado)
    k_med = pico_medido_forzado()
    if k_med is not None:
        esp["pico_medido_meseta"] = k_med
    modos = modos_verticales()
    k = esp["k_pico_1_m"]

    casos = {}
    for nombre, med in VELOCIDADES_MEDIDAS.items():
        u_sup = med["u_rms_superficie_m_s"]
        # el PIV mide la superficie; la clausura vive sobre el promedio vertical
        U = u_sup / c.sesgo_piv_superficie()
        ad = adimensionales(c, U, k)
        ad["u_rms_superficie_m_s"] = u_sup
        ad["verificado"] = med["verificado"]
        ad["procedencia"] = med["procedencia"]
        ad.update(distorsion(ad["delta"], modos))
        ad["t_hasta_distorsion_1pc_s"] = ventana_temporal(
            c, ad["distorsion_rms"], 0.01) if ad["distorsion_rms"] > 0 else 0.0
        casos[nombre] = ad

    # barrido en U, para leer el resultado sin depender de una sola cifra medida
    barrido = []
    for u_sup_mm in (0.25, 0.5, 1.0, 1.32, 2.0, 3.0, 4.0, 6.0):
        U = (u_sup_mm * 1e-3) / c.sesgo_piv_superficie()
        ad = adimensionales(c, U, k)
        ad["u_rms_superficie_mm_s"] = u_sup_mm
        ad.update(distorsion(ad["delta"], modos))
        barrido.append(ad)

    res = {
        "celda": {
            "nombre": c.nombre,
            "h_m": c.h,
            "procedencia_h": c.procedencia_h,
            "nu_m2_s": c.fluido.nu,
            "procedencia_nu": c.fluido.procedencia_nu,
            "procedencia_forzado": c.forzado.procedencia,
            "alfa_1_s": c.alfa(),
            "tau_alfa_s": c.tau_alfa(),
            "t_olvido_modal_s": c.t_olvido_modal(),
            "duracion_registro_s": c.camara.duracion(),
            "campo_de_vision_m": c.camara.campo_de_vision(),
        },
        "forzado": esp,
        "modos_verticales": modos,
        "casos_medidos": casos,
        "barrido_en_U": barrido,
    }

    os.makedirs(os.path.dirname(SALIDA), exist_ok=True)
    with open(SALIDA, "w") as f:
        json.dump(res, f, indent=1, ensure_ascii=False)

    # ---- informe corto por pantalla ------------------------------------------
    print("celda %s   h = %.1f mm   nu = %.1e m^2/s" % (c.nombre, c.h * 1e3,
                                                        c.fluido.nu))
    print("  alfa = %.5f 1/s   tau = %.2f s   registro = %.1f s"
          % (c.alfa(), c.tau_alfa(), c.camara.duracion()))
    print("  olvido del 2do modo vertical: %.2f s" % c.t_olvido_modal())
    print()
    print("FORZADO (medido sobre el patrón, no heredado)")
    print("  k dominante   = %.3f 1/m   (forma cerrada pi*sqrt(2)/a = %.3f, err %.1e)"
          % (esp["k_pico_1_m"], esp["k_pico_cerrado_1_m"],
             esp["error_relativo_vs_cerrado"]))
    print("  long. de onda = %.3f cm, orientación (%d,%d)*pi/a  -> %s"
          % (esp["longitud_de_onda_m"] * 100, esp["orientacion_pico"][0],
             esp["orientacion_pico"][1],
             "DIAGONAL" if 0 not in esp["orientacion_pico"] else "axial"))
    print("  energía en el fundamental: %.1f %%"
          % (100 * esp["fraccion_energia_en_fundamental"]))
    print("  longitud de onda equivalente al l = 2 cm heredado: k = 2 pi / l = %.1f 1/m"
          % (2 * math.pi / 0.02))
    if k_med is not None:
        print("  pico medido en la meseta forzada (mediana de %d espectros, t < 10 s): "
              "%.1f 1/m, bins de %.1f  -> cociente contra la red: %.3f"
              % (k_med["n"], k_med["k_1_m"], k_med["ancho_bin_1_m"],
                 k_med["k_1_m"] / esp["k_pico_1_m"]))
    print()
    print("MODOS VERTICALES (cuadratura contra forma cerrada)")
    print("  beta = %.12f  vs pi^2/8 = %.12f" % (modos["beta_medido"],
                                                 modos["beta_cerrado"]))
    print("  u(h)/<u> = %.12f  vs pi/2 = %.12f" % (modos["sesgo_superficie_medido"],
                                                   modos["sesgo_superficie_cerrado"]))
    print("  c_1 = %.9f  vs -8/(15 pi) = %.9f" %
          (modos["modos"][1]["acoplamiento_c_m"], modos["c1_cerrado"]))
    print("  c_m:", " ".join("%d:%+.4f" % (m["m"], m["acoplamiento_c_m"])
                             for m in modos["modos"][1:]))
    print()
    print("SOBRE LA CELDA   (eps = k h = %.4f, eps^2 = %.4f)" % (k * c.h,
                                                                 (k * c.h) ** 2))
    print("  la viscosidad horizontal sube la tasa lineal un %.1f %% sobre alfa"
          % (100 * (4 * (k * c.h) ** 2 / math.pi ** 2)))
    print()
    print("  u_sup   U=<u>    Re_h    delta   |a1/a0|  distorsión")
    print("  [mm/s]  [mm/s]                             rms")
    for b in barrido:
        print("  %5.2f   %5.3f  %6.2f  %6.2f   %6.4f   %6.4f"
              % (b["u_rms_superficie_mm_s"],
                 b["U_promedio_vertical_m_s"] * 1e3, b["Re_h"], b["delta"],
                 b["a1_sobre_a0"], b["distorsion_rms"]))
    print()
    for nombre, ad in casos.items():
        print("  %-20s delta = %5.2f  distorsión = %.3f  %s"
              % (nombre, ad["delta"], ad["distorsion_rms"],
                 "" if ad["verificado"] else "(velocidad NO verificada)"))
    print()
    print("escrito: %s" % os.path.relpath(SALIDA, RAIZ))
    return 0


if __name__ == "__main__":
    sys.exit(main())
