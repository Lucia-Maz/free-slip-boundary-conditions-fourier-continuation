#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Qué modelo de superficie libre hace falta en esta celda, con números.

Acompaña a teoria/superficie_libre_v_estrella_y_p.md. Tres preguntas, y las tres se
contestan con cuentas y no con criterio:

  1. ¿Cuánto se deforma la superficie? Si la deformación es de micrones sobre 6 mm, la
     superficie PLANA que ya está implementada es la correcta y la rama
     Neumann-Dirichlet de la presión no hace falta.
  2. ¿Gravedad o tensión superficial? Decide cuál es la fuerza restitutiva y, por lo
     tanto, qué término manda en la condición de tensión normal.
  3. ¿Cuánto cambia alfa una película superficial? Es la salvedad más fuerte del
     trabajo. El modelo de Boussinesq-Scriven da una familia de un parámetro que
     interpola entre tope libre y tope rígido, y con ella la salvedad deja de ser
     cualitativa.

Uso:
    /home/lucia/miniforge3/envs/piv-dt/bin/python numerico/fase4/superficie_libre_escalas.py
"""

import json
import math
import os
import sys

import numpy as np
from scipy.optimize import brentq

RAIZ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(RAIZ, "codigo"))
from celda import CAMPANA_02_06_25 as CELDA          # noqa: E402

SALIDA = os.path.join(RAIZ, "salidas", "tablas", "superficie_libre_escalas.json")

G = 9.81            # m/s^2
RHO = 1000.0        # kg/m^3, agua; el electrolito es más denso, ver salvedades
TENSION = 0.072     # N/m, agua-aire limpia; una interfaz contaminada tiene menos


def deformacion(U, k, h, rho=RHO, sigma=TENSION, g=G):
    """Amplitud de la deformación de la superficie ante fluctuaciones de presión.

    La presión dinámica de un flujo de velocidad U es del orden de rho U^2. Contra ella
    la superficie se restituye con rho g eta (gravedad) más sigma k^2 eta (tensión
    superficial), de modo que

        eta ~ rho U^2 / (rho g + sigma k^2) .

    El número de Froude U/sqrt(g h) da lo mismo cuando manda la gravedad: eta/h ~ Fr^2.
    """
    restitucion = rho * g + sigma * k ** 2
    eta = rho * U ** 2 / restitucion
    return {
        "eta_m": eta,
        "eta_sobre_h": eta / h,
        "froude": U / math.sqrt(g * h),
        "restitucion_gravedad_Pa_por_m": rho * g,
        "restitucion_capilar_Pa_por_m": sigma * k ** 2,
        "manda": "gravedad" if rho * g > sigma * k ** 2 else "tensión superficial",
        "longitud_capilar_m": math.sqrt(sigma / (rho * g)),
        "k_capilar_1_m": math.sqrt(rho * g / sigma),
    }


def alfa_con_pelicula(beta, h, nu):
    """alfa cuando el tope tiene una película con viscosidad superficial.

    Modelo de Boussinesq-Scriven en su forma más simple: el balance de tensión
    tangencial en la superficie deja de ser mu du/dz = 0 y pasa a ser

        mu du/dz|_h = - mu_s k^2 u(h) ,

    es decir una condición de Robin. Con u = sin(q z) —que ya cumple u(0)=0— la
    condición de arriba da

        (q h) cot(q h) = -beta ,   beta = mu_s k^2 h / mu ,

    y la tasa del modo más lento es lambda = nu q^2. Los dos límites conocidos salen
    solos: beta -> 0 da q h = pi/2 (tope libre, alfa = pi^2 nu/4h^2) y beta -> infinito
    da q h = pi (canal, alfa = pi^2 nu/h^2). El cociente entre ellos es 4, que es el
    resultado central de la Fase 1.
    """
    if beta < 0:
        raise ValueError("beta no puede ser negativo")
    lo, hi = math.pi / 2 + 1e-12, math.pi - 1e-12
    if beta == 0.0:
        x = math.pi / 2
    else:
        f = lambda x: x / math.tan(x) + beta          # noqa: E731
        x = brentq(f, lo, hi, xtol=1e-14, rtol=1e-15)
    return nu * (x / h) ** 2, x


def main():
    c = CELDA
    h, nu = c.h, c.fluido.nu
    k = c.forzado.k_fundamental()
    mu = RHO * nu

    # --- 1 y 2: deformación de la superficie ----------------------------------
    casos = {}
    for nombre, u_sup in (("forzado_2.15A", 1.3214710285395503e-3),
                          ("inicio_decaimiento_medido", 3.63e-3)):
        U = u_sup / c.sesgo_piv_superficie()
        casos[nombre] = deformacion(U, k, h)
        casos[nombre]["U_promedio_vertical_m_s"] = U
        casos[nombre]["u_rms_superficie_m_s"] = u_sup

    # --- 3: la película -------------------------------------------------------
    a_libre, _ = alfa_con_pelicula(0.0, h, nu)
    a_rigido = nu * (math.pi / h) ** 2
    filas = []
    for beta in (0.0, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 1000.0):
        a, x = alfa_con_pelicula(beta, h, nu)
        mu_s = beta * mu / (k ** 2 * h)
        filas.append({"beta": beta, "qh_sobre_pi": x / math.pi, "alfa_1_s": a,
                      "alfa_sobre_alfa_libre": a / a_libre,
                      "mu_s_N_s_por_m": mu_s})

    # beta que hace falta para mover alfa un 10 %
    def desvio(beta, objetivo):
        return alfa_con_pelicula(beta, h, nu)[0] / a_libre - objetivo
    beta_10 = brentq(lambda b: desvio(b, 1.10), 1e-6, 1e4, xtol=1e-10)
    beta_2x = brentq(lambda b: desvio(b, 2.00), 1e-6, 1e6, xtol=1e-10)

    res = {
        "celda": {"h_m": h, "nu_m2_s": nu, "k_forzado_1_m": k, "mu_Pa_s": mu,
                  "rho_kg_m3": RHO, "tension_N_m": TENSION,
                  "procedencia_rho_y_tension": (
                      "agua limpia a ~20 C. NO son del electrolito real (KNO3 16 % m/m, "
                      "más denso) ni de una interfaz con partículas flotando, que tiene "
                      "menos tensión. Supuesto declarado, no medido.")},
        "deformacion": casos,
        "pelicula": {
            "alfa_tope_libre_1_s": a_libre,
            "alfa_canal_1_s": a_rigido,
            "cociente": a_rigido / a_libre,
            "beta_para_10pc": beta_10,
            "mu_s_para_10pc_N_s_por_m": beta_10 * mu / (k ** 2 * h),
            "beta_para_duplicar": beta_2x,
            "mu_s_para_duplicar_N_s_por_m": beta_2x * mu / (k ** 2 * h),
            "tabla": filas,
        },
    }
    os.makedirs(os.path.dirname(SALIDA), exist_ok=True)
    with open(SALIDA, "w") as f:
        json.dump(res, f, indent=1, ensure_ascii=False)

    print("DEFORMACIÓN DE LA SUPERFICIE   (h = %.1f mm, k = %.1f 1/m)" % (h * 1e3, k))
    print("  longitud capilar = %.2f mm  ->  k_capilar = %.0f 1/m"
          % (res["deformacion"]["forzado_2.15A"]["longitud_capilar_m"] * 1e3,
             res["deformacion"]["forzado_2.15A"]["k_capilar_1_m"]))
    for nombre, d in casos.items():
        print("  %-26s U=%.3f mm/s  Fr=%.4f  eta=%.2f um  eta/h=%.1e  manda: %s"
              % (nombre, 1e3 * d["U_promedio_vertical_m_s"], d["froude"],
                 1e6 * d["eta_m"], d["eta_sobre_h"], d["manda"]))
    print()
    print("PELÍCULA SUPERFICIAL  (beta = mu_s k^2 h / mu)")
    print("  alfa tope libre = %.5f 1/s     alfa canal = %.5f 1/s   cociente %.3f"
          % (a_libre, a_rigido, a_rigido / a_libre))
    print("  beta      qh/pi   alfa [1/s]  alfa/alfa_libre   mu_s [N s/m]")
    for f in filas:
        print("  %8.3g  %.4f  %10.5f  %12.3f   %.3e"
              % (f["beta"], f["qh_sobre_pi"], f["alfa_1_s"],
                 f["alfa_sobre_alfa_libre"], f["mu_s_N_s_por_m"]))
    print("  +10 %% en alfa con beta = %.4f  ->  mu_s = %.2e N s/m"
          % (beta_10, res["pelicula"]["mu_s_para_10pc_N_s_por_m"]))
    print("  duplicar alfa con beta = %.3f  ->  mu_s = %.2e N s/m"
          % (beta_2x, res["pelicula"]["mu_s_para_duplicar_N_s_por_m"]))
    print("\nescrito: %s" % os.path.relpath(SALIDA, RAIZ))
    return 0


if __name__ == "__main__":
    sys.exit(main())
