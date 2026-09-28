#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Qué cambia la superficie deformable en la celda real, y cuánto cuesta integrarla.

Dos cuentas, ambas sobre el modelo lineal que implementa SPECTER ([D-49]):

1. Qué modos mueve la superficie. Un modo horizontal k se separa en una parte TOROIDAL
   (velocidad horizontal perpendicular a k, w = 0, p = 0: el flujo cuasi-2D del
   experimento) y una POLOIDAL (velocidad en el plano de k y z, con w). Con w = 0 y p = 0
   las tres condiciones de la superficie deformable se reducen a las de la plana, así que
   a orden lineal la parte toroidal decae EXACTAMENTE con nu((pi/2h)^2 + k^2), la tasa de
   la Fase 1 [V-07], con o sin deformación (derivación de esta sesión). La superficie
   sólo entra en la parte poloidal, que son ondas de gravedad-capilaridad amortiguadas;
   acá se calcula su frecuencia y su amortiguamiento para la celda con la referencia MAC
   de la puerta de la Fase 5 (importada, no copiada), con Richardson, en variables
   adimensionales (longitud h, tiempo h^2/nu). El acople con el flujo vortical es no
   lineal, eta ~ U^2/g, que es la estimación de [D-28].

2. El paso de tiempo que imponen las ondas de gravedad-capilaridad, integradas en forma
   explícita, contra el viscoso vertical, para las mallas de producción. Estimación de
   estabilidad lineal de esta sesión (no validada contra bibliografía): con ORD=2 el eje
   imaginario es marginalmente inestable, |G|^2 = 1 + (w dt)^4/4, y lo estabiliza sólo el
   amortiguamiento viscoso de la onda, que acá se toma como 2 nu k^2 (el de superficie).
   NO es conservador en todo el rango: para kh >~ 3 el amortiguamiento real es 2-4 % menor
   (lo mostró el verificador de [V-20]); con el autovalor exacto el límite RK2 queda apenas
   por encima del de esta fórmula en 64^2 y 128^2, así que la cota vale igual.

Los parámetros de la celda salen de codigo/celda.py; g, rho y sigma, de
numerico/fase4/superficie_libre_escalas.py, con su procedencia (agua limpia).

Salida: salidas/tablas/superficie_deformable_escalas.json y la figura
salidas/figuras/superficie_deformable_escalas.png.
"""

import json
import math
import os
import sys

import numpy as np
from scipy.linalg import eig

RAIZ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(RAIZ, "codigo"))
sys.path.insert(0, os.path.join(RAIZ, "verificacion"))
sys.path.insert(0, os.path.join(RAIZ, "numerico", "fase4"))
from celda import CAMPANA_02_06_25 as CELDA                 # noqa: E402
import superficie_libre_escalas as ESC                      # noqa: E402
import test_aceptacion_fase5 as PUERTA                      # noqa: E402

SALIDA = os.path.join(RAIZ, "salidas", "tablas", "superficie_deformable_escalas.json")
FIGURA = os.path.join(RAIZ, "salidas", "figuras", "superficie_deformable_escalas.png")


def onda_poloidal(K, G, n):
    """Autovalor poloidal menos amortiguado (una onda: parte imaginaria positiva),
    adimensional (nu = h = 1)."""
    J, _, _ = PUERTA._mac(n, K, 1.0, 1.0, G)
    lam = eig(J, right=False)
    lam = lam[lam.imag > 0]
    return complex(lam[np.argmax(lam.real)])


def main():
    h, nu = CELDA.h, CELDA.fluido.nu
    g, gam = ESC.G, ESC.TENSION / ESC.RHO
    k_forz = CELDA.forzado.k_fundamental()
    ks = np.array(sorted(set(np.round(np.geomspace(40.0, 600.0, 13), 1)) | {k_forz}))

    filas = []
    for k in ks:
        K = k * h
        G = (g + gam * k * k) * h ** 3 / nu ** 2
        l = {n: onda_poloidal(K, G, n) for n in (50, 100, 200)}
        r1 = (4 * l[100] - l[50]) / 3
        r2 = (4 * l[200] - l[100]) / 3
        lam = (16 * r2 - r1) / 15
        esc = nu / h ** 2
        toroidal = (math.pi ** 2 / 4 + K * K) * esc
        filas.append({
            "k_1_m": float(k), "kh": K, "G_adim": G,
            "onda_frecuencia_1_s": lam.imag * esc,
            "onda_amortiguamiento_1_s": -lam.real * esc,
            "onda_error_richardson_relativo": abs(r2 - r1) / abs(lam),
            "omega_inviscida_1_s": math.sqrt((g * k + gam * k ** 3) * math.tanh(k * h)),
            "lambda_toroidal_1_s": toroidal,
        })
        f = filas[-1]
        print("k = %6.1f 1/m  onda: omega %.2f 1/s, amortig. %.3f 1/s (inviscida %.2f); "
              "toroidal %.4f 1/s; Richardson %.1e"
              % (k, f["onda_frecuencia_1_s"], f["onda_amortiguamiento_1_s"],
                 f["omega_inviscida_1_s"], toroidal, f["onda_error_richardson_relativo"]))

    # 2. Paso de tiempo en las mallas de producción. La caja horizontal es el campo del
    #    PIV (como en verificacion/decaimiento_banda_ancha.py) y el k máximo retenido
    #    es el de la regla de 2/3.
    fov = CELDA.camara.campo_de_vision()
    mallas = []
    for nxy in (32, 64, 128):
        for nz_fis in (39, 103):
            kmax = (2.0 / 3.0) * math.pi * nxy / fov
            wmax = math.sqrt((g * kmax + gam * kmax ** 3) * math.tanh(kmax * h))
            dz = h / (nz_fis - 1)
            dt_visc = 2.0 / (nu * (math.pi / dz) ** 2)
            amort = 2 * nu * kmax ** 2
            dt_rk2 = (8.0 * amort / wmax ** 4) ** (1.0 / 3.0)
            dt_rk4 = 2.0 * math.sqrt(2.0) / wmax
            mallas.append({"nxy": nxy, "nz_fisico": nz_fis, "kmax_1_m": kmax,
                           "omega_max_1_s": wmax, "dt_viscoso_s": dt_visc,
                           "dt_ondas_rk2_s": dt_rk2, "dt_ondas_rk4_s": dt_rk4,
                           "cociente_rk2": dt_visc / dt_rk2})
            print("%3d^2 x %3d: dt viscoso %.2e s, ondas RK2 %.2e s (x%.1f), RK4 %.2e s"
                  % (nxy, nz_fis, dt_visc, dt_rk2, dt_visc / dt_rk2, dt_rk4))

    res = {
        "descripcion": __doc__.split("\n\n")[0].strip(),
        "celda": {"h_m": h, "nu_m2_s": nu, "k_forzado_1_m": k_forz,
                  "campo_de_vision_m": fov},
        "g_m_s2": g, "sigma_sobre_rho_m3_s2": gam,
        "procedencia_g_rho_sigma": "numerico/fase4/superficie_libre_escalas.py",
        "tasa_vortical": filas,
        "paso_de_tiempo": mallas,
        "salvedades": [
            "Modelo lineal en la frontera; la deformación no lineal forzada, "
            "eta ~ U^2/g, es la de [D-28] y no entra acá.",
            "El límite de dt con ORD=2 es una estimación de estabilidad lineal de esta "
            "sesión, con amortiguamiento 2 nu k^2 (no conservador para kh > ~3, por 2-4 %); "
            "no está validado con bibliografía.",
        ],
    }
    os.makedirs(os.path.dirname(SALIDA), exist_ok=True)
    with open(SALIDA, "w") as f:
        json.dump(res, f, indent=2, ensure_ascii=False)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    k = np.array([f["k_1_m"] for f in filas])
    ax[0].loglog(k, [f["onda_frecuencia_1_s"] for f in filas], "o-", color="C0",
                 label="onda poloidal: frecuencia")
    ax[0].loglog(k, [f["omega_inviscida_1_s"] for f in filas], ":", color="C0",
                 label="onda, límite inviscido (control)")
    ax[0].loglog(k, [f["onda_amortiguamiento_1_s"] for f in filas], "s-", color="C1",
                 label="onda poloidal: amortiguamiento")
    ax[0].loglog(k, [f["lambda_toroidal_1_s"] for f in filas], "^-", color="C2",
                 label="modo toroidal (vortical): tasa")
    ax[0].axvline(k_forz, color="k", lw=0.8, ls="--")
    ax[0].text(k_forz * 1.05, 3.0, "k forzado", fontsize=8)
    ax[0].set_xlabel("k [1/m]")
    ax[0].set_ylabel("[1/s]")
    ax[0].set_title("Celda real: lo que mueve la superficie", fontsize=10)
    ax[0].legend(fontsize=7)
    for nz_fis, mk in ((39, "s"), (103, "o")):
        sel = [m for m in mallas if m["nz_fisico"] == nz_fis]
        n = [m["nxy"] for m in sel]
        ax[1].loglog(n, [m["dt_viscoso_s"] for m in sel], mk + "-", color="C0",
                     label="viscoso, nz=%d" % nz_fis)
    sel = [m for m in mallas if m["nz_fisico"] == 103]
    n = [m["nxy"] for m in sel]
    ax[1].loglog(n, [m["dt_ondas_rk2_s"] for m in sel], "^-", color="C3",
                 label="ondas, ORD=2")
    ax[1].loglog(n, [m["dt_ondas_rk4_s"] for m in sel], "v-", color="C2",
                 label="ondas, ORD=4")
    ax[1].set_xlabel("puntos por lado")
    ax[1].set_ylabel("dt máximo [s]")
    ax[1].set_title("Paso de tiempo: caja = campo del PIV", fontsize=10)
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    os.makedirs(os.path.dirname(FIGURA), exist_ok=True)
    fig.savefig(FIGURA, dpi=150)
    print("escrito", SALIDA, "y", FIGURA)


if __name__ == "__main__":
    main()
