#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Chequeo algebraico de las fórmulas de orden N de la superficie deformable (Fase 6).

Qué verifica: que las fórmulas truncadas que implementa SPECTER para la condición
cinemática K, la tensión tangencial T_alpha y la tensión normal Phi, a orden N en la
amplitud, coinciden EXACTAMENTE (a redondeo) con el desarrollo de Taylor de las condiciones
exactas de una superficie gráfica z = h + eta, evaluadas en la superficie desplazada.

Cómo, sin álgebra simbólica: en un punto de la superficie, los datos locales son eta, la
pendiente s = grad(eta), las derivadas segundas de eta (para la curvatura) y las derivadas
normales f_m = d^m f/dz^m en z = h de cada campo (u, v, w, p y sus derivadas horizontales).
Con todo escalado por epsilon (campos ~ eta ~ epsilon: conteo de amplitud), la condición
exacta F(epsilon) es analítica en epsilon, y sus coeficientes de Taylor salen exactos de la
integral de Cauchy sobre un círculo (una FFT). La fórmula truncada de orden N tiene que ser
el polinomio de Taylor de grado N de F, coeficiente por coeficiente.

Las condiciones exactas son las de Wang, Tice & Kim, ARMA 212, 1-92 (2014),
DOI 10.1007/s00205-013-0700-2, §1.1 ecs. (1.4)-(1.8), leídas en el original: normal hacia
arriba N = (-s, 1), tangentes t_alpha = e_alpha + s_alpha e_z, curvatura
kappa = div(s / sqrt(1+|s|^2)), y p la presión cinemática menos la hidrostática del reposo.
Las fórmulas truncadas son de esta sesión ([D-51]).

Uso:
    /home/lucia/miniforge3/envs/piv-dt/bin/python numerico/fase6/chequeo_expansion.py
"""

import json
import math
import os
import sys

import numpy as np

RAIZ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SALIDA = os.path.join(RAIZ, "salidas", "tablas", "fase6_chequeo_expansion.json")

MTOP = 14                 # derivadas normales disponibles (la serie "exacta")
G, GAM, NU = 1.3, 0.7, 0.11
C_KAPPA = [1.0, -0.5, 3.0 / 8.0, -5.0 / 16.0]   # binom(-1/2, j): (1+x)^(-1/2)


def datos(rng):
    """Datos locales aleatorios en un punto de la superficie. Las derivadas normales
    decaen como 1/m! para que la serie de Taylor converja en el círculo de Cauchy."""
    d = {"eta": rng.normal(), "sx": rng.normal(), "sy": rng.normal(),
         "exx": rng.normal(), "exy": rng.normal(), "eyy": rng.normal()}
    for campo in ("u", "v", "w", "p"):
        d[campo] = rng.normal(size=MTOP + 2)
    for campo in ("u", "v", "w"):
        for der in ("x", "y"):
            d[campo + der] = rng.normal(size=MTOP + 2)
    return d


def tr(f, eta, M):
    """Tr_M[f] = sum_{m=0}^{M} eta^m/m! f_m (0 si M < 0)."""
    return sum(eta**m / math.factorial(m) * f[m] for m in range(0, M + 1)) if M >= 0 else 0.0


def componentes_S(d, m):
    """d^m/dz^m de las componentes de S = grad v + grad v^T en z = h."""
    return {"xx": 2 * d["ux"][m], "yy": 2 * d["vy"][m],
            "xy": d["uy"][m] + d["vx"][m],
            "xz": d["u"][m + 1] + d["wx"][m], "yz": d["v"][m + 1] + d["wy"][m],
            "zz": 2 * d["w"][m + 1]}


def trS(d, comp, eta, M):
    return sum(eta**m / math.factorial(m) * componentes_S(d, m)[comp]
               for m in range(0, M + 1)) if M >= 0 else 0.0


def exactas(d, e):
    """Condiciones exactas en z = h + eta, con todo escalado por e (amplitud)."""
    eta, sx, sy = e * d["eta"], e * d["sx"], e * d["sy"]
    exx, exy, eyy = e * d["exx"], e * d["exy"], e * d["eyy"]
    M = MTOP - 1
    # campos escalados por e, evaluados en la superficie desplazada (serie completa)
    bar = lambda f: e * tr(f, eta, M)
    barS = lambda c: e * trS(d, c, eta, M)
    s = (sx, sy)
    S = {c: barS(c) for c in ("xx", "yy", "xy", "xz", "yz", "zz")}
    Sab = [[S["xx"], S["xy"]], [S["xy"], S["yy"]]]
    Saz = [S["xz"], S["yz"]]
    K = bar(d["w"]) - bar(d["u"]) * sx - bar(d["v"]) * sy
    T = [Saz[a] - sum(Sab[a][b] * s[b] for b in range(2)) + s[a] * S["zz"]
         - s[a] * sum(s[b] * Saz[b] for b in range(2)) for a in range(2)]
    q2 = sx**2 + sy**2
    B = (S["zz"] - 2 * (sx * Saz[0] + sy * Saz[1])
         + sum(s[a] * s[b] * Sab[a][b] for a in range(2) for b in range(2))) / (1 + q2)
    kappa = (exx * (1 + sy**2) + eyy * (1 + sx**2) - 2 * sx * sy * exy) / (1 + q2)**1.5
    Phi = bar(d["p"]) - G * eta + GAM * kappa - NU * B
    return np.array([K, T[0], T[1], Phi])


def truncadas(d, e, N):
    """Las fórmulas de orden N que implementa SPECTER (D-51)."""
    eta, sx, sy = e * d["eta"], e * d["sx"], e * d["sy"]
    exx, exy, eyy = e * d["exx"], e * d["exy"], e * d["eyy"]
    s = (sx, sy)
    f = lambda nombre, M: e * tr(d[nombre], eta, M)
    fS = lambda c, M: e * trS(d, c, eta, M)
    K = f("w", N - 1) - s[0] * f("u", N - 2) - s[1] * f("v", N - 2)
    ab = {(0, 0): "xx", (0, 1): "xy", (1, 0): "xy", (1, 1): "yy"}
    az = {0: "xz", 1: "yz"}
    T = []
    for a in range(2):
        t = fS(az[a], N - 1)
        t -= sum(s[b] * fS(ab[(a, b)], N - 2) for b in range(2))
        t += s[a] * fS("zz", N - 2)
        t -= s[a] * sum(s[b] * fS(az[b], N - 3) for b in range(2))
        T.append(t)
    q2 = sx**2 + sy**2
    B = 0.0
    j = 0
    while 2 * j <= N - 1:
        B += (-q2)**j * (fS("zz", N - 1 - 2 * j)
                          - 2 * sum(s[a] * fS(az[a], N - 2 - 2 * j) for a in range(2))
                          + sum(s[a] * s[b] * fS(ab[(a, b)], N - 3 - 2 * j)
                                for a in range(2) for b in range(2)))
        j += 1
    # curvatura: div(s * sum_j c_j |s|^{2j}), 2j+1 <= N, desarrollada en el punto
    # d/dx(|s|^2) = 2(sx exx + sy exy), d/dy(|s|^2) = 2(sx exy + sy eyy)
    kappa = 0.0
    j = 0
    while 2 * j + 1 <= N:
        cj = C_KAPPA[j]
        divs = exx + eyy
        grad_q2 = (2 * (sx * exx + sy * exy), 2 * (sx * exy + sy * eyy))
        term = q2**j * divs
        if j > 0:
            term += j * q2**(j - 1) * (sx * grad_q2[0] + sy * grad_q2[1])
        kappa += cj * term
        j += 1
    Phi = f("p", N - 1) - G * eta + GAM * kappa - NU * B
    return np.array([K, T[0], T[1], Phi])


def coef_taylor(fun, grado, r=0.25, Q=128):
    """Coeficientes de Taylor 0..grado de una función analítica de e, por Cauchy."""
    ek = r * np.exp(2j * np.pi * np.arange(Q) / Q)
    vals = np.array([fun(x) for x in ek])                 # (Q, 4)
    return np.array([np.mean(vals * ek[:, None]**(-n), axis=0) for n in range(grado + 1)])


def main():
    rng = np.random.default_rng(20260928)
    peor = {}
    for caso in range(12):
        d = datos(rng)
        cex = coef_taylor(lambda x: exactas(d, x), 6)
        for N in range(1, 6):
            ctr = coef_taylor(lambda x: truncadas(d, x, N), 6)
            # grados 0..N tienen que coincidir; los de grado > N de la truncada, ser...
            # no nulos en general (la truncada es un polinomio de grado > N por los
            # productos), así que sólo se comparan 0..N
            escala = np.max(np.abs(cex[:N + 1])) + 1e-300
            err = float(np.max(np.abs(ctr[:N + 1] - cex[:N + 1])) / escala)
            peor[N] = max(peor.get(N, 0.0), err)
    # control negativo: una fórmula con un signo cambiado en la pendiente tangencial
    # tiene que fallar desde el orden 2
    d = datos(rng)
    cex = coef_taylor(lambda x: exactas(d, x), 3)

    def mala(x):
        v = truncadas(d, x, 2)
        eta, sx = x * d["eta"], x * d["sx"]
        v[1] += 2 * sx * x * 2 * d["ux"][0]      # + 2 s_x S_xx: invierte el término -S_xb s_b
        return v
    cmala = coef_taylor(mala, 3)
    err_malo = float(np.max(np.abs(cmala[:3] - cex[:3])) / np.max(np.abs(cex[:3])))
    res = {"peor_error_relativo_por_orden": {str(k): v for k, v in peor.items()},
           "control_negativo_signo_pendiente_orden2": err_malo}
    print(json.dumps(res, indent=1))
    ok = all(v < 1e-10 for v in peor.values()) and err_malo > 1e-3
    os.makedirs(os.path.dirname(SALIDA), exist_ok=True)
    res["pasa"] = ok
    with open(SALIDA, "w") as fsal:
        json.dump(res, fsal, indent=1)
    print("PASA" if ok else "FALLA")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
