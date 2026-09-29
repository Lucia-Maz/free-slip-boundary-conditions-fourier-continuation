#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Figura de los campos de superficie de la corrida con superficie deformable de orden 2
(numerico/fase6/campos_superficie.py), en unidades físicas. No corre SPECTER.

  (a) vorticidad vertical en la superficie, ω_z(x, y, z = h), en el último instante
  (b) elevación η(x, y), en μm, en el mismo instante
  (c) η_rms(t) del diagnóstico de orden, con la estimación cuasi-estática de su nivel

La estimación cuasi-estática es una derivación de esta sesión, no validada contra
bibliografía, y sólo da orden de magnitud y forma:
- la presión en la superficie sale de la ecuación de Poisson 2D con la velocidad de
  superficie, ∇²(p/ρ) = −∂_i u_j ∂_j u_i;
- el balance normal estático da (g − (σ/ρ)∇²)η = p/ρ.
El problema real es 3D (∂_zz p no es cero con kh ~ 1) y la superficie arranca plana, así
que η lleva además ondas de gravedad-capilaridad.

Uso: $PY numerico/fase6/figura_campos_superficie.py [--nxy 64] [--memorias 4] [--prueba]
     -> salidas/figuras/f_fs6_superficie.png (y .pdf)
Por omisión lee la corrida de 4 memorias de Sakura ([V-23]): es la figura 4 de la entrega
final. La de la laptop es --memorias 1.5.
"""

import argparse
import json
import os
import sys

import numpy as np

RAIZ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def cuasi_estatica(u, v, g, s_rho, L):
    """eta [m] de (g - s ∇²) eta = p/rho con ∇² p/rho = -du_i/dx_j du_j/dx_i, en 2D."""
    n = u.shape[0]
    k = np.fft.fftfreq(n, d=L / n) * 2 * np.pi          # 1/m, caja de lado L
    KX, KY = np.meshgrid(k, k, indexing="ij")
    d = lambda f, K: np.real(np.fft.ifft2(1j * K * np.fft.fft2(f)))
    ux, uy, vx, vy = d(u, KX), d(u, KY), d(v, KX), d(v, KY)
    rhs = -(ux * ux + 2 * uy * vx + vy * vy)
    K2 = KX ** 2 + KY ** 2
    P = np.zeros_like(K2, dtype=complex)
    m = K2 > 0
    P[m] = -np.fft.fft2(rhs)[m] / K2[m]
    E = np.zeros_like(P)
    E[m] = P[m] / (g + s_rho * K2[m])
    return np.real(np.fft.ifft2(E))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nxy", type=int, default=64)
    ap.add_argument("--memorias", type=float, default=4.0)
    ap.add_argument("--prueba", action="store_true")
    args = ap.parse_args()
    suf = "_prueba" if args.prueba else ""
    nombre = "fase6_campos_superficie_n%d_m%g%s" % (args.nxy, 0.2 if args.prueba else args.memorias, suf)
    J = json.load(open(os.path.join(RAIZ, "salidas", "tablas", nombre + ".json")))
    C = np.load(os.path.join(RAIZ, "salidas", "campos", nombre + ".npz"))
    cfg = J["config"]
    L = 2 * np.pi * cfg["L_c_m"]                         # lado de la caja [m]
    g, s = cfg["g_m_s2"], cfg["sigma_rho_m3_s2"]
    j = len(C["t_s"]) - 1
    eta, wz = C["eta_m"][j], C["omega_z_1_s"][j]
    qs = [cuasi_estatica(C["u_m_s"][i], C["v_m_s"][i], g, s, L) for i in range(len(C["t_s"]))]
    corr = float(np.corrcoef(eta.ravel(), qs[j].ravel())[0, 1])
    cs = [float(np.corrcoef(C["eta_m"][i].ravel(), qs[i].ravel())[0, 1]) for i in range(1, j + 1)]

    import matplotlib
    matplotlib.use("Agg")
    sys.path.insert(0, os.path.expanduser("~/Desktop/PhD/Num/Codigos"))
    try:
        import estetica  # noqa: F401  (rcParams compartidos)
        import matplotlib.pyplot as plt
        plt.rcParams["font.serif"] = ["Times New Roman", "Nimbus Roman", "Liberation Serif"]
        plt.rcParams["font.family"] = "serif"
        ancho = estetica.width_in
    except ImportError:
        import matplotlib.pyplot as plt
        ancho = 7.0
    fig, ax = plt.subplots(1, 3, figsize=(ancho, 0.36 * ancho),
                           gridspec_kw={"width_ratios": [1, 1, 1.15]})
    ext = [0, L * 100, 0, L * 100]                       # cm
    a = ax[0]
    vm = float(np.max(np.abs(wz)))
    im0 = a.imshow(wz.T, origin="lower", extent=ext, cmap="coolwarm", vmin=-vm, vmax=vm)
    a.set_title(r"(a) $\omega_z$ en la superficie [s$^{-1}$]", fontsize=8)
    fig.colorbar(im0, ax=a, fraction=0.046, pad=0.03)
    a = ax[1]
    ve = float(np.max(np.abs(eta))) * 1e6
    im1 = a.imshow(eta.T * 1e6, origin="lower", extent=ext, cmap="BrBG", vmin=-ve, vmax=ve)
    a.set_title(r"(b) $\eta$ [$\mu$m], $t=%.1f$ s" % C["t_s"][j], fontsize=8)
    fig.colorbar(im1, ax=a, fraction=0.046, pad=0.03)
    for a in ax[:2]:
        a.set_xlabel("x [cm]")
        a.set_aspect("equal")
    ax[0].set_ylabel("y [cm]")
    a = ax[2]
    S = J["serie_diagnostico"]
    a.plot(S["t_s"], np.asarray(S["eta_rms_m"]) * 1e6, "-", color="C0", lw=1.0,
           label=r"$\eta_{\rm rms}$, SPECTER ($N=2$)")
    qr = [np.sqrt(np.mean(q ** 2)) * 1e6 for q in qs]
    a.plot(C["t_s"], qr, "o--", color="C3", ms=3, lw=0.8, label=r"cuasi-estática 2D, $p/\rho g$ (§)")
    a.set_ylim(0, 1.3 * max(max(qr), 1e6 * max(S["eta_rms_m"])))          # lugar para la leyenda
    a.set_xlabel("t [s]")
    a.set_ylabel(r"$\eta_{\rm rms}$ [$\mu$m]")
    a.set_title(r"(c) corr. con la cuasi-estática: %.2f–%.2f" % (min(cs), max(cs)), fontsize=8)
    a.legend(fontsize=7)
    fig.savefig(os.path.join(RAIZ, "salidas", "figuras", "f_fs6_superficie%s.png" % suf), dpi=200)
    fig.savefig(os.path.join(RAIZ, "salidas", "figuras", "f_fs6_superficie%s.pdf" % suf))
    print("correlación eta / cuasi-estática en t = %.2f s: %.3f (rango %.2f-%.2f); eta_rms %.3e m; qs_rms %.3e m" % (
        C["t_s"][j], corr, min(cs), max(cs), float(np.sqrt(np.mean(eta ** 2))), float(np.sqrt(np.mean(qs[j] ** 2)))))
    print("-> salidas/figuras/f_fs6_superficie%s.png" % suf)


if __name__ == "__main__":
    main()
