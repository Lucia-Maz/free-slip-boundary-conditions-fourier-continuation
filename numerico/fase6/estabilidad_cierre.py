#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Estabilidad del cierre de la condición tangencial de orden N en SPECTER, en una réplica 1D
del paso de Runge-Kutta ([H-23], [H-24] de DECISIONES.md). Derivación y réplica de esta
sesión, no validadas contra bibliografía.

La réplica es la velocidad tangencial de una columna con las mismas piezas que SPECTER:
  - derivadas FC-Gram espectrales (tablas A25-5, Q5, Q1n5 del árbol de trabajo);
  - fondo Dirichlet u = 0 (no deslizante) y tope Neumann reconstruido con
    neumann_reconstruct (pesos Q5 · Q1n5, como fcgram_mod.f90);
  - RK2 de SPECTER (o = 2, 1; u^(o) = u^0 + dt/o nu D2 u^(o+1));
  - el dato de Neumann de orden N transferido con eta constante:
        N = 2:  u_z(Lz) = -eta u_zz(Lz)
        N = 3:  u_z(Lz) = -eta u_zz(Lz) - eta^2/2 u_zzz(Lz)
Cierres comparados para u_zz, u_zzz en el tope:
  'espectral'  traza espectral de la columna reconstruida, resuelta exactamente (implícita):
               es el cierre de la implementación final;
  'interior q' extrapolación polinomial desde q nodos interiores (sin el de borde): el
               cierre que congelaba intent_fase6.txt y resultó inestable.
Salida: radio espectral del paso completo contra r = eta/dz, la respuesta beta_m de la traza
al dato de Neumann (la que el Fortran mide en fs_gen_setup y debe coincidir), y la tasa del
modo espurio de N = 2 contra la del continuo, nu/eta^2.

  salidas/tablas/fase6_estabilidad_cierre.json
  salidas/figuras/f_fs6_estabilidad.png  (y .pdf)

Liviano: ~40 MB, unos segundos.
"""

import json
import math
import os
import sys

import numpy as np

RAIZ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TABLAS = os.path.join(RAIZ, "numerico", "SPECTER-trabajo", "tables")
SALIDA = os.path.join(RAIZ, "salidas", "tablas", "fase6_estabilidad_cierre.json")
FIG = os.path.join(RAIZ, "salidas", "figuras", "f_fs6_estabilidad")

C, D = 25, 5
NU, DT = 1e-2, 1e-3                  # los de la puerta de la Fase 6


def tablas():
    A = np.fromfile(os.path.join(TABLAS, "A%d-%d.dat" % (C, D)), dtype="<f8").reshape(
        (C, D), order="F")
    Q = np.fromfile(os.path.join(TABLAS, "Q%d.dat" % D), dtype="<f8").reshape((D, D), order="F")
    raw = np.fromfile(os.path.join(TABLAS, "Q1n%d.dat" % D), dtype="<f8")
    dxp, Qn = raw[0], raw[1:].reshape((D, D), order="F")
    return A @ Q.T, Q[D - 1, :] @ Qn.T, dxp


def fornberg(z0, xs, mmax):
    """Pesos de diferencias finitas (Fornberg 1988) para las derivadas 0..mmax en z0."""
    n = len(xs)
    c = np.zeros((n, mmax + 1))
    c1, c4 = 1.0, xs[0] - z0
    c[0, 0] = 1.0
    for i in range(1, n):
        mn, c2, c5, c4 = min(i, mmax), 1.0, c4, xs[i] - z0
        for j in range(i):
            c3 = xs[i] - xs[j]
            c2 *= c3
            if j == i - 1:
                for k in range(mn, 0, -1):
                    c[i, k] = c1 * (k * c[i - 1, k - 1] - c5 * c[i - 1, k]) / c2
                c[i, 0] = -c1 * c5 * c[i - 1, 0] / c2
            for k in range(mn, 0, -1):
                c[j, k] = (c4 * c[j, k] - k * c[j, k - 1]) / c3
            c[j, 0] = c4 * c[j, 0] / c3
        c1 = c2
    return c


class Columna:
    """Una columna z de SPECTER: nz puntos extendidos, nz - C físicos, Lz = 1."""

    def __init__(self, nz, Lz=1.0):
        DIR, neu, dxp = tablas()
        self.nz, self.n, self.top = nz, nz - C, nz - C - 1
        self.dz = Lz / (nz - C - 1)
        self.z = self.dz * np.arange(nz)
        kz = np.fft.fftfreq(nz, 1.0 / nz) * 2 * np.pi / (self.dz * nz)
        self.nn = neu.copy()
        self.nn[D - 1] *= self.dz / dxp
        E = np.zeros((nz, self.n))
        E[:self.n, :] = np.eye(self.n)
        for ii in range(C):
            for jj in range(D):
                E[nz - C + ii, nz - C - D + jj] += DIR[ii, jj]
                E[nz - C + ii, D - 1 - jj] += DIR[C - 1 - ii, jj]
        F = np.fft.fft(np.eye(nz), axis=0)
        Fi = np.fft.ifft(np.eye(nz), axis=0)
        self.D2 = np.real(Fi @ np.diag(-kz**2) @ F @ E)[:self.n, :]
        # trazas en el tope, con el término de Nyquist omitido en las derivadas impares
        # (la convención de fs_gen_setup)
        fase = np.exp(1j * kz * self.z[self.top]) / nz
        self.tr = {}
        for m in range(4):
            w = (1j * kz)**m * fase
            if m % 2:
                w[nz // 2] = 0.0
            self.tr[m] = np.real(w @ F @ E)

    def beta(self):
        """Respuesta de la traza d^m u/dz^m a un dato de Neumann unitario (columna nula)."""
        u = np.zeros(self.n)
        u[self.top] = self.nn[D - 1]
        return [float(self.tr[m] @ u) for m in range(4)]

    def paso(self, eta, N, cierre, q=8):
        """Matriz del paso RK2 completo sobre los valores físicos de la columna."""
        n, top = self.n, self.top
        W = fornberg(self.z[top], self.z[top - np.arange(1, q + 1)], 3) if cierre == "interior" else None
        f3 = eta**2 / 2 if N >= 3 else 0.0

        def bc(u):
            u = u.copy()
            u[0] = 0.0
            a, b = self.nn[:D - 1] @ u[top - D + 1:top], self.nn[D - 1]
            if cierre == "interior":
                ui = u[top - np.arange(1, q + 1)]
                g = -eta * (W[:, 2] @ ui) - f3 * (W[:, 3] @ ui)
                u[top] = a + b * g
            else:
                ut = u.copy()
                ut[top] = 0.0
                c0, c1 = self.tr[2] @ ut, self.tr[2][top]
                e0, e1 = self.tr[3] @ ut, self.tr[3][top]
                u[top] = (a - b * (eta * c0 + f3 * e0)) / (1 + b * (eta * c1 + f3 * e1))
            return u

        M = np.zeros((n, n))
        for k in range(n):
            u0 = np.zeros(n)
            u0[k] = 1.0
            u = u0.copy()
            for o in (2, 1):
                u = bc(u0 + DT / o * NU * (self.D2 @ u))
            M[:, k] = u
        return M


def main():
    res = {"descripcion": __doc__.strip().splitlines()[0], "nu": NU, "dt": DT}
    # 1. beta_m, independientes de la resolución en unidades de dz^(m-1)
    res["beta"] = {}
    for nz in (64, 128, 256):
        c = Columna(nz)
        b = c.beta()
        res["beta"][str(nz)] = [b[m] * c.dz**(m - 1) for m in range(4)]
        print("nz %3d  beta_m dz^(m-1) = %s" % (nz, " ".join("%.4f" % x for x in res["beta"][str(nz)])))
    col = Columna(128)
    b2, b3 = res["beta"]["128"][2], res["beta"]["128"][3]
    res["singular_N2"] = -1.0 / b2
    res["min_1_mas_sigma_N3"] = 1.0 - b2**2 / (2.0 * b3)
    print("N=2: 1 + sigma = 0 en eta/dz = %.3f;  N=3: min(1 + sigma) = %.3f" % (
        res["singular_N2"], res["min_1_mas_sigma_N3"]))
    # 2. radio espectral contra r = eta/dz
    rs = np.round(np.concatenate([np.linspace(-4, -1, 13), np.linspace(-0.9, 0.9, 19),
                                  np.linspace(1, 4, 13)]), 4)
    casos = [("espectral", 0), ("interior", 8), ("interior", 5), ("interior", 3)]
    res["r"] = rs.tolist()
    res["radio"] = {}
    for N in (2, 3):
        for cierre, q in casos:
            clave = "N%d_%s%s" % (N, cierre, "" if cierre == "espectral" else "_q%d" % q)
            rho = [float(max(abs(np.linalg.eigvals(col.paso(r * col.dz, N, cierre, q)))))
                   for r in rs]
            res["radio"][clave] = rho
            inest = [r for r, x in zip(rs, rho) if x > 1 + 1e-9]
            print("%-22s inestable en eta/dz = %s" % (
                clave, "ninguno" if not inest else "%.2f..%.2f (%d de %d)" % (
                    min(inest), max(inest), len(inest), len(rs))))
    # 3. el modo espurio de N = 2: tasa medida contra nu/eta^2 del continuo
    res["modo_espurio_N2"] = []
    for r in (-1.0, -2.0, -3.0, -4.0):
        rho = res["radio"]["N2_espectral"][list(rs).index(r)]
        eta = r * col.dz
        res["modo_espurio_N2"].append({"r": r, "rho": rho, "continuo": math.exp(NU * DT / eta**2)})
        print("N=2, eta/dz = %+.0f: rho = %.5f, exp(nu dt/eta^2) = %.5f" % (
            r, rho, math.exp(NU * DT / eta**2)))
    # 4. el umbral de N = 2 depende de D = nu dt/dz^2, y N = 3 cerca del límite explícito
    res["umbral_vs_D"] = umbrales()
    os.makedirs(os.path.dirname(SALIDA), exist_ok=True)
    json.dump(res, open(SALIDA, "w"), indent=1)
    figura(res)
    print("->", SALIDA)


def _rho(col, r, N, dt):
    global DT
    viejo, DT = DT, dt
    try:
        return float(max(abs(np.linalg.eigvals(col.paso(r * col.dz, N, "espectral")))))
    finally:
        DT = viejo


def umbrales():
    """Bisecciones en la réplica ([H-24], [V-22]): el límite explícito D_max (eta = 0), el
    umbral de N = 2 en eta/dz para varios D, y el máximo de rho a N = 3 cerca de D_max."""
    out = []
    for nz in (128, 256):
        col = Columna(nz)
        dz2 = col.dz**2
        # límite explícito en D, con eta = 0
        lo, hi = 0.05, 0.5
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            if _rho(col, 0.0, 2, mid * dz2 / NU) > 1 + 1e-9:
                hi = mid
            else:
                lo = mid
        dmax = lo
        fila = {"nz": nz, "D_max": dmax, "N2": [], "N3": []}
        print("nz %d: límite explícito D_max = %.4f" % (nz, dmax), flush=True)
        for frac in (0.1, 0.5, 0.65, 0.85):
            D = frac * dmax
            dt = D * dz2 / NU
            lo, hi = -0.7159, -0.3
            for _ in range(40):
                mid = 0.5 * (lo + hi)
                if _rho(col, mid, 2, dt) > 1 + 1e-9:
                    lo = mid
                else:
                    hi = mid
            r0 = 0.5 * (lo + hi)
            b2 = col.beta()[2] * col.dz
            fila["N2"].append({"D": D, "umbral": r0, "uno_mas_sigma": 1 + b2 * r0,
                               "guardia": 0.05 + D, "r_guardia": -(0.95 - D) / b2})
            print("  N=2, D = %.4f: umbral eta/dz = %.4f (1+sigma = %.3f); la guardia corta en %.4f" % (
                D, r0, 1 + b2 * r0, -(0.95 - D) / b2), flush=True)
        for frac in (0.90, 0.95, 0.97, 0.99, 0.995):
            dt = frac * dmax * dz2 / NU
            rs = np.linspace(-2.0, 0.0, 41)
            rho = [_rho(col, r, 3, dt) for r in rs]
            k = int(np.argmax(rho))
            fila["N3"].append({"D": frac * dmax, "fraccion": frac, "rho_max": rho[k], "r": float(rs[k])})
            print("  N=3, D = %.4f (%.3f D_max): max rho = %.6f en eta/dz = %.2f" % (
                frac * dmax, frac, rho[k], rs[k]), flush=True)
        out.append(fila)
    return out


def figura(res):
    import matplotlib
    matplotlib.use("Agg")
    sys.path.insert(0, os.path.expanduser("~/Desktop/PhD/Num/Codigos"))
    try:
        import estetica  # noqa: F401  (rcParams compartidos, ver la memoria de estética)
        import matplotlib.pyplot as plt
        plt.rcParams["font.serif"] = ["Times New Roman", "Nimbus Roman", "Liberation Serif"]
        plt.rcParams["font.family"] = "serif"
        ancho = estetica.width_in
    except ImportError:
        import matplotlib.pyplot as plt
        ancho = 7.0
    r = np.array(res["r"])
    fig, ax = plt.subplots(1, 2, figsize=(ancho, 0.42 * ancho), sharey=True)
    estilos = {"espectral": ("k", "-", "traza espectral, implícita"),
               "interior_q8": ("C3", "--", "extrapolación interior, q = 8"),
               "interior_q5": ("C1", "-.", "extrapolación interior, q = 5"),
               "interior_q3": ("C2", ":", "extrapolación interior, q = 3")}
    for i, N in enumerate((2, 3)):
        a = ax[i]
        for clave, (col, ls, lab) in estilos.items():
            rho = np.array(res["radio"]["N%d_%s" % (N, clave)])
            a.semilogy(r, np.maximum(rho - 1, 1e-6), ls, color=col, lw=1.3, label=lab)
        if N == 2:
            a.axvline(res["singular_N2"], color="0.5", lw=0.8)
            a.text(res["singular_N2"] - 0.1, 3e2, r"$1+\sigma=0$", ha="right", fontsize=8,
                   color="0.3")
            cont = [(m["r"], m["continuo"] - 1) for m in res["modo_espurio_N2"]]
            h, = a.plot(*zip(*cont), "o", mfc="none", color="k", ms=5,
                        label=r"modo espurio del continuo: $e^{\nu\,dt/\eta^2}-1$")
            a.legend(handles=[h], fontsize=7, loc="lower left")
        a.set_title("N = %d" % N)
        a.set_xlabel(r"$\eta/\Delta z$")
        a.axhline(1e-6, color="0.8", lw=0.6)
        a.set_ylim(5e-7, 1e3)
    ax[0].set_ylabel(r"$\rho-1$ del paso RK2  (piso: estable)")
    ax[1].legend(fontsize=7, loc="upper left")
    fig.savefig(FIG + ".png", dpi=200)
    fig.savefig(FIG + ".pdf")
    print("->", FIG + ".png")


if __name__ == "__main__":
    main()
