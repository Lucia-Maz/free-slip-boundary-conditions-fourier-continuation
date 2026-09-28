#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Referencia de geometría exacta para la puerta de la Fase 6: ondas viscosas no lineales en
2D con la superficie libre donde realmente está.

Es independiente de SPECTER y no usa ninguna fórmula de transferencia: el dominio
0 < z < H(x,t) = h + eta(x,t) se mapea a 0 < zeta = z/H < 1, y las condiciones de la
superficie se imponen EXACTAS en zeta = 1 (Wang, Tice & Kim, ARMA 212 (2014),
DOI 10.1007/s00205-013-0700-2, §1.1 ecs. (1.4)-(1.8), leídas en el original).

Formulación (derivación de esta sesión, [D-51]):
  - función corriente psi(x,z,t): u = psi_z, w = -psi_x (incompresible por construcción);
    vorticidad omega = w_x - u_z = -lap(psi);
  - volumen: omega_t + u omega_x + w omega_z = nu lap(omega), con omega_t|_z = -lap(psi_t|_z)
    y psi_t|_z = Psi_t - zeta (eta_t / H) Psi_zeta (derivada temporal a z fijo);
  - fondo zeta = 0: psi = 0, psi_z = 0 (no deslizante);
  - superficie zeta = 1, con q = P/rho + g z la presión dinámica total:
      tensión tangencial   (1 - eta_x^2) S_xz + eta_x (S_zz - S_xx) = 0
      tensión normal, derivada a lo largo de la superficie y combinada con el momento
      tangencial (la presión se elimina):
        [Du/Dt - nu lap u] + eta_x [Dw/Dt - nu lap w] + d/dx[g eta - gam kappa + nu B] = 0
      con kappa = d/dx(eta_x / sqrt(1+eta_x^2)) y B = n.S.n;
      cinemática  eta_t = w - u eta_x = -d/dx[psi(x, H(x))].
Discretización: Fourier en x, colocación de Chebyshev (Gauss-Lobatto) en zeta con
reemplazo de filas en los bordes, e integración en el tiempo con Radau IIA de 3 etapas
(orden 5, rígidamente precisa, apta para DAE de índice 1), Newton simplificado con
jacobiano por diferencias finitas.

Autoverificación (en la puerta, chequeo H6): el balance de energía exacto
  d/dt [ int 1/2 |v|^2 dA + g/2 int eta^2 + gam int (sqrt(1+eta_x^2) - 1) ] = -(nu/2) int S:S dA
cierra, el límite lineal coincide con la referencia MAC de la Fase 5, y converge en Nx,
Nz y dt.
"""

import numpy as np
from scipy.linalg import lu_factor, lu_solve


def cheb(n):
    """Matriz de derivada de Chebyshev (Trefethen) en x_j = cos(pi j/n), j = 0..n."""
    x = np.cos(np.pi * np.arange(n + 1) / n)
    c = np.ones(n + 1)
    c[0] = c[-1] = 2.0
    c *= (-1.0) ** np.arange(n + 1)
    X = np.tile(x, (n + 1, 1)).T
    D = np.outer(c, 1.0 / c) / (X - X.T + np.eye(n + 1))
    D -= np.diag(D.sum(axis=1))
    return x, D


def clenshaw_curtis(n):
    """Pesos de Clenshaw-Curtis en [-1, 1] para los nodos cos(pi j/n)."""
    th = np.pi * np.arange(n + 1) / n
    w = np.zeros(n + 1)
    v = np.ones(n - 1)
    if n % 2 == 0:
        w[0] = w[n] = 1.0 / (n * n - 1)
        for k in range(1, n // 2):
            v -= 2.0 * np.cos(2 * k * th[1:-1]) / (4 * k * k - 1)
        v -= np.cos(n * th[1:-1]) / (n * n - 1)
    else:
        w[0] = w[n] = 1.0 / n**2
        for k in range(1, (n - 1) // 2 + 1):
            v -= 2.0 * np.cos(2 * k * th[1:-1]) / (4 * k * k - 1)
    w[1:-1] = 2.0 * v / n
    return w


# Radau IIA, 3 etapas
_S6 = np.sqrt(6.0)
RADAU_A = np.array([[(88 - 7 * _S6) / 360, (296 - 169 * _S6) / 1800, (-2 + 3 * _S6) / 225],
                    [(296 + 169 * _S6) / 1800, (88 + 7 * _S6) / 360, (-2 - 3 * _S6) / 225],
                    [(16 - _S6) / 36, (16 + _S6) / 36, 1.0 / 9]])


class ReferenciaExacta:
    """Ondas viscosas en 2D con superficie libre exacta. Estado y = (Psi, eta)."""

    def __init__(self, Nx=16, Nz=24, h=1.0, nu=0.01, g=1.0, gam=0.0,
                 signo_pendiente_tangencial=1.0):
        self.Nx, self.Nz, self.h, self.nu, self.g, self.gam = Nx, Nz, h, nu, g, gam
        # control negativo de la puerta: invierte el término de pendiente de la tensión
        # tangencial (+1 es la física; -1 es un modelo deliberadamente equivocado)
        self.sgn_t = signo_pendiente_tangencial
        self.x = 2 * np.pi * np.arange(Nx) / Nx
        xc, Dc = cheb(Nz)
        self.zeta = (1.0 + xc) / 2.0                    # zeta_0 = 1 (superficie)
        self.Dz = 2.0 * Dc                               # d/dzeta
        self.wz = clenshaw_curtis(Nz) / 2.0              # pesos en [0, 1]
        k = np.fft.fftfreq(Nx, 1.0 / Nx)
        k[Nx // 2] = 0.0                                 # sin Nyquist
        self.ik = 1j * k
        self.nps = Nx * (Nz + 1)
        self.n = self.nps + Nx

    # --- operadores -----------------------------------------------------------------
    def dxs(self, F):
        """d/dx a zeta fijo (espectral, eje 0)."""
        return np.real(np.fft.ifft(self.ik[:, None] * np.fft.fft(F, axis=0), axis=0)) \
            if F.ndim == 2 else np.real(np.fft.ifft(self.ik * np.fft.fft(F)))

    def geom(self, eta):
        H = self.h + eta
        Hx = self.dxs(eta)
        return H, Hx

    def partir(self, y):
        return y[:self.nps].reshape(self.Nx, self.Nz + 1), y[self.nps:]

    def residuo(self, y, yp):
        nu, g, gam = self.nu, self.g, self.gam
        Psi, eta = self.partir(y)
        Pt, et = self.partir(yp)
        H, Hx = self.geom(eta)
        Hc, Hxc = H[:, None], Hx[:, None]
        zc = self.zeta[None, :]
        Dz = self.Dz

        def dzeta(F):
            return F @ Dz.T

        def dz(F):
            return dzeta(F) / Hc

        def dx(F):
            return self.dxs(F) - zc * (Hxc / Hc) * dzeta(F)

        def lap(F):
            return dx(dx(F)) + dz(dz(F))

        u = dz(Psi)
        w = -dx(Psi)
        om = dx(w) - dz(u)
        phi = Pt - zc * (et[:, None] / Hc) * dzeta(Psi)          # psi_t a z fijo
        R = np.zeros_like(Psi)
        vort = -lap(phi) + u * dx(om) + w * dz(om) - nu * lap(om)
        R[:, 2:self.Nz - 1] = vort[:, 2:self.Nz - 1]
        # superficie (zeta = 1, columna 0)
        ex = Hx
        Sxz = dz(u) + dx(w)
        Sxx = 2 * dx(u)
        Szz = 2 * dz(w)
        R[:, 0] = (1 - ex**2) * Sxz[:, 0] + self.sgn_t * ex * (Szz[:, 0] - Sxx[:, 0])
        ut = dz(phi)
        wt = -dx(phi)
        Mu = ut + u * dx(u) + w * dz(u) - nu * lap(u)
        Mw = wt + u * dx(w) + w * dz(w) - nu * lap(w)
        kappa = self.dxs(ex / np.sqrt(1 + ex**2))
        B = (Szz[:, 0] - 2 * ex * Sxz[:, 0] + ex**2 * Sxx[:, 0]) / (1 + ex**2)
        R[:, 1] = Mu[:, 0] + ex * Mw[:, 0] + self.dxs(g * eta - gam * kappa + nu * B)
        # fondo (zeta = 0, columna Nz): psi = 0, psi_zeta = 0
        R[:, self.Nz] = Psi[:, self.Nz]
        R[:, self.Nz - 1] = dzeta(Psi)[:, self.Nz]
        Reta = et + self.dxs(Psi[:, 0])
        return np.concatenate([R.ravel(), Reta])

    # --- diagnósticos ---------------------------------------------------------------
    def energia(self, y):
        """(E, D): energía total y disipación (nu/2) int S:S dA."""
        Psi, eta = self.partir(y)
        H, Hx = self.geom(eta)
        Hc, Hxc = H[:, None], Hx[:, None]
        zc = self.zeta[None, :]

        def dzeta(F):
            return F @ self.Dz.T

        def dz(F):
            return dzeta(F) / Hc

        def dx(F):
            return self.dxs(F) - zc * (Hxc / Hc) * dzeta(F)

        u, w = dz(Psi), -dx(Psi)
        dA = (2 * np.pi / self.Nx) * Hc * self.wz[None, :]
        KE = 0.5 * np.sum((u**2 + w**2) * dA)
        PE = 0.5 * self.g * np.sum(eta**2) * 2 * np.pi / self.Nx
        SE = self.gam * np.sum(np.sqrt(1 + Hx**2) - 1) * 2 * np.pi / self.Nx
        diss = self.nu * np.sum((2 * dx(u)**2 + 2 * dz(w)**2 + (dz(u) + dx(w))**2) * dA)
        return KE + PE + SE, diss

    # --- integración ----------------------------------------------------------------
    def jacobianos(self, y, yp):
        n = self.n
        F0 = self.residuo(y, yp)
        Jyp = np.empty((n, n))
        Jy = np.empty((n, n))
        # F es lineal en yp: la columna es exacta (salvo redondeo) con paso 1
        for k in range(n):
            e = np.zeros(n)
            e[k] = 1.0
            Jyp[:, k] = self.residuo(y, yp + e) - F0
            dk = 1e-7 * max(1.0, abs(y[k]))
            yk = y.copy()
            yk[k] += dk
            Jy[:, k] = (self.residuo(yk, yp) - F0) / dk
        return Jy, Jyp

    def _newton(self, y, K, paso, lu, tol, maxit, salida):
        """Newton simplificado sobre las derivadas de etapa K (se modifica en el lugar).
        Converge si |dK| <= tol |K|, o si se estancó en el piso de redondeo del sistema
        lineal (|dK| deja de bajar y ya es <= 1e-9 |K|): el piso depende del
        condicionamiento de la colocación, que crece con Nz."""
        n = self.n
        previo = np.inf
        for it in range(maxit):
            Kr = K.reshape(3, n)
            Y = y[None, :] + paso * (RADAU_A @ Kr)
            G = np.concatenate([self.residuo(Y[i], Kr[i]) for i in range(3)])
            dK = lu_solve(lu, -G)
            K += dK
            salida["newton"] += 1
            escala = max(1.0, np.max(np.abs(K)))
            ndk = np.max(np.abs(dK))
            if ndk <= tol * escala:
                return True
            if it >= 2 and ndk >= 0.3 * previo and ndk <= 1e-9 * escala:
                return True
            previo = ndk
        return False

    def integrar(self, eps, tiempos, modo=1, dt_max=0.02, dt0=1e-4, tol=1e-11,
                 eta0=None, energia=False):
        """Integra desde el reposo con eta = eps cos(modo x) (o eta0) y devuelve el
        coeficiente de Fourier normalizado de eta en los modos 1 y 2 (amplitud de
        coseno: 2 Re eta_k) en cada tiempo pedido, y opcionalmente E y D."""
        n = self.n
        y = np.zeros(n)
        y[self.nps:] = eps * np.cos(modo * self.x) if eta0 is None else eta0
        yp = np.zeros(n)
        t = 0.0
        dt = dt0
        tiempos = list(tiempos)
        salida = {"t": [], "eta1": [], "eta2": [], "eta3": [], "media": [], "E": [],
                  "D": [], "newton": 0, "factorizaciones": 0}

        def registrar():
            ehat = np.fft.fft(y[self.nps:]) / self.Nx
            salida["t"].append(t)
            salida["eta1"].append(2 * ehat[modo])
            salida["eta2"].append(2 * ehat[2 * modo])
            salida["eta3"].append(2 * ehat[3 * modo])
            salida["media"].append(ehat[0].real)
            if energia:
                E, D = self.energia(y)
                salida["E"].append(E)
                salida["D"].append(D)

        if tiempos and abs(tiempos[0]) < 1e-14:
            registrar()
            tiempos.pop(0)
        lu = None
        dt_lu = None
        while tiempos:
            paso = min(dt, tiempos[0] - t)
            if lu is None or abs(paso - dt_lu) > 1e-15 * max(1.0, paso):
                Jy, Jyp = self.jacobianos(y, yp)
                big = np.kron(np.eye(3), Jyp) + paso * np.kron(RADAU_A, Jy)
                lu = lu_factor(big)
                dt_lu = paso
                salida["factorizaciones"] += 1
            K = np.tile(yp, 3)
            ok = self._newton(y, K, paso, lu, tol, 25, salida)
            if not ok:
                # jacobiano viejo: rehacer en el estado actual y repetir el paso
                Jy, Jyp = self.jacobianos(y, yp)
                big = np.kron(np.eye(3), Jyp) + paso * np.kron(RADAU_A, Jy)
                lu = lu_factor(big)
                dt_lu = paso
                salida["factorizaciones"] += 1
                K = np.tile(yp, 3)
                ok = self._newton(y, K, paso, lu, tol, 40, salida)
                if not ok:
                    raise RuntimeError("Newton no converge en t = %g" % t)
            Kr = K.reshape(3, n)
            y = y + paso * (RADAU_A[2] @ Kr)
            yp = Kr[2].copy()
            t += paso
            if abs(t - tiempos[0]) < 1e-12:
                t = tiempos[0]
                registrar()
                tiempos.pop(0)
            dt = min(2 * dt, dt_max)
        for k in ("t", "E", "D", "eta1", "eta2", "eta3", "media"):
            salida[k] = np.array(salida[k])
        return salida
