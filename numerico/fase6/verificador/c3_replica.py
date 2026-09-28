#!/usr/bin/env python
"""C3 (verifier): fidelity of numerico/fase6/estabilidad_cierre.py and its stability map.

Part A: a literal transliteration (1-based loops) of the Fortran routines
  load_dirichlet_tables / load_neumann_tables (fcgram_mod.f90), neumann_reconstruct boun=6 ord=1,
  the continuation loop of fftp1d_real_to_complex_z (fftp.fpp), the kz grid of specter.fpp,
  fs_ph and the probe of fs_gen_setup (fsorder.f90),
compared number by number with the replica's Columna (nn, E, tr, beta()).
Part B: my own RK2 step (built from the Part-A operators, not from the replica's paso) for the
  implicit spectral closure, compared with the replica at a few r; then break-tests (fine r scan
  near the N=2 threshold, other dt, nz, nu, perturbed Neumann weight).
The replica module is imported read-only; its outputs are NOT written (main() is not called here).
"""
import math
import os
import sys

import numpy as np

sys.dont_write_bytecode = True
REPO = "/home/lucia/GW_AI/proyecto-final-piv"
sys.path.insert(0, os.path.join(REPO, "numerico", "fase6"))
import estabilidad_cierre as EC  # noqa: E402

TAB = os.path.join(REPO, "numerico", "SPECTER-trabajo", "tables")
C, d = 25, 5


def fortran_tables(dx):
    # READ(10) A with A(C,d): column-major
    A = np.fromfile(os.path.join(TAB, "A%d-%d.dat" % (C, d)), dtype="<f8").reshape((C, d), order="F")
    Qf = np.fromfile(os.path.join(TAB, "Q%d.dat" % d), dtype="<f8").reshape((d, d), order="F")
    raw = np.fromfile(os.path.join(TAB, "Q1n%d.dat" % d), dtype="<f8")
    dxp, Qn = raw[0], raw[1:1 + d * d].reshape((d, d), order="F")
    assert raw.size == 1 + d * d
    # load_dirichlet_tables: Q = TRANSPOSE(Q); dir = MATMUL(A, Q)
    dirm = A @ Qf.T
    # load_neumann_tables: neu = MATMUL(Q(d,:), TRANSPOSE(Qn)); neu(d) *= dx/dxp   (Q NOT transposed here)
    neu = Qf[d - 1, :] @ Qn.T
    neu = neu.copy()
    neu[d - 1] *= dx / dxp
    return dirm, neu, dxp


class Fortran:
    """1-based transliteration; arrays are numpy with index shift f[k-1]."""

    def __init__(self, nz, Lz=1.0):
        self.nz, self.Cz = nz, C
        rmp = Lz / (nz - C - 1)
        self.dz = rmp
        self.z = np.array([rmp * (k - 1) for k in range(1, nz + 1)])
        Dkz = 2.0 * math.pi / (rmp * nz)
        kz = np.zeros(nz)
        for k in range(1, nz // 2 + 1):
            kz[k - 1] = (k - 1)
            kz[k + nz // 2 - 1] = (k - nz // 2 - 1)
        self.kz = kz * Dkz
        self.dir, self.neu, _ = fortran_tables(self.z[1] - self.z[0])
        top = nz - C
        self.top = top
        ztop = self.z[top - 1]
        self.ph = np.zeros((nz, 4), dtype=complex)
        for k in range(1, nz + 1):
            for m in range(4):
                self.ph[k - 1, m] = (1j * self.kz[k - 1])**m * np.exp(1j * self.kz[k - 1] * ztop) / nz
                if k == nz // 2 + 1 and m % 2 == 1:
                    self.ph[k - 1, m] = 0.0

    def neumann_top(self, f):
        f = f.copy()
        nz, Cz, dd = self.nz, self.Cz, d
        i0 = nz - Cz
        acc = self.neu[dd - 1] * f[i0 - 1]
        for k in range(1, dd):
            acc = acc + self.neu[k - 1] * f[nz - Cz - dd + k - 1]
        f[i0 - 1] = acc
        return f

    def continue_(self, f):
        f = f.copy()
        n, Cc, dd = self.nz, self.Cz, d
        for ii in range(1, Cc + 1):
            f[n - Cc + ii - 1] = self.dir[ii - 1, 0] * f[n - Cc - dd + 1 - 1] + \
                self.dir[Cc - ii + 1 - 1, 0] * f[dd - 1]
        for jj in range(2, dd + 1):
            for ii in range(1, Cc + 1):
                f[n - Cc + ii - 1] = f[n - Cc + ii - 1] + self.dir[ii - 1, jj - 1] * f[n - Cc - dd + jj - 1] + \
                    self.dir[Cc - ii + 1 - 1, jj - 1] * f[dd - jj + 1 - 1]
        return f

    def fwd(self, f):          # fftp1d_real_to_complex_z: continuation + FFT (sign -1)
        return np.fft.fft(self.continue_(f))

    def trace(self, fhat, m):
        s = 0j
        for k in range(self.nz):
            s += fhat[k] * self.ph[k, m]
        return s

    def beta(self):
        pb = np.zeros(self.nz, dtype=complex)
        pb[self.top - 1] = 1.0
        pb = self.neumann_top(pb)
        pbh = self.fwd(pb)
        return [self.trace(pbh, m).real for m in range(4)]


def part_A():
    print("Part A: literal Fortran transliteration vs replica (estabilidad_cierre.Columna)")
    out = {}
    for nz in (64, 128, 256):
        Fo = Fortran(nz)
        col = EC.Columna(nz)
        n = nz - C
        # Neumann weights
        dnn = np.max(np.abs(Fo.neu - col.nn))
        # continuation matrix: apply my loop to unit vectors
        E = np.zeros((nz, n))
        for j in range(n):
            e = np.zeros(nz)
            e[j] = 1.0
            E[:, j] = Fo.continue_(e)
        Fh = np.fft.fft(E, axis=0)
        D2 = np.real(np.fft.ifft(-(Fo.kz**2)[:, None] * Fh, axis=0))[:n, :]
        dE = np.max(np.abs(D2 - col.D2)) / np.max(np.abs(col.D2))
        # trace functionals on the physical values
        dtr = max(np.max(np.abs(np.real(Fo.ph[:, m] @ Fh) - col.tr[m])) / np.max(np.abs(col.tr[m])) for m in range(4))
        bF = Fo.beta()
        bR = col.beta()
        sc = [Fo.dz**(m - 1) for m in range(4)]
        print("  nz=%3d: max|neu diff| %.1e, rel max|D2 diff| %.1e, rel max|trace-weight diff| %.1e" % (nz, dnn, dE, dtr))
        print("         beta_m dz^(m-1): Fortran-literal %s | replica %s" % (
            " ".join("%.4f" % (b * s) for b, s in zip(bF, sc)), " ".join("%.4f" % (b * s) for b, s in zip(bR, sc))))
        out[nz] = [b * s for b, s in zip(bF, sc)]
        if nz == 128:
            print("         neu weights (neu(1..d-1), neu(d)/dz):", np.round(Fo.neu[:-1], 5), round(Fo.neu[-1] / Fo.dz, 5))
            print("         max|dir| = %.3g, max row sum |dir| = %.3g" % (np.max(np.abs(Fo.dir)), np.max(np.sum(np.abs(Fo.dir), axis=1))))
    return out


class MiPaso:
    """My own RK2 step for the column (bottom u=0, top implicit spectral closure), from Part A ops."""

    def __init__(self, nz, nu, dt):
        self.F = Fortran(nz)
        self.nz, self.n = nz, nz - C
        self.nu, self.dt = nu, dt
        F = self.F
        n = self.n
        # matrices on physical values
        E = np.zeros((nz, n))
        for j in range(n):
            e = np.zeros(nz)
            e[j] = 1.0
            E[:, j] = F.continue_(e)
        Fh = np.fft.fft(E, axis=0)
        self.D2 = np.real(np.fft.ifft(-(F.kz**2)[:, None] * Fh, axis=0))[:n, :]
        self.T = [np.real(F.ph[:, m] @ Fh) for m in range(4)]
        self.top = F.top - 1   # 0-based

    def bc(self, u, eta, N):
        u = u.copy()
        u[0] = 0.0
        t = self.top
        f3 = eta**2 / 2 if N >= 3 else 0.0
        # u_top = a + b*g, g = -eta*T2(u) - f3*T3(u), T linear in u (incl. top node)
        a = sum(self.F.neu[k - 1] * u[t - d + k] for k in range(1, d))
        b = self.F.neu[d - 1]
        ut = u.copy()
        ut[t] = 0.0
        c0, c1 = self.T[2] @ ut, self.T[2][t]
        e0, e1 = self.T[3] @ ut, self.T[3][t]
        u[t] = (a - b * (eta * c0 + f3 * e0)) / (1.0 + b * (eta * c1 + f3 * e1))
        return u

    def matriz(self, eta, N):
        n = self.n
        M = np.zeros((n, n))
        for k in range(n):
            u0 = np.zeros(n)
            u0[k] = 1.0
            u = u0.copy()
            for o in (2, 1):
                u = self.bc(u0 + self.dt / o * self.nu * (self.D2 @ u), eta, N)
            M[:, k] = u
        return M

    def rho(self, r, N):
        return float(np.max(np.abs(np.linalg.eigvals(self.matriz(r * self.F.dz, N)))))


def part_B():
    print("Part B: my RK2 step vs the replica's paso (implicit spectral closure)")
    col = EC.Columna(128)
    mp = MiPaso(128, EC.NU, EC.DT)
    for N in (2, 3):
        for r in (-4.0, -2.0, -0.8, -0.7, -0.6, 0.5, 1.0, 4.0):
            a = mp.rho(r, N)
            b = float(np.max(np.abs(np.linalg.eigvals(col.paso(r * col.dz, N, "espectral")))))
            print("  N=%d r=%+.2f  rho mine %.8f  replica %.8f" % (N, r, a, b))
    # interior-extrapolation closure (replica only) at the quoted points of [H-23]
    for r in (1.0, 4.0):
        b = float(np.max(np.abs(np.linalg.eigvals(col.paso(r * col.dz, 2, "interior", 8)))))
        print("  replica q=8 N=2 r=%+.1f  rho %.3f" % (r, b))


def scan(mp, N, rs):
    return [(r, mp.rho(r, N)) for r in rs]


def part_C():
    print("Part C: break-tests with my step")
    fine = np.round(np.arange(-0.74, -0.39, 0.02), 3)
    for nz, nu, dt in ((128, 1e-2, 1e-3), (128, 1e-2, 1e-4), (128, 1e-2, 2e-3), (64, 1e-2, 1e-3),
                       (256, 1e-2, 1e-3), (128, 1e-3, 1e-3), (128, 1e-1, 1e-4)):
        mp = MiPaso(nz, nu, dt)
        rho0 = mp.rho(0.0, 2)
        res = scan(mp, 2, fine)
        unst = [r for r, x in res if x > 1 + 1e-9]
        print("  nz=%3d nu=%g dt=%g: rho(eta=0)=%.8f ; N=2 unstable r in fine scan: %s" % (
            nz, nu, dt, rho0, ("%.2f..%.2f" % (min(unst), max(unst))) if unst else "none"))
        print("      rho at r=-0.70,-0.66,-0.60,-0.50:", " ".join(
            "%.6f" % x for r, x in res if abs(r - (-0.70)) < 1e-9 or abs(r + 0.66) < 1e-9 or abs(r + 0.6) < 1e-9 or abs(r + 0.5) < 1e-9))
        wide = np.round(np.concatenate([np.arange(-8, -1.0, 0.5), np.arange(-1.0, 1.01, 0.1), np.arange(1.5, 9.01, 0.5)]), 3)
        res3 = scan(mp, 3, wide)
        un3 = [r for r, x in res3 if x > 1 + 1e-9]
        print("      N=3 over r in [-8, 9] (%d values): unstable %s ; max rho %.8f" % (
            len(wide), un3 if un3 else "none", max(x for _, x in res3)))
    # perturbed Neumann datum weight: threshold should move with beta_2
    for fac in (0.9, 1.1):
        mp = MiPaso(128, 1e-2, 1e-3)
        mp.F.neu = mp.F.neu.copy()
        mp.F.neu[d - 1] *= fac
        b2 = mp.F.neu[d - 1] * mp.T[2][mp.top] * mp.F.dz
        res = scan(mp, 2, np.round(np.arange(-0.90, -0.39, 0.02), 3))
        unst = [r for r, x in res if x > 1 + 1e-9]
        print("  neu(d) x %.1f: beta2*dz = %.4f, -1/beta2 = %.3f dz; N=2 unstable r: %s" % (
            fac, b2, -1 / b2, ("%.2f..%.2f" % (min(unst), max(unst))) if unst else "none"))


if __name__ == "__main__":
    part_A()
    part_B()
    part_C()
