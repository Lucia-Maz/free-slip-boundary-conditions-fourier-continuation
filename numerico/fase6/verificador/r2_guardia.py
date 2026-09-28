#!/usr/bin/env python
"""Round 2, item 1a (verifier): does the new guard sthr = 0.05 + nu dt/dz^2 (N = 2) cover every
unstable eta of the replica step, for every D = nu dt/dz^2 at which the eta = 0 step is stable?

Step: my round-1 construction (c3_replica.MiPaso: literal Fortran z-pipeline, bottom u = 0,
implicit spectral top closure, SPECTER RK2), here vectorized over the columns of the identity
and cross-checked against MiPaso.rho. The scaled problem depends on (D, r = eta/dz) only.
"""
import math
import sys
import time

import numpy as np

sys.dont_write_bytecode = True
import c3_replica as V

B2, B3 = 1.3967, 1.1515          # beta_2 dz, beta_3 dz^2 (nz = 64/128/256)


class Paso:
    def __init__(self, nz, D, nu=1e-2):
        self.mp = V.MiPaso(nz, nu, 1.0)
        self.dz = self.mp.F.dz
        self.nu = nu
        self.dt = D * self.dz**2 / nu
        self.D = D
        self.n = self.mp.n
        mp = self.mp
        self.t = mp.top
        self.neu = mp.F.neu
        self.T2, self.T3 = mp.T[2], mp.T[3]
        self.b2 = self.neu[-1] * self.T2[self.t] * self.dz
        self.b3 = self.neu[-1] * self.T3[self.t] * self.dz**2

    def bc(self, U, eta, N):
        U = U.copy()
        U[0, :] = 0.0
        t, d = self.t, V.d
        f3 = eta**2 / 2 if N >= 3 else 0.0
        a = self.neu[:d - 1] @ U[t - d + 1:t, :]
        b = self.neu[d - 1]
        Ut = U.copy()
        Ut[t, :] = 0.0
        c0, c1 = self.T2 @ Ut, self.T2[t]
        e0, e1 = self.T3 @ Ut, self.T3[t]
        U[t, :] = (a - b * (eta * c0 + f3 * e0)) / (1.0 + b * (eta * c1 + f3 * e1))
        return U

    def rho(self, r, N):
        eta = r * self.dz
        U0 = np.eye(self.n)
        U = U0
        for o in (2, 1):
            U = self.bc(U0 + self.dt / o * self.nu * (self.mp.D2 @ U), eta, N)
        return float(np.max(np.abs(np.linalg.eigvals(U))))


def dmax(nz):
    lo, hi = 0.15, 0.25
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if Paso(nz, mid).rho(0.0, 2) <= 1 + 1e-12:
            lo = mid
        else:
            hi = mid
    return lo


def unstable_set(P, N, rs, tol=1e-9):
    return [r for r in rs if P.rho(r, N) > 1 + tol]


def bisect_onset(P, lo, hi, N=2, it=30):
    """lo unstable, hi stable"""
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        if P.rho(mid, N) > 1 + 1e-9:
            lo = mid
        else:
            hi = mid
    return lo


def main():
    t0 = time.time()
    # cross-check with the round-1 step
    old = V.MiPaso(128, 1e-2, 1e-3)
    new = Paso(128, 1e-2 * 1e-3 / old.F.dz**2)
    for r in (-0.8, -0.667, -0.5, 1.0, 4.0):
        for N in (2, 3):
            a, b = old.rho(r, N), new.rho(r, N)
            assert abs(a - b) < 1e-10 * max(a, 1), (r, N, a, b)
    print("vectorized step == round-1 step (10 points, N=2,3)   [%.0f s]" % (time.time() - t0))
    for nz in (64, 128, 256):
        print("nz=%d: largest D stable at eta=0: D_max = %.6f  (2/pi^2 = %.6f)" % (nz, dmax(nz), 2 / math.pi**2))
    Dm = dmax(128)
    rs_neg = np.round(np.arange(-0.715, -0.0, 0.005), 4)
    rs_pos = np.round(np.concatenate([np.arange(0.05, 1.0, 0.05), np.arange(1.0, 12.01, 0.5)]), 4)
    rs3 = np.round(np.concatenate([np.arange(-12, -2.0, 0.5), np.arange(-2.0, 2.0, 0.05), np.arange(2.0, 12.01, 0.5)]), 4)
    print("\n   D      D/Dmax  sthr    r_guard   N=2 unstable r (neg side)    onset 1+sig   margin   N=2 pos side   N=3 unstable")
    for frac in (0.05, 0.1, 0.25, 0.5, 0.7, 0.85, 0.9, 0.95, 0.98, 0.99, 0.995, 0.999):
        D = frac * Dm
        P = Paso(128, D)
        sthr = 0.05 + D
        r_guard = (sthr - 1.0) / P.b2          # guard fires for r <= r_guard
        un = unstable_set(P, 2, rs_neg)
        if un:
            top = max(un)
            # refine the upper edge of the unstable set
            hi = min([r for r in rs_neg if r > top] or [0.0])
            edge = bisect_onset(P, top, hi)
        else:
            edge = float("nan")
        cont = (not un) or (sorted(un) == [r for r in rs_neg if r <= max(un)])
        pos = unstable_set(P, 2, rs_pos)
        un3 = unstable_set(P, 3, rs3)
        s_on = 1 + P.b2 * edge
        print(" %.5f  %.3f  %.4f  %+.4f   %s%s   %.4f      %+.4f   %s   %s" % (
            D, frac, sthr, r_guard,
            ("[%.3f .. %.4f]" % (min(un), edge)) if un else "none", "" if cont else " (NOT contiguous!)",
            s_on, sthr - s_on, (str(pos) if pos else "none"), ("%d pts, %.2f..%.2f" % (len(un3), min(un3), max(un3))) if un3 else "none"))
        sys.stdout.flush()
    # spot checks at other nz near the explicit limit
    print()
    for nz in (64, 256):
        Dmn = dmax(nz)
        for frac in (0.5, 0.95, 0.99):
            P = Paso(nz, frac * Dmn)
            un = unstable_set(P, 2, rs_neg)
            edge = bisect_onset(P, max(un), min([r for r in rs_neg if r > max(un)] or [0.0])) if un else float("nan")
            un3 = unstable_set(P, 3, rs3)
            pos = unstable_set(P, 2, rs_pos)
            print(" nz=%d D=%.5f (%.2f Dmax): N=2 onset r=%.4f, 1+sig=%.4f, sthr=%.4f, margin %+.4f ; N=2 pos unstable: %s ; N=3 unstable: %s" % (
                nz, frac * Dmn, frac, edge, 1 + P.b2 * edge, 0.05 + frac * Dmn, 0.05 + frac * Dmn - (1 + P.b2 * edge),
                pos if pos else "none", ("%d pts" % len(un3)) if un3 else "none"))
    print("[%.0f s]" % (time.time() - t0))


if __name__ == "__main__":
    main()
