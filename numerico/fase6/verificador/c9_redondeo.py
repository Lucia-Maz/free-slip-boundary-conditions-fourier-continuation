#!/usr/bin/env python
"""C9 (verifier): roundoff floor of the m = 3 spectral trace, measured in the Fortran-literal
z-pipeline (neumann_reconstruct top -> FC continuation loop -> FFT -> sum_k fhat_k fs_ph(k,m)).

The column is a wave-like tangential-velocity profile, calibrated to the author's numbers
(max|fhat| = 4.6e-5 in the (kz,ky,kx) normalization). The Neumann datum is perturbed by
+-h in random directions and the second difference T(g+h) - 2T(g) + T(g-h) (pure roundoff:
T is linear) is recorded. Then: which operation carries the noise (continuation in extended
precision vs float64), and what that gives for Dx through eta^2/2 and eta.
"""
import math
import sys
import numpy as np

sys.dont_write_bytecode = True
import c3_replica as V

rng = np.random.default_rng(7)
nz = 128
Fo = V.Fortran(nz)
n = nz - V.C
dz = Fo.dz
z = Fo.z
top = Fo.top - 1
k = math.sqrt(2.0)
eps = np.finfo(float).eps


def column(mismatch):
    zz = z[:n]
    prof = np.cosh(k * zz) / np.cosh(k) - np.exp(-zz / 0.1) * np.cos(zz / 0.1) / np.cosh(k)
    f = np.zeros(nz, dtype=complex)
    f[:n] = prof * (1.0 + 0.3j)
    g = (k * np.tanh(k) + mismatch) * (1.0 + 0.3j)    # Neumann datum (outward d/dz at the top)
    return f, g


def continue_ld(f):
    """continuation loop in extended precision, result rounded to float64 complex."""
    fr = f.real.astype(np.longdouble); fi = f.imag.astype(np.longdouble)
    dirl = Fo.dir.astype(np.longdouble)
    out_r, out_i = fr.copy(), fi.copy()
    Cc, dd = V.C, V.d
    for ii in range(1, Cc + 1):
        sr = np.longdouble(0); si = np.longdouble(0)
        for jj in range(1, dd + 1):
            sr += dirl[ii - 1, jj - 1] * fr[nz - Cc - dd + jj - 1] + dirl[Cc - ii, jj - 1] * fr[dd - jj]
            si += dirl[ii - 1, jj - 1] * fi[nz - Cc - dd + jj - 1] + dirl[Cc - ii, jj - 1] * fi[dd - jj]
        out_r[nz - Cc + ii - 1] = sr
        out_i[nz - Cc + ii - 1] = si
    return (out_r.astype(float) + 1j * out_i.astype(float))


def pipeline(f, g, scale, extended_cont=False):
    u = f.copy() * scale
    u[top] = g * scale
    u = Fo.neumann_top(u)
    uc = continue_ld(u) if extended_cont else Fo.continue_(u)
    uh = np.fft.fft(uc)
    return uh, [Fo.trace(uh, m) for m in range(4)]


def noise(f, g, scale, h_rel, extended_cont=False, trials=40):
    base = pipeline(f, g, scale, extended_cont)[1]
    out = {2: [], 3: []}
    for _ in range(trials):
        dg = h_rel * abs(g) * np.exp(2j * np.pi * rng.random())
        tp = pipeline(f, g + dg, scale, extended_cont)[1]
        tm = pipeline(f, g - dg, scale, extended_cont)[1]
        for m in (2, 3):
            out[m].append(abs(tp[m] - 2 * base[m] + tm[m]))
    return {m: float(np.sqrt(np.mean(np.square(out[m])))) for m in (2, 3)}


def main():
    kmax = math.pi / dz
    for mismatch in (0.0, 1.0, 5.0):
        f, g = column(mismatch)
        uh, _ = pipeline(f, g, 1.0)
        scale = 4.6e-5 / np.max(np.abs(uh))
        uh, tr = pipeline(f, g, scale)
        sabs = float(np.sum(np.abs(uh * Fo.ph[:, 3])))
        fb = np.max(np.abs(f[:n] * scale))
        print("column (datum mismatch %.1f): max|fhat| = %.2e, max|f| (z,ky,kx) = %.2e, "
              "sum_k|fhat_k ph3| = %.3g, |trace3| = %.3g" % (mismatch, np.max(np.abs(uh)), fb, sabs, abs(tr[3])))
        for h in (1e-9, 1e-6, 1e-3):
            a = noise(f, g, scale, h)
            b = noise(f, g, scale, h, extended_cont=True)
            print("   h=%.0e|g|: rms 2nd diff  trace3 %.2e  trace2 %.2e  | continuation in extended precision: "
                  "trace3 %.2e trace2 %.2e" % (h, a[3], a[2], b[3], b[2]))
    print("reference numbers (nz=128, dz=%.5f, kmax=pi/dz=%.1f, kmax^3=%.3g):" % (dz, kmax, kmax**3))
    A = 4.6e-5
    s6 = np.sqrt(np.sum(np.abs(Fo.ph[:, 3])**2))
    s1 = np.sum(np.abs(Fo.ph[:, 3]))
    print("   record's estimate eps*max|fhat|*kmax^3         = %.2e" % (eps * A * kmax**3))
    print("   with the 1/nz of fs_ph, white per-coefficient: eps*A*sqrt(sum|ph3|^2) = %.2e ; "
          "coherent bound eps*A*sum|ph3| = %.2e" % (eps * A * s6, eps * A * s1))
    print("   continuation conditioning: max row sum |dir| = %.3g" % np.max(np.sum(np.abs(Fo.dir), axis=1)))
    eta = 0.04
    print("   => Dx noise ~ (eta^2/2)*noise(trace3) + eta*noise(trace2) with eta = %.2f (eta^2/2 = %.1e)" % (eta, eta**2 / 2))


if __name__ == "__main__":
    main()
