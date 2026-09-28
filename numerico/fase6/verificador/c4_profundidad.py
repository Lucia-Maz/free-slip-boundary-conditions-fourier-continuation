#!/usr/bin/env python
"""C4 (verifier): depth-shift degeneracy.

(a) exact series (Fractions): with T = tanh(kh), x = k eta0,
    R_N(x) = (T + y_N)/(1 + T y_N), y_N = [sum_{odd j<N} x^j/j!]/[sum_{even j<N} x^j/j!]
    (this is S_N/C_N with S_N, C_N the order-(N-1) traces of sinh, cosh), versus
    exact tanh(kh + x) = (T + tanh x)/(1 + T tanh x).  Coefficients of x^1..x^5.
(b) viscous linear waves, finite depth, no-slip bottom, 5x5 determinant with modes
    e^{k(z-h)}, e^{-kz} (potential), e^{m(z-h)}, e^{-mz} (rotational), m^2 = k^2 + lam/nu.
    The order-N transferred top conditions = each mode's top entries times P_{N-1}(mu eta0);
    the exact problem at depth h + eta0 = the same with e^{mu eta0}.
"""
from fractions import Fraction as Fr
import math
import numpy as np

D = 7


def smul(a, b):
    r = [Fr(0)] * (D + 1)
    for i, x in enumerate(a):
        for j in range(D + 1 - i):
            r[i + j] += x * b[j]
    return r


def sinv(a):   # 1/a, a[0] != 0
    r = [Fr(0)] * (D + 1)
    r[0] = 1 / a[0]
    for n in range(1, D + 1):
        r[n] = -sum(a[k] * r[n - k] for k in range(1, n + 1)) / a[0]
    return r


def tanh_series():
    sh = [Fr(1, math.factorial(j)) if j % 2 else Fr(0) for j in range(D + 1)]
    ch = [Fr(0) if j % 2 else Fr(1, math.factorial(j)) for j in range(D + 1)]
    return smul(sh, sinv(ch))


def R_of(y, T):
    num = [T + (y[0] if True else 0)] + y[1:]
    num = [T * (1 if j == 0 else 0) + y[j] for j in range(D + 1)]
    den = [(1 if j == 0 else 0) + T * y[j] for j in range(D + 1)]
    return smul(num, sinv(den))


def parte_a():
    th = tanh_series()
    print("(a) tanh x series:", [str(c) for c in th[:6]])
    for T in (Fr(3, 7), Fr(8884, 10000), Fr(-2, 5)):
        ex = R_of(th, T)
        for N in (1, 2, 3, 4):
            odd = [Fr(1, math.factorial(j)) if (j % 2 == 1 and j <= N - 1) else Fr(0) for j in range(D + 1)]
            even = [Fr(1, math.factorial(j)) if (j % 2 == 0 and j <= N - 1) else Fr(0) for j in range(D + 1)]
            yN = smul(odd, sinv(even))
            RN = R_of(yN, T)
            diff = [RN[j] - ex[j] for j in range(D + 1)]
            first = next((j for j in range(D + 1) if diff[j] != 0), None)
            print("   T=%s N=%d: transferred - exact, first nonzero x^%s coeff = %s ; (1-T^2) * %s" % (
                T, N, first, diff[first], diff[first] / (1 - T * T)))
    # closed form of the numerator of S2/C2 - tanh(kh+x): x cosh x - sinh x
    ch = [Fr(0) if j % 2 else Fr(1, math.factorial(j)) for j in range(D + 1)]
    sh = [Fr(1, math.factorial(j)) if j % 2 else Fr(0) for j in range(D + 1)]
    xch = [Fr(0)] + ch[:-1]
    print("   x cosh x - sinh x =", [str(a - b) for a, b in zip(xch, sh)][:7], "(so the x^3 coefficient is 1/3, not 1/2)")


NU, G, GAM, H = 0.01, 1.0, 0.2, 1.0


def Pm(M, x):
    return sum(x**j / math.factorial(j) for j in range(M + 1))


def mat(lam, k, eta0, N, exact=False, nu=NU):
    m = np.sqrt(k * k + lam / nu + 0j)
    if m.real < 0:
        m = -m
    ik = 1j * k
    # (type, exponent mu, value at z=h, value at z=0) of the normalized modes
    modes = [("p", k, 1.0, np.exp(-k * H)), ("p", -k, np.exp(-k * H), 1.0),
             ("r", m, 1.0, np.exp(-m * H)), ("r", -m, np.exp(-m * H), 1.0)]
    M = np.zeros((5, 5), dtype=complex)
    for c, (typ, mu, vh, vb) in enumerate(modes):
        fac = np.exp(mu * eta0) if exact else Pm(N - 1, mu * eta0)
        if typ == "p":
            M[0, c] = ik * vb; M[1, c] = mu * vb
            w, sxz, p, szz = mu * vh, 2 * ik * mu * vh, -lam * vh, 2 * mu * mu * vh
        else:
            M[0, c] = -mu * vb; M[1, c] = ik * vb
            w, sxz, p, szz = ik * vh, -(mu * mu + k * k) * vh, 0.0, 2 * ik * mu * vh
        M[2, c] = -fac * w
        M[3, c] = fac * sxz
        M[4, c] = fac * (p - nu * szz)
    M[2, 4] = lam
    M[4, 4] = -(G + GAM * k * k)
    return M


def root(k, eta0, N, exact=False, nu=NU, lam0=None):
    if lam0 is None:
        Hh = H + eta0
        w = math.sqrt((G * k + GAM * k**3) * math.tanh(k * Hh))
        lam0 = -1j * w - 2 * nu * k * k
    f = lambda l: np.linalg.det(mat(l, k, eta0, N, exact, nu))
    x0, x1 = lam0, lam0 * (1 + 1e-4)
    f0, f1 = f(x0), f(x1)
    for _ in range(100):
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        x0, f0, x1, f1 = x1, f1, x2, f(x2)
        if abs(x1 - x0) < 1e-15 * abs(x1):
            break
    return x1


def parte_b():
    k = math.sqrt(2.0)
    print("(b) viscous, nu=%g, g=%g, gam=%g, h=%g, k=sqrt2; lam = -i omega + ..." % (NU, G, GAM, H))
    l_ex0 = root(k, 0.0, 1, exact=True)
    print("   exact root at depth 1: lam = %.10f %+.10fi" % (l_ex0.real, l_ex0.imag))
    etas = [0.0025, 0.005, 0.01, 0.02, 0.04]
    for N in (1, 2, 3):
        errs, errs_inv = [], []
        for e in etas:
            le = root(k, e, N, exact=True)
            lN = root(k, e, N, exact=False, lam0=le)
            errs.append(abs(lN - le))
            # inviscid prediction: omega_N^2 = (gk+gam k^3) R_N(k eta0)
            T, x = math.tanh(k * H), k * e
            odd = sum(x**j / math.factorial(j) for j in range(1, N, 2))
            even = sum(x**j / math.factorial(j) for j in range(0, N, 2))
            y = odd / even
            RN = (T + y) / (1 + T * y)
            wN = math.sqrt((G * k + GAM * k**3) * RN)
            wE = math.sqrt((G * k + GAM * k**3) * math.tanh(k * (H + e)))
            errs_inv.append(abs(wN - wE))
        sl = [math.log(errs[i + 1] / errs[i]) / math.log(2) for i in range(len(etas) - 1)]
        sli = [math.log(errs_inv[i + 1] / errs_inv[i]) / math.log(2) for i in range(len(etas) - 1)]
        print("   N=%d |lam_N - lam_exact(h+eta0)|: %s ; slopes %s" % (N, " ".join("%.3e" % v for v in errs),
                                                                    " ".join("%.2f" % s for s in sl)))
        print("       inviscid |omega_N - omega_ex|: %s ; slopes %s" % (" ".join("%.3e" % v for v in errs_inv),
                                                                     " ".join("%.2f" % s for s in sli)))
    # viscosity dependence of the N=2 residual O(eta0^2) piece: vary nu
    for nu in (0.1, 0.01, 0.001):
        e1, e2 = 0.005, 0.01
        d = []
        for e in (e1, e2):
            le = root(k, e, 2, exact=True, nu=nu)
            lN = root(k, e, 2, exact=False, nu=nu, lam0=le)
            d.append(lN - le)
        m = np.sqrt(k * k + le / nu + 0j)
        print("   nu=%g: N=2 error at eta0=0.005, 0.01: %.3e, %.3e (slope %.2f); e^{-Re(m)h} = %.2e" % (
            nu, abs(d[0]), abs(d[1]), math.log(abs(d[1]) / abs(d[0])) / math.log(2), math.exp(-m.real * H)))


if __name__ == "__main__":
    parte_a()
    parte_b()
