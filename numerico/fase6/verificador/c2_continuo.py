#!/usr/bin/env python
"""C2 (verifier): ill-posedness of the order-2 transfer for eta < 0, in the continuum.

(i) scalar model u_t = nu u_zz on z < Lz with sum_{j<N} eta^j/j! d^{j+1}u/dz^{j+1} = 0 at Lz:
    modes e^{lam t + kap (z-Lz)}, Re kap > 0, lam = nu kap^2, kap * P_{N-1}(eta kap) = 0.
(ii) full linearized 2D problem in a deep layer (half-space below z = h), base at rest with a
    frozen elevation eta0 and frozen slope s0, perturbation e^{ikx + lam t}:
    potential mode phi = A e^{k(z-h)} (p = -lam phi) + rotational psi = B e^{m(z-h)},
    m^2 = k^2 + lam/nu, Re m > 0; unknowns (A, B, eta1); rows = the order-N kinematic,
    tangential (with its slope terms) and normal conditions linearized about (eta0, s0).
    Roots of det M(lam) = 0 near the predicted spurious modes (Newton), and the number of
    roots with Re lam > 0 inside a large rectangle (argument principle).
Derivation of the rows is mine (verifier), from the header formulas of fsorder.f90.
"""
import math
import numpy as np

NU, G, GAM = 0.01, 1.0, 0.2


def P(M, x):
    return sum(x**j / math.factorial(j) for j in range(M + 1)) if M >= 0 else 0.0 * x


def matriz(lam, k, eta0, s0, N):
    m = np.sqrt(k * k + lam / NU + 0j)
    if m.real < 0:
        m = -m
    Pk = lambda M: P(M, k * eta0)
    Pm = lambda M: P(M, m * eta0)
    ik = 1j * k
    # mode quantities at z = h (A: potential, B: rotational)
    u = (ik, -m); w = (k, ik)
    Sxz = (2j * k * k, -(m * m + k * k)); Sxx = (-2 * k * k, -2j * k * m); Szz = (2 * k * k, 2j * k * m)
    pr = (-lam, 0.0)
    fac = lambda M: (Pk(M), Pm(M))
    row_kin = [0j, 0j, lam]
    row_tan = [0j, 0j, 0j]
    row_nor = [0j, 0j, 0j]
    for c in range(2):
        F = lambda M: fac(M)[c]
        row_kin[c] = -F(N - 1) * w[c] + s0 * F(N - 2) * u[c]
        row_tan[c] = (F(N - 1) * Sxz[c] - s0 * F(N - 2) * Sxx[c] + s0 * F(N - 2) * Szz[c]
                      - s0 * s0 * F(N - 3) * Sxz[c])
        B = 0j
        j = 0
        while 2 * j <= N - 1:
            B += (-s0 * s0)**j * (F(N - 1 - 2 * j) * Szz[c] - 2 * s0 * F(N - 2 - 2 * j) * Sxz[c]
                                  + s0 * s0 * F(N - 3 - 2 * j) * Sxx[c])
            j += 1
        row_nor[c] = F(N - 1) * pr[c] - NU * B
    ckap = 1.0 if N <= 2 else (1.0 - 1.5 * s0 * s0)
    row_nor[2] = -(G + GAM * k * k * ckap)
    return np.array([row_kin, row_tan, row_nor], dtype=complex)


def det(lam, *a):
    return np.linalg.det(matriz(lam, *a))


def newton(lam0, *a, it=60):
    lam = complex(lam0)
    for _ in range(it):
        f = det(lam, *a)
        h = 1e-7 * max(abs(lam), 1.0)
        df = (det(lam + h, *a) - det(lam - h, *a)) / (2 * h)
        step = f / df
        lam -= step
        if abs(step) < 1e-13 * max(abs(lam), 1.0):
            break
    return lam, abs(det(lam, *a)) / max(abs(det(lam * 1.01, *a)), 1e-300)


def roots_right(k, eta0, s0, N, R, delta, n=40000):
    """number of zeros of det in the rectangle [delta, R] x [-R, R] (argument principle)."""
    t = np.linspace(0, 1, n, endpoint=False)
    pts = np.concatenate([delta + (R - delta) * t - 1j * R,          # bottom edge
                          R + 1j * (-R + 2 * R * t),                   # right edge
                          R - (R - delta) * t + 1j * R,                # top edge
                          delta + 1j * (R - 2 * R * t)])               # left edge
    vals = np.array([det(z, k, eta0, s0, N) for z in pts])
    ph = np.unwrap(np.angle(np.concatenate([vals, vals[:1]])))
    return (ph[-1] - ph[0]) / (2 * np.pi)


def main():
    print("(i) scalar model: roots x = eta*kap of P_{N-1}(x) = 0 and lam*eta^2/nu = x^2")
    for N in (2, 3, 4):
        c = [1.0 / math.factorial(j) for j in range(N)][::-1]
        for x in np.roots(c):
            for sgn in (-1, +1):
                eta = sgn * 1.0
                kap = x / eta
                if kap.real > 0:
                    print("  N=%d eta %s0: x = %s, kap*|eta| = %s, lam*eta^2/nu = %s (%s)" % (
                        N, "<" if sgn < 0 else ">", np.round(x, 4), np.round(kap, 4),
                        np.round(x * x, 4), "UNSTABLE" if (x * x).real > 1e-12 else "neutral/damped"))
    print("(ii) full 2D linear problem, nu=%g g=%g gam=%g" % (NU, G, GAM))
    for N in (2, 3):
        for k in (math.sqrt(2), 10.0):
            for eta0 in (-0.01, -0.04):
                for s0 in (0.0, 0.05):
                    if N == 2:
                        lam0 = NU * (1 / eta0**2 - k * k)
                        pred = lam0
                    else:
                        mm = (1 - 1j) / abs(eta0)
                        lam0 = NU * (mm * mm - k * k)
                        pred = lam0
                    lam, qual = newton(lam0, k, eta0, s0, N)
                    print("  N=%d k=%5.2f eta0=%+.2f s0=%.2f : root lam = %.6g%+.6gi ; "
                          "scalar-model prediction %.6g%+.6gi ; rel.dev %.2e" % (
                              N, k, eta0, s0, lam.real, lam.imag, pred.real, pred.imag,
                              abs(lam - pred) / abs(pred)))
    print("  zeros with Re lam > 0 (argument principle, rectangle [1e-3, R] x [-R, R]):")
    for N in (2, 3):
        for eta0 in (-0.04, -0.01, 0.01, 0.04):
            for s0 in (0.0, 0.05):
                R = 50 * NU / eta0**2
                nz = roots_right(math.sqrt(2), eta0, s0, N, R, 1e-3)
                print("    N=%d eta0=%+.2f s0=%.2f  R=%.0f : %.3f" % (N, eta0, s0, R, nz))


if __name__ == "__main__":
    main()
