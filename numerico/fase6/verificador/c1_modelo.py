#!/usr/bin/env python
"""C1 (verifier): exact rational power-series check of the order-N transferred conditions.

Independent of numerico/fase6/chequeo_expansion.py (which uses numerical Cauchy integrals):
here every quantity is a truncated power series in eps with Fraction coefficients, so the
comparison is EXACT (difference identically zero or not).

Exact conditions derived here from first principles for a graph surface z = h + eta:
  n ~ (-s, 1) (up), t_a = e_a + s_a e_z, S = grad v + grad v^T
  kinematic : eta_t = w - s.u_h                              (material surface)
  tangential: t_a . S . n = 0  ->  S_az - s_b S_ab + s_a S_zz - s_a s_b S_bz = 0
  normal    : p = g eta - gam*kappa + nu * n.S.n / |n|^2,  (Young-Laplace, p minus rest hydrostatics)
              kappa = div(s/sqrt(1+|s|^2)) = ((1+sy^2)exx - 2 sx sy exy + (1+sx^2)eyy)/(1+|s|^2)^{3/2}
All fields evaluated at the displaced point z = h + eta = full Taylor series in eta.
Counting A (amplitude): eta, its derivatives and all fields scale with eps.
Counting B (deformation only): eta and its derivatives scale with eps, fields O(1).
"""
from fractions import Fraction as Fr
import math, random

D = 6  # series degree kept


class Ser:
    __slots__ = ("c",)

    def __init__(self, c):
        c = list(c)[:D + 1]
        self.c = c + [Fr(0)] * (D + 1 - len(c))

    @staticmethod
    def const(a):
        return Ser([Fr(a)])

    @staticmethod
    def eps(a):  # a*eps
        return Ser([Fr(0), Fr(a)])

    def __add__(self, o):
        o = o if isinstance(o, Ser) else Ser.const(o)
        return Ser([a + b for a, b in zip(self.c, o.c)])

    __radd__ = __add__

    def __neg__(self):
        return Ser([-a for a in self.c])

    def __sub__(self, o):
        return self + (-(o if isinstance(o, Ser) else Ser.const(o)))

    def __rsub__(self, o):
        return Ser.const(o) - self

    def __mul__(self, o):
        if not isinstance(o, Ser):
            return Ser([a * Fr(o) for a in self.c])
        r = [Fr(0)] * (D + 1)
        for i, a in enumerate(self.c):
            if a == 0:
                continue
            for j in range(D + 1 - i):
                r[i + j] += a * o.c[j]
        return Ser(r)

    __rmul__ = __mul__

    def powbin(self, alpha):
        """(1 + x)^alpha for x = self with zero constant term."""
        assert self.c[0] == 0
        r, term, coef = Ser.const(1), Ser.const(1), Fr(1)
        for j in range(1, D + 1):
            coef = coef * (Fr(alpha) - (j - 1)) / j
            term = term * self
            r = r + term * coef
        return r


def rnd():
    return Fr(random.randint(-9, 9), random.randint(1, 7))


def datos():
    d = {k: rnd() for k in ("eta", "sx", "sy", "exx", "exy", "eyy")}
    for f in ("u", "v", "w", "p", "ux", "uy", "vx", "vy", "wx", "wy"):
        d[f] = [rnd() for _ in range(D + 3)]
    return d


def disp(d, f, m, eta, amp):
    """d^m f/dz^m at the displaced point, as a series (Taylor in eta, all orders <= D)."""
    r = Ser.const(0)
    pw = Ser.const(1)
    for j in range(D + 1):
        r = r + pw * (Fr(1, math.factorial(j)) * d[f][m + j])
        pw = pw * eta
    return r * amp


def exactas(d, amp_fields, g, gam, nu):
    eta = Ser.eps(d["eta"])
    sx, sy = Ser.eps(d["sx"]), Ser.eps(d["sy"])
    exx, exy, eyy = Ser.eps(d["exx"]), Ser.eps(d["exy"]), Ser.eps(d["eyy"])
    amp = Ser.eps(1) if amp_fields else Ser.const(1)
    F = lambda f, m=0: disp(d, f, m, eta, amp)
    Sxx, Syy, Sxy = 2 * F("ux"), 2 * F("vy"), F("uy") + F("vx")
    Sxz, Syz, Szz = F("u", 1) + F("wx"), F("v", 1) + F("wy"), 2 * F("w", 1)
    K = F("w") - sx * F("u") - sy * F("v")
    Tx = Sxz - (sx * Sxx + sy * Sxy) + sx * Szz - sx * (sx * Sxz + sy * Syz)
    Ty = Syz - (sx * Sxy + sy * Syy) + sy * Szz - sy * (sx * Sxz + sy * Syz)
    q = sx * sx + sy * sy
    nSn = (Szz - 2 * (sx * Sxz + sy * Syz) + sx * sx * Sxx + 2 * sx * sy * Sxy + sy * sy * Syy) * q.powbin(-1)
    kap = ((1 + sy * sy) * exx - 2 * sx * sy * exy + (1 + sx * sx) * eyy) * q.powbin(Fr(-3, 2))
    Phi = F("p") - g * eta + gam * kap - nu * nSn
    return [K, Tx, Ty, Phi]


def truncadas(d, N, amp_fields, g, gam, nu, flip=False):
    """The header formulas of fsorder.f90, implemented from the text."""
    eta = Ser.eps(d["eta"])
    sx, sy = Ser.eps(d["sx"]), Ser.eps(d["sy"])
    exx, exy, eyy = Ser.eps(d["exx"]), Ser.eps(d["exy"]), Ser.eps(d["eyy"])
    amp = Ser.eps(1) if amp_fields else Ser.const(1)

    def tr(vals, M):  # Tr_M[f] with vals(m) the m-th normal derivative at h
        r, pw = Ser.const(0), Ser.const(1)
        for m in range(M + 1):
            r = r + pw * (Fr(1, math.factorial(m))) * vals(m)
            pw = pw * eta
        return r

    f = lambda name, k=0: (lambda m: amp * Ser.const(d[name][m + k]))
    Sxx = lambda m: 2 * f("ux")(m)
    Syy = lambda m: 2 * f("vy")(m)
    Sxy = lambda m: f("uy")(m) + f("vx")(m)
    Sxz = lambda m: f("u", 1)(m) + f("wx")(m)
    Syz = lambda m: f("v", 1)(m) + f("wy")(m)
    Szz = lambda m: 2 * f("w", 1)(m)
    s = (sx, sy)
    Sab = [[Sxx, Sxy], [Sxy, Syy]]
    Saz = [Sxz, Syz]
    K = tr(f("w"), N - 1) - sx * tr(f("u"), N - 2) - sy * tr(f("v"), N - 2)
    T = []
    for a in range(2):
        t = tr(Saz[a], N - 1)
        sl = sum((s[b] * tr(Sab[a][b], N - 2) for b in range(2)), Ser.const(0))
        t = t + sl if flip else t - sl
        t = t + s[a] * tr(Szz, N - 2)
        t = t - s[a] * sum((s[b] * tr(Saz[b], N - 3) for b in range(2)), Ser.const(0))
        T.append(t)
    q = sx * sx + sy * sy
    B = Ser.const(0)
    j = 0
    while 2 * j <= N - 1:
        wk = tr(Szz, N - 1 - 2 * j) - 2 * (sx * tr(Sxz, N - 2 - 2 * j) + sy * tr(Syz, N - 2 - 2 * j)) \
             + sx * sx * tr(Sxx, N - 3 - 2 * j) + 2 * sx * sy * tr(Sxy, N - 3 - 2 * j) \
             + sy * sy * tr(Syy, N - 3 - 2 * j)
        qj = Ser.const(1)
        for _ in range(j):
            qj = qj * (-q)
        B = B + qj * wk
        j += 1
    # kappa_N = div(s * sum_{2j+1<=N} c_j q^j) = sum_j c_j (q^j lap(eta) + 2 j q^{j-1} s_a s_b e_ab)
    lap = exx + eyy
    ssE = sx * sx * exx + 2 * sx * sy * exy + sy * sy * eyy
    kap = Ser.const(0)
    j, cj = 0, Fr(1)
    while 2 * j + 1 <= N:
        qj = Ser.const(1)
        for _ in range(j):
            qj = qj * q
        term = qj * lap
        if j > 0:
            qj1 = Ser.const(1)
            for _ in range(j - 1):
                qj1 = qj1 * q
            term = term + 2 * j * qj1 * ssE
        kap = kap + cj * term
        j += 1
        cj = cj * (Fr(-1, 2) - (j - 1)) / j   # binom(-1/2, j)
    Phi = tr(f("p"), N - 1) - g * eta + gam * kap - nu * B
    return [K, T[0], T[1], Phi]


def first_mismatch(a, b):
    for n in range(D + 1):
        if a.c[n] != b.c[n]:
            return n
    return None


def main():
    random.seed(20260928)
    g, gam, nu = Fr(13, 10), Fr(7, 10), Fr(11, 100)
    nombres = ["K", "T_x", "T_y", "Phi"]
    print("Counting A (eta and fields ~ eps): first degree where formula N differs from exact")
    for dim in ("3D", "2D"):
        for caso in range(4):
            d = datos()
            if dim == "2D":
                d["sy"] = d["exy"] = d["eyy"] = Fr(0)
                for f in ("uy", "vx", "vy", "wy", "v"):
                    d[f] = [Fr(0)] * len(d[f])
            ex = exactas(d, True, g, gam, nu)
            for N in (1, 2, 3, 4):
                tr_ = truncadas(d, N, True, g, gam, nu)
                mm = [first_mismatch(ex[i], tr_[i]) for i in range(4)]
                ok = all(x is None or x > N for x in mm)
                if caso == 0 or not ok:
                    print("  %s case %d N=%d first mismatch degree %s -> %s" % (
                        dim, caso, N, dict(zip(nombres, mm)), "OK (exact through degree N)" if ok else "FAIL"))
                assert ok, "mismatch at degree <= N"
    # negative control: flip the sign of the -s_b S_ab slope term
    d = datos()
    ex = exactas(d, True, g, gam, nu)
    bad = truncadas(d, 2, True, g, gam, nu, flip=True)
    print("negative control (sign of s_b S_ab flipped), N=2: first mismatch degree",
          [first_mismatch(ex[i], bad[i]) for i in range(4)])
    # Counting B: fields O(1)
    print("Counting B (only eta ~ eps, fields O(1)): first mismatch degree")
    d = datos()
    ex = exactas(d, False, g, gam, nu)
    for N in (1, 2, 3):
        tr_ = truncadas(d, N, False, g, gam, nu)
        print("  N=%d:" % N, dict(zip(nombres, [first_mismatch(ex[i], tr_[i]) for i in range(4)])))
    # the Phi entry under B mixes curvature (degree N kept) and velocity terms (degree N-1 kept):
    d2 = dict(d)
    for f in ("u", "v", "w", "p", "ux", "uy", "vx", "vy", "wx", "wy"):
        d2[f] = [Fr(0)] * len(d[f])
    ex = exactas(d2, False, g, gam, nu)
    for N in (1, 2, 3):
        tr_ = truncadas(d2, N, False, g, gam, nu)
        print("  N=%d, fields = 0 (pure gravity+curvature part of Phi): first mismatch degree %s" % (
            N, first_mismatch(ex[3], tr_[3])))


if __name__ == "__main__":
    main()
