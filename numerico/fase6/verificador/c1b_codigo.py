"""C1b: pointwise transliteration of fs_functionals (gph(:,:,1..6) and the kappa flux) with exact
series, vs the exact conditions of c1_modelo.py (amplitude counting)."""
import sys, math, random
from fractions import Fraction as Fr
sys.dont_write_bytecode = True
import c1_modelo as M
Ser = M.Ser

def codigo(d, order, g, gam, nu):
    amp = Ser.eps(1)
    et = Ser.eps(d["eta"]); sx = Ser.eps(d["sx"]); sy = Ser.eps(d["sy"])
    exx, exy, eyy = Ser.eps(d["exx"]), Ser.eps(d["exy"]), Ser.eps(d["eyy"])
    nm = order
    q2 = sx * sx + sy * sy
    Uu = [amp * d["u"][m] for m in range(nm + 1)]; Vv = [amp * d["v"][m] for m in range(nm + 1)]
    Ww = [amp * d["w"][m] for m in range(nm + 1)]
    Sxx = [2 * amp * d["ux"][m] for m in range(nm + 1)]; Syy = [2 * amp * d["vy"][m] for m in range(nm + 1)]
    Sxy = [amp * (d["uy"][m] + d["vx"][m]) for m in range(nm + 1)]
    Sxz = [Ser.const(0)] * (nm + 1); Syz = [Ser.const(0)] * (nm + 1); Szz = [Ser.const(0)] * (nm + 1)
    for m in range(nm):                      # DO m = 0,nm-1
        Sxz[m] = Uu[m + 1] + amp * d["wx"][m]
        Syz[m] = Vv[m + 1] + amp * d["wy"][m]
        Szz[m] = 2 * Ww[m + 1]
    def trz(A, Mt):
        r = Ser.const(0)
        if Mt < 0: return r
        pw = Ser.const(1)
        for mm in range(0, min(Mt, nm) + 1):
            r = r + pw * A[mm]
            pw = pw * et * Fr(1, mm + 1)
        return r
    K = trz(Ww, order - 1) - sx * trz(Uu, order - 2) - sy * trz(Vv, order - 2)
    Tx = trz(Sxz, order - 1) - (sx * trz(Sxx, order - 2) + sy * trz(Sxy, order - 2)) \
         + sx * trz(Szz, order - 2) - sx * (sx * trz(Sxz, order - 3) + sy * trz(Syz, order - 3))
    Ty = trz(Syz, order - 1) - (sx * trz(Sxy, order - 2) + sy * trz(Syy, order - 2)) \
         + sy * trz(Szz, order - 2) - sy * (sx * trz(Sxz, order - 3) + sy * trz(Syz, order - 3))
    B = Ser.const(0); jj = 0
    while 2 * jj <= order - 1:
        wk = trz(Szz, order - 1 - 2 * jj) - 2 * (sx * trz(Sxz, order - 2 - 2 * jj) + sy * trz(Syz, order - 2 - 2 * jj)) \
             + sx * sx * trz(Sxx, order - 3 - 2 * jj) + 2 * sx * sy * trz(Sxy, order - 3 - 2 * jj) + sy * sy * trz(Syy, order - 3 - 2 * jj)
        mq = Ser.const(1)
        for _ in range(jj): mq = mq * (-q2)
        B = B + mq * wk; jj += 1
    B = nu * B
    ck = [Fr(1), Fr(-1, 2), Fr(3, 8), Fr(-5, 16)]
    if order <= 2:
        kap = exx + eyy                                   # (g + gam k^2) eta: kappa = lap(eta)
    else:                                                 # div(s * sum ck q2^jj), taken pointwise
        wk = Ser.const(0); dwk_x = Ser.const(0); dwk_y = Ser.const(0); jj = 0
        while 2 * jj + 1 <= order:
            qj = Ser.const(1)
            for _ in range(jj): qj = qj * q2
            wk = wk + ck[jj] * qj
            if jj > 0:
                qj1 = Ser.const(1)
                for _ in range(jj - 1): qj1 = qj1 * q2
                dwk_x = dwk_x + ck[jj] * jj * qj1 * 2 * (sx * exx + sy * exy)
                dwk_y = dwk_y + ck[jj] * jj * qj1 * 2 * (sx * exy + sy * eyy)
            jj += 1
        kap = (exx + eyy) * wk + sx * dwk_x + sy * dwk_y
    Ptr = trz([amp * d["p"][m] for m in range(nm + 1)], order - 1)   # Tr_{N-1}[p] (fs_ptransfer + d(Lz))
    E = g * et - gam * kap + B
    return [K, Tx, Ty, Ptr - E]

random.seed(5)
g, gam, nu = Fr(13, 10), Fr(7, 10), Fr(11, 100)
for caso in range(3):
    d = M.datos()
    ex = M.exactas(d, True, g, gam, nu)
    for N in (1, 2, 3):
        c = codigo(d, N, g, gam, nu)
        mm = [M.first_mismatch(ex[i], c[i]) for i in range(4)]
        print("case %d N=%d: first mismatch degree (K, Tx, Ty, Phi) = %s -> %s" % (
            caso, N, mm, "OK" if all(x is None or x > N for x in mm) else "FAIL"))
