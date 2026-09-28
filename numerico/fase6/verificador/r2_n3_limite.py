"""Round 2: N=3 instability near the explicit limit: smallest D/Dmax with an unstable eta, and rho."""
import sys, numpy as np
sys.dont_write_bytecode = True
import r2_guardia as G
for nz in (64, 128, 256):
    Dm = G.dmax(nz)
    rs = np.round(np.arange(-3.0, -0.49, 0.05), 3)
    def worst(frac):
        P = G.Paso(nz, frac * Dm)
        vals = [(P.rho(r, 3), r) for r in rs]
        return max(vals)
    lo, hi = 0.90, 0.9999          # lo: stable for all r (checked), hi: unstable?
    rho_hi, r_hi = worst(hi)
    if rho_hi <= 1 + 1e-9:
        print("nz=%d: N=3 stable for all r in [-3,-0.5] up to 0.9999 Dmax (max rho %.9f)" % (nz, rho_hi)); continue
    assert worst(lo)[0] <= 1 + 1e-9
    for _ in range(18):
        mid = 0.5 * (lo + hi)
        if worst(mid)[0] > 1 + 1e-9: hi = mid
        else: lo = mid
    P = G.Paso(nz, 0.995 * Dm)
    r995 = max((P.rho(r, 3), r) for r in rs)
    print("nz=%d Dmax=%.5f: N=3 first unstable at D/Dmax = %.4f (D = %.5f), worst r there %.2f; "
          "at 0.995 Dmax max rho = %.6f at r=%.2f (x%.3g per 1000 steps)" % (
              nz, Dm, hi, hi * Dm, worst(hi)[1], r995[0], r995[1], r995[0]**1000))
