"""Round 2: guard coverage if SPECTER is built with ORD = 3 or 4 (u^(o) = u0 + dt/o RHS, o = ORD..1)."""
import sys, numpy as np
sys.dont_write_bytecode = True
import r2_guardia as G

class PasoOrd(G.Paso):
    def __init__(self, nz, D, ordk):
        super().__init__(nz, D); self.ordk = ordk
    def rho(self, r, N):
        eta = r * self.dz; U0 = np.eye(self.n); U = U0
        for o in range(self.ordk, 0, -1):
            U = self.bc(U0 + self.dt / o * self.nu * (self.mp.D2 @ U), eta, N)
        return float(np.max(np.abs(np.linalg.eigvals(U))))

rs_neg = np.round(np.arange(-0.715, -0.0, 0.005), 4)
rs3 = np.round(np.concatenate([np.arange(-6, -2.0, 0.5), np.arange(-2.0, 2.0, 0.05), np.arange(2.0, 9.01, 1.0)]), 4)
for ordk in (2, 3, 4):
    lo, hi = 0.1, 0.35
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        if PasoOrd(128, mid, ordk).rho(0.0, 2) <= 1 + 1e-12: lo = mid
        else: hi = mid
    Dm = lo
    print("ORD=%d: D_max(eta=0) = %.4f" % (ordk, Dm))
    for frac in (0.25, 0.5, 0.75, 0.9, 0.97):
        P = PasoOrd(128, frac * Dm, ordk)
        un = [r for r in rs_neg if P.rho(r, 2) > 1 + 1e-9]
        edge = G.bisect_onset(P, max(un), min([r for r in rs_neg if r > max(un)] or [0.0])) if un else float('nan')
        s_on = 1 + P.b2 * edge
        sthr = 0.05 + frac * Dm
        un3 = [r for r in rs3 if P.rho(r, 3) > 1 + 1e-9]
        print("   D=%.4f (%.2f Dmax): N=2 onset 1+sig=%.4f vs sthr=%.4f -> %s ; N=3 unstable: %s" % (
            frac * Dm, frac, s_on, sthr, "covered, margin %.3f" % (sthr - s_on) if sthr > s_on else "NOT COVERED",
            ("%d pts %.2f..%.2f" % (len(un3), min(un3), max(un3))) if un3 else "none"))
        sys.stdout.flush()
