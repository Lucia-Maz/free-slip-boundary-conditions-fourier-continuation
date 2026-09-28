"""C3b: (1) the replica's own main() with outputs redirected to the verifier scratch;
(2) onset of the N=2 instability (implicit spectral closure) vs the diffusion number D = nu dt/dz^2,
    compared with the guard 1+sigma <= 0.05 (r <= -0.6802) and with -1/beta2 = -0.716;
(3) N=3 and q=8 checks at valid (stable at eta=0) dt for nz = 256."""
import os, sys, json
import numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, "/home/lucia/GW_AI/proyecto-final-piv/numerico/fase6")
import estabilidad_cierre as EC
import c3_replica as V
AQUI = os.path.dirname(os.path.abspath(__file__))
EC.SALIDA = os.path.join(AQUI, "replica_fase6_estabilidad_cierre.json")
EC.FIG = os.path.join(AQUI, "replica_f_fs6_estabilidad")
print("=== (1) replica main() (outputs -> verifier scratch) ===")
EC.main()
print("=== (2) N=2 onset vs D = nu dt/dz^2 (nz=128, nu=0.01, my step) ===")
dz = 1.0 / 102
for dt in (1e-4, 2.5e-4, 5e-4, 1e-3, 1.25e-3, 1.5e-3, 1.7e-3):
    mp = V.MiPaso(128, 1e-2, dt)
    r0 = mp.rho(0.0, 2)
    if r0 > 1 + 1e-9:
        print("  dt=%.2e D=%.3f: base step unstable (rho(0)=%.5f)" % (dt, 1e-2 * dt / dz**2, r0)); continue
    # bisection for the onset between -0.716 (unstable) and -0.40 (stable)
    lo, hi = -0.7159, -0.40
    assert mp.rho(lo, 2) > 1 + 1e-9 and mp.rho(hi, 2) <= 1 + 1e-9
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        if mp.rho(mid, 2) > 1 + 1e-9: lo = mid
        else: hi = mid
    s = 1 + 1.3967 * lo
    print("  dt=%.2e D=%.3f: rho(0)=%.8f ; N=2 unstable for r <= %.4f (1+sigma there = %.4f) ; guard fires r <= -0.6802 -> %s" % (
        dt, 1e-2 * dt / dz**2, r0, lo, s, "window (%.4f, %.4f] UNGUARDED" % (-0.6802, lo) if lo > -0.6802 else "covered"))
print("=== (3) nz=256 with dt=2.5e-4 (valid), and nz=64 ===")
for nz, dt in ((256, 2.5e-4), (256, 1e-4), (64, 1e-3)):
    mp = V.MiPaso(nz, 1e-2, dt)
    r0 = mp.rho(0.0, 2)
    lo, hi = -0.7159, -0.40
    if mp.rho(lo, 2) > 1 + 1e-9 and mp.rho(hi, 2) <= 1 + 1e-9:
        for _ in range(25):
            mid = 0.5 * (lo + hi)
            if mp.rho(mid, 2) > 1 + 1e-9: lo = mid
            else: hi = mid
    wide = np.round(np.concatenate([np.arange(-8, -1.0, 0.5), np.arange(-1.0, 1.01, 0.1), np.arange(1.5, 9.01, 0.5)]), 3)
    r3 = [mp.rho(r, 3) for r in wide]
    print("  nz=%d dt=%.1e D=%.3f rho(0)=%.8f: N=2 onset r=%.4f ; N=3 max rho over [-8,9] = %.8f (%s)" % (
        nz, dt, 1e-2 * dt * (nz - 26)**2, r0, lo, max(r3), "stable" if max(r3) <= 1 + 1e-9 else "UNSTABLE"))
