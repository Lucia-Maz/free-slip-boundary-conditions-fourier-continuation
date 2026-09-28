"""Round 2: m-th spectral trace at z=Lz for polynomials and for the H2 profile, nz up to 1024 (float64)."""
import math, sys
import numpy as np
sys.dont_write_bytecode = True
import c3_replica as V
tests = {
    "z^2":  (lambda z: z**2, [1, 2, 2, 0]),
    "z^3":  (lambda z: z**3, [1, 3, 6, 6]),
    "z^4":  (lambda z: z**4, [1, 4, 12, 24]),
    "z^5":  (lambda z: z**5, [1, 5, 20, 60]),
    "1+0.3z^2 (H2 profile A3)": (lambda z: 1 + 0.3 * z**2, [1.3, 0.6, 0.6, 0.0]),
    "sin(1.2z) (H2 profile A2)": (lambda z: np.sin(1.2 * z), [1.2**m * math.sin(1.2 + m * math.pi / 2) for m in range(4)]),
}
nzs = (64, 128, 256, 512, 1024)
for name, (f, ex) in tests.items():
    rows = {m: [] for m in (2, 3)}
    for nz in nzs:
        Fo = V.Fortran(nz)
        n = nz - V.C
        col = np.zeros(nz, complex); col[:n] = f(Fo.z[:n])
        fh = Fo.fwd(col)
        for m in (2, 3):
            t = float(np.real(fh @ Fo.ph[:, m]))
            rows[m].append(abs(t - ex[m]) / max(abs(ex[m]), 1.0))
    for m in (2, 3):
        print("%-28s m=%d  err (abs, or rel if |exact|>1), nz=64..1024: %s" % (name, m, " ".join("%.1e" % x for x in rows[m])))
