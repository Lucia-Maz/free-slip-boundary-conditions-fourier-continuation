"""Round 2: convergence of the spectral traces d^m f/dz^m at z = Lz (FC-Gram C=25, d=5), in the
literal Fortran pipeline, float64 (as SPECTER) and fully in extended precision (longdouble:
continuation, DFT by matrix, trace sum) to separate FC truncation from roundoff."""
import math, sys
import numpy as np
sys.dont_write_bytecode = True
import c3_replica as V
LD = np.longdouble
funcs = {
    "sin(1.2z)":  (lambda z: np.sin(1.2 * z), lambda m: 1.2**m * math.sin(1.2 + m * math.pi / 2)),
    "exp(0.6z)":  (lambda z: np.exp(0.6 * z), lambda m: 0.6**m * math.exp(0.6)),
    "cosh(1.1z)": (lambda z: np.cosh(1.1 * z), lambda m: 1.1**m * (math.cosh(1.1) if m % 2 == 0 else math.sinh(1.1))),
    "0.3cos(3z+1)": (lambda z: 0.3 * np.cos(3 * z + 1), lambda m: 0.3 * 3**m * math.cos(3 + 1 + m * math.pi / 2)),
}
nzs = (64, 128, 256, 512)

def ld_trace(Fo, vals, nz):
    n = nz - V.C
    f = np.zeros(nz, dtype=LD); f[:n] = vals
    dirl = Fo.dir.astype(LD); Cc, dd = V.C, V.d
    for ii in range(1, Cc + 1):
        s = LD(0)
        for jj in range(1, dd + 1):
            s += dirl[ii - 1, jj - 1] * f[nz - Cc - dd + jj - 1] + dirl[Cc - ii, jj - 1] * f[dd - jj]
        f[nz - Cc + ii - 1] = s
    k = np.arange(nz)
    kk = np.where(k < nz // 2, k, k - nz).astype(LD)
    dz = LD(1) / LD(nz - Cc - 1)
    Dkz = LD(2) * LD(np.pi) / (dz * nz)          # pi to double precision is enough here
    kz = kk * Dkz
    ang = LD(2) * LD(np.pi) * np.outer(k, k).astype(LD) / LD(nz)
    Fr = (np.cos(ang) @ f); Fi = -(np.sin(ang) @ f)            # forward DFT, sign -1
    ztop = dz * LD(nz - Cc - 1)
    out = []
    for m in range(4):
        # (i kz)^m e^{i kz ztop} / nz, Nyquist dropped for odd m
        mag = kz**m; ph = kz * ztop + LD(m) * LD(np.pi) / 2
        wr, wi = mag * np.cos(ph) / nz, mag * np.sin(ph) / nz
        if m % 2: wr[nz // 2] = 0; wi[nz // 2] = 0
        out.append(float(np.sum(Fr * wr - Fi * wi)))
    return out

for name, (f, df) in funcs.items():
    e64 = {m: [] for m in range(4)}; eld = {m: [] for m in range(4)}
    for nz in nzs:
        Fo = V.Fortran(nz)
        n = nz - V.C
        col = np.zeros(nz, complex); col[:n] = f(Fo.z[:n])
        fh = Fo.fwd(col)
        tl = ld_trace(Fo, f(Fo.z[:n].astype(LD)), nz)
        for m in range(4):
            ref = df(m)
            e64[m].append(abs(Fo.trace(fh, m).real - ref) / max(abs(ref), 1e-300))
            eld[m].append(abs(tl[m] - ref) / max(abs(ref), 1e-300))
    print(name)
    for m in (1, 2, 3):
        o = [math.log(eld[m][i] / eld[m][i + 1]) / math.log((nzs[i + 1] - 26) / (nzs[i] - 26)) if eld[m][i + 1] > 0 else float('nan') for i in range(3)]
        print("   m=%d float64 %s | extended %s ; orders (extended) %s" % (
            m, " ".join("%.1e" % x for x in e64[m]), " ".join("%.1e" % x for x in eld[m]), " ".join("%.2f" % x for x in o)))
