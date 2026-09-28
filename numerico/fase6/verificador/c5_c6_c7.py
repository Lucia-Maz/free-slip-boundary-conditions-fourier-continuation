#!/usr/bin/env python
"""C5, C6, C7 (verifier) toy checks with literal transliterations of the Fortran routines.

C5: fs_to_phys / fs_to_spec (padded grid 2nx x 2ny, Nyquist dropped). (i) round trip on V;
    (ii) L(e) = to_spec(to_phys(e)/(1+sigma)) on V: fixed points of the preconditioned map
    coincide with those of F iff L has no kernel; eigenvalues for 1+sigma>0 and for a sign change;
    (iii) sign check of G = (F - 2 sigma dx(Wn - W) + sigma Dx)/(1+sigma) on a model pass
    F = P[R + sigma*(dx(w* - 2W) - Dx)] (the modelled response of the code's datum).
C6: fs_anderson transliterated (history slots, differences, real Gram, 1e-13 trace regularization)
    vs a textbook type-II Anderson in R^{2n}; and complex-gamma variant, on an R-linear map.
C7: Nyquist convention for odd derivatives at the node z = Lz.
"""
import math
import numpy as np

rng = np.random.default_rng(12345)


# ---------------------------------------------------------------- C5
class Surf:
    def __init__(self, nx, ny, Lx=1.0, Ly=1.0):
        self.nx, self.ny = nx, ny
        self.nkx = nx // 2
        self.mx = 2 * nx
        self.my = 2 * ny if ny > 1 else 1
        kx = np.array([i for i in range(nx // 2)] + [i - nx // 2 for i in range(nx // 2)], float) / Lx
        ky = np.array([j for j in range(ny // 2)] + [j - ny // 2 for j in range(ny // 2)], float) / Ly if ny > 1 else np.zeros(1)
        self.kx, self.ky = kx, ky
        self.nloc = nx // 2 + 1          # ista:iend = 1..nx/2+1 on one rank
        xp = 2 * np.pi * Lx * np.arange(self.mx) / self.mx
        yp = 2 * np.pi * Ly * np.arange(self.my) / self.my
        self.ex = np.exp(1j * np.outer(xp, kx[:self.nkx]))          # (mx, nkx)
        self.ey = np.exp(1j * np.outer(yp, ky))                      # (my, ny)

    def to_phys(self, F):               # F: (ny, nloc)
        full = np.zeros((self.ny, self.nkx), complex)
        full[:, :] = F[:, :self.nkx]
        if self.ny > 1:
            full[self.ny // 2, :] = 0
        fac = 1.0 / (self.nx * self.ny)
        g = self.ey @ full                                           # (my, nkx)
        out = np.zeros((self.mx, self.my))
        for q in range(self.my):
            out[:, q] = np.real(g[q, 0])
            for i in range(1, self.nkx):
                out[:, q] += 2.0 * np.real(self.ex[:, i] * g[q, i])
        return out * fac

    def to_spec(self, f):               # f: (mx, my) -> (ny, nloc)
        fac = self.nx * self.ny / (self.mx * self.my)
        Fo = np.zeros((self.ny, self.nloc), complex)
        h = np.conj(self.ex).T @ f                                   # (nkx, my)
        for i in range(self.nkx):
            for j in range(self.ny):
                if self.ny > 1 and j == self.ny // 2:
                    continue
                Fo[j, i] = fac * np.sum(np.conj(self.ey[:, j]) * h[i, :])
        return Fo

    def random_V(self, scale=1.0):
        """random real band-limited field in V (retained modes), SPECTER normalization."""
        f = np.zeros((self.mx, self.my))
        for i in range(self.nkx):
            for j in range(self.ny):
                if self.ny > 1 and j == self.ny // 2:
                    continue
                a = rng.normal() + 1j * rng.normal()
                ph = np.exp(1j * (self.kx[i] * 2 * np.pi * np.arange(self.mx)[:, None] / self.mx
                                  + self.ky[j] * 2 * np.pi * np.arange(self.my)[None, :] / self.my))
                f += np.real(a * ph)
        return self.to_spec(f * scale)


def realvec(F, S):
    """real parameterization of V: Re, Im of the independent coefficients (kx>0 all ky; kx=0 ky>=0)."""
    v = []
    for i in range(S.nkx):
        for j in range(S.ny):
            if S.ny > 1 and j == S.ny // 2:
                continue
            if i == 0 and S.ky[j] < 0:
                continue
            v += [F[j, i].real] if (i == 0 and S.ky[j] == 0) else [F[j, i].real, F[j, i].imag]
    return np.array(v)


def fromreal(v, S):
    F = np.zeros((S.ny, S.nloc), complex)
    c = 0
    for i in range(S.nkx):
        for j in range(S.ny):
            if S.ny > 1 and j == S.ny // 2:
                continue
            if i == 0 and S.ky[j] < 0:
                continue
            if i == 0 and S.ky[j] == 0:
                F[j, i] = v[c]
                c += 1
            else:
                F[j, i] = v[c] + 1j * v[c + 1]
                c += 2
    for j in range(S.ny):                        # Hermitian completion of the kx = 0 plane
        if S.ny > 1 and j == S.ny // 2:
            continue
        if S.ky[j] < 0:
            jj = int(np.where(S.ky == -S.ky[j])[0][0])
            F[j, 0] = np.conj(F[jj, 0])
    F[0, 0] = F[0, 0].real
    return F


def c5():
    print("C5")
    for nx, ny in ((8, 8), (16, 1)):
        S = Surf(nx, ny)
        F = S.random_V()
        rt = np.max(np.abs(S.to_spec(S.to_phys(F)) - F)) / np.max(np.abs(F))
        eta = S.to_phys(S.random_V())
        eta *= 1.0 / np.max(np.abs(eta))                          # max|eta| = 1 (units of dz below)
        b2, b3 = 1.3967, 1.1515
        dim = realvec(F, S).size
        out = []
        for label, sig in (("N=3, r=eta/dz up to 4 (1+sigma>=0.153)", b2 * 4 * eta + 0.5 * b3 * (4 * eta)**2),
                           ("N=2, eta/dz in [-0.6,0.6] (1+sigma>0)", b2 * 0.6 * eta),
                           ("N=2, eta/dz in [-1,1]   (1+sigma changes sign)", b2 * 1.0 * eta)):
            L = np.zeros((dim, dim))
            for c in range(dim):
                e = np.zeros(dim)
                e[c] = 1.0
                L[:, c] = realvec(S.to_spec(S.to_phys(fromreal(e, S)) / (1.0 + sig)), S)
            ev = np.linalg.eigvals(L)
            out.append((label, float(np.min(1 + sig)), float(np.min(ev.real)), float(np.max(np.abs(ev.imag)))))
        print("  nx=%d ny=%d: round trip rel.err %.1e ; dim V (real) = %d" % (nx, ny, rt, dim))
        for lab, smin, evmin, evim in out:
            print("    %-48s min(1+sigma)=%+.3f  min Re eig(L)=%+.4f  max|Im eig|=%.1e" % (lab, smin, evmin, evim))
        # scan a scaling of a sign-changing sigma to find a singular L (spurious fixed point)
        if (nx, ny) == (16, 1):
            best = None
            for r in np.linspace(0.72, 3.0, 229):
                sig = b2 * r * eta
                L = np.zeros((dim, dim))
                for c in range(dim):
                    e = np.zeros(dim); e[c] = 1.0
                    L[:, c] = realvec(S.to_spec(S.to_phys(fromreal(e, S)) / (1.0 + sig)), S)
                sv = np.linalg.svd(L, compute_uv=False)
                cand = (sv[-1] / sv[0], r)
                best = cand if best is None or cand < best else best
            print("    sign-changing sigma: smallest sigma_min/sigma_max of L over the scan = %.2e at r=%.2f"
                  " (L nearly singular -> G can have fixed points that F does not)" % best)
        # (iii) sign check on a model pass with the code's datum convention
        wst = S.random_V(); W = S.random_V(); Dx = S.random_V(); R = S.random_V()
        sig = b2 * 0.5 * eta                                         # a legal sigma
        ikx = 1j * S.kx[None, :S.nloc] if S.nloc <= S.nx else None
        ikx = np.zeros((S.ny, S.nloc), complex)
        for i in range(S.nloc):
            ikx[:, i] = 1j * (S.kx[i] if i < S.nx else 0)
        A = S.random_V(0.3); Bk = 0.1                                # W_new = A + Bk*Dx (model of the B_k chain)
        def passF(Win, Din):
            g = S.to_phys(ikx * (wst - 2 * Win)) - S.to_phys(Din)    # datum of u* (physical)
            Fx = S.to_spec(S.to_phys(R) + sig * g)                   # Dx output: responds with +sigma*datum
            Wn = A + Bk * Din
            return Wn, Fx
        # exact fixed point: D = P[R + sig*(dx(w*-2W(D)) - D)], W = A + Bk D -> solve by brute force in V
        dimv = realvec(R, S).size
        def resid(v):
            D = fromreal(v, S)
            Wn, Fx = passF(A + Bk * D, D)
            return realvec(Fx - D, S)
        J = np.zeros((dimv, dimv)); r0 = resid(np.zeros(dimv))
        for c in range(dimv):
            e = np.zeros(dimv); e[c] = 1.0
            J[:, c] = resid(e) - r0
        Dstar = fromreal(np.linalg.solve(J, -r0), S)
        # one step of the code's G from an arbitrary (W, Dx)
        Wn, Fx = passF(W, Dx)
        dW = S.to_phys(ikx * (Wn - W))
        G = S.to_spec((S.to_phys(Fx) - 2 * sig * dW + sig * S.to_phys(Dx)) / (1 + sig))
        # iterate G (Picard) a few times; with the model's exact sigma it converges geometrically (Bk coupling only)
        Wc, Dc = W, Dx
        errs = []
        for it in range(8):
            Wn, Fx = passF(Wc, Dc)
            dW = S.to_phys(ikx * (Wn - Wc))
            Gc = S.to_spec((S.to_phys(Fx) - 2 * sig * dW + sig * S.to_phys(Dc)) / (1 + sig))
            Wc, Dc = Wn, Gc
            errs.append(np.max(np.abs(Dc - Dstar)) / np.max(np.abs(Dstar)))
        # same with the W-term sign flipped (negative control)
        Wc, Dc = W, Dx
        errb = []
        for it in range(8):
            Wn, Fx = passF(Wc, Dc)
            dW = S.to_phys(ikx * (Wn - Wc))
            Gc = S.to_spec((S.to_phys(Fx) + 2 * sig * dW + sig * S.to_phys(Dc)) / (1 + sig))
            Wc, Dc = Wn, Gc
            errb.append(np.max(np.abs(Dc - Dstar)) / np.max(np.abs(Dstar)))
        print("    model pass, Picard on the code's G: rel.err per pass %s" % " ".join("%.1e" % e for e in errs))
        print("    same with +2 sigma dx(dW) (sign flipped):          %s" % " ".join("%.1e" % e for e in errb))


# ---------------------------------------------------------------- C6
class AndersonFortran:
    """literal transliteration of fs_anderson (slot 1 = newest)."""

    def __init__(self, shape, mhist):
        self.m = mhist
        self.ahx = np.zeros(shape + (mhist,), complex)
        self.ahr = np.zeros(shape + (mhist,), complex)
        self.nhist = 0

    def step(self, xa, fxa, esc):
        r = (fxa - xa) / esc[None, :]
        if self.nhist >= 1:
            mk = min(self.nhist, self.m)
            dR = np.zeros(xa.shape + (mk,), complex)
            dF = np.zeros(xa.shape + (mk,), complex)
            for m in range(1, mk + 1):
                if m == 1:
                    dR[..., 0] = r - self.ahr[..., 0]
                    dF[..., 0] = fxa - self.ahx[..., 0]
                else:
                    dR[..., m - 1] = self.ahr[..., m - 2] - self.ahr[..., m - 1]
                    dF[..., m - 1] = self.ahx[..., m - 2] - self.ahx[..., m - 1]
            Gm = np.array([[np.real(np.sum(np.conj(dR[..., a]) * dR[..., c])) for c in range(mk)] for a in range(mk)])
            bv = np.array([np.real(np.sum(np.conj(dR[..., a]) * r)) for a in range(mk)])
            tr_ = np.trace(Gm)
            Gm = Gm + 1e-13 * tr_ * np.eye(mk)
            gam = np.linalg.solve(Gm, bv)
            xna = fxa.copy()
            for m in range(mk):
                xna = xna - gam[m] * dF[..., m]
        else:
            xna = fxa.copy()
        for m in range(self.m, 1, -1):
            self.ahx[..., m - 1] = self.ahx[..., m - 2]
            self.ahr[..., m - 1] = self.ahr[..., m - 2]
        self.ahx[..., 0] = fxa
        self.ahr[..., 0] = r
        self.nhist += 1
        return xna


def anderson_textbook(Fmap, x0, esc, mhist, nit, complex_gamma=False):
    """type-II Anderson: x_{k+1} = F(x_k) - dF gamma, gamma = argmin |r_k - dR gamma|, over the last
    min(k, m) differences; real gamma = LS in R^{2n} (or complex LS if complex_gamma)."""
    xs, Fs, Rs = [], [], []
    x = x0.copy()
    hist = []
    for k in range(nit):
        f = Fmap(x)
        r = (f - x) / esc[None, :]
        hist.append(np.max(np.abs(r)))
        Fs.append(f); Rs.append(r)
        mk = min(k, mhist)
        if mk == 0:
            x = f
            continue
        dR = np.stack([(Rs[-1 - j] - Rs[-2 - j]).ravel() for j in range(mk)], axis=1)
        dF = np.stack([(Fs[-1 - j] - Fs[-2 - j]).ravel() for j in range(mk)], axis=1)
        if complex_gamma:
            gam = np.linalg.lstsq(dR, r.ravel(), rcond=None)[0]
        else:
            A = np.vstack([dR.real, dR.imag]); b = np.concatenate([r.ravel().real, r.ravel().imag])
            gam = np.linalg.solve(A.T @ A + 1e-13 * np.trace(A.T @ A) * np.eye(mk), A.T @ b)
        x = (f.ravel() - dF @ gam).reshape(f.shape)
    return hist, x


def c6():
    print("C6")
    n = 40
    A = rng.normal(size=(n * 4, n * 4)) * 0.08 + 0.1j * rng.normal(size=(n * 4, n * 4)) * 0.08
    B = rng.normal(size=(n * 4, n * 4)) * 0.08
    c = rng.normal(size=n * 4) + 1j * rng.normal(size=n * 4)
    esc = np.array([1.0, 2.0, 2.0, 0.5])
    Fmap = lambda x: (A @ x.ravel() + B @ np.conj(x.ravel()) + c).reshape(x.shape)   # R-linear, not C-linear
    x0 = np.zeros((n, 4), complex)
    # Fortran transliteration
    AF = AndersonFortran((n, 4), 5)
    x = x0.copy()
    hF, xsF = [], []
    for k in range(25):
        f = Fmap(x)
        hF.append(np.max(np.abs((f - x) / esc[None, :])))
        x = AF.step(x, f, esc)
    hT, _ = anderson_textbook(Fmap, x0, esc, 5, 25)
    hC, _ = anderson_textbook(Fmap, x0, esc, 5, 25, complex_gamma=True)
    hP = []
    x = x0.copy()
    for k in range(25):
        f = Fmap(x); hP.append(np.max(np.abs((f - x) / esc[None, :]))); x = f
    print("  residual per pass (every 4th): fortran-literal  ", " ".join("%.1e" % v for v in hF[::4]))
    print("                                 textbook real    ", " ".join("%.1e" % v for v in hT[::4]))
    print("                                 textbook complex ", " ".join("%.1e" % v for v in hC[::4]))
    print("                                 Picard           ", " ".join("%.1e" % v for v in hP[::4]))
    print("  max rel. difference fortran-literal vs textbook-real residual histories: %.1e" % (
        max(abs(a - b) / max(a, 1e-300) for a, b in zip(hF[:18], hT[:18]))))


# ---------------------------------------------------------------- C7
def c7():
    print("C7")
    nz, C = 128, 25
    dz = 1.0 / (nz - C - 1)
    Dkz = 2 * np.pi / (dz * nz)
    kz = np.array([k for k in range(nz // 2)] + [k - nz // 2 for k in range(nz // 2)], float) * Dkz
    z = dz * np.arange(nz)
    top = nz - C - 1
    print("  z(top) = %.17g (Lz = 1), (top) = node index %d, z(top)/dz = %.15g" % (z[top], top, z[top] / dz))
    f = rng.normal(size=nz)                                   # real column
    Fh = np.fft.fft(f)
    # symmetric real interpolant: Nyquist term = (F_N/2 / nz) cos(pi z / dz)
    def interp_deriv(m, zz):
        s = 0.0 + 0j
        for k in range(nz):
            if k == nz // 2:
                c = Fh[k] / nz
                s += c * (np.pi / dz)**m * np.cos(np.pi * zz / dz + m * np.pi / 2)
            else:
                s += Fh[k] / nz * (1j * kz[k])**m * np.exp(1j * kz[k] * zz)
        return s
    for m in range(4):
        ph = (1j * kz)**m * np.exp(1j * kz * z[top]) / nz
        ph_drop = ph.copy()
        if m % 2:
            ph_drop[nz // 2] = 0
        a = interp_deriv(m, z[top])
        print("  m=%d: symmetric interpolant %+.6e%+.1ei | fs_ph (Nyquist dropped for odd m) %+.6e%+.1ei | kz=-nz/2 kept %+.6e%+.1ei" % (
            m, a.real, a.imag, (Fh @ ph_drop).real, (Fh @ ph_drop).imag, (Fh @ ph).real, (Fh @ ph).imag))
    # complex columns of a real 3D field: Hermitian pair (kx,ky) and (-kx,-ky) must give conjugate traces
    g = rng.normal(size=nz) + 1j * rng.normal(size=nz)       # column of mode (0, ky)
    gm = np.conj(g)                                           # column of mode (0, -ky)
    for m in (1, 3):
        ph = (1j * kz)**m * np.exp(1j * kz * z[top]) / nz
        pd = ph.copy(); pd[nz // 2] = 0
        t1, t2 = np.fft.fft(g) @ ph, np.fft.fft(gm) @ ph
        d1, d2 = np.fft.fft(g) @ pd, np.fft.fft(gm) @ pd
        print("  m=%d Hermitian pair: kept |t(-k) - conj t(k)| = %.2e ; dropped: %.2e" % (m, abs(t2 - np.conj(t1)), abs(d2 - np.conj(d1))))


if __name__ == "__main__":
    c5()
    c6()
    c7()
