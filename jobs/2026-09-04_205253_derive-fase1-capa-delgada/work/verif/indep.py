"""Verificacion independiente (verifier). Metodos distintos a los del worker y a los de la puerta."""
import importlib.util, json, re, sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
from scipy.linalg import solve_banded

JOB = Path("/home/lucia/GW_AI/proyecto-final-piv/jobs/2026-09-04_205253_derive-fase1-capa-delgada")
spec = importlib.util.spec_from_file_location("capa_v", JOB/"out"/"capa.py")
capa = importlib.util.module_from_spec(spec); sys.modules["capa_v"] = capa
spec.loader.exec_module(capa)

def hdr(s): print("\n== " + s + " " + "="*(60-len(s)))

# ---------------------------------------------------------------- 1. espectro
# Ruta independiente A: colocacion de Chebyshev (NO diferencias finitas).
def cheb(N):
    x = np.cos(np.pi*np.arange(N+1)/N)
    c = np.hstack([2., np.ones(N-1), 2.])*(-1)**np.arange(N+1)
    X = np.tile(x, (N+1,1)).T
    dX = X - X.T
    D = np.outer(c, 1/c)/(dX + np.eye(N+1))
    D -= np.diag(D.sum(axis=1))
    return D, x

def cheb_rates(bordes, nu=1.0, h=1.0, n=4, N=64):
    D, x = cheb(N)              # x en [1,-1], x[0]=1
    z = h*(x+1)/2               # z[0]=h (tope), z[N]=0 (fondo)
    Dz = D*(2.0/h)
    L = -nu*(Dz@Dz)
    # incognitas: quitar el nodo del fondo (Dirichlet). Tope: DD quita, DN impone fila.
    if bordes == "noslip-noslip":
        idx = np.arange(1, N)               # interiores
        A = L[np.ix_(idx, idx)]
    else:
        idx = np.arange(0, N)               # incluye el tope z=h
        A = L[np.ix_(idx, idx)].copy()
        A[0, :] = -Dz[0, idx]               # fila del tope: dZ/dz(h)=0  (problema generalizado)
        B = np.eye(len(idx)); B[0,0] = 0.0
        lam = np.linalg.eigvals(np.linalg.solve(A - 0*B, B))  # placeholder, se resuelve abajo
        # resolver el generalizado A v = lam B v via eig generalizado
        import scipy.linalg as sla
        w = sla.eig(A, B, right=False)
        w = w[np.isfinite(w)]
        w = np.sort(w.real[w.real > 1e-8])
        return w[:n]
    w = np.sort(np.linalg.eigvals(A).real)
    return w[:n]

hdr("1. espectros: Chebyshev vs capa.py")
for b in ("noslip-noslip", "noslip-libre"):
    got = np.asarray(capa.tasas_decaimiento(b, 1.0, 1.0, 4), float)
    ref = cheb_rates(b, N=64)
    print(f"  {b:16s} capa={np.array2string(got, precision=9)}")
    print(f"  {'':16s} cheb={np.array2string(ref, precision=9)}  errrel_max={np.max(np.abs(got-ref)/ref):.3e}")

# Ruta independiente B: raices de la ecuacion secular por brentq (sin usar la formula cerrada)
hdr("1b. ecuacion secular por brentq (sin formula cerrada)")
# DD: sin(q h)=0 ; DN: cos(q h)=0.  Buscamos las raices numericamente en ventanas.
def secular_roots(f, n, hstep=0.05, qmax=60.0):
    qs = np.arange(1e-6, qmax, hstep); vals = f(qs); roots=[]
    for i in range(len(qs)-1):
        if vals[i]==0 or vals[i]*vals[i+1] < 0:
            roots.append(brentq(f, qs[i], qs[i+1], xtol=1e-14, rtol=8.9e-16))
        if len(roots) >= n: break
    return np.array(roots)
rDD = secular_roots(lambda q: np.sin(q*1.0), 4)
rDN = secular_roots(lambda q: np.cos(q*1.0), 4)
print("  DD lambda =", np.array2string(rDD**2, precision=10), " capa:", np.array2string(np.asarray(capa.tasas_decaimiento("noslip-noslip",1,1,4)), precision=10))
print("  DN lambda =", np.array2string(rDN**2, precision=10), " capa:", np.array2string(np.asarray(capa.tasas_decaimiento("noslip-libre",1,1,4)), precision=10))
print("  errrel DD:", np.max(np.abs(rDD**2 - capa.tasas_decaimiento("noslip-noslip",1,1,4))/ (rDD**2)))
print("  errrel DN:", np.max(np.abs(rDN**2 - capa.tasas_decaimiento("noslip-libre",1,1,4))/ (rDN**2)))

# escalado en nu y h
hdr("1c. escalado en nu y h")
nu, h = 3.7, 0.013
for b in ("noslip-noslip","noslip-libre"):
    got = np.asarray(capa.tasas_decaimiento(b, nu, h, 5), float)
    base = np.asarray(capa.tasas_decaimiento(b, 1.0, 1.0, 5), float)
    print(f"  {b:16s} max|got - nu/h^2*base|/... = {np.max(np.abs(got - nu*base/h**2)/(nu*base/h**2)):.3e}")

# ---------------------------------------------------------------- 2. sesgo PIV
hdr("2. u(h)/<u> del modo lento DN, por cuadratura adaptativa (quad)")
I, err = quad(lambda z: np.sin(np.pi*z/2.0), 0.0, 1.0, epsabs=1e-14, epsrel=1e-14)
ratio = np.sin(np.pi/2)/I
print(f"  integral={I!r} (+-{err:.1e}), u(h)/<u>={ratio!r}")
print(f"  capa.cociente_superficie_promedio()={capa.cociente_superficie_promedio()!r}")
print(f"  |dif| = {abs(ratio-capa.cociente_superficie_promedio()):.3e};  pi/2 = {np.pi/2!r}")
# beta
b1,_ = quad(lambda z: (np.pi/2*np.sin(np.pi*z/2))**2, 0,1, epsabs=1e-14)
b2,_ = quad(lambda z: (np.pi/2*np.sin(np.pi*z))**2, 0,1, epsabs=1e-14)
print(f"  beta_DN=<F^2>={b1!r}  beta_DD=<F^2>={b2!r}  pi^2/8={np.pi**2/8!r}")
m1,_ = quad(lambda z: np.pi/2*np.sin(np.pi*z/2),0,1,epsabs=1e-14)
m2,_ = quad(lambda z: np.pi/2*np.sin(np.pi*z),0,1,epsabs=1e-14)
print(f"  <F_DN>={m1!r}  <F_DD>={m2!r}")
print(f"  beta_s = beta*(2/pi) = {np.pi**2/8*2/np.pi!r}  vs pi/4 = {np.pi/4!r}")
