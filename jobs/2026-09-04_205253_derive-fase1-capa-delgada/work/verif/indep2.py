import importlib.util, sys, warnings
from pathlib import Path
import numpy as np
from scipy.linalg import solve_banded
JOB = Path("/home/lucia/GW_AI/proyecto-final-piv/jobs/2026-09-04_205253_derive-fase1-capa-delgada")
spec = importlib.util.spec_from_file_location("capa_v", JOB/"out"/"capa.py")
capa = importlib.util.module_from_spec(spec); sys.modules["capa_v"]=capa; spec.loader.exec_module(capa)
def hdr(s): print("\n== " + s + " " + "="*max(0,(62-len(s))))

# ---- 3. alfa como tasa de decaimiento del PROMEDIO VERTICAL, via Crank-Nicolson
# Dato inicial generico (no un autovector): u0 = z*(algo) que cumple los bordes.
hdr("3. -alpha U por Crank-Nicolson desde dato inicial generico")
def cn_decay(bordes, nu=1.0, h=1.0, nz=2000, T=2.0, nt=4000):
    dz = h/nz
    m = nz-1 if bordes=="noslip-noslip" else nz
    z = np.arange(1, m+1)*dz
    # operador -nu d2/dz2 en banda
    lo = np.full(m, -1.0); di = np.full(m, 2.0); up = np.full(m, -1.0)
    if bordes=="noslip-libre": lo[-1] = -2.0   # subdiagonal en la ultima fila
    # dato inicial generico, cumple u(0)=0 (y u(h)=0 en DD)
    if bordes=="noslip-noslip": u = z*(h-z)*(1+3*z/h)
    else:                       u = z*(1+2*(z/h)**3) + 0.5*np.sin(7*np.pi*z/(2*h))
    dt = T/nt
    th = 0.5*dt*nu/dz**2
    # (I + th*A) u^{n+1} = (I - th*A) u^n, A = tridiag(lo,di,up)/1
    ab_l = np.zeros((3,m)); ab_l[0,1:] = th*up[:-1]; ab_l[1,:] = 1+th*di; ab_l[2,:-1] = th*lo[1:]
    means=[]; ts=[]
    for n in range(nt):
        r = (1-th*di)*u
        r[:-1] += -th*up[:-1]*u[1:]
        r[1:]  += -th*lo[1:]*u[:-1]
        u = solve_banded((1,1), ab_l, r)
        if n%50==49:
            zf = np.concatenate(([0.0], z)); uf = np.concatenate(([0.0], u))
            means.append(np.trapezoid(uf, zf)/h); ts.append((n+1)*dt)
    ts=np.array(ts); means=np.array(means)
    sel = (ts>0.5)&(np.abs(means)>1e-12)
    p = np.polyfit(ts[sel], np.log(np.abs(means[sel])), 1)
    return -p[0]
for b,exact in (("noslip-noslip", np.pi**2), ("noslip-libre", np.pi**2/4)):
    r = cn_decay(b)
    print(f"  {b:16s} tasa de decaimiento de <u> = {r:.9f}   capa.alfa = {capa.alfa(1.0,1.0,b):.9f}  err_rel={abs(r-capa.alfa(1.0,1.0,b))/capa.alfa(1.0,1.0,b):.2e}")

# ---- 3b. alfa desde el promedio del termino viscoso: (nu/h)[du/dz]_0^h con el perfil F
hdr("3b. alpha = -(nu/h)[dF/dz]_0^h   (ruta simbolica/numerica directa)")
import numpy as _np
eps=1e-7
for b, F, in (("noslip-noslip", lambda z: _np.pi/2*_np.sin(_np.pi*z)),
              ("noslip-libre",  lambda z: _np.pi/2*_np.sin(_np.pi*z/2))):
    d0 = (F(eps)-F(-eps))/(2*eps); dh = (F(1+eps)-F(1-eps))/(2*eps)
    alpha = -(1.0/1.0)*(dh-d0)
    print(f"  {b:16s} dF/dz(0)={d0:.9f} dF/dz(h)={dh:.9f}  alpha={alpha:.9f}  capa.alfa={capa.alfa(1,1,b):.9f}")

# ---- 4. Laplace DN: BVP tridiagonal independiente vs capa.coef_laplace_dn
hdr("4. Laplace DN: solve_banded del BVP phi''=k^2 phi vs capa (punto a punto)")
def bvp_dn(k, Lz, g0, g1, n=20001):
    z = np.linspace(0, Lz, n); dz = z[1]-z[0]
    lo=np.zeros(n); di=np.zeros(n); up=np.zeros(n); rhs=np.zeros(n, dtype=complex)
    di[0]=1.0; rhs[0]=g0                             # phi(0)=g0
    for_i = slice(1,n-1)
    lo[for_i]=1/dz**2; di[for_i]=-2/dz**2-k**2; up[for_i]=1/dz**2
    # phi'(Lz)=g1 a segundo orden: (3phi_n -4phi_{n-1}+phi_{n-2})/(2dz)=g1 -> usar banda 2? usa 1er orden mejorado
    di[-1]=1/dz; lo[-1]=-1/dz; rhs[-1]=g1 + 0.5*dz*k**2*0  # 1er orden
    ab=np.zeros((3,n),dtype=complex); ab[0,1:]=up[:-1]; ab[1,:]=di; ab[2,:-1]=lo[1:]
    return z, solve_banded((1,1), ab, rhs)
for (k,Lz,g0,g1) in [(0.0,1.0,0.7,-0.3),(1.0,1.0,1.0,0.0),(3.7,2.0,-0.4,2.1),(0.35,1.5,1.0,1.0),(8.0,1.0,0.3,-2.0)]:
    z,phi = bvp_dn(k,Lz,g0,g1)
    sol = capa.coef_laplace_dn(k,Lz,g0,g1)
    pa = np.asarray(sol.phi(z),dtype=complex)
    scale = max(np.max(np.abs(pa)),1.0)
    print(f"  k={k:6.3g} Lz={Lz:4.2g}: max|phi_capa - phi_FD|/scale = {np.max(np.abs(pa-phi))/scale:.3e}  (FD 1er orden en el borde)")

# ---- 4b. residuos de borde barriendo k, cruzando el umbral kLz=0.5
hdr("4b. residuo de borde |phi(0)-g0|,|dphi(Lz)-g1| barriendo k (umbral kLz=0.5)")
worst=0; wk=None
for Lz in (0.5,1.0,2.0,7.3):
    for k in np.concatenate([np.logspace(-14,3,400), np.linspace(0.4/Lz,0.6/Lz,201)]):
        s=capa.coef_laplace_dn(float(k),Lz,0.9,-1.7)
        e=max(abs(s.phi(0.0)-0.9)/max(1,0.9), abs(s.dphi(Lz)-(-1.7))/max(1,1.7))
        if e>worst: worst, wk = e,(k,Lz)
print(f"  residuo maximo sobre 2404 casos = {worst:.3e}  en k={wk[0]:.6g}, Lz={wk[1]}  (tol puerta 1e-10)")

# ---- 4c. intento de romper: k pequeno con Lz enorme, justo por encima del umbral
hdr("4c. INTENTO DE ROTURA: cancelacion con k pequeno y k*Lz apenas > 0.5")
for k,Lz in [(1e-10,6e9),(1e-8,6e7),(1e-6,6e5),(1e-4,6e3),(1e-2,60.0),(1e-10,5e9*0.999)]:
    s=capa.coef_laplace_dn(k,Lz,1.0,1.0)
    e0=abs(s.phi(0.0)-1.0); e1=abs(s.dphi(Lz)-1.0)
    print(f"  k={k:9.1e} Lz={Lz:9.3e} kLz={k*Lz:6.3f} rama={'exp' if k*Lz>=0.5 else 'reg'}: |phi(0)-g0|={e0:.3e} |dphi(Lz)-g1|={e1:.3e}")

# ---- 4d. sinc complejo / warnings
hdr("4d. np.sinc con argumento imaginario: warnings?")
with warnings.catch_warnings():
    warnings.simplefilter("error")
    try:
        s=capa.coef_laplace_dn(1e-3,1.0,1.0+2j,-0.3+0.8j)
        print("  ok, phi(0)=",s.phi(0.0)," dphi(Lz)=",s.dphi(1.0))
    except Exception as e:
        print("  EXCEPCION/warning:", type(e).__name__, e)
print("  np.sinc(1j*0.3/np.pi) =", np.sinc(1j*0.3/np.pi), " sinh(0.3)/0.3 =", np.sinh(0.3)/0.3)

# ---- 5. solver_1d: es realmente un integrador? comparar contra CN independiente
hdr("5. solver_1d vs Crank-Nicolson independiente (mismo dato inicial, mismo nz)")
for b in ("noslip-noslip","noslip-libre"):
    for nz in (48,96):
        z,u = capa.solver_1d(1.0,1.0,nz,0.3,b)
        # CN independiente sobre la misma malla
        dz=1.0/nz; m=nz-1 if b=="noslip-noslip" else nz
        zi=np.arange(1,m+1)*dz
        u0 = np.sin(np.pi*zi) if b=="noslip-noslip" else np.sin(0.5*np.pi*zi)
        u0=u0/np.max(np.abs(u0))
        lo=np.full(m,-1.0);di=np.full(m,2.0);up=np.full(m,-1.0)
        if b=="noslip-libre": lo[-1]=-2.0
        nt=20000; dt=0.3/nt; th=0.5*dt/dz**2
        ab=np.zeros((3,m)); ab[0,1:]=th*up[:-1]; ab[1,:]=1+th*di; ab[2,:-1]=th*lo[1:]
        uu=u0.copy()
        for _ in range(nt):
            r=(1-th*di)*uu; r[:-1]+=-th*up[:-1]*uu[1:]; r[1:]+=-th*lo[1:]*uu[:-1]
            uu=solve_banded((1,1),ab,r)
        ref = u[1:-1] if b=="noslip-noslip" else u[1:]
        print(f"  {b:16s} nz={nz:3d}: max|solver_1d - CN| = {np.max(np.abs(ref-uu)):.3e}; amp={np.max(np.abs(u)):.6f}")
