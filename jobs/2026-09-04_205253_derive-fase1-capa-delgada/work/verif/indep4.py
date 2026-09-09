import numpy as np
from scipy.integrate import quad
print("== coeficientes de expansion (eq:ddexp / eq:dnexp), NO cubiertos por checks.py ==")
h=1.0
def ZDD(n,z): return np.sin(n*np.pi*z/h)
def ZDN(m,z): return np.sin((m+0.5)*np.pi*z/h)
for nm,Z,idx in (("DD",ZDD,range(1,6)),("DN",ZDN,range(0,5))):
    idx=list(idx); G=np.zeros((len(idx),len(idx)))
    for a,i in enumerate(idx):
        for b,j in enumerate(idx):
            G[a,b]=quad(lambda z: Z(i,z)*Z(j,z),0,h,epsabs=1e-13)[0]
    print(f"  {nm}: max|G - (h/2)I| = {np.max(np.abs(G-0.5*h*np.eye(len(idx)))):.2e}  -> normalizacion 2/h correcta")
for nm,Z,idx,u0 in (("DD",ZDD,range(1,201), lambda z: z*(h-z)*(1+3*z/h)),
                    ("DN",ZDN,range(0,200), lambda z: z*(1+2*(z/h)**3))):
    z=np.linspace(0,h,4001); rec=np.zeros_like(z)
    for i in idx:
        c=2/h*quad(lambda t: u0(t)*Z(i,t),0,h,epsabs=1e-13,limit=400)[0]
        rec+=c*Z(i,z)
    print(f"  {nm}: max|reconstruccion(200 modos) - u0| (interior) = {np.max(np.abs(rec-u0(z))[10:-10]):.2e}")
print("\n== gaps espectrales ==")
print("  DN lam1-lam0 =", np.pi**2*(2.25-0.25), " 2pi^2 =", 2*np.pi**2)
print("  DD lam2-lam1 =", np.pi**2*(4-1), " 3pi^2 =", 3*np.pi**2)
print("\n== validez: Re_l*eps^2 == Re_h*eps == tau_nu/tau_H ==")
U,l,hh,nu=0.03,0.02,6e-3,1e-6
eps=hh/l; Rel=U*l/nu; Reh=U*hh/nu; tnu=hh**2/nu; tH=l/U
print(f"  eps={eps:.4g} Re_l={Rel:.4g} Re_h={Reh:.4g}")
print(f"  tau_nu/tau_H={tnu/tH:.6g}  Re_l*eps^2={Rel*eps**2:.6g}  Re_h*eps={Reh*eps:.6g}")
