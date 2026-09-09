import importlib.util, sys
from pathlib import Path
import numpy as np
JOB=Path("/home/lucia/GW_AI/proyecto-final-piv/jobs/2026-09-04_205253_derive-fase1-capa-delgada")
spec=importlib.util.spec_from_file_location("capa_v", JOB/"out"/"capa.py")
capa=importlib.util.module_from_spec(spec); sys.modules["capa_v"]=capa; spec.loader.exec_module(capa)

print("== barrido dirigido: rama exponencial con Lz grande (cancelacion C1*r + C2) ==")
worst=(0,None)
rows=[]
for Lz in (1e2,1e4,1e6,1e8,1e10):
    for kLz in (0.5,0.5000001,0.51,0.6,0.8,1.0,1.5,2.0,3.0,5.0,10.0,20.0):
        k=kLz/Lz
        for g0,g1 in ((1.0,1.0),(1.0,-1.0),(1e-3,1.0),(1.0,7.3),(0.5,-1.0)):
            s=capa.coef_laplace_dn(k,Lz,g0,g1)
            e0=abs(s.phi(0.0)-g0)/max(1.0,abs(g0)); e1=abs(s.dphi(Lz)-g1)/max(1.0,abs(g1))
            e=max(e0,e1)
            if e>worst[0]: worst=(e,(k,Lz,kLz,g0,g1,e0,e1))
print(f"  peor residuo relativo (norma de la puerta) = {worst[0]:.3e}")
if worst[1]:
    k,Lz,kLz,g0,g1,e0,e1=worst[1]
    print(f"    k={k:.6g} Lz={Lz:.6g} kLz={kLz} g0={g0} g1={g1}  e0={e0:.3e} e1={e1:.3e}")
print(f"  tolerancia de la puerta: 1e-10 -> {'SUPERADA (fallaria)' if worst[0]>1e-10 else 'dentro de tolerancia'}")

print("\n== overflow/underflow: kLz extremos ==")
for k,Lz in ((350.0,2.0),(1e4,1.0),(1e5,10.0),(1e8,1.0),(700.0,1.0)):
    s=capa.coef_laplace_dn(k,Lz,0.5,-1.0)
    print(f"  k={k:9.3g} Lz={Lz:5.3g} kLz={k*Lz:10.4g}: phi(0)={s.phi(0.0):.17g} dphi(Lz)={s.dphi(Lz):.17g} C1={s.C1:.4g} C2={s.C2:.4g} finito={np.isfinite([s.phi(0.0),s.dphi(Lz)]).all()}")

print("\n== consistencia phi/dphi: dphi vs derivada numerica de phi ==")
for k,Lz,g0,g1 in ((0.0,1.0,0.7,-0.3),(0.2,1.0,1.0,1.0),(3.7,2.0,-0.4,2.1),(50.0,1.0,1.0,1.0)):
    s=capa.coef_laplace_dn(k,Lz,g0,g1)
    z=np.linspace(1e-3*Lz,Lz-1e-3*Lz,7); dz=1e-6*Lz
    num=(np.asarray(s.phi(z+dz))-np.asarray(s.phi(z-dz)))/(2*dz)
    ana=np.asarray(s.dphi(z))
    # tambien: phi'' - k^2 phi = 0
    d2=(np.asarray(s.phi(z+dz))-2*np.asarray(s.phi(z))+np.asarray(s.phi(z-dz)))/dz**2
    sc=max(np.max(np.abs(ana)),1.0)
    print(f"  k={k:5.3g}: max|dphi-num|/sc={np.max(np.abs(ana-num))/sc:.2e}   max|phi''-k^2phi|={np.max(np.abs(d2-k**2*np.asarray(s.phi(z)))):.2e}")

print("\n== rama regular vs rama exponencial: continuidad en kLz=0.5 ==")
for Lz in (1.0,2.0):
    for k in (0.4999999/Lz, 0.5/Lz, 0.5000001/Lz):
        s=capa.coef_laplace_dn(k,Lz,0.9,-1.7)
        z=np.linspace(0,Lz,5)
        print(f"  Lz={Lz} kLz={k*Lz:.7f} rama={'reg' if k*Lz<0.5 else 'exp'}: phi={np.array2string(np.asarray(s.phi(z)),precision=12)}")

print("\n== tasas_decaimiento: n grande y overflow ==")
t=capa.tasas_decaimiento("noslip-libre",1.0,1.0,2000)
print("  n=2000 ok, ultima =", t[-1], " creciente:", bool(np.all(np.diff(t)>0)))
try:
    t=capa.tasas_decaimiento("noslip-libre",1e300,1e-300,3); print("  nu=1e300,h=1e-300 ->", t)
except Exception as e: print("  nu=1e300,h=1e-300 ->", type(e).__name__, e)

print("\n== solver_1d: t_final=0 devuelve el dato inicial? y bordes ==")
for b in ("noslip-noslip","noslip-libre"):
    z,u=capa.solver_1d(1.0,1.0,20,0.0,b)
    print(f"  {b:16s} t=0: amp={np.max(np.abs(u)):.6f} u[0]={u[0]:.1e} u[-1]={u[-1]:.6f} len(z)={len(z)}")
    z,u=capa.solver_1d(1.0,1.0,20,0.0,b)
    ana = np.sin(np.pi*z) if b=="noslip-noslip" else np.sin(np.pi*z/2)
    print(f"  {'':16s} max|u(0)-sin| = {np.max(np.abs(u-ana)):.2e}")
