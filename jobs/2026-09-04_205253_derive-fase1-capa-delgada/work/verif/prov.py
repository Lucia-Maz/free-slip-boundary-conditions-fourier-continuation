import json,re
from pathlib import Path
JOB=Path("/home/lucia/GW_AI/proyecto-final-piv/jobs/2026-09-04_205253_derive-fase1-capa-delgada")
tex=(JOB/"out"/"notes.tex").read_text()
prov=json.loads((JOB/"out"/"provenance.json").read_text())
tags=re.findall(r"\\src\{([^}]*)\}", tex)
tags=[t for t in tags if "#1" not in t]      # descarta la definicion del macro
keys=[]
for t in tags:
    keys += [x.strip() for x in t.split(",") if x.strip()]
print(f"ocurrencias de \\src: {len(tags)}   claves citadas: {len(keys)}   claves unicas: {len(set(keys))}")
print(f"entradas en provenance.json: {len(prov)}")
falt=[k for k in keys if k not in prov]
sobra=[k for k in prov if k not in set(keys)]
print("claves citadas sin entrada:", falt or "ninguna")
print("entradas nunca citadas   :", sobra or "ninguna")
for k,v in prov.items():
    for f in ("statement","type","reproduce"):
        if f not in v or not str(v[f]).strip(): print("  INCOMPLETA:",k,f)
print("tipos:", {t: sum(1 for v in prov.values() if v['type']==t) for t in sorted({v['type'] for v in prov.values()})})

# reproduce: claves de tipo check apuntan a algo existente?
import subprocess
for k,v in prov.items():
    r=v["reproduce"]
    if r.startswith("out/checks.py"):
        fn=r.split("::")[1] if "::" in r else None
        src=(JOB/"out"/"checks.py").read_text()
        if fn and f"def {fn}" not in src: print("  REPRODUCE ROTO:",k,r)
    elif r.startswith("out/") or r.startswith("verificacion/") or r.startswith("numerico/"):
        p = JOB/r.split(",")[0] if r.startswith("out/") else Path("/home/lucia/GW_AI/proyecto-final-piv")/r.split(",")[0].split(" ")[0]
        if not p.exists(): print("  REPRODUCE ROTO (archivo):",k,p)

# ---- contrastar TODOS los decimales de notes.tex contra checks.py --values
print("\n== decimales de notes.tex vs checks.py --values ==")
out=subprocess.run(["/home/lucia/miniforge3/envs/piv-dt/bin/python",str(JOB/"out"/"checks.py"),"--values"],capture_output=True,text=True)
vals=json.loads(out.stdout)
flat={}
for k,v in vals.items():
    if isinstance(v,list):
        for i,x in enumerate(v): flat[f"{k}[{i}]"]=x
    else: flat[k]=v
# numeros literales que aparecen en la nota (>=6 digitos significativos)
lits=set(re.findall(r"(?<![\w.])(\d+\.\d{6,})", tex))
print("literales de >=6 decimales en notes.tex:", len(lits))
import math
for s in sorted(lits, key=float):
    x=float(s)
    hit=[k for k,v in flat.items() if isinstance(v,(int,float)) and (abs(v-x)<=5*10**-(len(s.split('.')[1])) )]
    # tambien derivados triviales
    extra={"pi/2":math.pi/2,"pi^2/8":math.pi**2/8,"pi/4":math.pi/4,"pi^2/4":math.pi**2/4,"pi^2":math.pi**2,
           "2pi^2":2*math.pi**2,"3pi^2":3*math.pi**2,"100(pi/2-1)":100*(math.pi/2-1),
           "alfa_h6":math.pi**2/4/6e-3**2,"alfa_h6_1e6":1e-6*math.pi**2/4/6e-3**2,"tau":1/(1e-6*math.pi**2/4/6e-3**2)}
    hit2=[k for k,v in extra.items() if abs(v-x)<=5*10**-(len(s.split('.')[1]))]
    print(f"  {s:22s} -> checks:{hit if hit else '-'}  analitico:{hit2 if hit2 else '-'}")
