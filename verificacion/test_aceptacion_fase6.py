#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Puerta de aceptación de la Fase 6: superficie libre deformable de orden N (N = 1, 2, 3).

Se escribe ANTES de la implementación y arranca en rojo. Convenciones congeladas en
verificacion/intent_fase6.txt. Ninguna comprobación contiene la fórmula de transferencia
que verifica: las referencias son (a) la geometría EXACTA evaluada directamente en la
superficie desplazada, (b) la referencia MAC lineal de la Fase 5 a otra profundidad, y
(c) un resolvedor 2D independiente de ondas viscosas no lineales con la superficie donde
realmente está (verificacion/referencia_fase6.py, congelado junto con esta puerta).

  H1  parser: fsorder 2 se acepta; 0 y 4 abortan con "fsorder"
  H2  trazas de d^m u/dz^m en la superficie (m = 0..3) contra las derivadas exactas
  H3  escalera algebraica: K, T_x, T_y y Phi de orden N contra la geometría exacta, con
      error ~ eps^N; la curvatura, ~ eps^2 (N <= 2) y eps^4 (N = 3)
  H4  regresión: fsorder = 1 sigue el camino de la Fase 5, y el camino general a N = 1
      (fsgen = 1) coincide con él
  H5  escalera dinámica de profundidad: con elevación media eta0, el modelo de orden N a
      profundidad Lz reproduce la onda lineal a profundidad Lz + eta0 (MAC) con error
      ~ eta0^N, descontado el piso de la malla
  H6  la referencia de geometría exacta se verifica a sí misma: límite lineal contra MAC,
      balance de energía no lineal, convergencia en Nz
  H7  escalera no lineal: ondas de amplitud eps contra la geometría exacta; el error del
      estado (armónicos 1-3) cae como eps^N
  H8  orden temporal con N = 2
  H9  la iteración converge sin tocar el tope, el volumen se conserva, y las corridas
      llegan a max|eta| > dz (el régimen donde la transferencia con trazas espectrales
      no convergería)
  H10 poder discriminante: una referencia con el término de pendiente tangencial con el
      signo cambiado queda lejos de SPECTER N = 2

Uso:
    /home/lucia/miniforge3/envs/piv-dt/bin/python verificacion/test_aceptacion_fase6.py
Variables de entorno: las de la Fase 2 (SPECTER_TOOLCHAIN, SPECTER_SCRATCH, SPECTER_KEEP).
"""

import math
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
from scipy.integrate import simpson

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import test_aceptacion_fase2 as F2                    # noqa: E402  (sólo lectura)
import test_aceptacion_fase5 as G5                    # noqa: E402  (sólo lectura)
import referencia_fase6 as REF                        # noqa: E402  (congelado con la puerta)

NU, GRAV, TENS = 1.0e-2, 1.0, 0.2
AMP_LIN = 1.0e-6
EPS_ONDA = (0.01, 0.02, 0.04)          # H7
ETA0 = (0.01, 0.02, 0.04)               # H5
EPS_JET = (0.2, 0.1, 0.05, 0.025)       # H3
DT, T_FIN, CSTEP = 1.0e-3, 4.0, 40      # H5, H7: salida cada 0,04
MALLA_3D = (8, 8, 128)
MALLA_2D = (16, 1, 128)
MALLA_TESTS = (16, 16, 128)


# ---------------------------------------------------------------------------------
# Arnés
# ---------------------------------------------------------------------------------

def construir(destino, target, malla):
    """Compila con la condición inicial de la Fase 5 y la malla pedida."""
    viejo = (F2.NX, F2.NY, F2.NZ, F2.NZFIS)
    F2.NX, F2.NY, F2.NZ = malla
    F2.NZFIS = F2.NZ - F2.CZ
    try:
        return G5.construir(destino, target)
    finally:
        F2.NX, F2.NY, F2.NZ, F2.NZFIS = viejo


PLANTILLA = G5.PLANTILLA.replace("fscoup = {fscoup}\n", """fscoup = {fscoup}
fsorder = {fsorder}
fsmean = {fsmean}
fstol = {fstol}
fsmaxit = {fsmaxit}
fsgen = {fsgen}
""")
assert "fsorder = {fsorder}" in PLANTILLA


def correr(bindir, etiqueta, *, fsorder=1, fsmean=0.0, fstol=1e-10, fsmaxit=30, fsgen=0,
           bcsta="noslip", bcend="freesurface", dt=DT, step=100, cstep=CSTEP, u0=0.0,
           nu=NU, vparam0=0.5, vparam1=0.0, fsgrav=GRAV, fstens=TENS, fsamp=AMP_LIN,
           fskx=1, fsky=0, fscoup=2, nproc=1, binario="HD", esperar_exito=True):
    dirrun = os.path.join(os.path.dirname(bindir), "corridas", etiqueta)
    os.makedirs(os.path.join(dirrun, "out"), exist_ok=True)
    os.makedirs(os.path.join(dirrun, "in"), exist_ok=True)
    with open(os.path.join(dirrun, "parameter.inp"), "w") as f:
        f.write(PLANTILLA.format(lz=1.0, dt=dt, step=step, cstep=cstep, u0=u0, nu=nu,
                                 vparam0=vparam0, vparam1=vparam1, bcsta=bcsta,
                                 bcend=bcend, fsgrav=fsgrav, fstens=fstens, fsamp=fsamp,
                                 fskx=fskx, fsky=fsky, fscoup=fscoup, fsorder=fsorder,
                                 fsmean=fsmean, fstol=fstol, fsmaxit=fsmaxit,
                                 fsgen=fsgen))
    enlace = os.path.join(os.path.dirname(dirrun), "tables")
    if not os.path.exists(enlace):
        os.symlink(os.path.join(os.path.dirname(bindir), "tables"), enlace)
    proc = subprocess.run(
        [os.path.join(G5.TOOLCHAIN, "bin", "mpirun"), "-np", str(nproc),
         "--oversubscribe", os.path.join(bindir, binario)],
        cwd=dirrun, env=F2._entorno(), capture_output=True, text=True, timeout=7200)
    if esperar_exito and proc.returncode != 0:
        raise RuntimeError("la corrida %r falló (rc=%d):\n%s" %
                           (etiqueta, proc.returncode, (proc.stdout + proc.stderr)[-2000:]))
    return proc.returncode, proc.stdout + proc.stderr, dirrun


def diag_orden(dirrun):
    ruta = os.path.join(dirrun, "freesurface_order_diagnostic.txt")
    if not os.path.exists(ruta):
        raise RuntimeError("no se escribió freesurface_order_diagnostic.txt en %s" % dirrun)
    d = np.atleast_2d(np.loadtxt(ruta))
    if d.shape[1] != 13:
        raise RuntimeError("freesurface_order_diagnostic.txt tiene %d columnas, no 13"
                           % d.shape[1])
    return d


def serie_modos(dirrun, dt=DT, cstep=CSTEP):
    d = diag_orden(dirrun)
    t_esp = np.arange(len(d)) * cstep * dt
    if np.max(np.abs(d[:, 0] - t_esp)) > 1e-9 * max(1.0, t_esp[-1]):
        raise RuntimeError("tiempos del diagnóstico de orden inesperados")
    return t_esp, d[:, 1] + 1j * d[:, 2], d[:, 3] + 1j * d[:, 4], d[:, 5] + 1j * d[:, 6], d


def mac_eta(K, G, L, paso, pasos):
    """Referencia MAC de la Fase 5 a profundidad L, con Richardson de dos niveles."""
    r = {n: G5._trayectoria(n, K, G, paso, pasos, L=L)[0][-1] for n in (100, 200, 400)}
    r1 = (4 * r[200] - r[100]) / 3
    r2 = (4 * r[400] - r[200]) / 3
    return ((16 * r2 - r1) / 15).real


def ordenes(errores):
    e = np.asarray(errores, dtype=float)
    return [float(np.log2(e[i + 1] / e[i])) if e[i] > 0 and e[i + 1] > 0 else float("nan")
            for i in range(len(e) - 1)]


# ---------------------------------------------------------------------------------
# Campos analíticos de H2/H3 (congelados en intent_fase6.txt) y sus derivadas exactas
# ---------------------------------------------------------------------------------
# Un término es c * cos(kx x + ky y + fase) * Z(z); Z es ('exp', a), ('cos', q, psi) o
# ('poly', coeficientes crecientes). Las derivadas son exactas.

def _dZ(Z):
    if Z[0] == "exp":
        return ("exp", Z[1]), Z[1]
    if Z[0] == "cos":
        return ("cos", Z[1], Z[2] + math.pi / 2), Z[1]
    c = Z[1]
    return ("poly", [i * c[i] for i in range(1, len(c))] or [0.0]), 1.0


def _Zval(Z, z):
    if Z[0] == "exp":
        return np.exp(Z[1] * z)
    if Z[0] == "cos":
        return np.cos(Z[1] * z + Z[2])
    return sum(c * z**i for i, c in enumerate(Z[1]))


def derivar(campo, dx=0, dy=0, dz=0):
    out = []
    for (c, kx, ky, fase, Z) in campo:
        for _ in range(dx):
            c, fase = c * kx, fase + math.pi / 2
        for _ in range(dy):
            c, fase = c * ky, fase + math.pi / 2
        for _ in range(dz):
            Z, fac = _dZ(Z)
            c = c * fac
        out.append((c, kx, ky, fase, Z))
    return out


def evaluar(campo, x, y, z):
    return sum(c * np.cos(kx * x + ky * y + fase) * _Zval(Z, z)
               for (c, kx, ky, fase, Z) in campo)


S2 = -math.pi / 2          # sin(t) = cos(t - pi/2)
U = [(-0.3, 2, 1, S2, ("poly", [1.0, 0.0, 0.3])), (0.48, 1, 1, S2, ("cos", 1.2, S2))]
V = [(0.3, 1, -1, 0.0, ("exp", 0.6)), (0.6, 2, 1, S2, ("poly", [1.0, 0.0, 0.3]))]
W = [(0.4, 1, 1, 0.0, ("cos", 1.2, 0.0)), (-0.5, 1, -1, S2, ("exp", 0.6))]
P = [(0.3, 1, 0, 0.0, ("exp", 1.1)), (0.3, 1, 0, 0.0, ("exp", -1.1)),
     (0.3, -1, 1, S2, ("exp", 0.4))]
GJ, GAMJ, NUJ = 1.3, 0.7, 0.11


def eta_jet(eps, x, y):
    e = eps * (np.cos(x) + 0.5 * np.sin(x + 2 * y))
    ex = eps * (-np.sin(x) + 0.5 * np.cos(x + 2 * y))
    ey = eps * (np.cos(x + 2 * y))
    exx = eps * (-np.cos(x) - 0.5 * np.sin(x + 2 * y))
    exy = eps * (-np.sin(x + 2 * y))
    eyy = eps * (-2.0 * np.sin(x + 2 * y))
    return e, ex, ey, exx, exy, eyy


def funcionales_exactos(eps, nf=64):
    """K, T_x, T_y, Phi, kappa EXACTOS en z = 1 + eta, en una malla fina."""
    x = 2 * np.pi * np.arange(nf) / nf
    X, Y = np.meshgrid(x, x, indexing="ij")
    e, ex, ey, exx, exy, eyy = eta_jet(eps, X, Y)
    Z = 1.0 + e
    ev = lambda f, a=0, b=0, c=0: evaluar(derivar(f, a, b, c), X, Y, Z)
    u, v, w, p = ev(U), ev(V), ev(W), ev(P)
    Sxx, Syy = 2 * ev(U, 1), 2 * ev(V, 0, 1)
    Sxy = ev(U, 0, 1) + ev(V, 1)
    Sxz, Syz = ev(U, 0, 0, 1) + ev(W, 1), ev(V, 0, 0, 1) + ev(W, 0, 1)
    Szz = 2 * ev(W, 0, 0, 1)
    s = (ex, ey)
    Sab = [[Sxx, Sxy], [Sxy, Syy]]
    Saz = [Sxz, Syz]
    K = w - u * ex - v * ey
    T = [Saz[a] - sum(Sab[a][b] * s[b] for b in range(2)) + s[a] * Szz
         - s[a] * sum(s[b] * Saz[b] for b in range(2)) for a in range(2)]
    q2 = ex**2 + ey**2
    B = (Szz - 2 * (ex * Sxz + ey * Syz)
         + sum(s[a] * s[b] * Sab[a][b] for a in range(2) for b in range(2))) / (1 + q2)
    kappa = (exx * (1 + ey**2) + eyy * (1 + ex**2) - 2 * ex * ey * exy) / (1 + q2)**1.5
    Phi = p - GJ * e + GAMJ * kappa - NUJ * B
    return {1: K, 2: T[0], 3: T[1], 4: Phi, 5: kappa}


def coef_specter(f, nx, ny):
    """Coeficiente de Fourier en la normalización de SPECTER (factor nx*ny)."""
    nf = f.shape[0]
    return np.fft.fft2(f) / nf**2 * nx * ny


# ---------------------------------------------------------------------------------
# Chequeos
# ---------------------------------------------------------------------------------

def h1_parser(ctx):
    b = ctx["b3"]
    rc, sal, _ = correr(b, "h1_orden2", fsorder=2, step=10, cstep=5, esperar_exito=False)
    if rc != 0:
        return False, "fsorder = 2 no corre:\n" + sal[-600:]
    for mal in (0, 4):
        rc, sal, _ = correr(b, "h1_orden%d" % mal, fsorder=mal, step=10, cstep=5,
                            esperar_exito=False)
        if rc == 0 or "fsorder" not in sal:
            return False, "fsorder = %d no aborta con 'fsorder' (rc=%d)" % (mal, rc)
    return True, "acepta fsorder = 2; rechaza 0 y 4"


def _correr_tests(ctx):
    if "tests_dir" not in ctx:
        _, sal, d = correr(ctx["bt"], "h23_tests", binario="TESTS", bcend="noslip",
                           step=10, cstep=5)
        ctx["tests_dir"], ctx["tests_sal"] = d, sal
    return ctx["tests_dir"]


def h2_trazas(ctx):
    d = _correr_tests(ctx)
    ruta = os.path.join(d, "fs_trazas.txt")
    if not os.path.exists(ruta):
        return False, "no se escribió fs_trazas.txt\n" + ctx["tests_sal"][-600:]
    a = np.loadtxt(ruta)
    nx, ny, _ = MALLA_TESTS
    nf = 64
    x = 2 * np.pi * np.arange(nf) / nf
    X, Y = np.meshgrid(x, x, indexing="ij")
    umbral = {0: 1e-9, 1: 1e-9, 2: 1e-6, 3: 1e-4}
    detalle, ok = [], True
    for m in range(4):
        ex = coef_specter(evaluar(derivar(U, 0, 0, m), X, Y, np.ones_like(X)), nx, ny)
        filas = a[a[:, 0].astype(int) == m]
        esc = np.max(np.abs(ex))
        err = 0.0
        for fila in filas:
            i, j = int(fila[1]), int(fila[2])
            ky = j - 1 if j - 1 <= ny // 2 else j - 1 - ny
            ref = ex[i - 1, ky % nf]
            err = max(err, abs(fila[3] + 1j * fila[4] - ref) / esc)
        ok = ok and len(filas) > 0 and err < umbral[m]
        detalle.append("m=%d: %.1e (umbral %.0e)" % (m, err, umbral[m]))
    return ok, "; ".join(detalle)


def h3_escalera_algebraica(ctx):
    d = _correr_tests(ctx)
    ruta = os.path.join(d, "fs_jet.txt")
    if not os.path.exists(ruta):
        return False, "no se escribió fs_jet.txt\n" + ctx["tests_sal"][-600:]
    a = np.loadtxt(ruta)
    nx, ny, _ = MALLA_TESTS
    exactos = {eps: {k: coef_specter(f, nx, ny) for k, f in funcionales_exactos(eps).items()}
               for eps in EPS_JET}
    nf = 64
    detalle, ok = [], True
    for N in (1, 2, 3):
        for fid, nombre in ((1, "K"), (2, "Tx"), (3, "Ty"), (4, "Phi"), (5, "kappa")):
            errs = []
            for eps in EPS_JET:
                sel = a[(a[:, 0].astype(int) == N) & (np.abs(a[:, 1] - eps) < 1e-12)
                        & (a[:, 2].astype(int) == fid)]
                if len(sel) == 0:
                    return False, "faltan filas N=%d eps=%g id=%d en fs_jet.txt" % (N, eps, fid)
                ex = exactos[eps][fid]
                esc = np.max(np.abs(ex))
                e = 0.0
                for fila in sel:
                    i, j = int(fila[3]), int(fila[4])
                    ky = j - 1 if j - 1 <= ny // 2 else j - 1 - ny
                    e = max(e, abs(fila[5] + 1j * fila[6] - ex[i - 1, ky % nf]) / esc)
                errs.append(e)
            # eps decrece en EPS_JET: el orden es log2(e(eps)/e(eps/2))
            p = [-q for q in ordenes(errs)]
            if fid == 5:
                esperado = 2 if N <= 2 else 4
            else:
                esperado = N
            med = float(np.median(p))
            bien = (esperado - 0.3) <= med <= (esperado + 0.6)
            ok = ok and bien
            detalle.append("N=%d %s: orden %.2f (esperado %d; errores %.1e..%.1e)%s"
                           % (N, nombre, med, esperado, errs[0], errs[-1],
                              "" if bien else " <- NO"))
    return ok, "; ".join(detalle)


def h4_regresion(ctx):
    b = ctx["b3"]
    pasos = int(round(2.0 / DT))
    kw = dict(step=pasos, fskx=1, fsky=1, fstens=TENS)
    _, _, d1 = correr(b, "h4_camino_f5", fsorder=1, fsgen=0, **kw)
    _, _, d2 = correr(b, "h4_camino_general", fsorder=1, fsgen=1, **kw)
    assert G5.AMP == AMP_LIN          # G5.serie_eta normaliza por G5.AMP
    t, r1, _ = G5.serie_eta(d1, DT, CSTEP)
    _, r2, _, _, _ = serie_modos(d2)
    n = min(len(r1), len(r2))
    dif = float(np.max(np.abs(r1[:n] - r2[:n] / AMP_LIN)))
    ref = mac_eta(math.sqrt(2.0), GRAV + TENS * 2.0, 1.0, CSTEP * DT, len(t) - 1)
    err = float(np.max(np.abs(r1.real - ref)))
    ok = dif < 1e-9 and err < 1e-4
    return ok, ("camino de la Fase 5 contra MAC %.1e; camino general N=1 contra el de la "
                "Fase 5 %.1e" % (err, dif))


def h5_escalera_profundidad(ctx):
    b = ctx["b3"]
    pasos = int(round(T_FIN / DT)) + 1      # SPECTER escribe al empezar el paso
    K = math.sqrt(2.0)
    G = GRAV + TENS * K * K
    kw = dict(step=pasos, fskx=1, fsky=1, fstens=TENS)
    _, _, d0 = correr(b, "h5_piso", fsorder=1, **kw)
    t, r0, _ = G5.serie_eta(d0, DT, CSTEP)
    piso = r0.real - mac_eta(K, G, 1.0, CSTEP * DT, len(t) - 1)   # r0 ya normalizada
    detalle, ok = [], True
    errores = {}
    for N in (1, 2, 3):
        e = []
        for eta0 in ETA0:
            _, _, d = correr(b, "h5_N%d_%g" % (N, eta0), fsorder=N, fsmean=eta0,
                             fsgen=1, **kw)
            _, r, _, _, dd = serie_modos(d)
            ref = mac_eta(K, G, 1.0 + eta0, CSTEP * DT, len(t) - 1)
            e.append(float(np.max(np.abs((r.real / AMP_LIN - ref) - piso))))
            ctx.setdefault("h9_diags", []).append(("h5 N=%d eta0=%g" % (N, eta0), dd))
        errores[N] = e
        p = ordenes(e)
        bien = all((N - 0.4) <= q <= (N + 0.6) for q in p)
        ok = ok and bien
        detalle.append("N=%d: errores %s, órdenes %s%s" % (
            N, " ".join("%.1e" % x for x in e), " ".join("%.2f" % q for q in p),
            "" if bien else " <- NO"))
    ok = ok and errores[3][-1] < errores[2][-1] < errores[1][-1]
    return ok, "; ".join(detalle) + ("; piso de malla %.1e" % float(np.max(np.abs(piso))))


def _referencia(ctx, eps, variante=1.0, Nz=40):
    clave = ("ref", eps, variante, Nz)
    if clave not in ctx:
        tiempos = np.arange(0, int(round(T_FIN / (CSTEP * DT))) + 1) * CSTEP * DT
        R = REF.ReferenciaExacta(Nx=MALLA_2D[0], Nz=Nz, nu=NU, g=GRAV, gam=TENS,
                                 signo_pendiente_tangencial=variante)
        ctx[clave] = R.integrar(eps, tiempos, dt_max=0.02, energia=True)
    return ctx[clave]


def h6_referencia(ctx):
    lin = _referencia(ctx, AMP_LIN)
    mac = mac_eta(1.0, GRAV + TENS, 1.0, CSTEP * DT, len(lin["t"]) - 1)
    e_lin = float(np.max(np.abs(lin["eta1"].real / AMP_LIN - mac)))
    nl = _referencia(ctx, EPS_ONDA[-1])
    E, D, t = nl["E"], nl["D"], nl["t"]
    bal = abs(E[0] - E[-1] - simpson(D, x=t)) / E[0]
    nl32 = _referencia(ctx, EPS_ONDA[-1], Nz=32)
    conv = max(float(np.max(np.abs(nl["eta%d" % k] - nl32["eta%d" % k]))) / EPS_ONDA[-1]
               for k in (1, 2, 3))
    ok = e_lin < 3e-8 and bal < 2e-6 and conv < 2e-7
    return ok, ("límite lineal contra MAC %.1e; balance de energía (eps=%g) %.1e; "
                "Nz 32 contra 40 %.1e" % (e_lin, EPS_ONDA[-1], bal, conv))


def _error_estado(S, R, piso1, eps):
    """max_t del error de los armónicos 1-3 de eta, relativo a eps; al modo 1 se le
    descuenta el piso lineal de la malla (piso1 * eps)."""
    t, s1, s2, s3 = S
    n = len(t)
    if len(R["eta1"]) != n or (np.ndim(piso1) and len(piso1) != n):
        raise RuntimeError("series de largo distinto: SPECTER %d, referencia %d"
                           % (n, len(R["eta1"])))
    e = np.sqrt(np.abs(s1 - R["eta1"] - piso1 * eps)**2 + np.abs(s2 - R["eta2"])**2
                + np.abs(s3 - R["eta3"])**2)
    return float(np.max(e)) / eps


def h7_escalera_no_lineal(ctx):
    b = ctx["b2"]
    pasos = int(round(T_FIN / DT)) + 1      # SPECTER escribe al empezar el paso
    kw = dict(step=pasos, fskx=1, fsky=0, fstens=TENS)
    _, _, d0 = correr(b, "h7_piso", fsorder=1, fsgen=1, fsamp=AMP_LIN, **kw)
    S0 = serie_modos(d0)
    R0 = _referencia(ctx, AMP_LIN)
    piso1 = (S0[1] - R0["eta1"]) / AMP_LIN
    detalle, ok = [], True
    E = {}
    for N in (1, 2, 3):
        e = []
        for eps in EPS_ONDA:
            _, _, d = correr(b, "h7_N%d_%g" % (N, eps), fsorder=N, fsgen=1, fsamp=eps, **kw)
            S = serie_modos(d)
            ctx["h7_%d_%g" % (N, eps)] = S
            ctx.setdefault("h9_diags", []).append(("h7 N=%d eps=%g" % (N, eps), S[4]))
            e.append(_error_estado(S[:4], _referencia(ctx, eps), piso1, eps))
        E[N] = e
        p = ordenes(e)
        bien = (N - 0.4) <= p[-1] <= (N + 0.6) and p[0] >= N - 0.6
        ok = ok and bien
        detalle.append("N=%d: errores %s, órdenes %s%s" % (
            N, " ".join("%.1e" % x for x in e), " ".join("%.2f" % q for q in p),
            "" if bien else " <- NO"))
    ctx["h7_E"] = E
    jer = E[2][-1] < 0.3 * E[1][-1] and E[3][-1] < 0.3 * E[2][-1]
    ok = ok and jer
    return ok, "; ".join(detalle) + ("" if jer else "; la jerarquía N=1>N=2>N=3 no se cumple")


def h8_orden_temporal(ctx):
    b = ctx["b2"]
    series = []
    for dt, cs in ((1e-3, 40), (5e-4, 80), (2.5e-4, 160)):
        _, _, d = correr(b, "h8_%g" % dt, fsorder=2, fsgen=1, fsamp=0.02, dt=dt,
                         step=int(round(2.0 / dt)), cstep=cs, fskx=1, fsky=0)
        t, s1, s2, s3, _ = serie_modos(d, dt, cs)
        series.append(np.concatenate([s1, s2, s3]))
    n = min(len(s) for s in series)
    e1 = float(np.max(np.abs(series[0][:n] - series[1][:n])))
    e2 = float(np.max(np.abs(series[1][:n] - series[2][:n])))
    p = math.log2(e1 / e2) if e2 > 0 else float("inf")
    return 1.6 < p < 2.6, "diferencias %.2e, %.2e: orden %.2f" % (e1, e2, p)


def h9_iteracion(ctx):
    diags = ctx.get("h9_diags", [])
    if not diags:
        return False, "no hay corridas de H5/H7 para revisar"
    peor_pasadas, peor_res, peor_media, max_rel = 0, 0.0, 0.0, 0.0
    for nombre, d in diags:
        peor_pasadas = max(peor_pasadas, int(np.max(d[:, 10])))
        peor_res = max(peor_res, float(np.max(d[:, 11])))
        escala = max(float(np.max(np.abs(d[:, 1]))), float(np.max(np.abs(d[:, 7]))), 1e-300)
        peor_media = max(peor_media, float(np.max(np.abs(d[:, 7] - d[0, 7]))) / escala)
        max_rel = max(max_rel, float(np.max(d[:, 12])))
    ok = peor_pasadas < 30 and peor_res <= 1e-9 and peor_media < 1e-10 and max_rel > 1.0
    return ok, ("pasadas máximas %d (tope 30); residuo final máximo %.1e; deriva de la "
                "media %.1e; max|eta|/dz alcanzado %.2f" % (peor_pasadas, peor_res,
                                                           peor_media, max_rel))


def h10_discriminante(ctx):
    eps = EPS_ONDA[-1]
    S = ctx.get("h7_2_%g" % eps)
    if S is None or "h7_E" not in ctx:
        return False, "H7 no produjo la corrida N=2"
    mala = _referencia(ctx, eps, variante=-1.0)
    buena = _referencia(ctx, eps)
    dev_mala = _error_estado(S[:4], mala, 0.0, eps)
    dev_buena = _error_estado(S[:4], buena, 0.0, eps)
    ok = dev_mala > 10 * max(dev_buena, ctx["h7_E"][2][-1])
    return ok, ("SPECTER N=2 contra la referencia con el signo de la pendiente tangencial "
                "cambiado: %.1e; contra la correcta: %.1e" % (dev_mala, dev_buena))


CHEQUEOS = [
    ("H1", "parser de fsorder", h1_parser),
    ("H2", "trazas d^m u/dz^m en la superficie", h2_trazas),
    ("H3", "escalera algebraica contra la geometría exacta", h3_escalera_algebraica),
    ("H4", "regresión del orden 1", h4_regresion),
    ("H5", "escalera de profundidad contra MAC", h5_escalera_profundidad),
    ("H6", "la referencia de geometría exacta", h6_referencia),
    ("H7", "escalera no lineal contra la geometría exacta", h7_escalera_no_lineal),
    ("H8", "orden temporal con N = 2", h8_orden_temporal),
    ("H9", "iteración, volumen y régimen eta > dz", h9_iteracion),
    ("H10", "poder discriminante", h10_discriminante),
]


def main():
    if not os.path.isdir(G5.ARBOL):
        print("[ERROR] no existe el árbol de trabajo %s" % G5.ARBOL)
        return 2
    scratch = os.environ.get("SPECTER_SCRATCH")
    temporal = scratch is None
    if temporal:
        scratch = tempfile.mkdtemp(prefix="fase6_")
    os.makedirs(scratch, exist_ok=True)
    print("=" * 78)
    print("Puerta de aceptación — Fase 6: superficie deformable de orden N")
    print("  árbol   : %s" % G5.ARBOL)
    print("  scratch : %s" % scratch)
    print("=" * 78, flush=True)
    ctx = {}
    try:
        print("\n[build] 3D %s, 2D %s, tests %s ..." % (MALLA_3D, MALLA_2D, MALLA_TESTS),
              flush=True)
        ctx["b3"] = construir(os.path.join(scratch, "m3d"), "main", MALLA_3D)
        ctx["b2"] = construir(os.path.join(scratch, "m2d"), "main", MALLA_2D)
        ctx["bt"] = construir(os.path.join(scratch, "tests"), "tests", MALLA_TESTS)
    except Exception as exc:                                  # noqa: BLE001
        print("\n[ERROR] %s\n\nRESULTADO: 0/%d — no se pudo compilar" % (exc, len(CHEQUEOS)))
        return 1
    resultados = []
    for clave, titulo, fn in CHEQUEOS:
        print("\n--- %s  %s" % (clave, titulo), flush=True)
        try:
            ok, msg = fn(ctx)
        except Exception as exc:                              # noqa: BLE001
            ok, msg = False, "excepción: %s" % exc
        print("    %s  %s" % ("PASA" if ok else "FALLA", msg), flush=True)
        resultados.append((clave, ok))
    n = sum(ok for _, ok in resultados)
    print("\n" + "=" * 78)
    print("RESULTADO: %d/%d  (%s)" % (n, len(resultados), " ".join(
        "%s:%s" % (c, "ok" if ok else "NO") for c, ok in resultados)))
    print("=" * 78)
    if temporal and os.environ.get("SPECTER_KEEP") != "1":
        shutil.rmtree(scratch, ignore_errors=True)
    else:
        print("scratch conservado en %s" % scratch)
    return 0 if n == len(resultados) else 1


if __name__ == "__main__":
    sys.exit(main())
