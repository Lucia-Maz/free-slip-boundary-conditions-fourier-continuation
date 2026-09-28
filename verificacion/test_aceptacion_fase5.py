#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Puerta de aceptación de la Fase 5: superficie libre DEFORMABLE (linealizada) en SPECTER.

Se escribe ANTES de la implementación y arranca en rojo. Mismo criterio que las Fases 1
y 2: la puerta no contiene la fórmula analítica del resultado que verifica. La
referencia es un problema de Stokes lineal por modo horizontal, discretizado acá mismo
en diferencias finitas sobre una malla escalonada (MAC), integrado en el tiempo con la
exponencial de la matriz, y extrapolado con Richardson. Antes de usarla se verifica a sí
misma (F3): su balance de energía cierra y converge.

Los ocho chequeos:

  F1  el parser acepta "freesurface" arriba, la rechaza abajo, y sigue abortando
      ante una cadena desconocida
  F2  la rama Neumann(z=0) - Dirichlet(z=Lz) de laplace_z (test unitario)
  F3  la referencia es confiable: balance de energía con gravedad y capilaridad, y
      convergencia de Richardson
  F4  la trayectoria eta(t) de SPECTER coincide con la referencia en dos casos (modo
      en x solo, gravedad; modo oblicuo, gravedad + capilaridad), con 1 y 2 procesos
  F5  la tensión tangencial COMPLETA se anula en la superficie a nivel de
      reconstrucción y sin depender de dt; el acoplamiento rezagado (control
      negativo) deja un residuo que escala con dt
  F6  orden temporal observado ~2 (ORD=2)
  F7  el modo de corte horizontal (k=0) decae con la tasa de la Fase 2: la superficie
      deformable no altera la física que no la mueve
  F8  poder discriminante: cuatro modelos físicamente equivocados (sin la tensión
      normal viscosa, sin capilaridad, capilaridad con el signo cambiado, y la
      condición plana de freeslip reusada) quedan todos lejos de lo que mide SPECTER

Uso:
    /home/lucia/miniforge3/envs/piv-dt/bin/python verificacion/test_aceptacion_fase5.py

Mismas variables de entorno que la Fase 2 (SPECTER_TOOLCHAIN, SPECTER_SCRATCH,
SPECTER_KEEP). Convenciones congeladas en verificacion/intent_fase5.txt.
"""

import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
from scipy.integrate import simpson
from scipy.linalg import expm, solve

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import test_aceptacion_fase2 as F2                     # noqa: E402  (sólo lectura)

# --------------------------------------------------------------------------------
# Configuración
# --------------------------------------------------------------------------------

ARBOL = F2.ARBOL
TOOLCHAIN = F2.TOOLCHAIN
NX, NY, NZ, CZ, OZ, ORD = F2.NX, F2.NY, F2.NZ, F2.CZ, F2.OZ, F2.ORD

LZ = 1.0
NU = 1.0e-2
GRAV = 1.0
AMP = 1.0e-6          # amplitud de eta: el término no lineal del volumen queda en ~1e-6

CADENA = "freesurface"

# Casos de F4: (etiqueta, fskx, fsky, sigma/rho)
CASOS = [("A_x_gravedad", 1, 0, 0.0),
         ("B_oblicuo_capilar", 1, 1, 0.2)]
DT_F4, T_F4, CSTEP_F4 = 5.0e-4, 10.0, 50

# Tolerancia de F4, sobre eta(t)/A. Presupuesto de error estimado ANTES de correr:
# tiempo (RK2, (omega dt)^2 omega T / 6) ~ 3e-7; espacio (FC-Gram de orden 5, capa
# límite de fondo con ~15 puntos) ~ 1e-6; término no lineal del volumen ~ A k ~ 1e-6;
# referencia (Richardson, medido en F3) ~ 1e-7. Suma ~ 3e-6. El umbral es ~30 veces eso.
# Que no sea laxo lo demuestra F8: el más cercano de los cuatro modelos equivocados queda
# a 2.5e-2, 250 veces el umbral.
TOL_F4 = 1.0e-4


# --------------------------------------------------------------------------------
# Referencia: Stokes lineal por modo horizontal, diferencias finitas MAC
# --------------------------------------------------------------------------------

def _mac(n, K, L, nu, G, visc=True, tension="completa"):
    """Devuelve J con dy/dt = J y, y = (w_1..w_n, eta), y la matriz U: y -> u.

    Malla escalonada: w en las caras z_j = j h (w_0 = 0 en el fondo), u y p en los
    centros. Continuidad i K u + w' = 0 en cada centro; momento x en los centros; momento
    z en las caras. p se elimina con el momento x. Fondo: u = 0 (fantasma u_0 = -u_1).
    Superficie: tensión tangencial u' + i K w = 0 (fantasma), tensión normal
    p = G eta + 2 nu w' con w' = -i K u, cinemática eta' = w_n. Segundo orden en h.

    visc=False quita el término 2 nu w' de la tensión normal; tension="plana" reemplaza
    la tensión tangencial por u' = 0 (la condición de la superficie plana). Son las
    variantes equivocadas de F8.
    """
    h = L / n
    m = n + 1
    E = np.eye(m)
    U = np.zeros((n, m), complex)
    for j in range(n):
        U[j, j] = 1j / (K * h)
        if j > 0:
            U[j, j - 1] = -1j / (K * h)
    u_fondo = -U[0]
    u_tope = U[n - 1].copy()
    if tension == "completa":
        u_tope[n - 1] += -1j * K * h
    Lu = np.zeros((n, m), complex)
    for j in range(n):
        izq = U[j - 1] if j > 0 else u_fondo
        der = U[j + 1] if j < n - 1 else u_tope
        Lu[j] = (der - 2 * U[j] + izq) / h**2 - K * K * U[j]
    P0 = nu * Lu / (1j * K)            # p = P0 y + P1 dy/dt
    P1 = -U / (1j * K)
    M = np.zeros((m, m), complex)
    A = np.zeros((m, m), complex)
    for j in range(n - 1):
        wm = E[j - 1] if j > 0 else np.zeros(m)
        lap = (E[j + 1] - 2 * E[j] + wm) / h**2 - K * K * E[j]
        M[j] = E[j] + (P1[j + 1] - P1[j]) / h
        A[j] = -(P0[j + 1] - P0[j]) / h + nu * lap
    # cara superior: w_n' = -p'(Lz) + nu (w'' - K^2 w)(Lz); p' con un esquema
    # unilateral de segundo orden (8 pL - 9 p_n + p_{n-1}) / (3h); w'' = -i K u'
    uL = U[n - 1].copy()
    if tension == "completa":
        uL[n - 1] += -1j * K * h / 2
        wzz = -K * K * E[n - 1]
    else:
        wzz = 0 * E[n - 1]
    pL = (-2j * nu * K * uL) if visc else np.zeros(m, complex)
    pL = pL.astype(complex)
    pL[n] += G
    M[n - 1] = E[n - 1] + (-9 * P1[n - 1] + P1[n - 2]) / (3 * h)
    A[n - 1] = (-(8 * pL - 9 * P0[n - 1] + P0[n - 2]) / (3 * h)
                + nu * (wzz - K * K * E[n - 1]))
    M[n] = E[n]
    A[n] = E[n - 1]
    return solve(M, A), U, h


def _trayectoria(n, K, G, paso, pasos, nu=NU, L=LZ, **kw):
    """y(t_k), t_k = k*paso, desde el reposo con eta = 1."""
    J, U, h = _mac(n, K, L, nu, G, **kw)
    P = expm(J * paso)
    y = np.zeros(n + 1, complex)
    y[-1] = 1.0
    Y = np.empty((n + 1, pasos + 1), complex)
    Y[:, 0] = y
    for k in range(pasos):
        y = P @ y
        Y[:, k + 1] = y
    return Y, U, h


def referencia_eta(K, G, paso, pasos, **kw):
    """eta(t)/eta(0) con Richardson de dos niveles sobre n = 100, 200, 400. Devuelve
    (ref, estimación de su error)."""
    r = {n: _trayectoria(n, K, G, paso, pasos, **kw)[0][-1] for n in (100, 200, 400)}
    r1 = (4 * r[200] - r[100]) / 3
    r2 = (4 * r[400] - r[200]) / 3
    ref = (16 * r2 - r1) / 15
    return ref, float(np.max(np.abs(r2 - r1)))


def balance_energia(n, K, G, paso, pasos):
    """|E(0) - E(T) - integral de la disipación| / E(0) para la trayectoria discreta.
    E incluye la cinética, (G/2)|eta|^2 = gravedad + capilaridad linealizada; la
    disipación es (nu/2) S:S con la tensión de la pared y la de la superficie (nula)."""
    Y, U, h = _trayectoria(n, K, G, paso, pasos)
    w, eta, u = Y[:n], Y[n], U @ Y
    peso = np.ones(n)
    peso[-1] = 0.5
    E = 0.5 * h * (np.sum(abs(u)**2, 0) + np.sum(peso[:, None] * abs(w)**2, 0)) \
        + 0.5 * G * abs(eta)**2
    wf = np.vstack([np.zeros(w.shape[1]), w])
    wz = (wf[1:] - wf[:-1]) / h
    centros = 2 * K * K * abs(u)**2 + 2 * abs(wz)**2
    caras = abs((u[1:] - u[:-1]) / h + 1j * K * w[:-1])**2
    fondo = abs(2 * u[0] / h)**2
    disip = NU * h * (np.sum(centros, 0) + np.sum(caras, 0) + 0.5 * fondo)
    return abs(E[0] - E[-1] - simpson(disip, dx=paso)) / E[0]


# --------------------------------------------------------------------------------
# Arnés
# --------------------------------------------------------------------------------

IC_FS = r"""
! Condición inicial de la puerta de la Fase 5.
!   vparam1 = 0 : fluido en reposo; la elevación inicial la pone &freesurface.
!   vparam1 = 1 : modo solenoidal Psi = cos(x) sin(pi z / Lz), el de V3 de la Fase 2.
!   vparam1 = 2 : modo de corte puro vx = u0 sin(vparam0 pi z / Lz), sin estructura
!                 horizontal (k = 0).
      vx = 0; vy = 0; vz = 0

      IF (ista .eq. 1) THEN
         IF ( int(vparam1) .eq. 1 ) THEN
            DO k = 1,nz-Cz
               vx(k,1,2) = u0*nx*ny*(pi/Lz)*COS(pi*z(k)/Lz)
               vz(k,1,2) = -im*kx(2)*u0*nx*ny*SIN(pi*z(k)/Lz)
            ENDDO
         ELSEIF ( int(vparam1) .eq. 2 ) THEN
            DO k = 1,nz-Cz
               vx(k,1,1) = u0*nx*ny*SIN(vparam0*pi*z(k)/Lz)
            ENDDO
         ENDIF
      ENDIF

      CALL fftp1d_real_to_complex_z(planfc,vx,MPI_COMM_WORLD)
      CALL fftp1d_real_to_complex_z(planfc,vy,MPI_COMM_WORLD)
      CALL fftp1d_real_to_complex_z(planfc,vz,MPI_COMM_WORLD)
"""


def construir(destino, target="main"):
    if os.path.exists(destino):
        shutil.rmtree(destino)
    shutil.copytree(ARBOL, destino,
                    ignore=shutil.ignore_patterns(".git", "*.o", "*.mod"))
    F2._parchear_makefile_in(os.path.join(destino, "src", "Makefile.in"))
    with open(os.path.join(destino, "src", "initialv.f90"), "w") as f:
        f.write(IC_FS)
    src = os.path.join(destino, "src")
    env = F2._entorno()
    env["PWD"] = src                    # GHOME = $(PWD), ver CLAUDE.md
    log = subprocess.run(["make", target], cwd=src, env=env,
                         capture_output=True, text=True)
    binario = os.path.join(destino, "bin", "HD" if target == "main" else "TESTS")
    if log.returncode != 0 or not os.path.exists(binario):
        raise RuntimeError("la compilación de '%s' falló:\n%s"
                           % (target, (log.stdout + log.stderr)[-3000:]))
    return os.path.join(destino, "bin")


PLANTILLA = """&status
idir = "in/"
odir = "out/"
tdir = "../tables/"
stat = 0
mult = 1
bench = 0
outs = 0
/
&boxparams
Lx = 1.00
Ly = 1.00
Lz = {lz}
/
&parameter
dt = {dt}
step = {step}
tstep = 1000000
cstep = {cstep}
seed = 1000
/
&velocity
f0 = 0.0
u0 = {u0}
kdn = 1.00
kup = 2.00
nu = {nu}
fparam0 = 0.0
vparam0 = {vparam0}
vparam1 = {vparam1}
/
&velbound
vbczsta = "{bcsta}"
vbczend = "{bcend}"
vxzsta  = 0.0
vyzsta  = 0.0
vxzend  = 0.0
vyzend  = 0.0
/
&freesurface
fsgrav = {fsgrav}
fstens = {fstens}
fsamp  = {fsamp}
fskx   = {fskx}
fsky   = {fsky}
fscoup = {fscoup}
/
"""


def correr(bindir, etiqueta, *, bcsta="noslip", bcend=CADENA, dt=DT_F4, step=100,
           cstep=CSTEP_F4, u0=0.0, nu=NU, vparam0=0.5, vparam1=0.0, fsgrav=GRAV,
           fstens=0.0, fsamp=AMP, fskx=1, fsky=0, fscoup=2, nproc=1, binario="HD",
           esperar_exito=True):
    dirrun = os.path.join(os.path.dirname(bindir), "corridas", etiqueta)
    os.makedirs(os.path.join(dirrun, "out"), exist_ok=True)
    os.makedirs(os.path.join(dirrun, "in"), exist_ok=True)
    with open(os.path.join(dirrun, "parameter.inp"), "w") as f:
        f.write(PLANTILLA.format(lz=LZ, dt=dt, step=step, cstep=cstep, u0=u0, nu=nu,
                                 vparam0=vparam0, vparam1=vparam1, bcsta=bcsta,
                                 bcend=bcend, fsgrav=fsgrav, fstens=fstens,
                                 fsamp=fsamp, fskx=fskx, fsky=fsky, fscoup=fscoup))
    enlace = os.path.join(os.path.dirname(dirrun), "tables")
    if not os.path.exists(enlace):
        os.symlink(os.path.join(os.path.dirname(bindir), "tables"), enlace)
    proc = subprocess.run(
        [os.path.join(TOOLCHAIN, "bin", "mpirun"), "-np", str(nproc),
         "--oversubscribe", os.path.join(bindir, binario)],
        cwd=dirrun, env=F2._entorno(), capture_output=True, text=True, timeout=3600)
    if esperar_exito and proc.returncode != 0:
        raise RuntimeError("la corrida %r falló (rc=%d):\n%s" %
                           (etiqueta, proc.returncode,
                            (proc.stdout + proc.stderr)[-2000:]))
    return proc.returncode, proc.stdout + proc.stderr, dirrun


def diag_fs(dirrun):
    ruta = os.path.join(dirrun, "freesurface_diagnostic.txt")
    if not os.path.exists(ruta):
        raise RuntimeError("no se escribió freesurface_diagnostic.txt en %s" % dirrun)
    d = np.atleast_2d(np.loadtxt(ruta))
    if d.shape[1] != 9:
        raise RuntimeError("freesurface_diagnostic.txt tiene %d columnas, no 9"
                           % d.shape[1])
    return d


def serie_eta(dirrun, dt, cstep):
    """(t, eta_k0(t)/A complejo), verificando que los tiempos son los esperados."""
    d = diag_fs(dirrun)
    t_esp = np.arange(len(d)) * cstep * dt
    if np.max(np.abs(d[:, 0] - t_esp)) > 1e-9 * max(1.0, t_esp[-1]):
        raise RuntimeError("los tiempos del diagnóstico no son k*cstep*dt")
    return t_esp, (d[:, 1] + 1j * d[:, 2]) / AMP, d


# --------------------------------------------------------------------------------
# Los ocho chequeos
# --------------------------------------------------------------------------------

def f1_parser(ctx):
    b = ctx["bindir"]
    rc, sal, _ = correr(b, "f1_acepta", step=10, cstep=5, esperar_exito=False)
    if rc != 0 or "Unknown boundary" in sal or "Unsupported" in sal:
        return False, "con vbczend=%r la corrida falla:\n%s" % (CADENA, sal[-600:])
    rc, sal, _ = correr(b, "f1_abajo", bcsta=CADENA, bcend="noslip", step=10,
                        cstep=5, esperar_exito=False)
    if rc == 0 or "only at z=Lz" not in sal:
        return False, ("la superficie deformable en z=0 no aborta con el mensaje "
                       "'only at z=Lz' (rc=%d):\n%s" % (rc, sal[-600:]))
    rc, sal, _ = correr(b, "f1_rechaza", bcend="superficie_blanda", step=10,
                        cstep=5, esperar_exito=False)
    if rc == 0 or "Unknown boundary" not in sal:
        return False, "una cadena desconocida no abortó con 'Unknown boundary'"
    return True, "acepta %r arriba, la rechaza abajo, y rechaza lo desconocido" % CADENA


def f2_laplace_neumann_dirichlet(ctx):
    _, sal, dirrun = correr(ctx["bindir_tests"], "f2_laplace", binario="TESTS",
                            bcend="noslip", step=10, cstep=5)
    ruta = os.path.join(dirrun, "laplace_neudir.txt")
    if not os.path.exists(ruta):
        return False, ("no se escribió laplace_neudir.txt (falta el test o la rama)."
                       "\n%s" % sal[-800:])
    d = np.loadtxt(ruta)
    peor_neu = peor_dir = peor_cerr = 0.0
    modos = 0
    for clave in sorted(set(map(tuple, d[:, :2].astype(int)))):
        m = (d[:, 0].astype(int) == clave[0]) & (d[:, 1].astype(int) == clave[1])
        k, z = d[m, 2][0], d[m, 3]
        a = d[m, 4] + 1j * d[m, 5]
        b = d[m, 6] + 1j * d[m, 7]
        g0 = d[m, 8][0] + 1j * d[m, 9][0]
        g1 = d[m, 10][0] + 1j * d[m, 11][0]
        largo = LZ if k == 0.0 else min(LZ, 1.0 / k)
        escala = max(abs(g1), abs(g0) * largo, 1e-30)
        modos += 1
        peor_neu = max(peor_neu, abs(b[0] - g0) / max(abs(g0), 1e-30))
        peor_dir = max(peor_dir, abs(a[-1] - g1) / escala)
        # forma cerrada en la base cosh/sinh: phi = A cosh(kz) + B sinh(kz)
        if k == 0.0:
            ref = g0 * z + (g1 - g0 * LZ)
        elif k * LZ < 30.0:
            B = g0 / k
            A = (g1 - B * np.sinh(k * LZ)) / np.cosh(k * LZ)
            ref = A * np.cosh(k * z) + B * np.sinh(k * z)
        else:
            ref = None
        if ref is not None:
            peor_cerr = max(peor_cerr, np.max(np.abs(a - ref)) / escala)
    if modos == 0:
        return False, "laplace_neudir.txt vacío"
    ok = peor_neu < 1e-12 and peor_dir < 1e-12 and peor_cerr < 1e-10
    return ok, ("%d modos; phi'(0)-g0 = %.2e, phi(Lz)-g1 = %.2e, contra cosh/sinh "
                "= %.2e" % (modos, peor_neu, peor_dir, peor_cerr))


def f3_referencia(ctx):
    detalle, ok = [], True
    for nombre, kx, ky, tens in CASOS:
        K = float(np.hypot(kx, ky))
        G = GRAV + tens * K * K
        b100 = balance_energia(100, K, G, 0.005, 2000)
        b400 = balance_energia(400, K, G, 0.005, 2000)
        _, err = referencia_eta(K, G, CSTEP_F4 * DT_F4, int(round(T_F4 / (CSTEP_F4 *
                                                                           DT_F4))))
        conv = b100 / b400
        caso_ok = b400 < 3e-6 and conv > 8.0 and err < 1e-6
        ok = ok and caso_ok
        detalle.append("%s: balance n=100 %.2e, n=400 %.2e (cociente %.1f, esperado "
                       "~16), Richardson %.2e" % (nombre, b100, b400, conv, err))
    return ok, "; ".join(detalle)


def _f4_caso(ctx, nombre, kx, ky, tens, nproc=1):
    pasos = int(round(T_F4 / DT_F4))
    _, _, dirrun = correr(ctx["bindir"], "f4_%s_np%d" % (nombre, nproc), step=pasos,
                          fskx=kx, fsky=ky, fstens=tens, nproc=nproc)
    t, r, d = serie_eta(dirrun, DT_F4, CSTEP_F4)
    return t, r, d


def f4_trayectorias(ctx):
    detalle, ok = [], True
    for nombre, kx, ky, tens in CASOS:
        K = float(np.hypot(kx, ky))
        G = GRAV + tens * K * K
        t, r, d = _f4_caso(ctx, nombre, kx, ky, tens)
        ref, err_ref = referencia_eta(K, G, CSTEP_F4 * DT_F4, len(t) - 1)
        ctx["f4_" + nombre] = (t, r, ref)
        err = float(np.max(np.abs(r.real - ref.real)))
        imag = float(np.max(np.abs(r.imag)))
        media = float(np.max(np.abs(d[:, 3]))) / AMP
        caso_ok = err < TOL_F4 and imag < TOL_F4 and media < 1e-12
        ok = ok and caso_ok
        detalle.append("%s: max|eta-ref|/A = %.2e (ref %.1e), Im %.1e, media %.1e"
                       % (nombre, err, err_ref, imag, media))
    # MPI: mismo caso con dos procesos
    nombre, kx, ky, tens = CASOS[0]
    _, r2, _ = _f4_caso(ctx, nombre, kx, ky, tens, nproc=2)
    r1 = ctx["f4_" + nombre][1]
    dif = float(np.max(np.abs(r2 - r1))) if len(r2) == len(r1) else np.inf
    ok = ok and dif < 1e-10
    detalle.append("np=2 contra np=1: %.1e" % dif)
    return ok, "; ".join(detalle)


def _residuo_tension(dirrun):
    d = diag_fs(dirrun)
    mitad = d[len(d) // 2:]
    return float(np.max(mitad[:, 5])), float(np.max(mitad[:, 6]))


def f5_tension(ctx):
    b = ctx["bindir"]
    med = {}
    for coup in (2, 1):
        for etiqueta, dt, pasos in (("dt", 1.0e-3, 200), ("dt2", 5.0e-4, 400)):
            _, _, dirrun = correr(b, "f5_c%d_%s" % (coup, etiqueta), dt=dt, step=pasos,
                                  cstep=20, fscoup=coup)
            sup, fondo = _residuo_tension(dirrun)
            med[(coup, etiqueta)] = sup / fondo if fondo > 0 else np.inf
    e1, e2 = med[(2, "dt")], med[(2, "dt2")]
    n1, n2 = med[(1, "dt")], med[(1, "dt2")]
    razon_e = e1 / e2 if e2 > 0 else 1.0
    razon_n = n1 / n2 if n2 > 0 else 0.0
    exacto_ok = e1 < 1e-9 and 0.2 < razon_e < 5.0
    control_ok = n1 > 100 * e1 and 2.5 < razon_n < 6.0
    msg = ("exacto: %.2e (razón dt %.2f); rezagado: %.2e (razón dt %.2f, esperado ~4 "
           "por ser O(dt) en amplitud)" % (e1, razon_e, n1, razon_n))
    if not exacto_ok:
        msg += " -> la tensión no se anula o depende de dt"
    if not control_ok:
        msg += " -> el control negativo no se comporta como O(dt): el chequeo no mide"
    return exacto_ok and control_ok, msg


def f6_orden_temporal(ctx):
    nombre, kx, ky, tens = CASOS[0]
    series = []
    for dt, cstep in ((1.0e-3, 40), (5.0e-4, 80), (2.5e-4, 160)):
        pasos = int(round(4.0 / dt))
        _, _, dirrun = correr(ctx["bindir"], "f6_%g" % dt, dt=dt, step=pasos,
                              cstep=cstep, fskx=kx, fsky=ky, fstens=tens)
        series.append(serie_eta(dirrun, dt, cstep)[1])
    n = min(len(s) for s in series)
    e1 = float(np.max(np.abs(series[0][:n] - series[1][:n])))
    e2 = float(np.max(np.abs(series[1][:n] - series[2][:n])))
    p = np.log2(e1 / e2) if e2 > 0 else np.inf
    return 1.6 < p < 2.6, ("diferencias sucesivas %.2e, %.2e: orden observado %.2f "
                           "(ORD=%d)" % (e1, e2, p, ORD))


def f7_corte_horizontal(ctx):
    ref = F2.tasas_referencia("neumann", 1)[0]
    _, _, dirrun = correr(ctx["bindir"], "f7_corte", dt=F2.DT_BASE, step=F2.PASOS,
                          cstep=F2.CSTEP, u0=1.0, vparam0=0.5, vparam1=2.0, fsamp=0.0)
    lam, discrep, disp = F2.tasa_de_corrida(dirrun)
    eta_rms = float(np.max(np.abs(diag_fs(dirrun)[:, 4])))
    err = abs(lam - ref) / ref
    ok = err < F2.TOL_TASA and discrep < 1e-6 and disp < 1e-9 and eta_rms < 1e-30
    return ok, ("lambda %.9e contra %.9e (%.2e); energía vs enstrofía %.1e; rms(eta) "
                "%.1e" % (lam, ref, err, discrep, eta_rms))


def f8_poder_discriminante(ctx):
    detalle, ok = [], True
    for nombre, kx, ky, tens in CASOS:
        clave = "f4_" + nombre
        if clave not in ctx:
            t, r, _ = _f4_caso(ctx, nombre, kx, ky, tens)
        else:
            t, r, _ = ctx[clave]
        K = float(np.hypot(kx, ky))
        variantes = [("sin 2nu dw/dz", dict(G=GRAV + tens * K * K, visc=False)),
                     ("freeslip plano reusado", dict(G=GRAV + tens * K * K,
                                                     tension="plana"))]
        if tens > 0:
            variantes += [("sin capilaridad", dict(G=GRAV)),
                          ("capilaridad con signo cambiado",
                           dict(G=GRAV - tens * K * K))]
        for etiqueta, kw in variantes:
            G = kw.pop("G")
            var = _trayectoria(200, K, G, CSTEP_F4 * DT_F4, len(t) - 1, **kw)[0][-1]
            dev = float(np.max(np.abs(r.real - var.real)))
            ok = ok and dev > 10 * TOL_F4
            detalle.append("%s / %s: %.2e" % (nombre, etiqueta, dev))
    return ok, ("desvío de SPECTER respecto de cada modelo equivocado (tiene que ser "
                "> %.0e): " % (10 * TOL_F4)) + "; ".join(detalle)


CHEQUEOS = [
    ("F1", "parser: 'freesurface' arriba sí, abajo no, desconocido no", f1_parser),
    ("F2", "rama Neumann-Dirichlet de laplace_z (test unitario)",
     f2_laplace_neumann_dirichlet),
    ("F3", "la referencia: balance de energía y Richardson", f3_referencia),
    ("F4", "trayectorias eta(t) contra la referencia, 1 y 2 procesos",
     f4_trayectorias),
    ("F5", "tensión tangencial completa nula; control negativo O(dt)", f5_tension),
    ("F6", "orden temporal observado", f6_orden_temporal),
    ("F7", "el modo k=0 decae como en la Fase 2", f7_corte_horizontal),
    ("F8", "poder discriminante contra cuatro modelos equivocados",
     f8_poder_discriminante),
]


def main():
    if not os.path.isdir(ARBOL):
        print("[ERROR] no existe el árbol de trabajo %s" % ARBOL)
        return 2
    if not os.path.exists(os.path.join(TOOLCHAIN, "bin", "mpif90")):
        print("[ERROR] no hay mpif90 en %s. Poné SPECTER_TOOLCHAIN." % TOOLCHAIN)
        return 2
    scratch = os.environ.get("SPECTER_SCRATCH")
    temporal = scratch is None
    if temporal:
        scratch = tempfile.mkdtemp(prefix="fase5_")
    os.makedirs(scratch, exist_ok=True)
    print("=" * 78)
    print("Puerta de aceptación — Fase 5: superficie libre deformable (lineal)")
    print("  árbol   : %s" % ARBOL)
    print("  scratch : %s" % scratch)
    print("  malla   : %dx%dx%d (Cz=%d), ORD=%d, nu=%g, g=%g, A=%g" %
          (NX, NY, NZ, CZ, ORD, NU, GRAV, AMP))
    print("=" * 78)
    ctx = {}
    try:
        print("\n[build] compilando el solver HD ...")
        ctx["bindir"] = construir(os.path.join(scratch, "main"), "main")
        print("[build] compilando los tests ...")
        ctx["bindir_tests"] = construir(os.path.join(scratch, "tests"), "tests")
    except Exception as exc:                            # noqa: BLE001
        print("\n[ERROR] %s\n\nRESULTADO: 0/%d — no se pudo compilar"
              % (exc, len(CHEQUEOS)))
        return 1
    resultados = []
    for clave, titulo, fn in CHEQUEOS:
        print("\n--- %s  %s" % (clave, titulo), flush=True)
        try:
            ok, msg = fn(ctx)
        except Exception as exc:                        # noqa: BLE001
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
