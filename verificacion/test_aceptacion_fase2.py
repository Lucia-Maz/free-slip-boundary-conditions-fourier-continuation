#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Puerta de aceptación de la Fase 2: la condición de superficie libre en SPECTER.

Se escribe ANTES de la implementación y por lo tanto arranca en rojo. El criterio de
diseño es el mismo de la Fase 1: la puerta no contiene ninguna fórmula analítica del
resultado que quiere verificar. Las tasas de decaimiento de referencia se obtienen
armando acá mismo el operador de diferencias finitas y sacándole los autovalores, con
extrapolación de Richardson; el factor 4 se MIDE corriendo el mismo binario con las dos
condiciones de contorno.

Los seis chequeos:

  V1  el parser reconoce "freeslip" y sigue abortando ante una cadena desconocida
  V2  la rama Dirichlet(z=0)-Neumann(z=Lz) de laplace_z, contra una forma cerrada
      independiente y contra sus propias propiedades definitorias
  V3  el residuo de tensión tangencial en la superficie libre es de nivel de
      reconstrucción y NO escala con dt (la predicción fuerte de la derivación)
  V4  el espectro de decaimiento coincide con el del problema de Stokes vertical con
      fondo no-deslizante y tope libre, para los tres modos más lentos
  V5  el cociente entre la fricción del canal y la de la capa con superficie libre,
      medido con el mismo binario, y la simetría de la cara libre abajo / arriba
  V6  poder discriminante: los criterios de V3 y V4 aplicados a la corrida
      no-deslizante tienen que FALLAR

Uso:
    /home/lucia/miniforge3/envs/piv-dt/bin/python verificacion/test_aceptacion_fase2.py

Variables de entorno opcionales:
    SPECTER_TOOLCHAIN   prefijo con bin/mpif90 y lib/libfftw3 (default: env conda
                        "specter" de miniforge)
    SPECTER_SCRATCH     directorio de compilación y corridas (default: un temporal)
    SPECTER_KEEP=1      no borrar el scratch al terminar
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

# --------------------------------------------------------------------------------
# Configuración
# --------------------------------------------------------------------------------

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ARBOL = os.path.join(RAIZ, "numerico", "SPECTER-trabajo")

TOOLCHAIN = os.environ.get(
    "SPECTER_TOOLCHAIN", "/home/lucia/miniforge3/envs/specter"
)

# Resolución de los tests. NZ es el total (físico + continuación), de modo que el
# dominio físico tiene NZ-CZ = 103 puntos y dz = Lz/(NZ-CZ-1).
NX, NY, NZ = 8, 8, 128
CZ, OZ = 25, 5
ORD = 2                      # orden del Runge-Kutta
NZFIS = NZ - CZ              # puntos físicos en z

LZ = 1.0                     # altura de la capa en las corridas de test
NU = 1.0e-2                  # viscosidad
U0 = 1.0                     # amplitud inicial

# Marcas de las corridas de decaimiento
DT_BASE = 2.0e-4
PASOS = 12000                # t final = 2.4
CSTEP = 20

CADENA_LIBRE = "freeslip"    # la cadena nueva, congelada en intent_fase2.txt

# Tolerancia de las tasas de decaimiento. NO es un número elegido a ojo: se midió
# corriendo el camino no-deslizante de upstream, que este trabajo no modifica, sobre
# esta misma malla y contra esta misma referencia de diferencias finitas. Los tres
# modos más lentos dieron errores relativos de 3.1e-9, 4.8e-8 y 2.4e-7. El umbral es
# unas 40 veces el peor de esos, de modo que un error de implementación tiene que ser
# grosero para pasar, y la malla no lo hace fallar sola.
TOL_TASA = 1.0e-5


# --------------------------------------------------------------------------------
# Referencias construidas acá, sin escribir la fórmula del resultado
# --------------------------------------------------------------------------------

def _operador_stokes(n, h, borde_arriba):
    """Matriz de d2/dz2 en (0,h] con u(0)=0 y, en z=h, o bien u(h)=0 ('dirichlet')
    o bien u'(h)=0 ('neumann'), en diferencias finitas de segundo orden."""
    dz = h / n
    if borde_arriba == "dirichlet":
        m = n - 1                      # incógnitas u_1..u_{n-1}
    elif borde_arriba == "neumann":
        m = n                          # incógnitas u_1..u_n
    else:
        raise ValueError(borde_arriba)
    a = np.zeros((m, m))
    for i in range(m):
        a[i, i] = -2.0
        if i > 0:
            a[i, i - 1] = 1.0
        if i < m - 1:
            a[i, i + 1] = 1.0
    if borde_arriba == "neumann":
        # punto fantasma u_{n+1} = u_{n-1} a segundo orden
        a[m - 1, m - 2] = 2.0
    return a / dz**2


def tasas_referencia(borde_arriba, cuantas, h=LZ, nu=NU, n=400):
    """Tasas de decaimiento lambda del problema de Stokes vertical, sacadas de los
    autovalores del operador de arriba, con Richardson entre n y 2n para matar el
    error O(dz^2) de la referencia."""
    def crudas(nn):
        val = np.linalg.eigvals(_operador_stokes(nn, h, borde_arriba))
        val = np.sort(-nu * np.real(val))
        return val[:cuantas]
    l1 = crudas(n)
    l2 = crudas(2 * n)
    return (4.0 * l2 - l1) / 3.0        # extrapolación de Richardson, orden 2


def laplace_referencia(k, lz, g0, g1, z):
    """Solución de phi'' - k^2 phi = 0 con phi(0)=g0 y phi'(Lz)=g1, escrita en la base
    cosh/sinh. Es una ruta algebraica distinta de la que usa el código (que normaliza
    con exp(-2 k Lz)), así que sirve de control independiente. Sólo es utilizable en
    doble precisión para k*Lz moderado; para k*Lz grande desborda y el chequeo se hace
    contra las propiedades definitorias."""
    if k == 0.0:
        return g0 + g1 * z
    ch, sh = np.cosh(k * lz), np.sinh(k * lz)
    b = (g1 / k - g0 * sh) / ch
    return g0 * np.cosh(k * z) + b * np.sinh(k * z)


# --------------------------------------------------------------------------------
# Arnés: compilar y correr SPECTER
# --------------------------------------------------------------------------------

def _entorno():
    env = dict(os.environ)
    env["PATH"] = os.path.join(TOOLCHAIN, "bin") + os.pathsep + env.get("PATH", "")
    env["LD_LIBRARY_PATH"] = (
        os.path.join(TOOLCHAIN, "lib") + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    )
    # OpenMPI se queja en contenedores y en máquinas de un solo nodo
    env["OMPI_MCA_plm_rsh_agent"] = ""
    env["OMPI_MCA_btl_vader_single_copy_mechanism"] = "none"
    env["OMPI_ALLOW_RUN_AS_ROOT"] = "1"
    env["OMPI_ALLOW_RUN_AS_ROOT_CONFIRM"] = "1"
    return env


def _parchear_makefile_in(ruta, solver="HD"):
    with open(ruta) as f:
        txt = f.read()

    reemplazos = [
        (r"^NX\s*=.*$",        "NX       = %d" % NX),
        (r"^NY\s*=.*$",        "NY       = %d" % NY),
        (r"^NZ\s*=.*$",        "NZ       = %d" % NZ),
        (r"^CZ\s*=.*$",        "CZ       = %d" % CZ),
        (r"^OZ\s*=.*$",        "OZ       = %d" % OZ),
        (r"^ORD\s*=.*$",       "ORD      = %d" % ORD),
        (r"^SOLVER\s*=.*$",    "SOLVER   = %s" % solver),
        (r"^PRECISION=.*$",    "PRECISION= DOUBLE"),
        (r"^FFTWDIR\s*=.*$",   "FFTWDIR = %s" % TOOLCHAIN),
        (r"^FC_GNU\s*=.*$",    "FC_GNU        = mpif90"),
        (r"^FFLAGS_GNU\s*=.*$",
         "FFLAGS_GNU    = -O2 -w -fallow-argument-mismatch"),
        (r"^FPSPEC_GNU\s*=.*$",
         "FPSPEC_GNU    = -O2 -w -fallow-argument-mismatch"),
        (r"^TARGET_GNU\s*=.*$", "TARGET_GNU    ="),
        (r"^MPIINC_OPEN\s*=.*$", "MPIINC_OPEN   ="),
        (r"^MPILIB_OPEN\s*=.*$", "MPILIB_OPEN   ="),
    ]
    for patron, nuevo in reemplazos:
        txt, n = re.subn(patron, nuevo, txt, count=1, flags=re.M)
        if n == 0:
            raise RuntimeError("no se pudo parchear %r en Makefile.in" % patron)

    # MPILD_OPEN ocupa varias líneas con continuaciones; se borra entero.
    txt = re.sub(r"^MPILD_OPEN\s*=(.*\\\n)*.*$", "MPILD_OPEN    =", txt, count=1,
                 flags=re.M)
    with open(ruta, "w") as f:
        f.write(txt)


IC_SHEAR = r"""
! Condición inicial de la puerta de aceptación de la Fase 2.
!   vparam1 = 0 : modo de corte puro vx = u0 sin(vparam0 pi z / Lz), sin estructura
!                 en x ni en y. Es autómodo del par (rígido abajo, libre arriba) para
!                 vparam0 = m + 1/2, y del canal para vparam0 entero. El término no
!                 lineal es exactamente nulo, así que la evolución es el problema de
!                 Stokes vertical y nada más.
!   vparam1 = 2 : el mismo modo para el par ESPEJO (libre abajo, rígido arriba),
!                 vx = u0 cos(vparam0 pi z / Lz), que cumple vx'(0)=0 y vx(Lz)=0.
!   vparam1 = 1 : modo solenoidal con dependencia en x, para que la condición de
!                 contorno se ejercite con kx distinto de cero.
      vx = 0; vy = 0; vz = 0

      IF (ista .eq. 1) THEN
         IF ( int(vparam1) .eq. 0 ) THEN
            DO k = 1,nz-Cz
               vx(k,1,1) = u0*nx*ny*SIN(vparam0*pi*z(k)/Lz)
            ENDDO
         ELSEIF ( int(vparam1) .eq. 2 ) THEN
            DO k = 1,nz-Cz
               vx(k,1,1) = u0*nx*ny*COS(vparam0*pi*z(k)/Lz)
            ENDDO
         ELSE
            ! v = (dPsi/dz, 0, -dPsi/dx) con Psi = a cos(kx x) sin(pi z/Lz):
            ! solenoidal por construcción, y con vz nulo en las dos caras.
            DO k = 1,nz-Cz
               vx(k,1,2) = u0*nx*ny*(pi/Lz)*COS(pi*z(k)/Lz)
               vz(k,1,2) = -im*kx(2)*u0*nx*ny*SIN(pi*z(k)/Lz)
            ENDDO
         ENDIF
      ENDIF

      CALL fftp1d_real_to_complex_z(planfc,vx,MPI_COMM_WORLD)
      CALL fftp1d_real_to_complex_z(planfc,vy,MPI_COMM_WORLD)
      CALL fftp1d_real_to_complex_z(planfc,vz,MPI_COMM_WORLD)
"""


def construir(destino, target="main"):
    """Copia el árbol de trabajo a `destino`, lo configura y compila. Devuelve la ruta
    del directorio bin/ donde quedó el binario."""
    if os.path.exists(destino):
        shutil.rmtree(destino)
    shutil.copytree(ARBOL, destino,
                    ignore=shutil.ignore_patterns(".git", "*.o", "*.mod"))

    _parchear_makefile_in(os.path.join(destino, "src", "Makefile.in"))
    with open(os.path.join(destino, "src", "initialv.f90"), "w") as f:
        f.write(IC_SHEAR)

    src = os.path.join(destino, "src")
    env = _entorno()
    # El Makefile de SPECTER usa GHOME = $(PWD), que make toma de la variable de
    # entorno PWD y no del cwd del proceso. Sin esto compila contra el árbol de quien
    # lo lanzó.
    env["PWD"] = src
    log = subprocess.run(["make", target], cwd=src, env=env,
                         capture_output=True, text=True)
    binario = os.path.join(destino, "bin", "HD" if target == "main" else "TESTS")
    if log.returncode != 0 or not os.path.exists(binario):
        cola = (log.stdout + log.stderr)[-3000:]
        raise RuntimeError("la compilación de '%s' falló:\n%s" % (target, cola))
    return os.path.join(destino, "bin")


PLANTILLA_INP = """&status
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
"""


def correr(bindir, etiqueta, *, bcsta="noslip", bcend="noslip", vparam0=1.0,
           vparam1=0.0, dt=DT_BASE, step=PASOS, cstep=CSTEP, lz=LZ, nu=NU, u0=U0,
           binario="HD", esperar_exito=True):
    """Corre el binario en un directorio propio y devuelve (returncode, stdout, dir)."""
    dirrun = os.path.join(os.path.dirname(bindir), "corridas", etiqueta)
    os.makedirs(os.path.join(dirrun, "out"), exist_ok=True)
    os.makedirs(os.path.join(dirrun, "in"), exist_ok=True)
    with open(os.path.join(dirrun, "parameter.inp"), "w") as f:
        f.write(PLANTILLA_INP.format(lz=lz, dt=dt, step=step, cstep=cstep, u0=u0,
                                     nu=nu, vparam0=vparam0, vparam1=vparam1,
                                     bcsta=bcsta, bcend=bcend))
    # tdir apunta a ../tables/ relativo al directorio de corrida
    enlace = os.path.join(os.path.dirname(dirrun), "tables")
    if not os.path.exists(enlace):
        os.symlink(os.path.join(os.path.dirname(bindir), "tables"), enlace)

    proc = subprocess.run(
        [os.path.join(TOOLCHAIN, "bin", "mpirun"), "-np", "1",
         "--oversubscribe", os.path.join(bindir, binario)],
        cwd=dirrun, env=_entorno(), capture_output=True, text=True, timeout=1800)
    if esperar_exito and proc.returncode != 0:
        raise RuntimeError("la corrida %r falló (rc=%d):\n%s" %
                           (etiqueta, proc.returncode,
                            (proc.stdout + proc.stderr)[-2000:]))
    return proc.returncode, proc.stdout + proc.stderr, dirrun


def leer_balance(dirrun):
    """balance.txt: tiempo, <v^2>, <omega^2>, inyección."""
    d = np.loadtxt(os.path.join(dirrun, "balance.txt"))
    return d[:, 0], d[:, 1], d[:, 2]


def tasa_de_ajuste(t, serie, desde=0.5):
    """lambda a partir de la pendiente de log(serie), descartando el arranque, y la
    dispersión del residuo, que dice si es una exponencial pura.

    Se ajusta la serie y no se usa nu*<w^2>/<v^2>: SPECTER promedia sobre el dominio
    EXTENDIDO por la continuación FC-Gram, no sobre el físico, así que <v^2> y <w^2>
    traen un factor de normalización que no se cancela entre ellos. En un cociente a
    lo largo del tiempo sí se cancela, y por eso la pendiente es exacta."""
    i0 = int(desde * len(t))
    coef = np.polyfit(t[i0:], np.log(serie[i0:]), 1)
    lam = -coef[0] / 2.0
    disp = float(np.std(np.log(serie[i0:]) - np.polyval(coef, t[i0:])))
    return lam, disp


def tasa_de_corrida(dirrun):
    """Devuelve (lambda, discrepancia entre energía y enstrofía, dispersión). Las dos
    series decaen como exp(-2*lambda*t) con la misma lambda, así que compararlas es un
    control de exponencial pura independiente del ajuste."""
    t, e, w = leer_balance(dirrun)
    lam_e, disp_e = tasa_de_ajuste(t, e)
    lam_w, _ = tasa_de_ajuste(t, w)
    return lam_e, abs(lam_e - lam_w) / lam_e, disp_e


def leer_diagnostico(dirrun, nombre):
    ruta = os.path.join(dirrun, nombre)
    if not os.path.exists(ruta):
        raise RuntimeError("no se escribió %s en %s" % (nombre, dirrun))
    return np.loadtxt(ruta)


# --------------------------------------------------------------------------------
# Los seis chequeos
# --------------------------------------------------------------------------------

def v1_parser(ctx):
    """El parser reconoce la cadena nueva y sigue rechazando lo desconocido."""
    bindir = ctx["bindir"]
    rc, salida, _ = correr(bindir, "v1_acepta", bcend=CADENA_LIBRE, step=10,
                           cstep=5, esperar_exito=False)
    if rc != 0:
        return False, ("con vbczend=%r el código falló (rc=%d). Últimas líneas:\n%s"
                       % (CADENA_LIBRE, rc, salida[-600:]))
    if "Unknown boundary" in salida or "Unsupported boundary" in salida:
        return False, "con vbczend=%r el código todavía protesta:\n%s" % (
            CADENA_LIBRE, salida[-600:])

    rc, salida, _ = correr(bindir, "v1_rechaza", bcend="superficie_libre_ideal",
                           step=10, cstep=5, esperar_exito=False)
    if rc == 0:
        return False, ("una cadena desconocida ('superficie_libre_ideal') fue "
                       "aceptada; el parser dejó de proteger")
    if "Unknown boundary" not in salida:
        return False, ("una cadena desconocida abortó, pero sin el mensaje "
                       "'Unknown boundary'; salida:\n%s" % salida[-600:])
    return True, "acepta %r y aborta ante una cadena desconocida" % CADENA_LIBRE


def v2_laplace_dirichlet_neumann(ctx):
    """La rama Dirichlet(z=0) - Neumann(z=Lz) de laplace_z.

    Salvedad declarada: el camino de superficie libre NO ejercita esta rama (ver
    intent_fase2.txt). Es un test unitario de la subrutina, no del camino físico.
    """
    bindir = ctx["bindir_tests"]
    _, salida, dirrun = correr(bindir, "v2_laplace", binario="TESTS", step=10,
                               cstep=5)
    ruta = os.path.join(dirrun, "laplace_dirneu.txt")
    if not os.path.exists(ruta):
        return False, ("no se escribió laplace_dirneu.txt: falta el test "
                       "src/tests/laplace_dirneu.f90 o la rama no está implementada."
                       "\nSalida:\n%s" % salida[-800:])
    d = np.loadtxt(ruta)
    # columnas: i, j, khom, z, Re(a), Im(a), Re(b), Im(b), Re(g0), Im(g0), Re(g1), Im(g1)
    peor_cerrada, peor_dir, peor_neu = 0.0, 0.0, 0.0
    modos = 0
    for clave in sorted(set(map(tuple, d[:, :2].astype(int)))):
        m = (d[:, 0].astype(int) == clave[0]) & (d[:, 1].astype(int) == clave[1])
        k = d[m, 2][0]
        z = d[m, 3]
        a = d[m, 4] + 1j * d[m, 5]
        b = d[m, 6] + 1j * d[m, 7]
        g0 = d[m, 8][0] + 1j * d[m, 9][0]
        g1 = d[m, 10][0] + 1j * d[m, 11][0]
        # La parte de la solución que fija el dato de Neumann tiene amplitud
        # |g1| por la longitud característica: Lz para k=0, y 1/k cuando el modo
        # decae más rápido que el dominio.
        long_car = LZ if k == 0.0 else min(LZ, 1.0 / k)
        escala = max(abs(g0), abs(g1) * long_car, 1e-30)
        modos += 1
        # propiedades definitorias, que valen para cualquier k
        peor_dir = max(peor_dir, abs(a[0] - g0) / escala)
        peor_neu = max(peor_neu, abs(b[-1] - g1) / max(abs(g1), 1e-30))
        # forma cerrada independiente, sólo donde no desborda
        if k * LZ < 30.0:
            ref = laplace_referencia(k, LZ, g0, g1, z)
            peor_cerrada = max(peor_cerrada,
                               np.max(np.abs(a - ref)) / escala)
    if modos == 0:
        return False, "laplace_dirneu.txt está vacío"
    ok = peor_dir < 1e-12 and peor_neu < 1e-12 and peor_cerrada < 1e-10
    msg = ("%d modos; phi(0)-g0 = %.2e, phi'(Lz)-g1 = %.2e, "
           "contra la forma cosh/sinh = %.2e" %
           (modos, peor_dir, peor_neu, peor_cerrada))
    return ok, msg


def _residuos_frontera(dirrun):
    """Del freeslip_diagnostic.txt: (residuo de tensión en z=0, en z=Lz), tomando el
    máximo sobre la segunda mitad de la corrida."""
    d = leer_diagnostico(dirrun, "freeslip_diagnostic.txt")
    # columnas: t, <|div v|^2>, <|d vt/dz|^2>|z=0, |z=Lz, <|vn|^2>|z=0, |z=Lz
    mitad = d[len(d) // 2:]
    return float(np.max(mitad[:, 2])), float(np.max(mitad[:, 3]))


def v3_tension_nula(ctx):
    """El residuo de tensión tangencial en la cara libre, y su dependencia con dt.

    Dos cosas distintas se chequean acá.

    (a) Magnitud. El residuo se normaliza contra el de la cara RÍGIDA de la MISMA
        corrida, donde la tensión tangencial es la física del problema y no un error.
        Así el criterio es adimensional y no depende de la amplitud ni de la malla.
        No se pide precisión de máquina: la condición se impone reconstruyendo el
        valor de borde con las tablas FC-Gram, que son de orden finito. La vara es la
        que el propio código alcanza para el no-deslizamiento, que sobre esta malla
        deja <|v_t|^2> = 1.0e-13 en la pared contra un campo de orden 10.

    (b) Dependencia con dt. La derivación dice que la condición se impone sobre la
        vorticidad tangencial y que la proyección la conserva EXACTAMENTE, así que el
        residuo NO debe escalar con dt. El término que acopla v*_z en el dato de borde
        es de orden dt: si su signo estuviera mal, el residuo pasaría a ser de orden
        dt y esta parte lo detecta. Es la parte falsable de la derivación.

    Se hace en las dos caras, porque la convención de signo de la derivada normal
    saliente es distinta arriba y abajo.
    """
    bindir = ctx["bindir"]
    detalle, ok = [], True
    for cara, kw in (("arriba", dict(bcsta="noslip", bcend=CADENA_LIBRE)),
                     ("abajo", dict(bcsta=CADENA_LIBRE, bcend="noslip"))):
        med = {}
        for nombre, dt, pasos in (("dt", DT_BASE, 400), ("dt2", DT_BASE / 2, 800)):
            _, _, dirrun = correr(bindir, "v3_%s_%s" % (cara, nombre), vparam1=1.0,
                                  dt=dt, step=pasos, cstep=50, **kw)
            bot, top = _residuos_frontera(dirrun)
            libre, rigida = (top, bot) if cara == "arriba" else (bot, top)
            med[nombre] = (libre, rigida)
        (l1, r1), (l2, _) = med["dt"], med["dt2"]
        rel = l1 / r1 if r1 > 0 else np.inf
        razon = l1 / l2 if l2 > 0 else 1.0
        chico = rel < 1e-9
        sin_dt = 0.2 < razon < 5.0
        ok = ok and chico and sin_dt
        detalle.append("cara %s: libre %.3e / rígida %.3e = %.2e, razón dt %.2f"
                       % (cara, l1, r1, rel, razon))
        if not chico:
            detalle.append("  -> el residuo no baja respecto de la cara rígida")
        if not sin_dt:
            detalle.append("  -> el residuo escala con dt: la condición NO se "
                           "conserva bajo la proyección")
    return ok, "; ".join(detalle)


def v4_espectro(ctx):
    """Las tres tasas más lentas del caso fondo no-deslizante / tope libre."""
    bindir = ctx["bindir"]
    ref = tasas_referencia("neumann", 3)
    medidas, errores = [], []
    for m in (0, 1, 2):
        _, _, dirrun = correr(bindir, "v4_m%d" % m, bcend=CADENA_LIBRE,
                              vparam0=m + 0.5)
        lam, discrep, disp = tasa_de_corrida(dirrun)
        if discrep > 1e-6 or disp > 1e-9:
            return False, ("modo m=%d: no decae como una exponencial pura "
                           "(energía vs enstrofía difieren en %.2e, dispersión del "
                           "ajuste %.2e). La condición inicial dejó de ser autómodo, "
                           "que es lo que pasa si la condición de contorno impuesta "
                           "no es la de tensión nula." % (m, discrep, disp))
        medidas.append(lam)
        errores.append(abs(lam - ref[m]) / ref[m])
    ctx["lambda_libre_m0"] = medidas[0]
    ok = max(errores) < TOL_TASA
    msg = "; ".join("m=%d: %.9e vs %.9e (%.2e)" % (m, medidas[m], ref[m], errores[m])
                    for m in range(3))
    return ok, msg


def v5_factor_cuatro(ctx):
    """El cociente entre la fricción del canal y la de la capa libre, medido; y la
    simetría entre poner la cara libre arriba o abajo."""
    bindir = ctx["bindir"]

    _, _, d_canal = correr(bindir, "v5_canal", bcsta="noslip", bcend="noslip",
                           vparam0=1.0)
    lam_canal, _, _ = tasa_de_corrida(d_canal)

    lam_libre = ctx.get("lambda_libre_m0")
    if lam_libre is None:
        _, _, d = correr(bindir, "v5_libre", bcend=CADENA_LIBRE, vparam0=0.5)
        lam_libre, _, _ = tasa_de_corrida(d)

    razon = lam_canal / lam_libre

    # Misma física con la cara libre abajo: verifica la convención de signo de la
    # derivada normal saliente en z=0, que es distinta de la de z=Lz.
    # El par espejo tiene su propio autómodo, cos(pi z / 2Lz), que cumple v'(0)=0 y
    # v(Lz)=0. Usar el seno de arriba daría una superposición y la tasa medida sería
    # la del transitorio, no la del modo.
    _, _, d_esp = correr(bindir, "v5_espejo", bcsta=CADENA_LIBRE, bcend="noslip",
                         vparam0=0.5, vparam1=2.0)
    lam_espejo, discrep_esp, disp_esp = tasa_de_corrida(d_esp)
    if discrep_esp > 1e-6 or disp_esp > 1e-9:
        return False, ("el caso espejo no decae como exponencial pura (energía vs "
                       "enstrofía %.2e, dispersión %.2e): la cara libre de abajo no "
                       "está imponiendo tensión nula" % (discrep_esp, disp_esp))
    dif_espejo = abs(lam_espejo - lam_libre) / lam_libre

    ok = abs(razon - 4.0) / 4.0 < TOL_TASA and dif_espejo < TOL_TASA
    msg = ("canal %.9e / libre %.9e = %.9f (esperado 4); "
           "cara libre abajo %.9e, difiere %.2e" %
           (lam_canal, lam_libre, razon, lam_espejo, dif_espejo))
    return ok, msg


def v6_poder_discriminante(ctx):
    """El caso que falla a propósito: los criterios de V3 y V4 aplicados a la corrida
    no-deslizante tienen que dar rojo. Si dieran verde, la puerta no estaría midiendo
    nada."""
    bindir = ctx["bindir"]

    # (a) el residuo de tensión no es chico cuando la cara es no-deslizante
    _, _, dirrun = correr(bindir, "v6_noslip_diag", bcend="noslip", vparam1=1.0,
                          dt=DT_BASE, step=400, cstep=50)
    bot_ns, top_ns = _residuos_frontera(dirrun)
    rel_ns = top_ns / bot_ns
    a_falla = rel_ns > 1e-9

    # (b) la tasa del canal no coincide con la del tope libre
    lam_ref = tasas_referencia("neumann", 1)[0]
    _, _, d_canal = correr(bindir, "v6_canal", bcend="noslip", vparam0=1.0)
    lam_canal, _, _ = tasa_de_corrida(d_canal)
    b_falla = abs(lam_canal - lam_ref) / lam_ref > TOL_TASA

    ok = a_falla and b_falla
    msg = ("no-deslizante: residuo relativo de tensión arriba %.2e (el criterio de "
           "V3 pide < 1e-9); lambda = %.9e contra la referencia libre %.9e" %
           (rel_ns, lam_canal, lam_ref))
    if not a_falla:
        return False, ("el diagnóstico de tensión no distingue no-deslizamiento de "
                       "superficie libre. " + msg)
    if not b_falla:
        return False, ("la tasa del canal coincide con la de la capa libre: el test "
                       "no mide la condición de contorno. " + msg)
    return True, msg


CHEQUEOS = [
    ("V1", "el parser acepta la cadena nueva y sigue rechazando lo desconocido",
     v1_parser),
    ("V2", "la rama Dirichlet-Neumann de laplace_z (test unitario)",
     v2_laplace_dirichlet_neumann),
    ("V3", "residuo de tensión tangencial en la cara libre, e independencia de dt",
     v3_tension_nula),
    ("V4", "el espectro lambda_m del problema fondo rígido / tope libre",
     v4_espectro),
    ("V5", "el factor 4 medido, y la simetría de la cara libre",
     v5_factor_cuatro),
    ("V6", "poder discriminante: la corrida no-deslizante tiene que fallar",
     v6_poder_discriminante),
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
        scratch = tempfile.mkdtemp(prefix="fase2_")
    os.makedirs(scratch, exist_ok=True)

    print("=" * 78)
    print("Puerta de aceptación — Fase 2: superficie libre en SPECTER")
    print("  árbol      : %s" % ARBOL)
    print("  toolchain  : %s" % TOOLCHAIN)
    print("  scratch    : %s" % scratch)
    print("  malla      : %dx%dx%d (Cz=%d, %d puntos físicos en z), ORD=%d" %
          (NX, NY, NZ, CZ, NZFIS, ORD))
    print("=" * 78)

    ctx = {}
    try:
        print("\n[build] compilando el solver HD ...")
        ctx["bindir"] = construir(os.path.join(scratch, "main"), target="main")
        print("[build] compilando los tests ...")
        ctx["bindir_tests"] = construir(os.path.join(scratch, "tests"),
                                        target="tests")
    except Exception as exc:                       # noqa: BLE001
        print("\n[ERROR] %s" % exc)
        print("\nRESULTADO: 0/6 — no se pudo compilar")
        return 1

    resultados = []
    for clave, titulo, fn in CHEQUEOS:
        print("\n--- %s  %s" % (clave, titulo))
        try:
            ok, msg = fn(ctx)
        except Exception as exc:                   # noqa: BLE001
            ok, msg = False, "excepción: %s" % exc
        print("    %s  %s" % ("PASA" if ok else "FALLA", msg))
        resultados.append((clave, ok))

    n = sum(1 for _, ok in resultados if ok)
    print("\n" + "=" * 78)
    print("RESULTADO: %d/%d  (%s)" %
          (n, len(resultados),
           " ".join("%s:%s" % (c, "ok" if ok else "NO") for c, ok in resultados)))
    print("=" * 78)

    if temporal and os.environ.get("SPECTER_KEEP") != "1":
        shutil.rmtree(scratch, ignore_errors=True)
    else:
        print("scratch conservado en %s" % scratch)
    return 0 if n == len(resultados) else 1


if __name__ == "__main__":
    sys.exit(main())
