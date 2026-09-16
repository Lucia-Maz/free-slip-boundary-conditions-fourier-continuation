#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Etapa 1 del barrido en delta: ¿se aparta la tasa de decaimiento de la predicción
lineal cuando el término no lineal deja de ser chico?

QUÉ CONTESTA
============
[P-01] dejó dos flecos: el factor geométrico O(1) de la estimación de la distorsión del
perfil vertical, que puse en 1 sin calcularlo, y si el sesgo sobre alfa es de orden delta
o delta^2. Los dos se cierran midiendo cuánto se aparta la tasa medida de la lineal a
medida que sube la amplitud.

CÓMO, Y POR QUÉ ASÍ
===================
La referencia lineal NO se escribe a mano: se MIDE, corriendo el mismo binario, la misma
malla y la misma condición inicial a amplitud minúscula, donde el término no lineal es
despreciable por construcción. Así la referencia incluye automáticamente el k discreto de
la malla, el sesgo del integrador y el filtro de dealiasing, y el cociente
lambda_eff/lambda_lin es un número puro: cualquier apartamiento es no linealidad y nada
más. Escribir la fórmula analítica habría mezclado el efecto buscado con los errores de
discretización.

Dos resoluciones horizontales, y no es opcional: [H-09] mostró que con superficie libre el
término no lineal se desestabiliza sobre mallas donde el canal aguanta. Un apartamiento
que NO converge con la resolución es numérico, no físico, y ese es el confundido que hay
que descartar.

La condición inicial es propia ([H-13]): celular y tridimensional, el modo fundamental
del par rígido-libre en z por una suma de dos modos horizontales que no es autofunción
del laplaciano, para que el término no lineal no se anule por simetría.

Dos casos, por [D-41]/[H-16]: `forzado` (eps = 1,78, el fundamental de la red de imanes,
delta hasta 30) y `ventana` (eps = 0,71, la ventana del ajuste experimental, delta hasta
15). La ventana de ajuste es tardía, a más de tres memorias modales, porque la pregunta
([P-02]) es sobre el régimen cuasi-estacionario y no sobre el transitorio.

Uso:
    /home/lucia/miniforge3/envs/piv-dt/bin/python verificacion/barrido_delta_etapa1.py --caso forzado
    /home/lucia/miniforge3/envs/piv-dt/bin/python verificacion/barrido_delta_etapa1.py --caso ventana
    ... --prueba   corrida de humo (una malla, dos amplitudes, pasos/20), JSON aparte
"""

import importlib.util
import json
import math
import os
import sys

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(AQUI)

_spec = importlib.util.spec_from_file_location(
    "fase2", os.path.join(AQUI, "test_aceptacion_fase2.py"))
fase2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fase2)

def ruta_salida(caso):
    return os.path.join(RAIZ, "salidas", "tablas",
                        "barrido_delta_etapa1_%s.json" % caso)


# --- geometría y física de las corridas ------------------------------------------
# Los números de onda de SPECTER son ENTEROS (specter.fpp:777-780, ky(j)=j-1), o sea
# que el Lx=1 del archivo de parámetros son 2*pi de largo físico y el modo fundamental
# es k=1. Verificado comparando <v^2> del código contra la fórmula analítica de la
# condición inicial: 4.750e-3 contra 4.750e-3 con u0=0.1. Por eso Lz NO lleva 2*pi:
# es la longitud literal.
#
# Dos casos, porque el experimento tiene dos regímenes con eps distinto ([D-41], [H-16]):
# el estado forzado y el arranque del decaimiento están en el fundamental de la red de
# imanes, k = 296 1/m, eps = k h = 1,78; y la ventana donde se ajusta alfa (t > 10 s)
# tiene la energía ya migrada a k ~ 60-130 1/m, eps ~ 0,35-0,8. El modo dominante de
# la condición inicial es el diagonal, |k| = sqrt(2) en unidades de código, así que
# Lz = eps / sqrt(2). (La versión anterior tenía Lz = 0,5 y decía "eps = 0,5" contando
# k = 1, el axial; con el diagonal es 0,71, que cae en la ventana del ajuste y por eso
# el caso se conserva tal cual, reinterpretado.)
#
# nu escala con Lz^2 para que el tiempo de código sea el mismo múltiplo de h^2/nu en
# los dos casos: adimensionalizando con h y h^2/nu los únicos parámetros son eps y
# delta, y nu sólo fija la unidad de tiempo. Con eso dz^2/nu, la memoria modal
# Lz^2/(2 pi^2 nu) = 1,27 y el margen de dt son idénticos entre casos.
LZ_REF, NU_REF = 0.5, 1.0e-2
CASOS = {
    # eps del fundamental de la red; delta en la ventana hasta ~30 (arranque: 27)
    "forzado": {"eps": 1.7772, "delta_objetivo": (0.3, 1.0, 3.0, 10.0, 20.0, 30.0)},
    # eps de la ventana del ajuste; delta hasta ~15 (ahí va de ~10 a ~1)
    "ventana": {"eps": 0.7071, "delta_objetivo": (0.3, 1.0, 3.0, 6.0, 10.0, 15.0)},
}
for _c in CASOS.values():
    _c["lz"] = _c["eps"] / math.sqrt(2.0)
    _c["nu"] = NU_REF * (_c["lz"] / LZ_REF) ** 2
NZ, CZ = 64, 25        # 39 puntos físicos en z: de sobra para sin(pi z / 2h)
DT = 5.0e-4            # ~1/6 del límite de estabilidad viscosa explícita (dz^2/nu no
                       # cambia entre casos, así que el margen tampoco)
RESOLUCIONES = (16, 32)
CAMPOS_EN_VENTANA = 6   # etapa 1b (--superficie): campos escritos dentro de la ventana
U0_REFERENCIA = 1.0e-4  # referencia lineal: delta ~ 1e-3
VPARAM1 = 1.0

# --- condición inicial propia -----------------------------------------------------
# La de la puerta NO sirve acá y conviene decir por qué: es psi ~ cos(kx x) sin(pi z/Lz)
# en el plano (x,z), que es autofunción del laplaciano 2D, de modo que J(psi, lap psi)
# se anula idénticamente y el término no lineal es CERO. Además, en un flujo 2D en
# (x,z) la continuidad integrada en z da d/dx <u> = -[w] = 0, o sea que el promedio
# vertical no puede tener estructura horizontal, que es justo la variable de la
# clausura. La prueba tiene que ser tridimensional.
#
# Se usa entonces  u_h = F(z) grad_perp(psi(x,y)),  w = 0,  con
#     F(z)  = sin(pi z / 2Lz)          el modo fundamental del par rígido-libre
#     psi   = cos(k x)cos(k y) + b cos(2k x)cos(k y)
# El campo es horizontalmente solenoidal por construcción (u = d psi/dy, v = -d psi/dx)
# y tiene w = 0 exacto, así que cumple impermeabilidad en las dos caras. La suma de dos
# modos con |k|^2 distinto (2k^2 y 5k^2) rompe la degeneración: psi ya no es autofunción
# del laplaciano horizontal y J(psi, lap psi) no se anula.
B_SEGUNDO_MODO = 0.6

IC_CELULAR = r"""
! Condicion inicial celular 3D para el barrido en delta.
!   u_h = F(z) grad_perp(psi),  w = 0
!   F(z) = sin(pi z / 2 Lz);  psi = cos(k x)cos(k y) + b cos(2k x)cos(k y)
      vx = 0; vy = 0; vz = 0

      IF (ista .eq. 1) THEN
         DO k = 1,nz-Cz
            rmt = SIN(0.5_GP*pi*z(k)/Lz)
            ! modo 1: kx(2), +-ky(2)
            rmq = 0.25_GP*u0*nx*ny
            vx(k,2,2)  = vx(k,2,2)  + im*ky(2)*rmq*rmt
            vy(k,2,2)  = vy(k,2,2)  - im*kx(2)*rmq*rmt
            vx(k,ny,2) = vx(k,ny,2) - im*ky(2)*rmq*rmt
            vy(k,ny,2) = vy(k,ny,2) - im*kx(2)*rmq*rmt
            ! modo 2: kx(3), +-ky(2), con amplitud relativa b
            rmq = 0.25_GP*BSEG*u0*nx*ny
            vx(k,2,3)  = vx(k,2,3)  + im*ky(2)*rmq*rmt
            vy(k,2,3)  = vy(k,2,3)  - im*kx(3)*rmq*rmt
            vx(k,ny,3) = vx(k,ny,3) - im*ky(2)*rmq*rmt
            vy(k,ny,3) = vy(k,ny,3) - im*kx(3)*rmq*rmt
         ENDDO
      ENDIF

      CALL fftp1d_real_to_complex_z(planfc,vx,MPI_COMM_WORLD)
      CALL fftp1d_real_to_complex_z(planfc,vy,MPI_COMM_WORLD)
      CALL fftp1d_real_to_complex_z(planfc,vz,MPI_COMM_WORLD)
"""

# --- amplitudes: de delta objetivo a u0 --------------------------------------------
# Para esta condición inicial <v^2> = (1/2)|grad_perp psi|^2_rms = 0,475 u0^2 (con
# b = 0,6; verificado contra el código, ver arriba), la velocidad de la clausura es el
# promedio vertical U = (2/pi) sqrt(2 <v^2>) = 0,62 u0, y el k del campo pesado por
# enstrofía es sqrt((2·0,5 + 5·0,45)/0,95) = 1,85. Entonces delta(t=0) = U k Lz^2/nu =
# 28,7 u0 en los dos casos (Lz^2/nu es el mismo). En la ventana la energía ya cayó
# exp(-2 lambda_lin t), y eso sí depende de eps; el factor se estima acá con la
# lambda_lin analítica y el delta REAL se calibra después desde el <v^2> medido en la
# ventana, así que estos u0 son un punto de partida y no un resultado.
def _u0_para_delta(caso, delta):
    c = CASOS[caso]
    v2_0 = 0.5 * (0.5 + 1.25 * B_SEGUNDO_MODO ** 2)
    U_por_u0 = (2.0 / math.pi) * math.sqrt(2.0 * v2_0)
    k_ef = math.sqrt((2.0 * 0.5 + 5.0 * 1.25 * B_SEGUNDO_MODO ** 2)
                     / (0.5 + 1.25 * B_SEGUNDO_MODO ** 2))
    delta_por_u0 = U_por_u0 * k_ef * c["lz"] ** 2 / c["nu"]
    lam = c["nu"] * (k_ef ** 2 + (math.pi / (2.0 * c["lz"])) ** 2)
    t_medio = 0.5 * (VENT[0] + VENT[1]) * T_FINAL
    return delta / (delta_por_u0 * math.exp(-lam * t_medio))


# --- la ventana de ajuste, y por qué es TARDÍA ---------------------------------------
# La versión anterior ajustaba temprano, en [0,25 T, 0,70 T] con T = 3, con el
# argumento de que sobre la cola el flujo ya sería lineal. Pero la pregunta que hay que
# contestar ([P-02], candidato iv) es sobre la ventana del ajuste experimental, t > 10 s
# después del corte: más de cinco memorias modales (1,82 s), o sea con el perfil
# vertical en balance cuasi-estacionario con el término no lineal. Una ventana a menos
# de una memoria modal del arranque mide el transitorio, que es otra cosa. Se ajusta
# entonces en [3,0; 4,25] memorias modales, y se corre hasta 4,7. La energía en la
# ventana está al 13 % (forzado) o al 34 % (ventana) de la inicial: lejos de minúscula,
# y delta se lee ahí, no en t = 0. Dentro de la ventana la energía cae ~25 %, así que
# delta se informa en el medio y con su rango.
# La condición inicial ES el modo fundamental del par rígido-libre, así que en el caso
# lineal decae como una exponencial pura desde t=0, y la referencia no tiene transitorio.
T_MEM = LZ_REF ** 2 / (2.0 * math.pi ** 2 * NU_REF)     # = 1,27, igual en los dos casos
T_FINAL = 4.7 * T_MEM
PASOS = int(round(T_FINAL / DT))
CSTEP = 30
VENT = (3.0 / 4.7, 4.25 / 4.7)


def diagnostico_campo(dirrun, indice, nx, ny, nzf, lz, nu):
    """k horizontal, k tridimensional y alfa_eff, MEDIDOS sobre el campo escrito.

    Por qué sobre el campo y no sobre balance.txt: el <omega^2> que SPECTER escribe ahí
    no reproduce el rotor del campo que el mismo código guarda. Sobre la condición
    inicial del caso forzado (Lz = 1,257) el campo da <omega^2>/<v^2> = 4,983 contra
    4,982 analítico, y balance.txt dice 2,12; con Lz = 0,5 coincidían (13,4 vs 13,3),
    que es por lo que no se había notado. La puerta ya advertía que ese cociente trae
    una normalización que no se cancela. No se investigó la causa dentro de SPECTER:
    se deja de usar ([H-17]).

    El campo se escribe en (nx, ny, nz-Cz) doble precisión, orden Fortran, sólo el
    dominio físico. Derivadas espectrales en x,y (k enteros) y diferencias finitas de
    segundo orden en z. Devuelve:
      k_h        sqrt(<omega_z^2>/<u_h^2>): el k horizontal pesado por energía, que es
                 el que entra en delta
      k_3d       sqrt(<|omega|^2>/<v^2>), para contrastar con balance.txt
      alfa_eff   nu <d_z u_h|_0 . U_h> / (h <|U_h|^2>) con U_h el promedio vertical:
                 la fricción de fondo SIN clausura, dividida por la de la clausura
                 nu pi^2/(4 h^2). Vale 1 exactamente para el perfil sin(pi z/2h).
    """
    ext = "%04d" % indice
    ruta = os.path.join(dirrun, "out")
    def leer(nombre):
        return np.fromfile(os.path.join(ruta, nombre), dtype="<f8").reshape(
            (nx, ny, nzf), order="F")
    u, v, w = leer("vx.%s.out" % ext), leer("vy.%s.out" % ext), leer("vz.%s.out" % ext)
    z = np.loadtxt(os.path.join(dirrun, "z.txt"))
    kx = np.fft.fftfreq(nx, d=1.0 / nx)
    ky = np.fft.fftfreq(ny, d=1.0 / ny)
    def dx(f):
        return np.real(np.fft.ifft(1j * kx[:, None, None] * np.fft.fft(f, axis=0), axis=0))
    def dy(f):
        return np.real(np.fft.ifft(1j * ky[None, :, None] * np.fft.fft(f, axis=1), axis=1))
    def dz(f):
        return np.gradient(f, z, axis=2, edge_order=2)
    oz = dx(v) - dy(u)
    e_h = np.mean(u ** 2 + v ** 2)
    e = e_h + np.mean(w ** 2)
    k_h = math.sqrt(np.mean(oz ** 2) / e_h)
    ox = dy(w) - dz(v)
    oy = dz(u) - dx(w)
    k_3d = math.sqrt(np.mean(ox ** 2 + oy ** 2 + oz ** 2) / e)
    del ox, oy, oz
    # alfa_eff: tensión en el fondo contra promedio vertical, proyectada sobre U_h
    dzu0 = dz(u)[:, :, 0]
    dzv0 = dz(v)[:, :, 0]
    Ub = np.trapezoid(u, z, axis=2) / lz
    Vb = np.trapezoid(v, z, axis=2) / lz
    alfa_eff = nu * np.mean(dzu0 * Ub + dzv0 * Vb) / (lz * np.mean(Ub ** 2 + Vb ** 2))
    alfa_clausura = nu * math.pi ** 2 / (4.0 * lz ** 2)
    return {"k_h": float(k_h), "k_3d": float(k_3d),
            "alfa_eff_sobre_alfa": float(alfa_eff / alfa_clausura),
            "v2_campo": float(e), "w2_sobre_v2": float(np.mean(w ** 2) / e),
            # lo que ve el PIV: rms de la velocidad horizontal en el plano de arriba,
            # contra el rms del promedio vertical y el del campo 3D
            "rms_superficie": float(math.sqrt(np.mean(u[:, :, -1] ** 2 + v[:, :, -1] ** 2))),
            "rms_promedio": float(math.sqrt(np.mean(Ub ** 2 + Vb ** 2))),
            "rms_3d": float(math.sqrt(e))}


def pendientes_superficie(dirrun, indices_t, nx, ny, nzf, lz, nu, t0, t1):
    """Pendientes de log(rms) de la superficie, del promedio vertical y del campo 3D,
    sobre los campos escritos dentro de la ventana [t0, t1] ([V-16], candidato (v)).

    El PIV mide u en la superficie, porque las partículas flotan, y el alfa medido es la
    pendiente de log u(h). La clausura predice la del promedio vertical. Si el perfil se
    distorsiona con delta, y delta cae en el tiempo, u(h)/<u> no es pi/2 ni constante, y
    las dos pendientes difieren. Acá se miden las dos sobre los mismos campos.

    `indices_t` es la lista de (índice de archivo, t) de los campos escritos. Devuelve
    las tres tasas (definidas como -d log(rms)/dt, comparables con lambda), la
    dispersión del ajuste de superficie, y u(h)/<u> al principio y al final de la
    ventana. La de 3D tiene que coincidir con la lambda de balance.txt: es el chequeo.
    """
    tt, sup, prom, tot = [], [], [], []
    for indice, t in indices_t:
        if t < t0 - 1e-9 or t > t1 + 1e-9:
            continue
        dg = diagnostico_campo(dirrun, indice, nx, ny, nzf, lz, nu)
        tt.append(t)
        sup.append(dg["rms_superficie"])
        prom.append(dg["rms_promedio"])
        tot.append(dg["rms_3d"])
    if len(tt) < 3:
        return None
    tt = np.array(tt)
    def tasa(serie):
        coef = np.polyfit(tt, np.log(serie), 1)
        return -float(coef[0]), float(np.std(np.log(serie) - np.polyval(coef, tt)))
    l_sup, d_sup = tasa(sup)
    l_prom, _ = tasa(prom)
    l_tot, _ = tasa(tot)
    return {"n_campos": len(tt), "t_campos": tt.tolist(),
            "lambda_superficie": l_sup, "lambda_promedio": l_prom,
            "lambda_3d": l_tot, "dispersion_superficie": d_sup,
            "sup_sobre_prom_inicio": sup[0] / prom[0],
            "sup_sobre_prom_fin": sup[-1] / prom[-1]}


def tasa_en_ventana(t, e, t0, t1):
    """lambda del ajuste de log(<v^2>) ~ -2*lambda*t en [t0, t1], y su dispersión."""
    m = (t >= t0) & (t <= t1) & (e > 0)
    if m.sum() < 5:
        return float("nan"), float("nan")
    coef = np.polyfit(t[m], np.log(e[m]), 1)
    disp = float(np.std(np.log(e[m]) - np.polyval(coef, t[m])))
    return -float(coef[0]) / 2.0, disp


def _compilar_liviano():
    """Baja la optimización de gfortran de -O2 a -O1 sólo para estas corridas.

    No es una decisión física: compilar SPECTER con -O2 es lo que hace picar la memoria
    en esta máquina, y ya mató dos corridas. Con -O1 el binario es algo más lento pero
    idéntico en resultados, y estas corridas son chicas. La puerta de la Fase 2 sigue
    compilando con -O2, porque no se la toca.
    """
    orig = fase2._parchear_makefile_in

    def envoltura(ruta, solver="HD"):
        orig(ruta, solver)
        with open(ruta) as f:
            txt = f.read()
        with open(ruta, "w") as f:
            f.write(txt.replace("-O2", "-O1"))

    fase2._parchear_makefile_in = envoltura


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--caso", choices=sorted(CASOS), default="forzado",
                    help="régimen del experimento a reproducir (eps); ver CASOS")
    ap.add_argument("--resolucion", type=int, default=None,
                    help="correr una sola malla, para no tener dos compilaciones "
                         "vivas en la misma sesión de memoria")
    ap.add_argument("--prueba", action="store_true",
                    help="humo: una malla, dos amplitudes, 1/20 de los pasos; escribe "
                         "a un JSON aparte y no toca el de producción")
    ap.add_argument("--superficie", action="store_true",
                    help="etapa 1b: escribe %d campos en la ventana y mide la pendiente "
                         "de u(h) contra la de <u>; JSON con sufijo _superficie"
                         % CAMPOS_EN_VENTANA)
    args = ap.parse_args()
    caso = CASOS[args.caso]
    LZ, NU = caso["lz"], caso["nu"]
    U0S = tuple(round(_u0_para_delta(args.caso, d), 4) for d in caso["delta_objetivo"])
    resoluciones = (args.resolucion,) if args.resolucion else RESOLUCIONES
    pasos = PASOS
    sufijo = "_superficie" if args.superficie else ""
    salida = ruta_salida(args.caso + sufijo)
    if args.prueba:
        resoluciones = resoluciones[:1]
        U0S = U0S[-2:]
        pasos = max(PASOS // 20, 2 * CSTEP)
        salida = ruta_salida(args.caso + sufijo + "_prueba")
    _compilar_liviano()

    scratch = os.environ.get("SPECTER_SCRATCH")
    temporal = scratch is None
    if temporal:
        import tempfile
        scratch = tempfile.mkdtemp(prefix="delta1_")
    print("scratch: %s" % scratch, flush=True)
    print("caso %s: eps=%.3f  Lz=%.4f  nu=%.3e  dt=%.1e  pasos=%d  t_final=%.3f  "
          "ventana=[%.3f, %.3f]  (t_mem=%.3f)"
          % (args.caso, caso["eps"], LZ, NU, DT, pasos, T_FINAL, VENT[0]*T_FINAL,
             VENT[1]*T_FINAL, T_MEM), flush=True)
    print("u0 del barrido (para delta objetivo %s): %s"
          % (caso["delta_objetivo"], U0S), flush=True)
    # Campos escritos: el archivo 0001 es t=0 y el k-ésimo sale a los (k-1)*tstep pasos
    # (specter.fpp:859, 1005). Se toca la plantilla en memoria, no la puerta.
    # Etapa 1: un solo campo, en el MEDIO de la ventana. Etapa 1b (--superficie):
    # CAMPOS_EN_VENTANA campos repartidos en la ventana, para ajustar pendientes.
    paso_medio = int(round(0.5 * (VENT[0] + VENT[1]) * pasos))
    if args.superficie:
        tstep = int(round((VENT[1] - VENT[0]) * pasos / (CAMPOS_EN_VENTANA - 1)))
    else:
        tstep = paso_medio
    assert "tstep = 1000000" in fase2.PLANTILLA_INP
    fase2.PLANTILLA_INP = fase2.PLANTILLA_INP.replace("tstep = 1000000",
                                                      "tstep = %d" % tstep)
    # (índice de archivo, t) de cada campo que se escribe, y el más cercano al medio
    indices_t = [(k + 1, k * tstep * DT) for k in range(pasos // tstep + 1)]
    indice_medio = min(indices_t, key=lambda it: abs(it[1] - paso_medio * DT))[0]
    print("campos cada %d pasos (t = %.3f); el del medio de la ventana es el %04d"
          % (tstep, tstep * DT, indice_medio), flush=True)

    # la malla vertical es la misma en todas las corridas
    fase2.NZ, fase2.CZ = NZ, CZ
    fase2.NZFIS = NZ - CZ
    # se reemplaza la condición inicial de la puerta por la celular 3D
    # rmt y rmq ya están declaradas REAL(KIND=GP) en el ámbito donde se
    # incluye initialv.f90 (specter.fpp:104-105), así que se usan tal cual.
    fase2.IC_SHEAR = IC_CELULAR.replace("BSEG", "%.4f_GP" % B_SEGUNDO_MODO)

    res = {"caso": args.caso, "eps": caso["eps"], "Lz": LZ, "nu": NU, "dt": DT,
           "pasos": pasos, "t_final": pasos * DT, "t_mem": T_MEM, "tstep": tstep,
           "superficie": bool(args.superficie),
           "ventana": [VENT[0] * pasos * DT, VENT[1] * pasos * DT],
           "delta_objetivo": list(caso["delta_objetivo"]),
           "nz": NZ, "cz": CZ, "ord": fase2.ORD, "u0s": list(U0S),
           "u0_referencia": U0_REFERENCIA, "resoluciones": list(RESOLUCIONES),
           "corridas": {}}

    # Se ANEXA sobre lo ya calculado y se guarda DESPUÉS DE CADA MALLA. Cada malla
    # cuesta varios minutos y esta máquina ya mató dos corridas por falta de memoria;
    # perder una malla entera porque el guardado estaba al final es un defecto.
    if os.path.exists(salida):
        try:
            with open(salida) as fh:
                previo = json.load(fh)
            if all(previo.get(c) == res[c] for c in
                   ("Lz", "nu", "dt", "pasos", "nz", "cz", "u0s", "tstep")):
                res["corridas"].update(previo.get("corridas", {}))
                res.setdefault("lambda_referencia", {}).update(
                    previo.get("lambda_referencia", {}))
                # mallas a medias: se guarda después de CADA corrida, porque el
                # vigilante de memoria de la laptop puede matar el proceso en
                # cualquier momento y una malla de 32x32 son ~40 min
                res["parciales"] = previo.get("parciales", {})
                print("anexando sobre: %s%s" % (
                    ", ".join(sorted(res["corridas"])),
                    "".join("; %s a medias (%d corridas)" % (k, len(v))
                            for k, v in res["parciales"].items())), flush=True)
            else:
                print("[aviso] el JSON previo es de otros parámetros; se ignora",
                      flush=True)
        except (ValueError, OSError) as exc:
            print("[aviso] no se pudo leer %s (%s)" % (salida, exc), flush=True)

    for nxy in resoluciones:
        if "n%d" % nxy in res["corridas"]:
            print("\n=== malla %dx%d ya calculada, se saltea" % (nxy, nxy), flush=True)
            continue
        fase2.NX = fase2.NY = nxy
        bindir = fase2.construir(os.path.join(scratch, "arbol_%d" % nxy))
        print("\n=== malla %dx%dx%d compilada" % (nxy, nxy, NZ), flush=True)

        clave = "n%d" % nxy
        filas = list(res.get("parciales", {}).get(clave, []))
        hechas = {f["u0"] for f in filas}
        for u0 in (U0_REFERENCIA,) + U0S:
            if u0 in hechas:
                print("   u0=%.1e ya calculada, se saltea" % u0, flush=True)
                continue
            etiqueta = "n%d_u%.0e" % (nxy, u0)
            _, _, dirrun = fase2.correr(
                bindir, etiqueta, bcsta="noslip", bcend=fase2.CADENA_LIBRE,
                vparam1=VPARAM1, dt=DT, step=pasos, cstep=CSTEP, lz=LZ, nu=NU, u0=u0)
            t, e, w = fase2.leer_balance(dirrun)

            t0, t1 = VENT[0] * pasos * DT, VENT[1] * pasos * DT
            lam, disp = tasa_en_ventana(t, e, t0, t1)
            # amplitud en el MEDIO de la ventana, de balance.txt; k y alfa_eff del
            # campo escrito ahí (ver diagnostico_campo: el <omega^2> de balance.txt
            # no sirve)
            tm = 0.5 * (t0 + t1)
            im = int(np.argmin(np.abs(t - tm)))
            i0 = int(np.argmin(np.abs(t - t0)))
            i1 = int(np.argmin(np.abs(t - t1)))
            k_bal = float(np.sqrt(w[im] / e[im]))
            dg = diagnostico_campo(dirrun, indice_medio, nxy, nxy, NZ - CZ, LZ, NU)
            k_h, k_ef = dg["k_h"], dg["k_3d"]
            # el campo y balance.txt tienen que contar la misma energía a ese instante
            t_medio = (indice_medio - 1) * tstep * DT
            im_c = int(np.argmin(np.abs(t - t_medio)))
            if abs(dg["v2_campo"] / e[im_c] - 1.0) > 2e-2:
                print("   [aviso] <v2> del campo %.4e vs balance %.4e en t=%.3f"
                      % (dg["v2_campo"], e[im_c], t[im_c]), flush=True)
            fila = {
                "u0": u0, "lambda": lam, "dispersion": disp,
                "k_efectivo": k_h, "k_3d_campo": k_ef, "k_3d_balance": k_bal,
                "alfa_eff_sobre_alfa": dg["alfa_eff_sobre_alfa"],
                "w2_sobre_v2": dg["w2_sobre_v2"],
                "sup_sobre_prom_medio": dg["rms_superficie"] / dg["rms_promedio"],
                "v2_inicial": float(e[0]), "v2_medio": float(e[im]),
                "v2_ventana": [float(e[i0]), float(e[i1])],
                "v2_final": float(e[-1]),
                "crece": bool(e[-1] > e[0]),
            }
            extra = ""
            if args.superficie:
                ps = pendientes_superficie(dirrun, indices_t, nxy, nxy, NZ - CZ, LZ, NU,
                                           t0, t1)
                if ps is not None:
                    fila.update(ps)
                    extra = ("\n            superficie: lambda_sup=%.6e  lambda_prom=%.6e  "
                             "lambda_3d=%.6e (balance %.6e)  u(h)/<u>: %.4f -> %.4f  "
                             "(%d campos)"
                             % (ps["lambda_superficie"], ps["lambda_promedio"],
                                ps["lambda_3d"], lam, ps["sup_sobre_prom_inicio"],
                                ps["sup_sobre_prom_fin"], ps["n_campos"]))
            filas.append(fila)
            print("   u0=%.1e  lambda=%.6e  disp=%.2e  k_h=%.2f  k_3d=%.2f (bal %.2f)  "
                  "alfa_eff/alfa=%.4f  <v2>: %.3e -> %.3e%s%s"
                  % (u0, lam, disp, k_h, k_ef, k_bal, dg["alfa_eff_sobre_alfa"],
                     e[0], e[-1],
                     "   [CRECE: inestable]" if e[-1] > e[0] else "", extra), flush=True)
            del t, e, w
            res.setdefault("parciales", {})[clave] = filas
            os.makedirs(os.path.dirname(salida), exist_ok=True)
            with open(salida, "w") as f:
                json.dump(res, f, indent=1)

        ref = filas[0]["lambda"]
        for f in filas:
            f["lambda_sobre_referencia"] = f["lambda"] / ref
            if "lambda_superficie" in f:
                # las dos contra la MISMA referencia lineal (la de balance.txt): en la
                # corrida lineal lambda_sup = lambda_prom = lambda, así que el cociente
                # sup/prom es lo que el PIV mediría de más o de menos
                f["lambda_superficie_sobre_referencia"] = f["lambda_superficie"] / ref
                f["lambda_promedio_sobre_referencia"] = f["lambda_promedio"] / ref
                f["lambda_superficie_sobre_promedio"] = (f["lambda_superficie"]
                                                         / f["lambda_promedio"])
            # delta CALIBRADO, con el k horizontal del propio campo y h = Lz. Para esta condición inicial <v^2> = <F^2>|grad_perp
            # psi|^2_rms = (1/2)|grad_perp psi|^2_rms, y la velocidad que entra en la
            # clausura es el PROMEDIO VERTICAL, U = <F>|grad_perp psi|_rms con
            # <F> = 2/pi. De ahí U = (2/pi) sqrt(2 <v^2>), y delta = U k h^2 / nu.
            def _delta(v2):
                U = (2.0 / math.pi) * math.sqrt(2.0 * max(v2, 0.0))
                return U * f["k_efectivo"] * LZ ** 2 / NU
            f["delta"] = _delta(f["v2_medio"])
            f["delta_rango"] = [_delta(v) for v in f["v2_ventana"]]
        res["corridas"][clave] = filas
        res.setdefault("lambda_referencia", {})[clave] = ref
        res.get("parciales", {}).pop(clave, None)
        os.makedirs(os.path.dirname(salida), exist_ok=True)
        with open(salida, "w") as f:
            json.dump(res, f, indent=1)
        print("   guardado parcial (%s)" % ", ".join(sorted(res["corridas"])),
              flush=True)

    os.makedirs(os.path.dirname(salida), exist_ok=True)
    with open(salida, "w") as f:
        json.dump(res, f, indent=1)

    if len(res["corridas"]) < len(resoluciones):
        print("\nfalta(n) la(s) malla(s): %s"
              % ", ".join(str(n) for n in resoluciones
                          if "n%d" % n not in res["corridas"]))
        print("escrito %s" % os.path.relpath(salida, RAIZ))
        return 0

    print("\nRESUMEN  lambda/lambda_lineal  (caso %s, eps=%.2f)" % (args.caso, caso["eps"]))
    print("  u0        delta    " + "".join("  %dx%d      " % (n, n) for n in resoluciones)
          + ("   | lambda_sup/lin  lambda_prom/lin  sup/prom   u(h)/<u> ini->fin"
             if args.superficie else ""))
    for i, u0 in enumerate(U0S):
        fila = "  %.2e " % u0
        for j, n in enumerate(resoluciones):
            f = res["corridas"]["n%d" % n][i + 1]
            if j == 0:
                fila += " %6.2f " % f["delta"]
            fila += "  %9.6f%s" % (f["lambda_sobre_referencia"],
                                   "*" if f["crece"] else " ")
        if args.superficie and "lambda_superficie_sobre_referencia" in f:
            fila += "   | %9.5f      %9.5f      %7.4f    %.4f -> %.4f" % (
                f["lambda_superficie_sobre_referencia"],
                f["lambda_promedio_sobre_referencia"],
                f["lambda_superficie_sobre_promedio"],
                f["sup_sobre_prom_inicio"], f["sup_sobre_prom_fin"])
        print(fila)
    print("  (* = la energía creció: corrida inestable, no usable)")
    print("\nescrito %s" % os.path.relpath(salida, RAIZ))

    if temporal and os.environ.get("SPECTER_KEEP") != "1":
        import shutil
        shutil.rmtree(scratch, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
