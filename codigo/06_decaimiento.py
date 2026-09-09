#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
El decaimiento medido, corregido por el piso de ruido de PIV.

Por qué hace falta corregir. Las series `.vec` archivadas de los cuatro decaimientos
(`Deca1..4_destino`) están armadas con pares de cuadros CONSECUTIVOS, Δt = 1/60 s, que
es el peor Δt posible: el desplazamiento verdadero al arrancar el decaimiento es de
~0,2 px por cuadro y el piso de ruido medido en el capítulo de instrumento es de
~0,17 px. Es decir que el ruido es del orden de la señal, y la velocidad leída
directamente de esos archivos está inflada. Además esas series tienen 66 % de vectores
reemplazados, contra 2,5 % del multipaso propio sobre los mismos cuadros.

Qué hace este script, en tres partes:

  A. Lee las series `.vec` completas (308 pares por registro, uno cada 10 cuadros) y
     saca u_rms(t) sobre los vectores VÁLIDOS. Da la forma de la curva, gratis, pero
     inflada por ruido.
  B. En unos pocos instantes por registro rearma los pares desde los TIF crudos con
     varias separaciones n y ajusta

         desp_rms(n)^2 = (n d1)^2 + sigma^2

     que separa la velocidad verdadera d1 del ruido sigma. Es el método de dos Δt del
     capítulo de instrumento, aplicado ahora dentro de una ventana corta del
     decaimiento: lícito porque el flujo cambia en ~15 s y la ventana dura <0,5 s.
  C. Con el campo de mejor relación señal-ruido de cada instante calcula el espectro
     de energía, para sacar de los datos la escala que fija delta en [P-01], en vez de
     heredarla del espaciado de los imanes.

Memoria: no se guarda ningún campo completo salvo el del espectro del instante en
curso. Los `.vec` se leen de a uno y sólo sobreviven escalares.

Uso:
    /home/lucia/miniforge3/envs/piv-dt/bin/python codigo/06_decaimiento.py
    /home/lucia/miniforge3/envs/piv-dt/bin/python codigo/06_decaimiento.py --rapido
"""

import argparse
import glob
import importlib.util
import json
import math
import os
import sys

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(AQUI)
sys.path.insert(0, AQUI)

import piv_comun as pc                                    # noqa: E402
from celda import CAMPANA_02_06_25 as CELDA               # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "barrido", os.path.join(AQUI, "02_barrido_dt.py"))
barrido = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(barrido)

SALIDA = os.path.join(RAIZ, "salidas", "tablas", "decaimiento.json")

# med_S00nn  ->  carpeta con la serie .vec ya procesada
REGISTROS = {
    "med_S0003": "Deca1_destino",
    "med_S0005": "Deca2_destino",
    "med_S0006": "Deca3_destino",
    "med_S0009": "Deca4_destino",
}

# Separaciones entre cuadros del ajuste de ruido. El rango tiene que cubrir tanto el
# arranque (rápido: n grande viola la regla del cuarto) como la cola (lenta: n chico
# no despega del ruido).
NS = (1, 2, 4, 8, 16, 32)
FRAC_VALIDA_MINIMA = 0.90


# --------------------------------------------------------------------------------
# A. La serie completa, de los .vec archivados
# --------------------------------------------------------------------------------

def serie_vec(destino, fps):
    """u_rms(t) sobre los vectores válidos de una serie .vec archivada.

    El paso temporal entre pares se deduce del propio conteo: los TIF presentes en la
    carpeta son pares consecutivos separados por un intervalo fijo.
    """
    fs = sorted(glob.glob(os.path.join(pc.RAIZ, destino, "*.vec")))
    if not fs:
        raise RuntimeError("sin .vec en %s" % destino)
    tifs = sorted(glob.glob(os.path.join(pc.RAIZ, destino, "*.tif")))
    paso = (len(tifs) // len(fs)) * 5 if len(tifs) else 10   # 2 cuadros por par
    paso = 10 if len(tifs) == 2 * len(fs) else paso
    t, u, frac = [], [], []
    for i, f in enumerate(fs):
        v = pc.leer_vec(f)
        t.append(i * paso / fps)
        u.append(pc.urms(v))
        frac.append(v["frac_valida"])
        del v
    return (np.array(t), np.array(u), np.array(frac), len(fs), paso)


# --------------------------------------------------------------------------------
# B. La corrección de ruido, rearmando pares desde los TIF
# --------------------------------------------------------------------------------

def ajuste_dos_dt(ns, desp2):
    """desp^2 = d1^2 n^2 + sigma^2, por mínimos cuadrados en (n^2, 1).

    Devuelve d1 (px por cuadro), sigma (px, suma de las dos componentes) y R^2. Un
    sigma^2 ajustado negativo se informa como tal y NO se recorta a cero: es la señal
    de que el modelo no describe esos puntos, y taparlo fabricaría una barra de error
    optimista. Ver [H-06], donde un sigma^2 negativo delató el reemplazo de vectores.
    """
    x = np.asarray(ns, float) ** 2
    y = np.asarray(desp2, float)
    A = np.vstack([x, np.ones_like(x)]).T
    (m, b), *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ np.array([m, b])
    ss = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - float(np.sum((y - pred) ** 2)) / ss if ss > 0 else float("nan")
    d1 = math.sqrt(m) if m > 0 else float("nan")
    sigma = math.sqrt(b) if b > 0 else -math.sqrt(-b)
    return d1, sigma, r2, float(m), float(b)


def calibrar_instante(tifs, i0, ns, escala, fps, con_campo=False):
    """Rearma el par (i0, i0+n) para cada n y ajusta el modelo de dos Δt."""
    a = pc.leer_tif(tifs[i0])
    puntos, campo_mejor, n_mejor = [], None, None
    for n in ns:
        if i0 + n >= len(tifs):
            continue
        b = pc.leer_tif(tifs[i0 + n])
        campo = barrido.piv_multipaso(a, b)
        m = barrido.metricas(campo)
        m["n"] = n
        puntos.append(m)
        # el campo para el espectro: el de mayor desplazamiento que siga siendo
        # válido y por debajo de la regla del cuarto (ventana final de 16 px)
        if (con_campo and m["frac_valida"] >= FRAC_VALIDA_MINIMA
                and m["desp_rms_px"] < 4.0):
            campo_mejor = {k: campo[k] for k in ("u", "v", "s2n")}
            n_mejor = n
        del b, campo
    del a

    usa = [p for p in puntos if p["frac_valida"] >= FRAC_VALIDA_MINIMA]
    if len(usa) < 3:
        return None, campo_mejor, n_mejor
    d1, sigma, r2, m, b = ajuste_dos_dt([p["n"] for p in usa],
                                        [p["desp_rms_px"] ** 2 for p in usa])

    # Estimación ROBUSTA, y es la que se usa. A n grande el ruido es despreciable por
    # construcción: con desplazamiento n*d1 de varios px y sigma ~ 0,2 px, la
    # corrección relativa es (sigma/(n d1))^2/2, que se informa abajo. No depende de
    # que el modelo de dos Δt valga en todo el rango, y en particular no depende de
    # los n chicos, donde el desplazamiento es de dos décimas de píxel y el
    # estimador subpíxel sesga hacia el entero (peak locking).
    p_max = max((p for p in usa if p["desp_rms_px"] < 4.0),
                key=lambda p: p["n"], default=None)
    if p_max is None:
        p_max = min(usa, key=lambda p: p["desp_rms_px"])
    d1_dir = p_max["desp_rms_px"] / p_max["n"]
    correccion = (abs(b) / p_max["desp_rms_px"] ** 2) / 2.0

    res = {
        "t_s": i0 / fps,
        "d1_px_por_cuadro": d1_dir,
        "n_de_d1": p_max["n"],
        "desp_en_n_px": p_max["desp_rms_px"],
        "correccion_relativa_por_ruido": correccion,
        "u_rms_superficie_m_s": d1_dir * fps / escala,
        "ajuste_dos_dt": {
            "d1_px_por_cuadro": d1,
            "sigma_px": sigma,
            "sigma_por_componente_px": sigma / math.sqrt(2.0) if sigma > 0 else sigma,
            "sigma2_negativo": b < 0,
            "R2": r2,
            "u_rms_superficie_m_s": d1 * fps / escala,
        },
        "n_usados": [p["n"] for p in usa],
        "n_excluidos": [p["n"] for p in puntos if p not in usa],
        "puntos": puntos,
    }
    return res, campo_mejor, n_mejor


# --------------------------------------------------------------------------------
# C. El espectro, para sacar k de los datos
# --------------------------------------------------------------------------------

def espectro(campo, paso_m, s2n_umbral):
    """Espectro de energía radial de un campo de desplazamiento.

    Los vectores inválidos se ponen en cero DESPUÉS de restar la media, que con >95 %
    de válidos es una perturbación chica. Se aplica una ventana de Hann en las dos
    direcciones porque el campo no es periódico y el salto de borde metería una ley de
    potencias espuria.

    Devuelve k [1/m] y E(k) en unidades arbitrarias, más el piso de ruido estimado
    como la mediana del cuarto superior en k, que es donde el espectro de PIV se
    aplana ([V-05], marco de Foucaut).
    """
    u = np.array(campo["u"], float)
    v = np.array(campo["v"], float)
    ok = (campo["s2n"] >= s2n_umbral) & np.isfinite(u) & np.isfinite(v)
    u -= u[ok].mean()
    v -= v[ok].mean()
    u[~ok] = 0.0
    v[~ok] = 0.0

    ny, nx = u.shape
    w = np.outer(np.hanning(ny), np.hanning(nx))
    fu = np.fft.fft2(u * w)
    fv = np.fft.fft2(v * w)
    pot = (np.abs(fu) ** 2 + np.abs(fv) ** 2) / (ny * nx) ** 2
    del u, v, fu, fv, w

    kx = 2.0 * np.pi * np.fft.fftfreq(nx, d=paso_m)
    ky = 2.0 * np.pi * np.fft.fftfreq(ny, d=paso_m)
    K = np.hypot(*np.meshgrid(kx, ky, indexing="xy"))
    dk = 2.0 * np.pi / (nx * paso_m)
    nb = int(K.max() / dk)
    idx = np.clip((K / dk).astype(int), 0, nb)
    E = np.bincount(idx.ravel(), weights=pot.ravel(), minlength=nb + 1)
    cuenta = np.bincount(idx.ravel(), minlength=nb + 1)
    k = (np.arange(nb + 1) + 0.5) * dk
    del K, idx, pot

    E = E[1:]            # se descarta el modo medio
    k = k[1:]
    cuenta = cuenta[1:]
    # El ruido de PIV es blanco en el campo de desplazamiento: aporta la MISMA potencia
    # por modo, así que en el espectro sumado por anillos aporta p0 * (modos del
    # anillo), que crece con k. Restar una constante en vez de p0*cuenta deja casi
    # todo el ruido adentro y corre el k pesado hacia arriba.
    cola = slice(int(0.75 * len(E)), None)
    p0 = float(np.median((E / np.maximum(cuenta, 1))[cola]))
    Elimpio = np.maximum(E - p0 * cuenta, 0.0)

    def pesado(w_):
        s = w_.sum()
        return float((k * w_).sum() / s) if s > 0 else float("nan")

    return {
        "k_1_m": k.tolist(),
        "E": E.tolist(),
        "modos_por_anillo": cuenta.tolist(),
        "piso_por_modo": p0,
        "k_pico_1_m": float(k[int(np.argmax(Elimpio))]),
        "k_pesado_1_m": pesado(Elimpio),
        "k_pesado_sin_restar_piso_1_m": pesado(E),
        "paso_grilla_m": paso_m,
    }


# --------------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rapido", action="store_true",
                    help="3 instantes por registro en vez de 8, para probar la cadena")
    ap.add_argument("--registros", nargs="*", default=sorted(REGISTROS))
    args = ap.parse_args()

    fps = CELDA.camara.fps
    escala = CELDA.camara.escala_px_por_m
    n_inst = 3 if args.rapido else 8
    # instantes repartidos sobre el registro, más densos al principio, que es donde
    # delta es peor y donde el perfil vertical todavía no relajó
    fracs = np.linspace(0.0, 0.85, n_inst) ** 1.5

    out = {"celda": {"h_m": CELDA.h, "nu_m2_s": CELDA.fluido.nu,
                     "alfa_1_s": CELDA.alfa(), "fps": fps,
                     "escala_px_por_m": escala,
                     "k_forzado_1_m": CELDA.forzado.k_fundamental()},
           "registros": {}}

    # Se ANEXA sobre lo que ya haya. Cada registro tarda ~12 min y el guardado es
    # parcial, así que una corrida interrumpida —por ejemplo porque la máquina se
    # quedó sin memoria, que ya pasó— no puede costar los registros que sí terminaron.
    if os.path.exists(SALIDA):
        try:
            with open(SALIDA) as fh:
                previo = json.load(fh)
            out["registros"].update(previo.get("registros", {}))
            print("anexando sobre %d registro(s) ya calculado(s): %s"
                  % (len(out["registros"]), ", ".join(sorted(out["registros"]))),
                  flush=True)
        except (ValueError, OSError) as exc:
            print("[aviso] no se pudo leer %s (%s); se empieza de cero"
                  % (SALIDA, exc), flush=True)

    for med in args.registros:
        destino = REGISTROS[med]
        print("=== %s (%s)" % (med, destino), flush=True)

        t, u, frac, npares, paso = serie_vec(destino, fps)
        print("    serie .vec: %d pares, cada %d cuadros, válidos %.1f%%–%.1f%%"
              % (npares, paso, 100 * frac.min(), 100 * frac.max()), flush=True)

        tifs = pc.tifs(os.path.join(pc.RAIZ, med))
        print("    %d TIF crudos" % len(tifs), flush=True)

        calib, espectros = [], []
        for f in fracs:
            i0 = int(f * (len(tifs) - max(NS) - 1))
            res, campo, n_esp = calibrar_instante(
                tifs, i0, NS, escala, fps, con_campo=True)
            if res is None:
                print("    t=%.2f s: sin puntos suficientes" % (i0 / fps), flush=True)
                continue
            aj = res["ajuste_dos_dt"]
            print("    t=%5.2f s  u=%.3f mm/s (n=%d, %.2f px, corr %.2f%%)  "
                  "| dos-Δt: u=%.3f mm/s sigma=%+.3f px R2=%.4f%s "
                  "| crudo n=1: %.3f mm/s"
                  % (res["t_s"], 1e3 * res["u_rms_superficie_m_s"], res["n_de_d1"],
                     res["desp_en_n_px"], 100 * res["correccion_relativa_por_ruido"],
                     1e3 * aj["u_rms_superficie_m_s"], aj["sigma_px"], aj["R2"],
                     " [sigma2<0]" if aj["sigma2_negativo"] else "",
                     1e3 * res["puntos"][0]["desp_rms_px"] * fps / escala),
                  flush=True)
            if campo is not None:
                e = espectro(campo, 4.0 / escala, barrido.S2N_UMBRAL)
                e["t_s"] = res["t_s"]
                e["n"] = n_esp
                espectros.append(e)
                print("        espectro (n=%d): k_pico=%.1f 1/m  k_pesado=%.1f 1/m"
                      % (n_esp, e["k_pico_1_m"], e["k_pesado_1_m"]), flush=True)
                del campo
            calib.append(res)

        out["registros"][med] = {
            "destino": destino,
            "serie_vec": {"t_s": t.tolist(), "u_rms_m_s": u.tolist(),
                          "frac_valida": frac.tolist(),
                          "paso_cuadros": paso, "dt_cuadros": 1},
            "calibracion": calib,
            "espectros": espectros,
        }

        os.makedirs(os.path.dirname(SALIDA), exist_ok=True)
        with open(SALIDA, "w") as fh:
            json.dump(out, fh, indent=1)
        print("    guardado parcial", flush=True)

    print("escrito %s" % os.path.relpath(SALIDA, RAIZ))
    return 0


if __name__ == "__main__":
    sys.exit(main())
