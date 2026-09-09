"""Calibra el estimador subpíxel con desplazamientos conocidos, sobre imágenes reales.

Por qué hace falta
------------------
El barrido de Δt (`02_barrido_dt.py`) ajusta el modelo

    U_med(n)² = U_verdadera² + ( σ / (Δt(n)·s) )²

y de ahí despeja σ. Ese modelo supone que **todo** lo que depende de n es ruido de
media cero. Es falso en este régimen. Los datos del barrido de ventanas muestran que
la velocidad medida *sube* con n entre gap 5 y gap 10 en las cuatro ventanas
(64/32/16/8 px), y el modelo exige que baje: el ajuste devuelve σ² < 0.

La causa es **peak locking**: con desplazamientos de una fracción de píxel, el
estimador subpíxel gaussiano tira el resultado hacia el entero más cercano, que acá
es el cero. Medido sobre los campos guardados en
`GW-AI-course/piv-s0008-forzado/outputs/sweep_multipass_s0008_20260902`, el exceso de
masa en desplazamiento nulo vale 2.39 cuando el desplazamiento medio es 0.172 px y
baja a 1.07 cuando es 0.565 px.

Eso deja los tres valores de σ del proyecto contaminados, en direcciones opuestas:

  0.169 px  barrido de Δt      — sus puntos de n chico mezclan ruido (que infla) con
                                 peak locking (que deprime); el σ ajustado absorbe los dos
  0.083 px  dos Δt sobre .vec  — además con los filtros de validación de la GUI encima
  0.030 px  desplazamiento nulo — el caso peor: con desplazamiento verdadero cero, el
                                 estimador tira todo al cero. Es una cota inferior deprimida

Qué hace este script
--------------------
Toma un cuadro real de la campaña, lo desplaza una cantidad **conocida y exacta**, y
mide qué devuelve el multipaso. Con el desplazamiento verdadero conocido, el sesgo y
la dispersión se separan sin suponer ningún modelo:

    sesgo(d)  = <u> - d          la curva de peak locking
    σ_u(d)    = sd(u)            dispersión con desplazamiento verdadero d
    σ_v(d)    = sd(v)            dispersión con desplazamiento verdadero 0

Se barre d de 0 a 1.5 px, que cubre un período completo del peak locking (es periódico
en 1 px) y todo el régimen del barrido de Δt, cuyo n=1 mide 0.201 px.

Lo que este σ NO incluye
------------------------
Un par sintético está perfectamente correlacionado: no hay movimiento fuera del plano,
ni cambio de iluminación entre cuadros, ni ruido de lectura independiente. Así que el σ
de acá es **sólo la contribución de la correlación y del estimador subpíxel**, y es una
cota inferior del σ que sufre un par real. Sirve para lo que se necesita: separar el
sesgo del ruido y decidir si el modelo del barrido es aplicable.

El control de interpolación
---------------------------
Desplazar por Fourier introduce su propio artefacto: la imagen no es periódica, así que
el corrimiento envuelve por los bordes y puede timbrar. Dos defensas:

  1. se descartan los vectores a menos de `MARGEN_PX` del borde;
  2. se repite un subconjunto de los corrimientos con interpolación por splines cúbicos,
     que tiene artefactos distintos. Si las dos dan la misma curva de sesgo, el sesgo es
     del estimador de PIV y no de la interpolación. Ese es el punto del control.

Verificación interna: a d = 0.0 y d = 1.0 px, ambos enteros, el estimador no tiene nada
que interpolar y debe devolver exactamente 0 y 1. Los dos están en el barrido.

Uso
---
    python 05_corrimiento_sintetico.py [--cuadros 2] [--carpeta med_S0008]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
from scipy.ndimage import fourier_shift, shift as shift_spline

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import piv_comun as pc  # noqa: E402
from importlib import import_module  # noqa: E402

bd = import_module("02_barrido_dt")

AQUI = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Corrimientos verdaderos, en píxeles. Cubren un período de peak locking (1 px) más
# medio período, para ver que la curva se repite y no es un artefacto de un tramo.
CORRIMIENTOS = [round(0.1 * i, 2) for i in range(16)]          # 0.0 .. 1.5
CORRIMIENTOS_CONTROL = [0.0, 0.2, 0.5, 0.8, 1.0]             # con splines; subconjunto de CORRIMIENTOS
MARGEN_PX = 64          # se descartan los vectores a menos de esto del borde


def desplazar_fourier(img: np.ndarray, dx: float) -> np.ndarray:
    """B(x) = A(x - dx): un rasgo en x de A aparece en x+dx de B, así que PIV mide +dx."""
    return np.fft.ifftn(fourier_shift(np.fft.fftn(img), (0.0, dx))).real


def desplazar_spline(img: np.ndarray, dx: float) -> np.ndarray:
    """Mismo corrimiento por splines cúbicos. Control de [la interpolación de] Fourier."""
    return shift_spline(img, (0.0, dx), order=3, mode="reflect")


def medir(a: np.ndarray, b: np.ndarray, d_verdadero: float) -> dict:
    """Corre el multipaso sobre el par y compara contra el desplazamiento conocido."""
    campo = bd.piv_multipaso(a, b)
    cx = campo["cx"]
    dentro = (cx >= MARGEN_PX) & (cx <= a.shape[0] - MARGEN_PX)
    sel = np.ix_(dentro, dentro)
    u, v, s2n = campo["u"][sel], campo["v"][sel], campo["s2n"][sel]
    ok = (s2n >= bd.S2N_UMBRAL) & np.isfinite(u) & np.isfinite(v)
    u, v = u[ok], v[ok]
    del campo
    return {
        "d_verdadero_px": d_verdadero,
        "u_medio_px": float(u.mean()),
        "sesgo_px": float(u.mean() - d_verdadero),
        "sigma_u_px": float(u.std()),
        "v_medio_px": float(v.mean()),
        "sigma_v_px": float(v.std()),
        "frac_valida": float(ok.mean()),
        "n_vectores": int(ok.sum()),
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--carpeta", default="med_S0008")
    p.add_argument("--cuadros", type=int, default=2,
                   help="cuántos cuadros distintos se usan como imagen base")
    args = p.parse_args()

    carpeta = os.path.join(pc.RAIZ, args.carpeta)
    archivos = pc.tifs(carpeta)
    cih = pc.leer_cih(carpeta)
    escala = pc.escala_px_por_m()
    tipo, desc, _ = pc.MEDICIONES[args.carpeta]

    print(f"# {args.carpeta}: {desc}")
    print(f"# {cih['ancho']}x{cih['alto']}, {cih['bits']} bits, {len(archivos)} TIF")
    print(f"# PIV: multipaso {'->'.join(map(str, bd.VENTANAS))} px, solape {bd.SOLAPE:.0%}, "
          f"subpixel {bd.SUBPIXEL}, s2n >= {bd.S2N_UMBRAL}")
    print(f"# margen descartado: {MARGEN_PX} px por borde")
    print(f"# NOTA: par sintetico, correlacion perfecta -> sigma es COTA INFERIOR\n")

    # Cuadros base, repartidos a lo largo del registro.
    indices = np.linspace(0, len(archivos) - 1, args.cuadros + 2).astype(int)[1:-1]
    t0 = time.time()

    print(f"{'d_real':>7}{'u_medio':>10}{'sesgo':>9}{'sigma_u':>9}"
          f"{'v_medio':>9}{'sigma_v':>9}{'valid':>8}{'t':>7}")
    print("-" * 68)

    fourier, control = [], []
    for d in CORRIMIENTOS:
        acum = []
        for i in indices:
            a = pc.leer_tif(archivos[i]).astype(np.float64)
            b = desplazar_fourier(a, d)
            acum.append(medir(a, b, d))
            del a, b
        fila = {k: float(np.mean([x[k] for x in acum])) for k in acum[0]}
        fila["sesgo_sd_px"] = float(np.std([x["sesgo_px"] for x in acum]))
        fila["n_cuadros"] = len(acum)
        fourier.append(fila)
        print(f"{d:>7.2f}{fila['u_medio_px']:>10.4f}{fila['sesgo_px']:>+9.4f}"
              f"{fila['sigma_u_px']:>9.4f}{fila['v_medio_px']:>+9.4f}"
              f"{fila['sigma_v_px']:>9.4f}{fila['frac_valida']:>8.1%}"
              f"{time.time()-t0:>6.0f}s", flush=True)

    print(f"\ncontrol con splines cubicos (misma cuenta, otra interpolacion)")
    print(f"{'d_real':>7}{'u_medio':>10}{'sesgo':>9}{'sigma_u':>9}{'dif_vs_fourier':>16}")
    print("-" * 51)
    for d in CORRIMIENTOS_CONTROL:
        a = pc.leer_tif(archivos[indices[0]]).astype(np.float64)
        b = desplazar_spline(a, d)
        fila = medir(a, b, d)
        del a, b
        par = next((f for f in fourier if abs(f["d_verdadero_px"] - d) < 1e-9), None)
        dif = fila["sesgo_px"] - par["sesgo_px"] if par else float("nan")
        fila["dif_sesgo_vs_fourier_px"] = dif
        control.append(fila)
        print(f"{d:>7.2f}{fila['u_medio_px']:>10.4f}{fila['sesgo_px']:>+9.4f}"
              f"{fila['sigma_u_px']:>9.4f}{dif:>+16.4f}", flush=True)

    # ---- lectura del resultado -------------------------------------------------
    enteros = [f for f in fourier if abs(f["d_verdadero_px"] % 1.0) < 1e-9]
    peor = max(fourier, key=lambda f: abs(f["sesgo_px"]))
    sig_ent = float(np.mean([f["sigma_u_px"] for f in enteros])) if enteros else float("nan")
    difs = [abs(c["dif_sesgo_vs_fourier_px"]) for c in control]

    print(f"\n{'='*68}")
    print(f"verificacion en enteros (d=0, 1): sesgo maximo "
          f"{max(abs(f['sesgo_px']) for f in enteros):.4f} px  -> deberia ser ~0")
    print(f"sesgo maximo del barrido:        {peor['sesgo_px']:+.4f} px en d={peor['d_verdadero_px']:.2f}")
    print(f"sigma_u en corrimiento entero:   {sig_ent:.4f} px")
    print(f"control splines: |dif| maxima    {max(difs):.4f} px  -> si es chica, el sesgo es del PIV")
    print(f"{'='*68}")
    print(f"comparar contra: barrido de dt 0.169 px | desp. nulo 0.030 px | Foucaut 0.087 px")

    salida = {
        "carpeta": args.carpeta, "descripcion": desc,
        "escala_px_por_m": escala, "fps": cih["fps"],
        "ventanas": list(bd.VENTANAS), "solape": bd.SOLAPE,
        "subpixel": bd.SUBPIXEL, "s2n_umbral": bd.S2N_UMBRAL,
        "margen_px": MARGEN_PX, "cuadros_base": indices.tolist(),
        "fourier": fourier, "control_splines": control,
        "sesgo_max_px": peor["sesgo_px"], "sesgo_max_en_d": peor["d_verdadero_px"],
        "sigma_u_en_entero_px": sig_ent,
        "control_dif_max_px": max(difs),
        "advertencia": ("par sintetico: correlacion perfecta, sin movimiento fuera de "
                        "plano ni cambio de iluminacion. sigma es cota inferior."),
    }
    destino = os.path.join(AQUI, "salidas", "tablas", "corrimiento_sintetico.json")
    with open(destino, "w", encoding="utf-8") as fh:
        json.dump(salida, fh, indent=1, ensure_ascii=False)
    print(f"\n-> {destino}")


if __name__ == "__main__":
    main()
