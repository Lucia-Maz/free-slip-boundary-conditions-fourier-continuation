"""Barrido de Δt sobre una corrida forzada, rearmando los pares desde los TIF crudos.

La idea
-------
`med_S0008` (Forzado 2.1 A) tiene 3072 cuadros consecutivos a 60 fps. Eso permite
construir pares con cualquier separación Δt = n/60 s sin medir nada nuevo ([D-05]).
El flujo es estacionario, así que la velocidad verdadera **no depende de n**. Todo lo
que dependa de n es error de medición, y eso es lo que el barrido aísla.

El modelo
---------
Si el ruido de PIV es un error de desplazamiento de desviación σ (en píxeles),
independiente de Δt porque es una propiedad de la correlación y no del flujo, la
velocidad medida cumple

    U_med(n)² = U_verdadera² + ( σ / (Δt(n) · s) )²        con Δt(n) = n/fps

es decir, una **recta** si se grafica U_med² contra 1/n². La ordenada al origen es la
velocidad verdadera y la pendiente da σ. A n grande la recta se rompe: aparece pérdida
de correlación y crecen los vectores inválidos. Los dos extremos del rango útil de Δt
salen de la misma figura.

El σ así obtenido es contrastable contra Foucaut, Carlier & Stanislas 2004
(DOI 10.1088/0957-0233/15/6/003), que mide 0.087 px para ventana de 16×16 px.

Uso
---
    python 02_barrido_dt.py [--carpeta med_S0008] [--pares 12] [--rapido]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
from scipy.ndimage import map_coordinates

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import piv_comun as pc  # noqa: E402

from openpiv import pyprocess  # noqa: E402

AQUI = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Configuración de PIV copiada de las sesiones de Lucía y dejada fija ([D-06]).
VENTANAS = (64, 32, 16)
SOLAPE = 0.75
SUBPIXEL = "gaussian"
S2N_METODO = "peak2peak"
S2N_UMBRAL = 1.05
# Correlación normalizada y resta de la media global. Sin esto el pico se apoya sobre
# un pedestal de continua y el peak2peak deja de discriminar: medido sobre el par
# (5, 10) de med_S0008, s2n mediano 1.17 y 79 % de vectores válidos sin normalizar,
# contra 2.75 y 96,5 % normalizando. Ver [D-13].
NORMALIZAR = True
RESTAR_MEDIA = True


def _grilla_de_centros(lado: int, W: int, paso: int) -> np.ndarray:
    n = (lado - W) // paso + 1
    return W / 2.0 + paso * np.arange(n)


def _deformar(img: np.ndarray, cx: np.ndarray, cy: np.ndarray,
              du: np.ndarray, dv: np.ndarray, signo: float) -> np.ndarray:
    """Desplaza la imagen por el campo predicho, interpolando el campo a cada píxel.

    Deformación simétrica (Scarano 2002, DOI 10.1088/0957-0233/13/1/201): se mueve
    cada imagen media predicción en sentidos opuestos, en vez de mover una entera.
    """
    lado = img.shape[0]
    ejes = np.arange(lado, dtype=np.float32)
    # predicción interpolada a la grilla de píxeles (bilineal por separabilidad)
    iy = np.interp(ejes, cy, np.arange(len(cy)))
    ix = np.interp(ejes, cx, np.arange(len(cx)))
    U = map_coordinates(du, np.meshgrid(iy, ix, indexing="ij"), order=1, mode="nearest")
    V = map_coordinates(dv, np.meshgrid(iy, ix, indexing="ij"), order=1, mode="nearest")
    Y, X = np.meshgrid(ejes, ejes, indexing="ij")
    return map_coordinates(img, [Y + signo * V / 2.0, X + signo * U / 2.0],
                           order=3, mode="nearest").astype(np.float32)


def piv_multipaso(a: np.ndarray, b: np.ndarray) -> dict:
    """Multipaso 64 → 32 → 16 con deformación simétrica de imagen.

    Se implementa acá en vez de llamar a `openpiv.windef` para no depender de una API
    que cambió entre versiones, y para que el único parámetro que varía en el barrido
    sea Δt. Las correlaciones las hace `pyprocess.extended_search_area_piv`, que es
    la parte estable de OpenPIV. Ver [D-08].

    Devuelve el desplazamiento en PÍXELES (dt=1), no en m/s: la conversión a
    velocidad se hace después, y así el ruido queda expresado en las unidades en las
    que es constante.
    """
    if RESTAR_MEDIA:
        a = a - a.mean()
        b = b - b.mean()
    lado = a.shape[0]
    u_ac = v_ac = None
    cx = cy = None
    for W in VENTANAS:
        paso = int(round(W * (1 - SOLAPE)))
        ncx = _grilla_de_centros(lado, W, paso)
        if u_ac is None:
            a_d, b_d = a, b
        else:
            a_d = _deformar(a, cx, cy, u_ac, v_ac, -1.0)
            b_d = _deformar(b, cx, cy, u_ac, v_ac, +1.0)
        u, v, s2n = pyprocess.extended_search_area_piv(
            a_d, b_d, window_size=W, overlap=W - paso, dt=1.0,
            search_area_size=W, subpixel_method=SUBPIXEL,
            sig2noise_method=S2N_METODO, correlation_method="circular",
            normalized_correlation=NORMALIZAR,
        )
        # NO se invierte el signo de v. OpenPIV devuelve v en la convención de
        # índices de fila (y hacia abajo) y el campo realimenta la deformación: con
        # el signo invertido la iteración diverge en vez de converger. Ver [D-14].
        if u_ac is not None:
            # llevar el acumulado a la grilla nueva y sumarle el residuo
            iy = np.interp(ncx, cy, np.arange(len(cy)))
            ix = np.interp(ncx, cx, np.arange(len(cx)))
            gy, gx = np.meshgrid(iy, ix, indexing="ij")
            u = map_coordinates(u_ac, [gy, gx], order=1, mode="nearest") + u
            v = map_coordinates(v_ac, [gy, gx], order=1, mode="nearest") + v
        u_ac, v_ac, cx, cy = u, v, ncx, ncx
        del a_d, b_d
    return {"u": u_ac, "v": v_ac, "s2n": s2n, "cx": cx, "cy": cy}


def metricas(campo: dict) -> dict:
    """Escalares por par. No se guarda ningún campo completo (memoria)."""
    ok = campo["s2n"] >= S2N_UMBRAL
    ok &= np.isfinite(campo["u"]) & np.isfinite(campo["v"])
    u, v = campo["u"][ok], campo["v"][ok]
    d = np.sqrt(u ** 2 + v ** 2)
    return {
        "frac_valida": float(ok.mean()),
        "desp_rms_px": float(np.sqrt(np.mean(u ** 2 + v ** 2))) if ok.any() else np.nan,
        "desp_mediano_px": float(np.median(d)) if ok.any() else np.nan,
        "s2n_mediano": float(np.median(campo["s2n"][ok])) if ok.any() else np.nan,
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--carpeta", default="med_S0008")
    p.add_argument("--pares", type=int, default=12)
    p.add_argument("--rapido", action="store_true", help="pocos n, para probar")
    args = p.parse_args()

    carpeta = os.path.join(pc.RAIZ, args.carpeta)
    cih = pc.leer_cih(carpeta)
    escala = pc.escala_px_por_m()
    archivos = pc.tifs(carpeta)
    tipo, desc, corriente = pc.MEDICIONES[args.carpeta]

    print(f"# {args.carpeta}: {desc}")
    print(f"# cámara: {cih['ancho']}x{cih['alto']}, {cih['bits']} bits, "
          f"{cih['fps']} fps, shutter {cih['shutter']}, {len(archivos)} TIF")
    print(f"# escala: {escala} px/m  ->  {1e3/escala:.4f} mm/px")
    print(f"# PIV: multipaso {'->'.join(map(str, VENTANAS))} px, solape {SOLAPE:.0%}, "
          f"subpixel {SUBPIXEL}, s2n {S2N_METODO} >= {S2N_UMBRAL}")
    assert tipo == "forzado", f"{args.carpeta} no es una corrida forzada; el modelo supone U estacionaria"

    ns = [1, 4, 16, 64] if args.rapido else [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128]
    filas = []
    t0 = time.time()
    for n in ns:
        arranques = np.linspace(0, len(archivos) - n - 1, args.pares).astype(int)
        m = []
        for i in arranques:
            a = pc.leer_tif(archivos[i])
            b = pc.leer_tif(archivos[i + n])
            campo = piv_multipaso(a, b)
            m.append(metricas(campo))
            del a, b, campo
        dt = n / cih["fps"]
        fila = {
            "n": n, "dt_s": dt,
            "desp_rms_px": float(np.mean([x["desp_rms_px"] for x in m])),
            "desp_rms_px_sd": float(np.std([x["desp_rms_px"] for x in m])),
            "frac_valida": float(np.mean([x["frac_valida"] for x in m])),
            "s2n_mediano": float(np.mean([x["s2n_mediano"] for x in m])),
            "n_pares": len(m),
        }
        fila["U_rms_m_s"] = fila["desp_rms_px"] / (dt * escala)
        filas.append(fila)
        print(f"n={n:4d}  dt={dt:6.3f}s  desp={fila['desp_rms_px']:7.3f}±"
              f"{fila['desp_rms_px_sd']:.3f} px  U={fila['U_rms_m_s']*1e3:7.3f} mm/s"
              f"  valid={fila['frac_valida']:5.1%}  s2n={fila['s2n_mediano']:.2f}"
              f"  [{time.time()-t0:.0f} s]", flush=True)

    salida = {
        "carpeta": args.carpeta, "descripcion": desc, "corriente_A": corriente,
        "fps": cih["fps"], "escala_px_por_m": escala, "n_tif": len(archivos),
        "ventanas": list(VENTANAS), "solape": SOLAPE, "s2n_umbral": S2N_UMBRAL,
        "pares_por_n": args.pares, "filas": filas,
        "generado": time.strftime("%Y-%m-%d %H:%M:%S"),
        "script": os.path.basename(__file__),
    }
    destino = os.path.join(AQUI, "salidas", "tablas", f"barrido_dt_{args.carpeta}.json")
    with open(destino, "w", encoding="utf-8") as fh:
        json.dump(salida, fh, indent=1, ensure_ascii=False)
    print(f"\n-> {destino}")


if __name__ == "__main__":
    main()
