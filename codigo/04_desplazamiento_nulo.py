"""Mide σ directamente, sobre registros de fluido en reposo con partículas.

Es el experimento de Foucaut, Carlier & Stanislas 2004 (DOI 10.1088/0957-0233/15/6/003,
sección 2): si el desplazamiento verdadero es nulo, todo lo que mide PIV es su propio
ruido. Ellos lo hacían en un acuario con agua quieta; acá hay dos registros que sirven
y que ya estaban grabados.

  `fondo_202409_1643`  campaña del 05-09-24, 3072 cuadros. El nombre y la estadística
                       de imagen dicen que es la celda sembrada en reposo: 3432 picos
                       locales y asimetría 4,86, iguales a los del forzado del mismo día,
                       pero con la mitad de diferencia entre cuadros consecutivos.

  `med_S0010`          campaña del 02-06-25, Quenching 2. El cuaderno dice "Guardo desde
                       antes de encenderlo (frames 500 a 3072)", así que el **principio**
                       del registro es la celda en reposo, con la óptica y el sembrado
                       exactos de todo lo demás que analizamos.

Sobre la separación en dos componentes: Foucaut informa σ por componente (0,087 px en x
y 0,050 px en y para ventana 16×16, su tabla 1), no el módulo. Acá se informan las dos
cosas, porque compararlas mal es un factor √2 fácil de cometer.

Uso
---
    python 04_desplazamiento_nulo.py
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import piv_comun as pc  # noqa: E402
from importlib import import_module  # noqa: E402

bd = import_module("02_barrido_dt")

AQUI = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAIZ_DISCO = os.path.dirname(pc.RAIZ)


def sigma_de(archivos: list[str], indices: list[int], n: int = 1) -> dict:
    """Estadística del desplazamiento medido sobre pares con separación n.

    El desplazamiento **medio** es deriva residual: el fluido nunca está perfectamente
    quieto, y una deriva uniforme no es ruido de correlación. σ se define como la
    dispersión alrededor de esa media, no el valor cuadrático medio crudo. Se informan
    los dos para que se vea cuánto pesa la deriva.
    """
    su, sv, mu, mv, val, s2n = [], [], [], [], [], []
    for i in indices:
        a = pc.leer_tif(archivos[i])
        b = pc.leer_tif(archivos[i + n])
        c = bd.piv_multipaso(a, b)
        ok = (c["s2n"] >= bd.S2N_UMBRAL) & np.isfinite(c["u"]) & np.isfinite(c["v"])
        u, v = c["u"][ok], c["v"][ok]
        mu.append(float(u.mean())); mv.append(float(v.mean()))
        su.append(float(u.std())); sv.append(float(v.std()))
        val.append(float(ok.mean())); s2n.append(float(np.median(c["s2n"][ok])))
        del a, b, c
    su, sv = np.array(su), np.array(sv)
    return {
        "n_pares": len(indices), "n_separacion": n,
        "deriva_media_px": [float(np.mean(mu)), float(np.mean(mv))],
        "sigma_u_px": float(su.mean()), "sigma_v_px": float(sv.mean()),
        "sigma_u_sd": float(su.std()), "sigma_v_sd": float(sv.std()),
        "sigma_modulo_px": float(np.sqrt(su.mean() ** 2 + sv.mean() ** 2)),
        "frac_valida": float(np.mean(val)), "s2n_mediano": float(np.mean(s2n)),
    }


def main() -> None:
    resultados = {}

    # ---- 1) el fondo de la campaña 05-09-24 ---------------------------------
    fondo = sorted(glob.glob(os.path.join(
        RAIZ_DISCO, "05-09-24", "fondo_202409_1643", "*.tif")))
    if fondo:
        idx = np.linspace(0, len(fondo) - 2, 12).astype(int).tolist()
        resultados["fondo_05-09-24"] = sigma_de(fondo, idx)
        resultados["fondo_05-09-24"]["n_cuadros"] = len(fondo)

    # ---- 2) el tramo en reposo de med_S0010 ---------------------------------
    # Primero hay que encontrar dónde se enciende la corriente: se recorre el principio
    # del registro y se busca el salto en el desplazamiento medido.
    s10 = pc.tifs(os.path.join(pc.RAIZ, "med_S0010"))
    perfil = []
    for i in range(0, 600, 40):
        a, b = pc.leer_tif(s10[i]), pc.leer_tif(s10[i + 1])
        c = bd.piv_multipaso(a, b)
        ok = (c["s2n"] >= bd.S2N_UMBRAL) & np.isfinite(c["u"])
        perfil.append((i, float(np.sqrt(np.mean(c["u"][ok] ** 2 + c["v"][ok] ** 2)))))
        del a, b, c
    resultados["perfil_encendido_med_S0010"] = perfil
    base = np.median([d for _, d in perfil[:5]])
    encendido = next((i for i, d in perfil if d > 2 * base), None)
    resultados["cuadro_encendido_estimado"] = encendido
    tope = (encendido if encendido else 400) - 40
    if tope > 40:
        idx = np.linspace(0, tope - 2, 10).astype(int).tolist()
        resultados["reposo_med_S0010"] = sigma_de(s10, idx)
        resultados["reposo_med_S0010"]["tramo_usado"] = [0, tope]

    destino = os.path.join(AQUI, "salidas", "tablas", "desplazamiento_nulo.json")
    with open(destino, "w", encoding="utf-8") as fh:
        json.dump(resultados, fh, indent=1, ensure_ascii=False)

    print("perfil de encendido de med_S0010 (cuadro, desplazamiento rms en px):")
    print("  " + "  ".join(f"{i}:{d:.2f}" for i, d in perfil))
    print(f"  encendido estimado en el cuadro {encendido}\n")
    for k, r in resultados.items():
        if not isinstance(r, dict) or "sigma_u_px" not in r:
            continue
        print(f"{k}")
        print(f"   deriva media     = ({r['deriva_media_px'][0]:+.3f}, "
              f"{r['deriva_media_px'][1]:+.3f}) px")
        print(f"   sigma_u          = {r['sigma_u_px']:.4f} ± {r['sigma_u_sd']:.4f} px")
        print(f"   sigma_v          = {r['sigma_v_px']:.4f} ± {r['sigma_v_sd']:.4f} px")
        print(f"   sigma modulo     = {r['sigma_modulo_px']:.4f} px")
        print(f"   validos {r['frac_valida']:.1%}   s2n {r['s2n_mediano']:.2f}   "
              f"{r['n_pares']} pares\n")
    print("Referencia: Foucaut 2004 tabla 1, ventana 16x16: sigma_u = 0.087, sigma_v = 0.050 px")
    print(f"\n-> {destino}")


if __name__ == "__main__":
    main()
