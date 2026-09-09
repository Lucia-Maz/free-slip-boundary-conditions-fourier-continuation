"""Inventario verificable de la campaña del 02-06-25.

No interpreta nada: lee lo que hay, lo cruza contra el cuaderno de laboratorio y
contra sí mismo, y reporta qué verificaciones pasan y cuáles no. Escribe
`salidas/tablas/inventario.json`.

Las verificaciones que corre:

  [V-01] el orden temporal de los archivos coincide con el orden del cuaderno
  [V-02] la ventana de PIV registrada en cada sesión coincide con la grilla de los .vec
  [V-03] el Δt de cada sesión coincide con la separación real entre los cuadros usados
  [V-04] la cantidad de TIF coincide con la cabecera .cih de la cámara
  [V-05] el prefijo de los .vec coincide con el de los TIF de la misma carpeta
"""
from __future__ import annotations

import glob
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import piv_comun as pc  # noqa: E402

AQUI = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def numeros_de_cuadro(carpeta: str, cuantos: int = 12) -> list[int]:
    ns = []
    for p in pc.tifs(carpeta)[:cuantos]:
        m = re.search(r"(\d{6})\.tif$", p)
        if m:
            ns.append(int(m.group(1)))
    return ns


def separacion_real(nums: list[int], secuencia: str) -> float | None:
    """Cuántos cuadros hay entre las dos imágenes de un par, según la numeración.

    `(1+2),(3+4)` toma pares disjuntos: la separación es nums[1]-nums[0].
    `(1+2),(2+3)` toma pares deslizantes: la separación es la misma, pero los campos
    consecutivos comparten un cuadro y no son independientes.
    """
    if len(nums) < 2:
        return None
    return float(nums[1] - nums[0])


def main() -> None:
    escala = pc.escala_px_por_m()
    inventario, chequeos = [], []

    # --- carpetas crudas -----------------------------------------------------
    for nombre, (tipo, desc, corriente) in pc.MEDICIONES.items():
        carpeta = os.path.join(pc.RAIZ, nombre)
        if not os.path.isdir(carpeta):
            continue
        cih = pc.leer_cih(carpeta)
        ts = pc.tifs(carpeta)
        entrada = {
            "carpeta": nombre, "tipo": tipo, "descripcion": desc,
            "corriente_A": corriente, "n_tif": len(ts),
            "cih_n_cuadros": cih["n_cuadros"], "fps": cih["fps"],
            "shutter": cih["shutter"], "bits": cih["bits"],
            "tamano": f"{cih['ancho']}x{cih['alto']}",
            "prefijo_tif": os.path.basename(ts[0])[:40] if ts else None,
            "primeros_cuadros": numeros_de_cuadro(carpeta),
            "n_vec": len(glob.glob(os.path.join(carpeta, "*.vec"))),
        }
        inventario.append(entrada)
        # [V-04]
        ok = len(ts) == cih["n_cuadros"]
        nota = ""
        if not ok and nombre == "med_S0010":
            # el cuaderno dice "Guardo desde antes de encenderlo (frames 500 a 3072)"
            ok = len(ts) == cih["n_cuadros"] - 500 + 1
            nota = "guardado parcial desde el cuadro 500, según el cuaderno"
        chequeos.append({"id": "V-04", "objeto": nombre,
                         "pregunta": "¿n de TIF coincide con la cabecera .cih?",
                         "esperado": cih["n_cuadros"], "obtenido": len(ts),
                         "pasa": bool(ok), "nota": nota})

    # [V-01] orden temporal contra orden del cuaderno
    marcas = []
    for e in inventario:
        m = re.search(r"_(\d{8})_(\d{6})_", e["prefijo_tif"] or "")
        marcas.append(int(m.group(1) + m.group(2)) if m else None)
    ordenado = all(a is not None and b is not None and a < b
                   for a, b in zip(marcas, marcas[1:]))
    chequeos.append({"id": "V-01", "objeto": "campaña completa",
                     "pregunta": "¿los timestamps siguen el orden del cuaderno?",
                     "esperado": "creciente", "obtenido": "creciente" if ordenado else "NO",
                     "pasa": bool(ordenado), "nota": "mapea carpeta -> medición"})

    # --- carpetas de salidas -------------------------------------------------
    for carpeta in sorted(glob.glob(os.path.join(pc.RAIZ, "*_destino*"))):
        nombre = os.path.basename(carpeta)
        vs = sorted(glob.glob(os.path.join(carpeta, "*.vec")))
        ts = pc.tifs(carpeta)
        if not vs:
            continue
        v = pc.leer_vec(vs[len(vs) // 2])
        candidatos = pc.ventana_desde_grilla(v["n_x"])
        pref_vec = os.path.basename(vs[0]).rsplit("_", 1)[0]
        pref_tif = os.path.basename(ts[0]).split("_2025")[0] if ts else None
        entrada = {
            "carpeta": nombre, "n_vec": len(vs), "n_tif": len(ts),
            "grilla": f"{v['n_x']}x{v['n_y']}",
            "ventana_solape_compatibles": [f"{W} px @ {f:.0%}" for W, f in candidatos],
            "frac_valida": round(v["frac_valida"], 4),
            "prefijo_vec": pref_vec, "prefijo_tif": pref_tif,
            "primeros_cuadros": numeros_de_cuadro(carpeta),
        }
        # sesión guardada dentro de la carpeta, si la hay
        ses = [p for p in glob.glob(os.path.join(carpeta, "Ses*ion_*"))
               if os.path.isfile(p)]
        if ses:
            s = pc.leer_sesion(ses[0])
            entrada["sesion"] = os.path.basename(ses[0])
            entrada["sesion_dt"] = s.get("dt")
            entrada["sesion_secuencia"] = s.get("sequence")
            entrada["sesion_ventana_final"] = (
                s.get("corr_window_3") if s.get("custom_windowing") else s.get("corr_window"))
            # [V-03] Δt registrado contra separación real de cuadros
            sep = separacion_real(entrada["primeros_cuadros"], s.get("sequence", ""))
            fps = 60.0
            dt_esperado = sep / fps if sep else None
            pasa = dt_esperado is not None and abs(s["dt"] - dt_esperado) < 0.002
            chequeos.append({"id": "V-03", "objeto": nombre,
                             "pregunta": "¿el dt de la sesión coincide con la separación real de cuadros?",
                             "esperado": dt_esperado, "obtenido": s.get("dt"),
                             "pasa": bool(pasa), "nota": f"secuencia {s.get('sequence')}"})
            # [V-02] ventana registrada contra grilla real
            W = entrada["sesion_ventana_final"]
            pasa2 = any(W == c[0] for c in candidatos)
            chequeos.append({"id": "V-02", "objeto": nombre,
                             "pregunta": "¿la ventana de la sesión es compatible con la grilla de los .vec?",
                             "esperado": entrada["ventana_solape_compatibles"],
                             "obtenido": f"{W} px", "pasa": bool(pasa2), "nota": ""})
        # [V-05] etiquetado
        coincide = bool(pref_tif and pref_vec and
                        pref_tif.lower().replace("_", "")[:6] in pref_vec.lower().replace("_", ""))
        chequeos.append({"id": "V-05", "objeto": nombre,
                         "pregunta": "¿el prefijo de los .vec coincide con el de los TIF?",
                         "esperado": pref_tif, "obtenido": pref_vec,
                         "pasa": coincide,
                         "nota": "" if coincide else "ETIQUETADO INCONSISTENTE"})
        inventario.append(entrada)

    salida = {"raiz": pc.RAIZ, "escala_px_por_m": escala,
              "mm_por_px": 1e3 / escala, "inventario": inventario,
              "chequeos": chequeos}
    destino = os.path.join(AQUI, "salidas", "tablas", "inventario.json")
    with open(destino, "w", encoding="utf-8") as fh:
        json.dump(salida, fh, indent=1, ensure_ascii=False)

    fallan = [c for c in chequeos if not c["pasa"]]
    print(f"{len(inventario)} carpetas inventariadas, {len(chequeos)} verificaciones, "
          f"{len(fallan)} fallan\n")
    for c in chequeos:
        print(f"  [{c['id']}] {'PASA ' if c['pasa'] else 'FALLA'} {c['objeto']:22s} {c['pregunta']}")
        if not c["pasa"]:
            print(f"           esperado={c['esperado']}  obtenido={c['obtenido']}  {c['nota']}")
    print(f"\n-> {destino}")


if __name__ == "__main__":
    main()
