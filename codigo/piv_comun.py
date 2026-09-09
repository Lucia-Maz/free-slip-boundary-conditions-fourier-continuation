"""Utilidades compartidas: leer la campaña del 02-06-25 sin modificarla.

Todo lo que sabe este módulo sobre el experimento sale de archivos que están en el
disco de Lucía, no de constantes escritas a mano:

  - la cabecera .cih de Photron da fps, tamaño de imagen y profundidad de bits;
  - los JSON de sesión de OpenPIV dan la escala px/m y la configuración de PIV;
  - los .vec dan la grilla realmente usada, que sirve para verificar lo anterior.

Ver DECISIONES.md [D-02]: la carpeta de origen se trata como sólo lectura.
"""
from __future__ import annotations

import glob
import json
import os
import re

import numpy as np

RAIZ = "/mnt/usb-Seagate_Expansion_HDD_00000000NACCVMS5-0:0-part1/Lucía Mazaira/02-06-25"

# Mapa carpeta -> medición, según la entrada del 2 de junio de 2025 del cuaderno de
# laboratorio. Los timestamps de los archivos siguen exactamente este orden, que es
# la verificación de que el mapa es correcto (ver [V-01]).
MEDICIONES = {
    "med_S0001": ("grid", "Grid — calibración espacial", None),
    "med_S0002": ("forzado", "Forzado 32 V, 1.9–2.0 A (CI 1)", 1.95),
    "med_S0003": ("decaimiento", "Decaimiento 1 desde CI 1", None),
    "med_S0004": ("forzado", "Forzado 16 V, 0.9–1.0 A (CI 2)", 0.95),
    "med_S0005": ("decaimiento", "Decaimiento 2 desde CI 2", None),
    "med_S0006": ("decaimiento", "Decaimiento 3 desde CI tipo 1", None),
    "med_S0007": ("quenching", "Quenching 1 (cambalache)", None),
    "med_S0008": ("forzado", "Forzado 2.1–2.2 A (CI 3)", 2.15),
    "med_S0009": ("decaimiento", "Decaimiento 4 desde CI 3", None),
    "med_S0010": ("quenching", "Quenching 2 (desde el reposo)", None),
    "med_S0011": ("quenching", "Quenching 3 (regrabación de med 7)", None),
}


def leer_cih(carpeta: str) -> dict:
    """Lee la cabecera de la cámara Photron. Devuelve fps, tamaño y bits."""
    cih = glob.glob(os.path.join(carpeta, "*.cih"))
    if not cih:
        raise FileNotFoundError(f"no hay .cih en {carpeta}")
    d = {}
    with open(cih[0], "rb") as fh:
        for linea in fh.read().decode("latin-1").splitlines():
            if ":" in linea:
                k, _, v = linea.partition(":")
                d[k.strip()] = v.strip()
    return {
        "archivo": os.path.basename(cih[0]),
        "fps": float(d["Record Rate(fps)"]),
        "shutter": d["Shutter Speed(s)"],
        "ancho": int(d["Image Width"]),
        "alto": int(d["Image Height"]),
        "bits": int(d["Color Bit"]),
        "n_cuadros": int(d["Total Frame"]),
        "fecha": d.get("Date", "?"),
    }


def leer_sesion(ruta: str) -> dict:
    """Lee un JSON de sesión de la GUI de OpenPIV."""
    with open(ruta, encoding="utf-8") as fh:
        return json.load(fh)


def escala_px_por_m(raiz: str = RAIZ) -> float:
    """La escala de calibración, tomada de las sesiones y verificada consistente.

    Se lee de todas las sesiones disponibles y se exige que coincidan; si no
    coincidieran habría que decidir cuál vale, y eso es una decisión científica.
    """
    valores = set()
    for ruta in glob.glob(os.path.join(raiz, "Sesion_*")) + glob.glob(
        os.path.join(raiz, "*", "Ses*ion_*")
    ):
        try:
            valores.add(float(leer_sesion(ruta)["scale"]))
        except Exception:
            continue
    if len(valores) != 1:
        raise ValueError(f"las sesiones no coinciden en la escala: {sorted(valores)}")
    return valores.pop()


def tifs(carpeta: str) -> list[str]:
    """Los TIF de una carpeta, ordenados por su número de cuadro."""
    ts = glob.glob(os.path.join(carpeta, "*.tif"))
    def numero(p):
        m = re.search(r"(\d{6})\.tif$", p)
        return int(m.group(1)) if m else -1
    return sorted(ts, key=numero)


def leer_tif(ruta: str) -> np.ndarray:
    """Un cuadro como float32.

    Se devuelve float32 y no float64 a propósito: son 4 MB por imagen en vez de 8,
    y la correlación no gana nada con la precisión extra (ver memoria de RAM).

    Se usa `tifffile` y no `PIL`: los TIFF de la Photron son de 10 bits guardados en
    uint16 y comprimidos, y PIL falla con `UnidentifiedImageError` sobre ellos. Ver
    [D-12].
    """
    import tifffile
    return tifffile.imread(ruta).astype(np.float32)


def leer_vec(ruta: str) -> dict:
    """Lee un .vec de OpenPIV. Columnas: x, y, vx, vy, val, sig2noise.

    `val` vale 1 en los vectores marcados inválidos/reemplazados y 0 en los que son
    medición. La convención se dedujo comparando fracciones entre procesamientos
    (0.1 % con ventana de 16 px contra 65 % con ventana de 8 px); la lectura opuesta
    daría 99.9 % de inválidos con la ventana buena, que es absurdo. Ver [V-02].
    """
    a = np.loadtxt(ruta)
    valido = a[:, 4] == 0
    return {
        "x": a[:, 0], "y": a[:, 1], "vx": a[:, 2], "vy": a[:, 3],
        "valido": valido, "s2n": a[:, 5],
        "frac_valida": float(valido.mean()),
        "n_x": len(np.unique(a[:, 0])), "n_y": len(np.unique(a[:, 1])),
    }


def urms(vec: dict) -> float:
    """Velocidad cuadrática media sobre los vectores válidos únicamente."""
    v = vec["valido"]
    return float(np.sqrt(np.mean(vec["vx"][v] ** 2 + vec["vy"][v] ** 2)))


def ventana_desde_grilla(n_vec: int, lado_img: int = 1024) -> list[tuple[int, float]]:
    """Qué combinaciones (ventana, solape) producen n_vec vectores por lado.

    Sirve para verificar la configuración registrada en un JSON contra la grilla que
    realmente tienen los .vec, que es el único registro que no se puede falsear.
    """
    salida = []
    for W in (8, 16, 32, 64, 128):
        for frac in (0.0, 0.25, 0.5, 0.75, 0.875):
            paso = int(round(W * (1 - frac)))
            if paso <= 0:
                continue
            if (lado_img - W) // paso + 1 == n_vec:
                salida.append((W, frac))
    return salida
