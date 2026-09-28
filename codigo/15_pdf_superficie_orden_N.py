#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compila informe/superficie_orden_N.pdf, el resumen corto de la superficie deformable de orden N
([D-51], [D-52], [H-23]–[H-26], [V-21], [V-22]), desde informe/superficie_orden_N.tex.

Los números salen de verificacion/logs/fase6_*.log y de salidas/tablas/fase6_*.json. Las
figuras se regeneran con numerico/fase6/estabilidad_cierre.py (f_fs6_estabilidad) y
numerico/fase6/escaleras.py figura (f_fs6_escaleras). El registro completo, para reproducir, es
numerico/fase6/REGISTRO_fase6.md.

Uso:
    /home/lucia/miniforge3/envs/piv-dt/bin/python codigo/15_pdf_superficie_orden_N.py

Dos pasadas de pdflatex, para las referencias cruzadas.
"""

import os
import subprocess
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(AQUI)
INFORME = os.path.join(RAIZ, "informe")
NOMBRE = "superficie_orden_N"


def main():
    for pasada in (1, 2):
        p = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
                            NOMBRE + ".tex"], cwd=INFORME, capture_output=True)
        if p.returncode != 0:
            # pdflatex escribe el log en latin-1 con los acentos del castellano
            with open(os.path.join(INFORME, NOMBRE + ".log"), encoding="latin-1") as f:
                lineas = f.read().splitlines()
            errores = [l for l in lineas if l.startswith("!")]
            print("pdflatex falló en la pasada %d:\n  " % pasada + "\n  ".join(errores[:5]))
            return 1
    with open(os.path.join(INFORME, NOMBRE + ".log"), encoding="latin-1") as f:
        log = f.read()
    paginas = log.rsplit("Output written on", 1)[-1].split("(")[1].split(" page")[0] \
        if "Output written on" in log else "?"
    print("escrito informe/%s.pdf (%s páginas); cajas desbordadas: %d"
          % (NOMBRE, paginas, log.count("Overfull")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
