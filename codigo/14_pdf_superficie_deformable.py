#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compila informe/superficie_deformable.pdf, el informe de la superficie libre deformable
([D-49], [V-20], [H-20], [H-21]), desde informe/superficie_deformable.tex y sus partes
informe/superficie_deformable_*.tex.

Los números salen de los logs de verificacion/logs/fase5_*.log y fase2_regresion_*.log, y de
salidas/tablas/superficie_deformable_{corridas,escalas}.json. Las figuras se regeneran con
numerico/fase5/figuras_superficie_deformable.py y superficie_deformable_escalas.py.

Uso:
    /home/lucia/miniforge3/envs/piv-dt/bin/python codigo/14_pdf_superficie_deformable.py

Dos pasadas de pdflatex, para el índice y las referencias cruzadas.
"""

import os
import subprocess
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(AQUI)
INFORME = os.path.join(RAIZ, "informe")


def main():
    for pasada in (1, 2):
        p = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
                            "superficie_deformable.tex"], cwd=INFORME, capture_output=True)
        if p.returncode != 0:
            # pdflatex escribe el log en latin-1 con los acentos del castellano
            with open(os.path.join(INFORME, "superficie_deformable.log"), encoding="latin-1") as f:
                lineas = f.read().splitlines()
            errores = [l for l in lineas if l.startswith("!")]
            print("pdflatex falló en la pasada %d:\n  " % pasada + "\n  ".join(errores[:5]))
            return 1
    with open(os.path.join(INFORME, "superficie_deformable.log"), encoding="latin-1") as f:
        log = f.read()
    paginas = log.rsplit("Output written on", 1)[-1].split("(")[1].split(" page")[0] \
        if "Output written on" in log else "?"
    print("escrito informe/superficie_deformable.pdf (%s páginas); cajas desbordadas: %d"
          % (paginas, log.count("Overfull")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
