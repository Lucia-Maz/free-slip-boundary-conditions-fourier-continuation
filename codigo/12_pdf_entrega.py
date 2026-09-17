#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compila el PDF de la entrega, informe/entrega.pdf, desde informe/entrega.tex.

El .tex es prosa escrita a mano: la versión explicada de salidas/entrega1_resumen.html,
con las derivaciones. Los números que contiene son los de salidas/tablas/*.json y de
DECISIONES.md, y cada uno cita su procedencia en el propio texto (apéndice B). Las figuras
son las de salidas/figuras/, que se regeneran con codigo/07_*.py y codigo/10_*.py.

Uso:
    /home/lucia/miniforge3/envs/piv-dt/bin/python codigo/12_pdf_entrega.py

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
                            "entrega.tex"], cwd=INFORME, capture_output=True)
        if p.returncode != 0:
            # pdflatex escribe el log en latin-1 con los acentos del castellano
            with open(os.path.join(INFORME, "entrega.log"), encoding="latin-1") as f:
                lineas = f.read().splitlines()
            errores = [l for l in lineas if l.startswith("!")]
            print("pdflatex falló en la pasada %d:\n  " % pasada + "\n  ".join(errores[:5]))
            return 1
    with open(os.path.join(INFORME, "entrega.log"), encoding="latin-1") as f:
        log = f.read()
    paginas = log.rsplit("Output written on", 1)[-1].split("(")[1].split(" page")[0] \
        if "Output written on" in log else "?"
    print("escrito informe/entrega.pdf (%s páginas); cajas desbordadas: %d"
          % (paginas, log.count("Overfull")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
