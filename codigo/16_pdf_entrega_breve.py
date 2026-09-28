#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compila informe/entrega_breve.pdf, el informe de la entrega final del proyecto entero
(≤ 5 páginas), desde informe/entrega_breve.tex.

Los números salen de DECISIONES.md y de salidas/tablas/ (cada uno con su etiqueta en el
texto). Las figuras: f7 (verificacion/barrido_delta_etapa1.py), f1 (codigo/07), y
f_fs6_superficie (numerico/fase6/campos_superficie.py y figura_campos_superficie.py).

Uso:
    /home/lucia/miniforge3/envs/piv-dt/bin/python codigo/16_pdf_entrega_breve.py

Dos pasadas de pdflatex, para las referencias cruzadas.
"""

import os
import subprocess
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(AQUI)
INFORME = os.path.join(RAIZ, "informe")
NOMBRE = "entrega_breve"


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
