#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compila el informe de la entrega final del proyecto entero (≤ 5 páginas):
informe/FINAL/entrega_final.tex -> informe/FINAL/entrega_final.pdf.

La carpeta informe/FINAL/ compila sola (para editarla a mano o en Overleaf). Este script
primero copia a informe/FINAL/figuras/ las figuras actuales de salidas/figuras/:
  f7_barrido_delta_{ventana,forzado}.pdf   verificacion/barrido_delta_etapa1.py ([V-16])
  f1_decaimiento_claro.png                 codigo/07_figuras_decaimiento.py ([V-13])
  f_fs6_escaleras.png                      numerico/fase6/escaleras.py figura ([H-25], [H-26])
  f_fs6_superficie.png                     numerico/fase6/figura_campos_superficie.py
Los números del texto salen de DECISIONES.md y de salidas/tablas/, cada uno con su etiqueta.

Uso:
    /home/lucia/miniforge3/envs/piv-dt/bin/python codigo/16_pdf_entrega_breve.py

Dos pasadas de pdflatex, para las referencias cruzadas.
"""

import os
import shutil
import subprocess
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(AQUI)
FINAL = os.path.join(RAIZ, "informe", "FINAL")
NOMBRE = "entrega_final"
FIGURAS = ("f7_barrido_delta_ventana.pdf", "f7_barrido_delta_forzado.pdf",
           "f1_decaimiento_claro.png", "f_fs6_escaleras.png", "f_fs6_superficie.png")


def main():
    os.makedirs(os.path.join(FINAL, "figuras"), exist_ok=True)
    for f in FIGURAS:
        shutil.copy2(os.path.join(RAIZ, "salidas", "figuras", f), os.path.join(FINAL, "figuras", f))
    for pasada in (1, 2):
        p = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
                            NOMBRE + ".tex"], cwd=FINAL, capture_output=True)
        if p.returncode != 0:
            # pdflatex escribe el log en latin-1 con los acentos del castellano
            with open(os.path.join(FINAL, NOMBRE + ".log"), encoding="latin-1") as f:
                lineas = f.read().splitlines()
            errores = [l for l in lineas if l.startswith("!")]
            print("pdflatex falló en la pasada %d:\n  " % pasada + "\n  ".join(errores[:5]))
            return 1
    with open(os.path.join(FINAL, NOMBRE + ".log"), encoding="latin-1") as f:
        log = f.read()
    paginas = log.rsplit("Output written on", 1)[-1].split("(")[1].split(" page")[0] \
        if "Output written on" in log else "?"
    for ext in (".aux", ".log", ".out"):
        try:
            os.remove(os.path.join(FINAL, NOMBRE + ext))
        except OSError:
            pass
    print("escrito informe/FINAL/%s.pdf (%s páginas); cajas desbordadas: %d"
          % (NOMBRE, paginas, log.count("Overfull")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
