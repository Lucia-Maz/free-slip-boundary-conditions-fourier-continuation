#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Arma el PDF del barrido en delta, etapa 1.

El veredicto no se escribe a mano: sale de un criterio fijado ACÁ, en el código, y
aplicado a los números del JSON. Así el informe no puede decir "no hay apartamiento"
porque a mí me parezca.

CRITERIO, declarado antes de mirar los resultados
=================================================
  - HAY apartamiento  si alguna corrida estable se aparta más que la barra del alfa
    medido en el experimento (7,3 %) Y ese apartamiento CONVERGE con la resolución
    (las dos mallas coinciden dentro de un tercio del apartamiento).
  - NO HAY apartamiento  si todas las corridas estables quedan dentro del 7,3 %.
  - NO CONCLUYENTE  si hay apartamiento pero no converge con la resolución: entonces es
    numérico, no físico, y hace falta una malla más fina.

Uso:
    /home/lucia/miniforge3/envs/piv-dt/bin/python codigo/11_pdf_barrido_delta.py
"""

import json
import math
import os
import subprocess
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(AQUI)
TABLAS = os.path.join(RAIZ, "salidas", "tablas")
FIGS = os.path.join(RAIZ, "salidas", "figuras")
DESTINO = os.path.join(RAIZ, "informe")

BARRA_EXPERIMENTAL = 0.073        # 0,00488 / 0,06687, de [V-13]
TOL_CONVERGENCIA = 1.0 / 3.0


def veredicto(d):
    """Aplica el criterio de arriba. Devuelve (clave, diccionario de evidencia)."""
    claves = sorted(d["corridas"])
    estables = {c: [f for f in d["corridas"][c][1:] if not f["crece"]] for c in claves}
    inestables = {c: [f for f in d["corridas"][c][1:] if f["crece"]] for c in claves}

    ap = {c: (max(abs(f["lambda_sobre_referencia"] - 1) for f in estables[c])
              if estables[c] else float("nan")) for c in claves}
    ap_max = max(ap.values())

    ev = {"apartamiento_por_malla": ap, "apartamiento_maximo": ap_max,
          "barra": BARRA_EXPERIMENTAL,
          "n_estables": {c: len(estables[c]) for c in claves},
          "n_inestables": {c: len(inestables[c]) for c in claves},
          "u0_inestables": {c: [f["u0"] for f in inestables[c]] for c in claves}}

    if ap_max <= BARRA_EXPERIMENTAL:
        return "no_hay", ev

    # ¿converge? se comparan las dos mallas en las amplitudes donde las dos son estables
    fina, gruesa = claves[-1], claves[0]
    pares = []
    for f in estables[fina]:
        g = next((x for x in estables[gruesa] if abs(x["u0"] - f["u0"]) < 1e-12), None)
        if g is not None:
            pares.append((f["u0"], f["lambda_sobre_referencia"],
                          g["lambda_sobre_referencia"]))
    ev["pares"] = pares
    if not pares:
        return "no_concluyente", ev
    peor = max(abs(a - b) for _, a, b in pares)
    ev["discrepancia_entre_mallas"] = peor
    if peor <= TOL_CONVERGENCIA * ap_max:
        return "hay", ev
    return "no_concluyente", ev


TEXTO = {
 "no_hay": r"""\textbf{No hay apartamiento medible.} Sobre todas las corridas estables, la
tasa de decaimiento coincide con la lineal dentro de %(ap).2f\,\%%, es decir por debajo de
la barra de error del $\alpha$ medido en el experimento (%(barra).1f\,\%%). En el rango de
$\delta$ cubierto, la clausura de un modo no se rompe de manera detectable con este
método.

Esto \emph{no} demuestra que el efecto sea nulo: demuestra que, si existe, es menor que la
incertidumbre con la que se lo compararía. Para el uso que tiene en este trabajo
—contrastar $\alpha$ calculado contra $\alpha$ medido— eso alcanza, y quiere decir que la
comparación de~[V-13] no está contaminada por no linealidad en este régimen. Lo que sigue
abierto en~[P-02] —que el $\alpha$ medido queda por debajo de $\lambda(k)$— no se explica
entonces por este candidato.""",
 "hay": r"""\textbf{Hay apartamiento, y converge con la resolución.} El máximo es
%(ap).2f\,\%%, por encima de la barra experimental (%(barra).1f\,\%%), y las dos mallas
coinciden dentro de %(conv).2f\,\%%, de modo que no es un artefacto de resolución. La
clausura de un modo se rompe en el rango de $\delta$ cubierto, y la comparación
de~[V-13] tiene que corregirse por este efecto o restringirse a la cola del registro.""",
 "no_concluyente": r"""\textbf{No concluyente.} Se observa un apartamiento de hasta
%(ap).2f\,\%%, por encima de la barra experimental, pero las dos mallas discrepan entre sí
en %(conv).2f\,\%%, que no es chico frente al efecto. Con estos datos no se puede separar
la no linealidad de la falta de resolución, que es exactamente lo que [H-09] advertía.
Hace falta una malla más fina antes de sacar conclusiones.""",
}


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--caso", default="forzado")
    args = ap.parse_args()
    with open(os.path.join(TABLAS, "barrido_delta_etapa1_%s.json" % args.caso)) as f:
        d = json.load(f)
    ver, ev = veredicto(d)

    claves = sorted(d["corridas"])
    filas = []
    for i, u0 in enumerate(d["u0s"]):
        celdas = []
        for c in claves:
            f = d["corridas"][c][i + 1]
            celdas.append((f["delta"], f["lambda_sobre_referencia"],
                           f.get("alfa_eff_sobre_alfa", float("nan")), f["crece"]))
        filas.append((u0, celdas))

    cuerpo = []
    for u0, celdas in filas:
        s = r"\num{%.2e}" % u0
        for delta, rat, aeff, crece in celdas:
            s += " & %.3g & %s & %s" % (
                delta, r"\emph{inestable}" if crece else "%.5f" % rat,
                "--" if crece else "%.4f" % aeff)
        cuerpo.append(s + r" \\")

    tex = r"""\documentclass[11pt,a4paper]{article}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage[spanish,es-noquoting]{babel}
\usepackage{amsmath}
\usepackage{booktabs}
\usepackage{graphicx}
\usepackage{siunitx}
\usepackage[margin=2.4cm]{geometry}
\usepackage[colorlinks=true,linkcolor=black,urlcolor=blue,citecolor=black]{hyperref}
\sisetup{output-exponent-marker=\text{e}}

\title{¿Se aparta la fricción de fondo calculada cuando el flujo deja de ser lento?\\
\large Barrido en $\delta$, etapa 1 --- caso \emph{%(caso)s}, $\epsilon = %(eps).2f$}
\author{Proyecto final --- capa delgada con superficie libre}
\date{%(fecha)s}

\begin{document}
\maketitle

\section{Qué está en juego}

El trabajo afirma que el coeficiente de fricción de fondo se puede \emph{calcular} en vez
de ajustarlo: $\alpha = \pi^2\nu/4h^2$ para fondo no deslizante y tope free-slip. Esa
afirmación descansa en una clausura de un solo modo vertical, y la pregunta abierta
[P-01] es hasta dónde vale. El parámetro que la gobierna es
\begin{equation}
  \delta = \frac{U k h^2}{\nu} = \mathrm{Re}_h\,\epsilon , \qquad \epsilon = k h ,
\end{equation}
que compara la advección con la difusión viscosa vertical. En la celda medida $\delta$ va
de $\sim 20$ al arranque del registro a $\sim 1$ al final, así que la pregunta no es
académica.

Dos cosas quedaron sin resolver analíticamente. La estimación de la distorsión del perfil
vertical lleva un factor geométrico $O(1)$ que no se calculó, y el orden del sesgo sobre
$\alpha$ depende de la estructura horizontal del campo: para un único modo celular la
realimentación se anula hasta $O(\delta^2)$ por selección de números de onda, pero un
campo de banda ancha no tiene esa protección. Las dos se cierran midiendo.

\section{El método, y por qué así}

Se integra la ecuación de Navier--Stokes incompresible con el esquema de partición de
Fontana \emph{et al.}~\cite{fontana2020}, con fondo no deslizante y tope free-slip,
y se mide la tasa de decaimiento $\lambda$ de $\langle v^2\rangle$ sobre una ventana
tardía, entre %(v0).2f y %(v1).2f memorias modales $h^2/(2\pi^2\nu)$ desde el arranque:
la pregunta [P-02] es sobre la ventana del ajuste experimental, a más de cinco memorias
del corte, con el perfil vertical en balance cuasi-estacionario con el término no lineal
([D-42]). Además, sobre el campo escrito en el medio de esa ventana se mide directamente
$\alpha_\text{eff} = \nu\,\langle \partial_z u_h|_0 \cdot \bar U_h\rangle / (h\,\langle|\bar U_h|^2\rangle)$,
la fricción de fondo \emph{sin} clausura, que vale exactamente $\alpha$ para el perfil
$\sin(\pi z/2h)$; y el $k$ horizontal pesado por energía, que es el que entra en $\delta$.

\paragraph{El caso.} $\epsilon = kh = %(eps).2f$ reproduce el régimen \emph{%(caso)s} del
experimento ([D-41]): con el paso real de la red de imanes, el estado forzado y el arranque
del decaimiento están en $\epsilon = 1{,}78$, y la ventana del ajuste de $\alpha$
($t > 10$~s), con la energía ya migrada a $k \approx 60$--$130$~m$^{-1}$, en
$\epsilon \approx 0{,}35$--$0{,}8$. El modo dominante de la condición inicial es el diagonal,
$|k| = \sqrt2$ en unidades de código, así que $L_z = \epsilon/\sqrt2$.

\paragraph{La referencia lineal se mide, no se escribe.} El punto de comparación sale de
correr \emph{el mismo binario, la misma malla y la misma condición inicial} a amplitud
$u_0 = \num{1e-4}$, dos órdenes por debajo de la corrida más chica del barrido, donde
el término no lineal es despreciable por construcción. Así la referencia incluye
automáticamente el $k$ discreto de la malla, el sesgo del integrador y el filtro de
\emph{dealiasing}, y el cociente $\lambda/\lambda_\text{lin}$ es un número puro:
cualquier apartamiento es no linealidad y nada más. Escribir la fórmula analítica habría
mezclado el efecto buscado con los errores de discretización.

\paragraph{Dos resoluciones horizontales, y no es opcional.} El hallazgo [H-09] de este
mismo proyecto mostró que con el tope free-slip el término no lineal se desestabiliza sobre
mallas donde el canal aguanta: a amplitud alta y malla chica la energía \emph{crece} sin
forzado, que es imposible. Un apartamiento que no converge al refinar es numérico, no
físico. Por eso todo el barrido se corre a $16\times16$ y a $32\times32$.

\paragraph{Parámetros.} Caja $2\pi\times2\pi\times%(lz).4f$, $\nu = \num{%(nu)s}$
(escalada con $L_z^2$ para que el tiempo de código sea el mismo múltiplo de $h^2/\nu$ en
los dos casos), $N_z = %(nz)d$ con $C_z = %(cz)d$ puntos de continuación, Runge--Kutta de
orden %(ord)d, $\Delta t = \num{%(dt)s}$ (un sexto del límite de estabilidad viscosa
explícita), %(pasos)d pasos hasta $t = %(tf).3f$, es decir %(tf_tmem).1f memorias
modales. Las amplitudes $u_0$ se eligieron para cubrir $\delta$ objetivo
%(delta_obj)s en la ventana; el $\delta$ de la tabla es el medido ahí.

\paragraph{El criterio, fijado antes de mirar los resultados.} Se declara que
\emph{hay} apartamiento si alguna corrida estable se aparta más que la barra de error del
$\alpha$ medido en el experimento, $%(barra).1f\,\%%$~[V-13], \emph{y} ese apartamiento
converge con la resolución (las dos mallas coinciden dentro de un tercio del efecto). Si
todas las corridas estables quedan dentro de esa barra, \emph{no hay}. Si hay apartamiento
pero no converge, el resultado es \emph{no concluyente}. El veredicto de la
sección~\ref{sec:res} lo produce el código, no el autor.

\section{Resultados}
\label{sec:res}

\begin{table}[h]
\centering
\begin{tabular}{l%(cols)s}
\toprule
& \multicolumn{%(ncol)d}{c}{$\delta$, $\lambda/\lambda_\text{lin}$ y $\alpha_\text{eff}/\alpha$ por malla}\\
\cmidrule(l){2-%(lastcol)d}
$u_0$ %(cab)s \\
\midrule
%(cuerpo)s
\bottomrule
\end{tabular}
\caption{Tasa de decaimiento medida, normalizada por la de la corrida lineal de
referencia de la misma malla. $\delta$ está a menos de una constante de forma $O(1)$;
lo que importa acá es el rango cubierto y la tendencia, no su valor absoluto.}
\end{table}

\begin{figure}[h]
\centering
\includegraphics[width=\linewidth]{%(fig)s}
\caption{Cociente entre la tasa medida y la lineal, contra $\delta$, para las dos mallas.
La banda gris es la barra de error del $\alpha$ medido en el experimento. Las cruces
marcan corridas en las que la energía creció, o sea numéricamente inestables, que no se
usan.}
\end{figure}

\section{Veredicto}

%(veredicto)s

\section{Lo que esto no contesta}

La etapa 1 mide la tasa de decaimiento de la energía \emph{total}, que mezcla tres cosas:
la fricción de fondo, la disipación viscosa horizontal y la transferencia no lineal entre
escalas. No aísla el corte en el fondo. La cantidad limpia es
\begin{equation}
  \alpha_\text{eff} \equiv \frac{\nu}{h}\,
  \frac{\partial u/\partial z\big|_{0}}{\langle u\rangle} ,
\end{equation}
que es exacta y no supone ninguna clausura, pero pedirle a SPECTER que la escriba requiere
agregar un diagnóstico en Fortran. Esa es la etapa 2, y sólo hace falta si esta etapa
muestra algo.

Otras salvedades: $\delta$ está a menos de una constante de forma $O(1)$, porque
$\langle v^2\rangle$ que informa el código está promediado sobre el dominio extendido por
la continuación FC--Gram y no sobre el físico; la condición inicial no cumple el no
deslizamiento del fondo, así que hay un transitorio que se descarta ajustando sobre el
último tercio; y la estructura horizontal es un único modo con dependencia en $x$, que es
justamente el caso donde la realimentación no lineal está más suprimida.

\begin{thebibliography}{9}
\bibitem{fontana2020} M.~Fontana, O.~P.~Bruno, P.~D.~Mininni y P.~Dmitruk,
\emph{Fourier continuation method for incompressible fluids with boundaries},
Comput. Phys. Commun. \textbf{256}, 107482 (2020).
\href{https://doi.org/10.1016/j.cpc.2020.107482}{doi:10.1016/j.cpc.2020.107482}.
Leído en el original en sus secciones 2, 4.1, 4.2 y 4.4.
\bibitem{orszag1986} S.~A.~Orszag, M.~Israeli y M.~O.~Deville,
\emph{Boundary conditions for incompressible flows},
J. Sci. Comput. \textbf{1}, 75 (1986).
\href{https://doi.org/10.1007/BF01061454}{doi:10.1007/BF01061454}.
Origen de las condiciones de contorno sobre $v^*$ y sobre $p$ que usa el solver.
Localizado y verificado en Crossref; \emph{no} leído en el original, es de acceso cerrado.
\bibitem{kio1991} G.~E.~Karniadakis, M.~Israeli y S.~A.~Orszag,
\emph{High-order splitting methods for the incompressible Navier--Stokes equations},
J. Comput. Phys. \textbf{97}, 414 (1991).
\href{https://doi.org/10.1016/0021-9991(91)90007-8}{doi:10.1016/0021-9991(91)90007-8}.
Su sección 2.3 se aplicó al esquema de este solver en~[V-15].
\end{thebibliography}

\end{document}
"""

    ncol = 3 * len(claves)
    cab = "".join(r" & \multicolumn{3}{c}{$%s\times%s$}" % (c[1:], c[1:])
                  for c in claves)
    sub = " ".join(r"& $\delta$ & $\lambda/\lambda_\text{lin}$ & $\alpha_\text{eff}/\alpha$"
                   for _ in claves)

    datos = {
        "lz": d["Lz"], "nu": "%.3e" % d["nu"], "nz": d["nz"], "cz": d["cz"],
        "ord": d["ord"], "dt": "%.0e" % d["dt"], "pasos": d["pasos"],
        "tf": d["t_final"], "barra": 100 * BARRA_EXPERIMENTAL,
        "caso": d.get("caso", args.caso), "eps": d.get("eps", float("nan")),
        "fecha": __import__("datetime").date.today().strftime("%d de %B de %Y"),
        "v0": d["ventana"][0] / d["t_mem"], "v1": d["ventana"][1] / d["t_mem"],
        "tf_tmem": d["t_final"] / d["t_mem"],
        "delta_obj": ", ".join("%g" % x for x in d.get("delta_objetivo", [])),
        "cols": "rrr" * len(claves), "ncol": ncol, "lastcol": ncol + 1,
        "cab": cab + r" \\ " + sub,
        "cuerpo": "\n".join(cuerpo),
        "fig": os.path.relpath(os.path.join(FIGS, "f7_barrido_delta_%s.pdf"
                                            % d.get("caso", args.caso)), DESTINO),
        "veredicto": TEXTO[ver] % {
            "ap": 100 * ev["apartamiento_maximo"],
            "barra": 100 * BARRA_EXPERIMENTAL,
            "conv": 100 * ev.get("discrepancia_entre_mallas", float("nan")),
        },
    }

    os.makedirs(DESTINO, exist_ok=True)
    ruta_tex = os.path.join(DESTINO, "barrido_delta_etapa1_%s.tex" % datos["caso"])
    with open(ruta_tex, "w") as f:
        f.write(tex % datos)

    for _ in range(2):
        p = subprocess.run(["pdflatex", "-interaction=nonstopmode",
                            "-halt-on-error", os.path.basename(ruta_tex)],
                           cwd=DESTINO, capture_output=True)
    pdf = ruta_tex.replace(".tex", ".pdf")
    if not os.path.exists(pdf):
        # pdflatex escribe el log en latin-1 con los acentos del castellano
        log = open(ruta_tex.replace(".tex", ".log"), "rb").read().decode(
            "utf-8", "replace")
        print(log[-3000:])
        raise RuntimeError("pdflatex falló")
    print("veredicto: %s" % ver)
    print(json.dumps(ev, indent=1, default=str))
    print("escrito %s (%.0f kB)" % (os.path.relpath(pdf, RAIZ),
                                    os.path.getsize(pdf) / 1e3))
    return 0


if __name__ == "__main__":
    sys.exit(main())
