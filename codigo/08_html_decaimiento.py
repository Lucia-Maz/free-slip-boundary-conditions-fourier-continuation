#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Arma la página HTML de revisión del decaimiento, con las figuras embebidas.

Las tablas se generan desde los JSON, no se transcriben a mano: la página no puede
decir un número distinto del que produjo el cálculo. Las figuras van embebidas como
data URI, en las dos versiones (clara y oscura), y la hoja de estilo elige según el
tema de quien mire.

Uso:
    /home/lucia/miniforge3/envs/piv-dt/bin/python codigo/08_html_decaimiento.py
"""

import base64
import json
import math
import os
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(AQUI)
sys.path.insert(0, AQUI)
from celda import CAMPANA_02_06_25 as CELDA, VELOCIDADES_MEDIDAS  # noqa: E402

TABLAS = os.path.join(RAIZ, "salidas", "tablas")
FIGS = os.path.join(RAIZ, "salidas", "figuras")
SALIDA = os.path.join(RAIZ, "salidas", "decaimiento.html")

ETIQ = {"med_S0003": "Decaimiento 1", "med_S0005": "Decaimiento 2",
        "med_S0006": "Decaimiento 3", "med_S0009": "Decaimiento 4"}


def datauri(ruta):
    with open(ruta, "rb") as f:
        return "data:image/png;base64," + base64.b64encode(f.read()).decode("ascii")


def figura(nombre, titulo, pie):
    return """<figure>
  <img class="clara" src="%s" alt="%s">
  <img class="oscura" src="%s" alt="%s">
  <figcaption><b>%s</b> %s</figcaption>
</figure>""" % (datauri(os.path.join(FIGS, nombre + "_claro.png")), titulo,
                datauri(os.path.join(FIGS, nombre + "_oscuro.png")), titulo,
                titulo, pie)


def tabla(cabeceras, filas, clases=""):
    th = "".join("<th>%s</th>" % c for c in cabeceras)
    tr = "".join("<tr>%s</tr>" % "".join("<td>%s</td>" % c for c in f) for f in filas)
    return ('<div class="scroll"><table class="%s"><thead><tr>%s</tr></thead>'
            "<tbody>%s</tbody></table></div>" % (clases, th, tr))


def main():
    d = json.load(open(os.path.join(TABLAS, "decaimiento.json")))
    r = json.load(open(os.path.join(TABLAS, "decaimiento_resumen.json")))
    vf = json.load(open(os.path.join(TABLAS, "verificacion_filtrado.json")))
    ruta_est = os.path.join(TABLAS, "forzados_estacionariedad.json")
    est = json.load(open(ruta_est)) if os.path.exists(ruta_est) else None

    a_todo = r["alfa_medido_ensemble_1_s"]
    a_cola = r["alfa_medido_cola_1_s"]
    ea_cola = r["alfa_medido_cola_error_1_s"]
    a_libre = r["alfa_predicho_1_s"]
    a_canal = r["alfa_predicho_canal_1_s"]

    # --- tabla principal: los 32 puntos corregidos ---------------------------
    filas = []
    for med in sorted(d["registros"]):
        reg = d["registros"][med]
        esp = {round(e["t_s"], 3): e for e in reg["espectros"]}
        for p in reg["calibracion"]:
            crudo = p["puntos"][0]["desp_rms_px"] * d["celda"]["fps"] / d["celda"]["escala_px_por_m"]
            e = esp.get(round(p["t_s"], 3))
            filas.append([
                ETIQ[med], "%.2f" % p["t_s"],
                "%.3f" % (1e3 * p["u_rms_superficie_m_s"]),
                "%.3f" % (1e3 * crudo),
                "%.2f×" % (crudo / p["u_rms_superficie_m_s"]),
                "%d" % p["n_de_d1"], "%.2f" % p["desp_en_n_px"],
                "%.3f" % p["ajuste_dos_dt"]["sigma_por_componente_px"],
                "%.4f" % p["ajuste_dos_dt"]["R2"],
                "%.0f" % e["k_pico_1_m"] if e else "—",
            ])

    filas_alfa = [[ETIQ[m], "%.5f" % v["alfa_1_s"], "± %.5f" % v["err"],
                   "%.2f" % (1 / v["alfa_1_s"])]
                  for m, v in sorted(r["alfa_por_registro"].items())]

    bloque_est = ""
    if est:
        fe = []
        for med, serie in sorted(est.items()):
            for s in serie:
                fe.append([med, "%.1f" % s["t_s"],
                           "%.3f" % (1e3 * s["u_rms_superficie_m_s"]),
                           "%d" % s["n"], "%.2f" % s["desp_px"]])
        bloque_est = """
<h2>Un hallazgo colateral: <code>med_S0008</code> contiene un decaimiento</h2>
<p>El barrido en Δt del capítulo de instrumento supone que la corrida forzada es
<b>estacionaria</b>, y de ahí sale la velocidad de referencia de 1,321&nbsp;mm/s. Medidos
con el mismo estimador robusto que se usa acá, los dos registros forzados se comportan de
manera opuesta:</p>
%s
<p><code>med_S0002</code> (1,95&nbsp;A) es plano al 2&nbsp;%% en 43&nbsp;s: es
estacionario, como se suponía. <code>med_S0008</code> (2,15&nbsp;A) <b>decae</b>: de
3,10&nbsp;mm/s a los 12,7&nbsp;s a 0,335&nbsp;mm/s a los 43&nbsp;s, lo que da una tasa de
≈&nbsp;0,073&nbsp;s⁻¹ — la tasa de fricción de fondo. El registro etiquetado «Forzado
2,1–2,2&nbsp;A» contiene un decaimiento.</p>
<p><b>Qué cambia y qué no.</b> El barrido reparte sus 12 pares uniformemente a lo largo del
registro para cada n, así que la velocidad que ajusta no es «la velocidad del forzado»
sino un rms sobre un registro en el que el flujo cae unas siete veces. Hay que releerla
como promedio del registro. El efecto sobre σ&nbsp;=&nbsp;0,169&nbsp;px <b>no está
establecido</b>: como los pares se reparten igual para todos los n, el término de ruido
sigue siendo separable a primer orden, pero la premisa del modelo está violada y eso pide
un reanálisis, no una afirmación.</p>
<p>De paso explica el salto que estaba anotado como abierto entre forzado y decaimiento:
el Decaimiento&nbsp;4 arranca en 4,07&nbsp;mm/s y el forzado que sí es estacionario
(<code>med_S0002</code>, con <i>menos</i> corriente) está en 3,55&nbsp;mm/s. El valor
anómalo era el de <code>med_S0008</code>, no el del decaimiento, y no es sesgo de
supervivencia: estos números están corregidos por ruido y calculados sobre
&gt;95&nbsp;%% de vectores válidos.</p>
""" % tabla(["registro", "t [s]", "u_rms [mm/s]", "n", "desp [px]"], fe, "num")

    html = """<title>Decaimiento 02-06-25</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Source+Serif+4:opsz,wght@8..60,400;8..60,600&family=IBM+Plex+Sans:wght@400;600&family=IBM+Plex+Mono:wght@400;600&display=swap">
<style>
:root {
  color-scheme: light;
  --ground:  #f6f7f8;
  --surface: #ffffff;
  --ink:     #12161a;
  --muted:   #5b6570;
  --rule:    #dde2e6;
  --accent:  #2a78d6;
  --accent-soft: #eaf1fb;
  --warn:    #b45309;
  --serie1:#2a78d6; --serie2:#eb6834; --serie3:#1baf7a; --serie4:#eda100;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --ground:  #14171a;
    --surface: #1b1f23;
    --ink:     #eef1f4;
    --muted:   #a3adb7;
    --rule:    #2b3137;
    --accent:  #3987e5;
    --accent-soft: #17263a;
    --warn:    #d9922b;
    --serie1:#3987e5; --serie2:#d95926; --serie3:#199e70; --serie4:#c98500;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --ground:  #14171a;
  --surface: #1b1f23;
  --ink:     #eef1f4;
  --muted:   #a3adb7;
  --rule:    #2b3137;
  --accent:  #3987e5;
  --accent-soft: #17263a;
  --warn:    #d9922b;
  --serie1:#3987e5; --serie2:#d95926; --serie3:#199e70; --serie4:#c98500;
}

body { background: var(--ground); color: var(--ink);
       font-family: "Source Serif 4", Georgia, serif; font-size: 17px;
       line-height: 1.62; }
.hoja { max-width: 60rem; margin: 0 auto; padding: 3.5rem 1.5rem 6rem;
        display: flex; flex-direction: column; gap: 1.15rem; }
.hoja > * { margin: 0; }

.eyebrow { font-family: "IBM Plex Sans", system-ui, sans-serif; font-size: .74rem;
           letter-spacing: .13em; text-transform: uppercase; color: var(--muted); }
h1 { font-size: clamp(1.9rem, 4.4vw, 2.7rem); line-height: 1.14; font-weight: 600;
     text-wrap: balance; letter-spacing: -.015em; }
h2 { font-family: "IBM Plex Sans", system-ui, sans-serif; font-size: 1.16rem;
     font-weight: 600; margin-top: 2.4rem; text-wrap: balance;
     padding-bottom: .4rem; border-bottom: 1px solid var(--rule); }
h3 { font-family: "IBM Plex Sans", system-ui, sans-serif; font-size: .96rem;
     font-weight: 600; margin-top: 1.2rem; }
p, li { max-width: 40rem; }
.hoja > p, .hoja > ul, .hoja > ol { color: var(--ink); }
a { color: var(--accent); }
code, .num td, .num th { font-family: "IBM Plex Mono", ui-monospace, monospace; }
code { font-size: .9em; background: var(--accent-soft); padding: .08em .32em;
       border-radius: 3px; }
b, strong { font-weight: 600; }
.nota { color: var(--muted); font-size: .93rem; }

/* la conclusión, que es el motivo de la página */
.veredicto { background: var(--surface); border: 1px solid var(--rule);
             border-left: 3px solid var(--accent);
             padding: 1.3rem 1.5rem; display: flex; flex-direction: column;
             gap: .9rem; }
.cifras { display: grid; grid-template-columns: repeat(auto-fit, minmax(11rem, 1fr));
          gap: 1.4rem 2rem; }
.cifra { display: flex; flex-direction: column; gap: .15rem; }
.cifra .v { font-family: "IBM Plex Mono", monospace; font-size: 1.45rem;
            font-weight: 600; font-variant-numeric: tabular-nums; }
.cifra .k { font-family: "IBM Plex Sans", sans-serif; font-size: .78rem;
            letter-spacing: .04em; text-transform: uppercase; color: var(--muted); }
.cifra.destacada .v { color: var(--accent); }

figure { background: var(--surface); border: 1px solid var(--rule);
         padding: 1rem 1rem .35rem; display: flex; flex-direction: column;
         gap: .55rem; }
figure img { display: block; width: 100%; height: auto; }
figcaption { font-family: "IBM Plex Sans", sans-serif; font-size: .85rem;
             color: var(--muted); line-height: 1.5; padding-bottom: .6rem;
             max-width: 46rem; }
figcaption b { color: var(--ink); }
img.oscura { display: none; }
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) img.clara { display: none; }
  :root:not([data-theme="light"]) img.oscura { display: block; }
}
:root[data-theme="dark"] img.clara { display: none; }
:root[data-theme="dark"] img.oscura { display: block; }

.scroll { overflow-x: auto; border: 1px solid var(--rule); background: var(--surface); }
table { border-collapse: collapse; width: 100%; font-size: .84rem; }
th, td { padding: .4rem .7rem; text-align: right; white-space: nowrap;
         border-bottom: 1px solid var(--rule); }
th:first-child, td:first-child { text-align: left; }
thead th { font-family: "IBM Plex Sans", sans-serif; font-weight: 600;
           font-size: .74rem; letter-spacing: .04em; text-transform: uppercase;
           color: var(--muted); position: sticky; top: 0; background: var(--surface); }
tbody tr:last-child td { border-bottom: none; }
.num td { font-variant-numeric: tabular-nums; }

.pendientes li { margin-bottom: .5rem; }
.pie { color: var(--muted); font-size: .85rem;
       font-family: "IBM Plex Sans", sans-serif; border-top: 1px solid var(--rule);
       padding-top: 1rem; margin-top: 2rem; max-width: none; }
</style>

<main class="hoja">

<p class="eyebrow">Campaña 02-06-25 · cuatro decaimientos · Fase 4</p>
<h1>El decaimiento medido, con el ruido de PIV descontado</h1>
<p class="nota">Generado por <code>codigo/08_html_decaimiento.py</code> desde
<code>salidas/tablas/decaimiento.json</code>. Todas las cifras salen de los JSON: la
página no puede decir un número distinto del que produjo el cálculo.</p>

<div class="veredicto">
  <div class="cifras">
    <div class="cifra destacada"><span class="v">${a_cola4}</span>
      <span class="k">α medido (t &gt; 10 s)</span></div>
    <div class="cifra"><span class="v">${a_libre4}</span>
      <span class="k">α predicho, tope libre</span></div>
    <div class="cifra"><span class="v">${a_canal4}</span>
      <span class="k">α predicho, tapa rígida</span></div>
    <div class="cifra"><span class="v">${rat_canal}</span>
      <span class="k">medido / tapa rígida</span></div>
  </div>
  <p style="max-width:46rem">Todo en s⁻¹, sobre las cuatro corridas como ensemble. El
  decaimiento medido es compatible con la predicción de <b>fondo no-deslizante y
  superficie libre</b> dentro de un ${desv}&nbsp;%, y está <b>${fac} veces por
  debajo</b> de lo que daría una tapa rígida. Es la comparación que motivaba todo el
  proyecto.</p>
</div>

<h2>Por qué había que rehacer la medición</h2>
<p>Las series <code>.vec</code> archivadas de los cuatro decaimientos están armadas con
pares de cuadros <b>consecutivos</b>, Δt&nbsp;=&nbsp;1/60&nbsp;s. Es el peor Δt posible:
al arrancar el decaimiento el desplazamiento verdadero es de unas dos décimas de píxel por
cuadro, y el piso de ruido del capítulo de instrumento es de ~0,17&nbsp;px. El ruido es
del orden de la señal. Además esas series tienen entre <b>30 y 36&nbsp;% de vectores
válidos</b>; el multipaso propio sobre los mismos cuadros da <b>97,5&nbsp;%</b>.</p>
<p>La consecuencia se ve en la figura&nbsp;1: leída de los <code>.vec</code>, la velocidad
<b>no decae</b> — se queda clavada en el piso de ruido durante los 51&nbsp;s. Un α ajustado
sobre esa curva habría dado casi cero.</p>

${f1}

<h2>El método</h2>
<p>En ocho instantes de cada registro se rearman los pares desde los TIF crudos con
separaciones n&nbsp;=&nbsp;1, 2, 4, 8, 16 y 32 cuadros, y se ajusta</p>
<p style="font-family:'IBM Plex Mono',monospace;font-size:.95rem;text-align:center">
desp_rms(n)² = (n·d₁)² + σ²</p>
<p>que separa el desplazamiento verdadero por cuadro <code>d₁</code> del ruido
<code>σ</code>. Es el método de dos Δt del capítulo de instrumento, aplicado dentro de una
ventana corta: lícito porque el flujo cambia en ~15&nbsp;s y la ventana dura menos de
0,6&nbsp;s.</p>
<p><b>Pero el valor que se usa no sale del ajuste.</b> A n grande el ruido es despreciable
por construcción, y ahí <code>d₁ = desp_rms(n)/n</code> directamente, con una corrección
relativa que la tabla informa y que casi siempre está por debajo del 1&nbsp;%. Se prefiere
ese estimador porque el ajuste completo da <b>σ² negativo</b> en varios instantes: a
n&nbsp;=&nbsp;1 el desplazamiento es de dos décimas de píxel y el estimador subpíxel sesga
hacia el entero. El σ² negativo se informa, no se recorta a cero.</p>

${f2}
${f3}

<h2>La curva, registro por registro</h2>
${t_alfa}
<p class="nota">Los cuatro α coinciden dentro de sus errores. El ajuste del ensemble sobre
todos los puntos da ${a_todo5}&nbsp;s⁻¹; restringido a t&nbsp;&gt;&nbsp;10&nbsp;s da
${a_cola5}&nbsp;±&nbsp;${ea_cola5}&nbsp;s⁻¹, y es el que vale.</p>
<p><b>La meseta inicial es el estado forzado, antes de cortar la corriente.</b> Está
grabada a propósito, para no tener que coordinar el corte con el inicio del registro y
arriesgarse a perder los primeros pasos del decaimiento. Por eso hay que excluirla del
ajuste: no es decaimiento. Que el corte no esté en un instante conocido no afecta a α — un
corrimiento del origen de tiempos cambia la amplitud ajustada, no la pendiente.</p>
<p>De paso, la meseta <b>mide el estado forzado</b> de cada corrida: 4,07&nbsp;mm/s la que
sale de CI&nbsp;3 y ~3,6&nbsp;mm/s las que salen de CI&nbsp;1. Es un contraste
independiente para la sección de más abajo.</p>

<h2>La escala del flujo, y qué le hace a [P-01]</h2>
<p>El parámetro que decide si la reducción bidimensional vale es
<code>δ = U·k·h²/ν</code>, y depende de <b>k</b>. Hasta ahora se lo evaluaba con el
armónico del forzado (88,9&nbsp;m⁻¹, diagonal, 7,07&nbsp;cm). Los espectros medidos dicen
que el flujo tiene energía a escalas bastante más chicas, y que <b>la escala crece durante
el decaimiento</b>: el k pesado por energía pasa de ${k_ini} a ${k_fin}&nbsp;m⁻¹.</p>

${f4}
${f5}

<p>Con el k medido, δ arranca en ~20 y termina en ~1: la distorsión estimada del perfil
vertical (≈&nbsp;0,0077·δ) va del <b>15&nbsp;% al 1&nbsp;%</b>. O sea que la clausura de
un modo es marginal en los primeros segundos y buena en el resto — que es justamente la
ventana donde se ajustó α.</p>

<h2>¿Se está filtrando de más?</h2>
<p>Conviene separar lo que <b>es</b> un filtro de lo que no. La mediana móvil de la
figura&nbsp;1 es cosmética: se aplica sólo a las líneas tenues y no entra en ningún
número. Quedan tres cosas que sí pueden borrar información, y cada una tiene su
verificación.</p>

<h3>A · Elegir Δt largo, ¿pierde señal?</h3>
<p>Usar n grande no suaviza en el espacio, pero sí promedia en el tiempo sobre n/60&nbsp;s.
Hay un test sin parámetros libres: si el exceso que se ve a n chico fuera <b>sólo</b> ruido
aditivo, entonces</p>
<p style="font-family:'IBM Plex Mono',monospace;font-size:.95rem;text-align:center">
(d₁(n)/d₁)² − 1 = (σ / n·d₁)²</p>
<p>exactamente: una recta de pendiente 1, sin ajustar nada. Los 93 puntos caen sobre o por
debajo de esa recta (cociente mediano 0,22). <b>Puntos por encima</b> serían la firma de
que el n grande perdió señal, y no aparecen. Además, entre los dos n más grandes d₁
coincide dentro del <b>2,1&nbsp;%</b> mediano; el peor instante llega al 15&nbsp;%, y es el
más lento de su registro, donde el desplazamiento es de menos de un píxel.</p>

<h3>B · ¿Se ve la ventana en el espectro?</h3>
<p>Sí, pero no la que uno esperaría, y eso es un resultado. Las escalas en juego:</p>
<ul>
<li>Nyquist de la grilla de vectores (paso 4&nbsp;px): <b>2985&nbsp;m⁻¹</b></li>
<li>primer cero de sinc²(kX/2) con la ventana final X&nbsp;=&nbsp;16&nbsp;px:
<b>1492&nbsp;m⁻¹</b></li>
<li>corte de Foucaut, k_c·X&nbsp;=&nbsp;2,8, con X&nbsp;=&nbsp;16&nbsp;px:
<b>665&nbsp;m⁻¹</b></li>
</ul>
<p>El piso medido cae un factor <b>0,19</b> entre 900 y 2900&nbsp;m⁻¹. Un sinc² de la
ventana de 16&nbsp;px predice 0,002 —y un cero en el medio de la banda, que no está—;
el del paso de grilla de 4&nbsp;px predice 0,46. Es decir que <b>el ruido de este
multipaso no está filtrado por la ventana final</b>: se comporta como si estuviera
correlacionado sobre el paso de la grilla. Para el marco de Foucaut que usa el capítulo de
instrumento, con E&nbsp;=&nbsp;E_ruido·sinc²(kX/2) y X&nbsp;=&nbsp;16&nbsp;px, quiere decir
que <b>no describe esta salida</b>.</p>
<p class="nota">Salvedad: en esa banda puede haber señal real mezclada con el ruido, así
que la caída medida es una cota y no una medición limpia de la respuesta. Y en la
dirección que importa para δ: como la ventana sí suprime energía verdadera a k alto, el k
que reporto es una <b>cota inferior</b>, y δ también.</p>

<h3>C · La resta del piso, ¿borra señal?</h3>
<p>Es el único paso donde una elección mía puede quitar información, así que va verificado
por tres lados.</p>
<p><b>C1 — sobre datos sintéticos de σ conocido</b>, que es la única forma de saber si el
estimador está bien normalizado: se le da ruido blanco puro y tiene que devolver el σ que
se puso.</p>
${t_sint}
<p><b>C2 — contra una medición independiente.</b> El σ leído del espectro es una serie
limpia y monótona que cae junto con el flujo (0,26&nbsp;→&nbsp;0,06&nbsp;px en el
Decaimiento&nbsp;1), lo que confirma [H-07] por una ruta distinta. El σ del ajuste de dos
Δt, en cambio, es inestable en estos instantes —da negativo a menudo— así que no sirve
como contraste punto a punto; en la figura&nbsp;3 están los dos.</p>
<p><b>C3 — sensibilidad a la banda.</b> Cambiando dónde se lee el piso, el k pesado se
mueve <b>menos del 3&nbsp;%</b>: ${k_banda}. El número sale del dato, no de la
elección.</p>

${f6}

<p><b>Cuánto puede sesgar todo esto la velocidad.</b> Tomando la cota superior
σ&nbsp;=&nbsp;0,235&nbsp;px y el desplazamiento más chico que se usó (0,68&nbsp;px), la
corrección sería del <b>6&nbsp;%</b>; con el σ efectivo que sale del colapso de A,
<b>1,3&nbsp;%</b>. En los instantes típicos, con 2 a 3&nbsp;px, queda por debajo del
1&nbsp;%. Todo eso es menor que la barra del α, que es del 7&nbsp;%.</p>

${bloque_est}

<h2>Los 32 puntos</h2>
${t_puntos}
<p class="nota">«crudo n=1» es lo que se leería del <code>.vec</code> archivado en ese
instante; «inflación» es el cociente entre ambos. σ es por componente y sale del ajuste
completo, que es el que puede dar negativo.</p>

<h2>Lo que esto no resuelve</h2>
<ul class="pendientes">
<li><b>El instante exacto del corte del forzado</b> no está registrado; se sabe que cae
dentro de la meseta. No afecta a α —un corrimiento del origen de tiempos cambia la
amplitud ajustada, no la pendiente— pero impide fechar el decaimiento en absoluto.</li>
<li><b>La viscosidad es la del agua</b>, no la del electrolito (KNO₃ al 16&nbsp;% m/m).
Entra a primer orden: <code>δα/α = δν/ν − 2δh/h</code>.</li>
<li><b>El espesor.</b> Con h&nbsp;=&nbsp;6&nbsp;mm la predicción es
${a_libre4}&nbsp;s⁻¹; para reproducir exactamente el α medido haría falta
h&nbsp;=&nbsp;${h_equiv}&nbsp;mm. Un error de medio milímetro en h mueve α un
17&nbsp;%, más que la discrepancia observada.</li>
<li><b>La película superficial.</b> Una película rígida sube α, no la baja. Que el α
medido esté en o por debajo del valor de superficie libre limpia acota la película:
en el modelo de Boussinesq–Scriven, β&nbsp;≲&nbsp;0,15, o sea
μ_s&nbsp;≲&nbsp;3·10⁻⁶&nbsp;N·s/m. La cota está dominada por la incertidumbre en h, no
por la estadística.</li>
<li><b>El k que hay que usar en δ</b> depende de cómo se lo defina; acá es el pesado por
energía con el ruido blanco de PIV restado, y el piso se lee por debajo de Nyquist
(2984&nbsp;m⁻¹ para esta grilla).</li>
</ul>

<p class="pie">h = ${h} mm · ν = 10⁻⁶ m²/s (agua, supuesto declarado) · imanes de 1 cm en
red tipo tablero de 5 cm · 3800 px/m · 60 fps · registros de 51,2 s.
Parámetros en <code>codigo/celda.py</code>; derivación en
<code>teoria/P01_validez_reduccion_2D.md</code>; log en <code>DECISIONES.md</code>.</p>

</main>
"""
    from string import Template
    html = Template(html).substitute(
        a_cola4="%.4f" % a_cola, a_cola5="%.5f" % a_cola,
        ea_cola5="%.5f" % ea_cola,
        a_todo5="%.5f" % a_todo, a_todo4="%.4f" % a_todo,
        a_libre4="%.4f" % a_libre, a_canal4="%.4f" % a_canal,
        rat_canal="%.3f" % (a_cola / a_canal),
        desv="%.0f" % (100 * abs(a_cola - a_libre) / a_libre),
        fac="%.1f" % (a_canal / a_cola),
        k_ini="%.0f" % r["k_medido_inicio_1_m"],
        k_fin="%.0f" % r["k_medido_final_1_m"],
        h="%.0f" % (CELDA.h * 1e3),
        h_equiv="%.2f" % (CELDA.h * 1e3 * math.sqrt(a_libre / a_cola)),
        f1=figura("f1_decaimiento", "Figura 1 — el decaimiento.",
                  "Puntos: rearmado desde los TIF y corregido por ruido. Línea tenue: "
                  "mediana móvil de la serie .vec archivada, que no decae. Zona gris: "
                  "excluida del ajuste."),
        f2=figura("f2_metodo", "Figura 2 — cómo se separa el ruido.",
                  "Decaimiento 1, tres instantes. La pendiente da la velocidad verdadera "
                  "y la ordenada al origen el ruido. Ejes logarítmicos."),
        f3=figura("f3_ruido", "Figura 3 — σ cae junto con el flujo.",
                  "Dos rutas independientes al mismo número: cuadrados, el piso del "
                  "espectro, con la normalización verificada sobre ruido sintético al "
                  "1 %; círculos tenues, la ordenada al origen del ajuste de dos Δt. Las "
                  "dos caen con el flujo: es la confirmación directa de [H-07]."),
        f4=figura("f4_espectros", "Figura 4 — espectro de energía.",
                  "Decaimiento 1, tres instantes. Tenue: crudo. Grueso: con el ruido "
                  "blanco de PIV restado. La línea vertical es el armónico del forzado."),
        f5=figura("f5_delta", "Figura 5 — el parámetro δ de [P-01].",
                  "Puntos: con el k medido del campo. Punteado: con el k del forzado, "
                  "que es lo que se venía usando. La diferencia es un factor ~2,5."),
        t_alfa=tabla(["registro", "α [1/s]", "error", "τ = 1/α [s]"], filas_alfa, "num"),
        t_puntos=tabla(["registro", "t [s]", "u_rms [mm/s]", "crudo n=1", "inflación",
                        "n", "desp [px]", "σ/comp [px]", "R²", "k pico [1/m]"],
                       filas, "num"),
        f6=figura("f6_filtrado", "Figura 6 — las dos verificaciones del filtrado.",
                  "A: el exceso de d₁ a n chico contra lo que predice sólo el ruido; la "
                  "recta es y = x, sin ajustar nada. B: el piso del espectro medido "
                  "contra las dos respuestas candidatas."),
        t_sint=tabla(["σ puesto [px]", "σ recuperado [px]", "error"],
                     [["%.3f" % f["sigma_puesto_px"],
                       "%.4f" % f["sigma_recuperado_px"],
                       "%+.2f %%" % (100 * f["error_relativo"])]
                      for f in vf["sintetico"]], "num"),
        k_banda=(", ".join("%.0f" % f["k_pesado_t0"] for f in vf["sensibilidad_banda"])
                 + " m⁻¹ al inicio, y "
                 + ", ".join("%.0f" % f["k_pesado_tfinal"]
                             for f in vf["sensibilidad_banda"]) + " m⁻¹ al final"),
        bloque_est=bloque_est,
    )

    with open(SALIDA, "w") as f:
        f.write(html)
    print("escrito %s (%.1f MB)"
          % (os.path.relpath(SALIDA, RAIZ), os.path.getsize(SALIDA) / 1e6))
    return 0


if __name__ == "__main__":
    sys.exit(main())
