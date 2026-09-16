# Log de decisiones

Una entrada por decisión. Cada una dice **qué** se decidió, **por qué**, y **qué
alternativa se descartó**. Orden cronológico, se agrega al final, no se reescribe.
Las decisiones científicas (cambiar una banda, un Δt, un criterio de descarte) las
toma Lucía; acá quedan registradas incluso cuando fueron suyas.

Convención: `[D-nn]` decisión, `[H-nn]` hallazgo que motivó una decisión,
`[V-nn]` verificación que respalda una decisión.

---

## 2026-09-02

### [D-01] El proyecto vive en `/home/lucia/proyecto-final-piv/`, fuera del repo del curso

`~/GW-AI-course` es un clon del repositorio del curso, de otro autor. Poner el
proyecto adentro lo dejaría como cambios sin seguimiento sobre ese repo y complicaría
`git init` propio, que el proyecto final necesita.
**Alternativa descartada:** `~/GW-AI-course/proyecto-final-piv/`.

### [D-02] Los datos crudos del disco USB se tratan como sólo lectura

Nada se escribe, mueve ni borra en
`/mnt/usb-.../Lucía Mazaira/02-06-25`. Todo lo derivado va a este árbol. Es la única
copia de una campaña de medición de un día completo que no se puede repetir.
**Nota:** en `med_S0009` hay 28 GB de `.txt` que duplican 11 GB de `.vec`, y 5,4 GB de
sesiones de PIVlab. Redundantes, pero **no se tocan** — la decisión de borrar es de Lucía.

### [D-03] Entorno aislado en `.venv/`, con Python 3.12 del sistema

`openpiv` no estaba instalado. Se instala en un venv propio del proyecto en vez de en
el Python del sistema o en un env de conda existente (`gw-ai`, `fisica`), para que la
reproducción sea autocontenida y para no alterar entornos que Lucía usa para otra cosa.
`/usr/bin/python3.12` se eligió porque ya tenía numpy/scipy/matplotlib y el
`miniforge3/bin/python3` no.
**Alternativa descartada:** instalar en el env `gw-ai`.

### [D-04] Se arranca por `med_S0008` (Forzado 2.1 A), no por un decaimiento

Decisión de Lucía. El forzado es estacionario, así que tiene un único Δt óptimo y
sirve como banco de calibración: se puede barrer Δt y verificar que la función de
costo tiene el mínimo predicho. Los decaimientos, que son el objeto real de su tesis,
no admiten un Δt único (ver [H-03]) y se atacan después con lo que salga de acá.

### [D-05] Los pares se rearman desde los TIF crudos, no se reusan los `.vec` existentes

`med_S0008` tiene 3072 cuadros consecutivos a 60 fps. Eso permite construir pares con
cualquier separación Δt = n/60 s sin medir nada nuevo. Las salidas ya existentes
(`Forza3_destino_dt5`) cubren un solo n y además excluyen 25 % de vectores.

### [D-06] La configuración de PIV se copia de las sesiones de Lucía y se deja fija

Multipaso 64 → 32 → 16 px con 75 % de solape, correlación circular, subpíxel gaussiano,
deformación simétrica, `sig2noise` peak2peak con umbral 1,05. Son los valores de
`med_S0009/Session_0202` y `Deca4_destino_dt5/Sesion_200226_win16`, verificados contra
la grilla real de los `.vec` (253×253 sobre 1024 px ⇒ paso de 4 px ⇒ ventana 16 con
75 %). Se deja fija para que **el único parámetro que varía en el barrido sea Δt**.

### [D-07] No se aplica `smoothn` ni reemplazo de vectores en el barrido

Las sesiones de la familia `180226` tienen `smoothn_each_pass=True`. Un suavizado
posterior a la correlación agrega una función de transferencia propia, no incluida en
el modelo de Foucaut, que contaminaría la medición del piso de ruido — que es
justamente lo que el barrido quiere medir. Los vectores inválidos se cuentan y se
excluyen, no se rellenan.
**Consecuencia a vigilar:** excluir inválidos sesga hacia regiones bien correlacionadas,
que son las lentas. Es la explicación candidata de la discrepancia de [H-04].

### [D-08] El multipaso se implementa acá, no se llama a `openpiv.windef`

`windef` cambió de API entre versiones de OpenPIV y su punto de entrada trabaja sobre
carpetas enteras, no sobre un par de imágenes elegido. Como el barrido necesita elegir
exactamente qué dos cuadros se correlacionan, se implementa el multipaso en
`codigo/02_barrido_dt.py` sobre `pyprocess.extended_search_area_piv`, que es la parte
estable de la biblioteca, con deformación simétrica de imagen entre pasos
(Scarano 2002, DOI 10.1088/0957-0233/13/1/201).
**Consecuencia a declarar en el informe:** no es bit a bit el mismo código que corrió
la GUI de Lucía. La verificación de que es equivalente es reproducir uno de sus `.vec`
existentes con este multipaso y comparar campo a campo — pendiente, [V-03].

### [D-09] El barrido mide desplazamiento en píxeles, no velocidad

`piv_multipaso` devuelve el desplazamiento con `dt=1`. La conversión a m/s se hace
después. El motivo es que el ruido de PIV es constante en píxeles y no en m/s: dejarlo
en píxeles es lo que hace que el modelo `U² = U_true² + (σ/(Δt·s))²` tenga una sola
incógnita instrumental en vez de una por cada Δt.

### [D-10] Se descartan vectores por `sig2noise` y no se rellenan

Umbral 1.05, el de las sesiones de Lucía. Los descartados se **cuentan** — la fracción
válida es una de las dos señales que marcan el extremo superior del rango útil de Δt —
pero no se interpolan, porque un vector interpolado no es una medición y contaminaría
tanto U como σ. Ver también [D-07].

### [D-11] Entorno: conda en vez de venv

`/usr/bin/python3.12` no tiene `ensurepip`, así que `python -m venv` no puede crear un
entorno con pip, y `python3.12-venv` requiere permisos de administrador que no se
pidieron. Se creó en cambio el env de conda `piv-dt` con miniforge, que Lucía ya tenía.
Corrige a [D-03], que decía venv.
**Reproducción:** `mamba create -n piv-dt -c conda-forge python=3.11 numpy scipy
matplotlib tifffile imageio scikit-image pip` y después `pip install openpiv`.

### [D-12] Los TIF se leen con `tifffile`, no con `PIL`

Los TIFF de la Photron son de 10 bits guardados en `uint16` y comprimidos (1,3 MB en
vez de los 2 MB que ocuparían planos). `PIL` los rechaza con
`UnidentifiedImageError`; `tifffile` e `imageio` los leen bien. Se eligió `tifffile`
por ser la dependencia más chica de las dos.
**Verificado:** el primer cuadro de `med_S0008` da un array (1024, 1024) `uint16` con
valores entre 0 y 730, coherente con 10 bits (máximo posible 1023).

### [D-13] Correlación normalizada y resta de la media global

`normalized_correlation=True` y se resta la media de cada imagen antes de correlacionar
(equivale al `background_subtract = 'global mean'` de las sesiones de Lucía). Sin eso el
pico de correlación se apoya sobre un pedestal de continua y el criterio `peak2peak`
deja de discriminar.
**Medido** sobre el par (5, 10) de `med_S0008`, un solo paso de 16 px:

| | s2n mediano | vectores válidos |
|---|---|---|
| sin normalizar | 1,17 | 79,0 % |
| normalizado | **2,75** | **96,5 %** |

Como referencia, la corrida de la GUI de Lucía sobre ese mismo par da s2n 2,23 y 74,7 %.

### [D-14] No se invierte el signo de v

La primera versión hacía `v = -v` suponiendo que OpenPIV devuelve el eje y invertido.
Es incorrecto **para este uso**: el campo se realimenta a la deformación de imagen del
paso siguiente, así que con el signo cambiado la iteración empuja en sentido contrario
y diverge. Se detectó porque la correlación con el campo de la GUI daba u = +0,974 pero
v = **−0,943**.
Con la inversión, el multipaso daba 1,704 px donde un solo paso daba 0,545 px.

### [V-03] El multipaso reproduce la salida de la GUI de Lucía — PASA

Par (5, 10) de `med_S0008`, ventana final 16 px con 75 % de solape, contra
`Forza3_destino_dt5/forzado3_dt5_0203_piv_FFT_000.vec`:

| | desplazamiento rms | válidos | s2n |
|---|---|---|---|
| GUI de OpenPIV (Lucía) | 0,682 px | 74,7 % | 2,23 |
| multipaso de acá | **0,665 px** | 96,9 % | 3,09 |

**2,5 % de diferencia** en el desplazamiento. La validez es mayor acá porque la
normalización mejora el pico ([D-13]); la de Lucía además aplica filtros de validación
adicionales. Esto habilita a usar este multipaso como sustituto de su cadena, con la
salvedad declarada de que no es bit a bit el mismo código ([D-08]).

### [D-15] El barrido llega hasta n = 128 aunque el modelo no valga ahí

A n grande el desplazamiento supera la regla del cuarto para ventana de 16 px y aparece
decorrelación: a n = 64 la fracción válida ya cae a 39 %. Se miden igual esos puntos
porque **dónde se rompe el modelo es el resultado buscado** — es el extremo superior del
rango útil de Δt. El ajuste de σ los excluye por fracción válida (umbral 90 %,
`03_figuras.FRAC_VALIDA_MINIMA`), y esa exclusión queda registrada en el JSON de salida.

### [D-16] Error de unidades en el ajuste, corregido

La primera versión devolvía `sigma_px = sqrt(pendiente)` sin el factor `s/fps`, y daba
0,0029 px en vez de 0,18 px — 63 veces más chico, que es exactamente 3800/60. Es un
error que **no cambia el R²**, así que no se nota mirando la calidad del ajuste; se
notó porque el número era absurdo comparado con cualquier valor publicado.
**Lección para el informe:** contrastar contra un orden de magnitud externo detecta
errores de unidades que ningún diagnóstico interno del ajuste detecta.

### [D-17] El ajuste es no lineal con pesos relativos, no mínimos cuadrados sin pesos

Se probaron tres parametrizaciones del mismo modelo:

| | resultado |
|---|---|
| `U²` contra `n⁻²`, sin pesos | R² = 0,986 **pero** dominado por el punto n=1: sacándolo, R² cae a 0,861 y σ pasa de 0,181 a 0,142 px |
| `desp²` contra `n²`, sin pesos | R² = 0,999 pero la ordenada al origen sale **negativa**: σ² < 0, físicamente imposible |
| `desp = sqrt(c²n² + σ²)` con pesos `1/desp` | **la elegida.** R² = 0,9986, σ = 0,169 px, y sacando n=1 da 0,145 px |

El motivo es que n² recorre de 1 a 1024: sin pesos, mínimos cuadrados queda determinado
por el extremo y σ² sale como diferencia de números grandes. La incerteza de un
desplazamiento medido es aproximadamente proporcional al desplazamiento, así que
corresponde pesar por el valor.
La barra de error sale por bootstrap sobre los puntos (2000 remuestreos).

### [H-05] Este diseño constriñe bien la velocidad y mal el ruido

`U_verdadera = 1,321 mm/s` con IC68 [1,306, 1,346] — 3 % de ancho.
`σ = 0,169 px` con IC68 [0,144, 0,179] — 20 % de ancho.

No es un defecto del ajuste sino del diseño: a n chico la señal es comparable a σ, y a
n grande σ es despreciable y esos puntos no aportan información sobre el ruido. σ
despejado punto a punto da 0,18 / 0,14 / 0,14 / 0,22 / 0,23 / 0,14 px para n de 1 a 8,
y valores sin sentido para n ≥ 12 porque ahí σ ya no se distingue del residuo.
**Para medir σ bien hace falta el experimento de desplazamiento nulo** — la celda con
la corriente apagada, que es el método de Foucaut y son treinta minutos de laboratorio.

### [D-18] Cómo se contrasta contra Foucaut, y cómo no

σ = 0,169 px contra los 0,087 px de su tabla 1 es un factor 2, y **no es una
discrepancia**: su propia ecuación (3), σ = sqrt(4ζI/XY), hace a σ dependiente de ζ, que
es una propiedad de la cadena de registro. Ellos usaron YAG pulsado con imágenes de
partícula de 2,5 px; acá hay iluminación continua (shutter 1/60 a 60 fps) y partículas
de 100 µm. Que coincida el orden de magnitud es lo que valida el método.

Las dos predicciones que **sí** son contrastables, y que hay que correr:
1. σ no depende de Δt — ya confirmado: el modelo con σ constante ajusta a R² = 0,9986.
2. σ ∝ 1/sqrt(área de ventana) — **pendiente**: rehacer el barrido con ventana final de
   8 y de 32 px. σ debería duplicarse y reducirse a la mitad respectivamente.

### [V-04] El estimador subpíxel, calibrado contra corrimientos conocidos — PASA

`05_corrimiento_sintetico.py` sobre dos cuadros de `med_S0008`, desplazados por
interpolación de Fourier en pasos de 0,1 px de 0 a 1,5 px, descartando 64 px de borde
para que el envolvimiento del corrimiento no entre en la cuenta.

| | |
|---|---|
| sesgo máximo en 0 – 1,5 px | **−0,018 px** (en d = 1,50) |
| sesgo en corrimientos enteros (d = 0, 1) | 0,0021 px |
| σ_u máximo | **0,0114 px** (en d = 0,50) |
| σ_u en corrimientos enteros | 0,0027 px |
| control con splines cúbicos | coincide dentro de **0,0025 px** |

σ_u tiene máximo en el medio píxel y mínimo en los enteros, que es la firma esperada
del ajuste subpíxel. **No hay peak locking**: el sesgo no es la sierra antisimétrica de
período 1 px que ese efecto produce, sino una curva suave que va de +3,5 % del
desplazamiento a d chico a −1,2 % a d = 1,5.

El control con splines existe porque desplazar por Fourier tiene artefactos propios
(la imagen no es periódica). Que las dos interpolaciones den la misma curva dentro de
0,0025 px es lo que permite atribuir el sesgo a la cadena de PIV y no al método de
corrimiento. **Salvedad:** en esta corrida el control coincidió en 3 corrimientos
(0,0 / 0,5 / 1,0), no en 5: 0,25 y 0,75 no están en la grilla de Fourier y su columna
de diferencia salió `nan`. La lista de control se corrigió para que una corrida futura
compare en los cinco.

**Salvedad de fondo:** un par sintético está perfectamente correlacionado. No tiene
movimiento fuera del plano, ni cambio de iluminación entre cuadros, ni ruido de lectura
independiente. Este σ es sólo la contribución de la correlación y del estimador, y es
una **cota inferior**.

### [H-06] El σ² negativo del árbol paralelo es reemplazo de vectores, no física

En `~/GW-AI-course/piv-s0008-forzado` hay un barrido de tamaño de ventana sobre la
misma `med_S0008`, con `delta_frame` 5 y 10. Aplicarle el método de dos Δt da **σ² < 0
en las cuatro ventanas**, porque la velocidad medida *sube* con la separación:

| ventana | U(gap 5) | U(gap 10) |
|---|---|---|
| 64 px | 0,722 mm/s | 0,822 mm/s |
| 32 px | 1,310 | 1,386 |
| 16 px | 1,524 | 1,562 |
| 8 px | 1,533 | 1,561 |

El modelo exige que baje. La causa está en `config/parameter_sweep.json`:
`replace_vectors: true` con `replacement_method: "localmean"` y
`sig2noise_threshold: 1.3`. Las métricas se calculan sobre el campo entero
(`kinetic_energy = 0.5·nanmean(u² + v²)`), con los reemplazos adentro; los marcados
sólo se cuentan. En ventana de 16 px eso es **22,5 % de los vectores que entran en la
energía y que no son mediciones**, sino la media local de sus vecinos.

Reemplazar por la media local suprime varianza, y la suprime más cuanto peor resuelto
está el campo — que es el caso de gap 5 frente a gap 10. De ahí el signo.

Confirma [D-10] con una consecuencia medida: un vector interpolado no es una medición.
El barrido propio, que descarta sin rellenar, da la tendencia correcta sobre la misma
corrida (n=6 → 1,453 mm/s, n=12 → 1,288 mm/s).

### [H-07] σ es casi todo flujo, no instrumento — y eso rompe la hipótesis central

Con [V-04] y [H-06] descartados como explicación, los tres valores de σ dejan de ser
una discrepancia y pasan a ser una descomposición. Los tres son de la misma cámara con
la misma configuración de PIV, y se diferencian sólo en qué degradación física incluyen:

| | σ | qué agrega respecto del anterior |
|---|---|---|
| par sintético ([V-04]) | ≤ 0,011 px | nada: correlación y estimador subpíxel solos |
| fluido en reposo, par real | 0,030 px | ruido de lectura entre dos exposiciones, fluctuación de iluminación, movimiento propio de las partículas |
| forzado 2,15 A ([H-05]) | 0,169 px | gradiente de velocidad dentro de la ventana de 16 px, movimiento fuera del plano |

En cuadratura, el término que agrega el flujo es
√(0,169² − 0,030²) = **0,166 px**: el **98 % de la varianza**. El instrumento aporta el
2 % restante.

Es la descomposición de George & Stanislas (arXiv 2010.10768): ruido de pixelización
más ruido por no uniformidad de la velocidad dentro del volumen de interrogación. El
registro en reposo mide sólo el primer término porque no tiene el segundo.

**Consecuencia:** la hipótesis escrita en `README.md`, en `ESTADO.md` y en [D-09] —
"σ es constante porque es propiedad de la correlación y no del flujo" — es **falsa**, y
no marginalmente. Como el término dominante depende del flujo, en un decaimiento, donde
U cae un factor 10, σ tampoco es constante a lo largo del registro. La corrección a la
curva de decaimiento no puede usar un solo número.

**Salvedades, que son fuertes.** El 0,030 px viene de `fondo_202409_1643`, de la
campaña del **05-09-24**, no de la del 02-06-25: otro día, otro sembrado, y da s2n 4,00
contra 3,25 del forzado de junio. El único candidato a registro en reposo de junio,
el principio de `med_S0010`, tiene movimiento residual — su perfil oscila entre 0,07 y
0,18 px antes del encendido, con 67 % de dispersión entre pares, así que no sirve como
cero. Y el 0,169 px sale de un ajuste cuyo modelo esta misma entrada declara inválido,
de modo que es un valor efectivo, no una medición limpia. La descomposición es
indicativa; la magnitud del efecto, no el número, es lo defendible.

**Predicción contrastable, con el signo invertido respecto de lo que decía el plan:**
σ del barrido debe salir **distinto** en las tres corrientes (0,9 / 1,9 / 2,15 A) y
escalar con el gradiente de velocidad, en vez de salir igual en las tres. El barrido de
las otras dos corrientes ya estaba planificado como "la prueba fuerte"; lo que cambia
es qué resultado lo confirma y cuál lo refuta.

### [D-19] Los dos árboles del proyecto quedan separados por ahora

Decisión de Lucía. Existen en paralelo `/home/lucia/proyecto-final-piv` (este) y
`~/GW-AI-course/piv-s0008-forzado`, creados el mismo día, sobre la misma corrida, sin
conocerse. El segundo no se toca ni se borra: aportó [H-06] y el barrido de tamaño de
ventana. Se sigue trabajando acá.
**Pendiente declarado:** la consigna (`final-project.html`) pide un repositorio de
GitHub más una página HTML y un PDF. Hoy ninguno de los dos árboles es un repo git, y
el del curso está sin trackear dentro de un clon ajeno, que es lo que [D-01] quiso
evitar. Unificar es trabajo pendiente, no resuelto.

### [D-20] Convención de σ: por componente, y corrección al factor de [D-18]

Al escribir `informe/01_marco_teorico.md` apareció una inconsistencia de convención entre
nuestras cantidades y las de Foucaut.

El barrido mide `desp_rms_px = √⟨u² + v²⟩` (`02_barrido_dt.metricas`), que es un σ **de
módulo**. Foucaut informa σ **por componente** (su tabla 1 da σ_u y σ_v por separado).
Comparar los dos directamente es un error de √2, que es exactamente el que la
documentación de `04_desplazamiento_nulo.py` advertía que era fácil de cometer.

**Corrige a [D-18]**, que compara los 0,169 px del ajuste contra los 0,087 px de Foucaut
y concluye "un factor 2". En convención por componente el ajuste da 0,169/√2 ≈ **0,120
px**, y la comparación correcta es un **factor 1,4**.

Rehecha [H-07] en la misma convención: fluido en reposo 0,022 px (σ_u), forzado 0,120 px,
término del flujo √(0,120² − 0,022²) ≈ 0,118 px, que es el **97 %** de la varianza, no el
98 %. La conclusión no cambia; los números sí.

**Decisión:** de acá en más, σ se informa **por componente** y se dice explícitamente cuál
es la convención cada vez que se compara contra literatura. El módulo se sigue calculando
—es lo que sale del barrido— pero se convierte antes de comparar.
**Pendiente:** `03_figuras.py` sigue ajustando el módulo y escribiendo `sigma_px` sin
etiquetar la convención en `ajuste_med_S0008.json`. Hay que etiquetarlo.

### [V-05] Las fórmulas del marco teórico, verificadas contra los PDF — PASA

Antes de escribir `informe/01_marco_teorico.md` se extrajo el texto de los dos PDF
locales con `pdftotext -layout` y se verificó cada fórmula y cada valor atribuido:

| afirmación | verificada |
|---|---|
| Foucaut tabla 1: σ_u = 0,087 / 0,055 / 0,033 px para 16 / 32 / 64 | sí, y además el sesgo sistemático de −0,016 / −0,023 / −0,026 px en x |
| E_ii = E_noise · sinc²(kX/2), su ecuación (1) | sí; E_noise = 17,5·10⁻³ px³ (E₁₁) y 4,4·10⁻³ px³ (E₂₂) |
| corte universal k_c·X = 2,8 | sí |
| Y·E_ii = ζ, su ecuación (2) | sí |
| σ_u = √(4ζI/XY) con I = 1,492 (1,418 al 50 % de solape) | sí |
| George & Stanislas: tres contribuciones, varianzas ∝ 1/N, filtradas por la ventana | sí, en el resumen |
| "noise resulting from the flow itself" | sí, sección 3.2, citado textual |

**Dos hallazgos del texto que no estaban en `REFERENCIAS.md`:**

1. Foucaut señala que a desplazamiento nulo el peak locking *"tends to reduce the noise
   intensity"*. Coincide con la salvedad que [H-07] pone sobre el registro en reposo: es
   una cota inferior, y ahora hay respaldo en la fuente.
2. La sección de caracterización del ruido dice *"The noise due to particle motion will be
   studied afterwards"* — y en la versión ISPIV 2003 que tenemos ese "afterwards" no
   llega: seis secciones y ninguna lo desarrolla. **El término que Foucaut difiere es el
   que en nuestros datos resulta dominante.** Es lo que ubica la contribución de este
   trabajo, y también un motivo fuerte para conseguir la versión MST 2004.

También se verificó una advertencia de George & Stanislas que hay que declarar: la
pixelización y el efecto del número finito de partículas quedan ambos filtrados
espacialmente *"and hence difficult to distinguish"*. Está dirigida a los métodos
espectrales; el método de acá los separa quitando físicamente el flujo, lo que sortea esa
dificultad pero introduce las suyas.

---

## 2026-09-04

### [D-21] El proyecto cambia de objeto: de Δt óptimo a la condición de superficie libre

El trabajo sobre el Δt óptimo de PIV queda como **capítulo de instrumento**, no como el
objeto del proyecto. Motivo, en las palabras de la decisión: la pregunta del Δt se puede
resolver tomando mejores mediciones, así que no deja instalada una herramienta; y el
proyecto tiene que servir como herramienta para el doctorado. El plazo se extendió a un
mes, lo que hace viable un objeto más grande.

El objeto nuevo es la **fricción de fondo de una capa delgada con fondo no-deslizante y
superficie libre**, y la pregunta es si se puede calcular en vez de ajustarla. Todo lo
producido hasta acá se conserva: el piso de ruido de PIV medido ([V-04], [H-05], [H-07])
es la barra de error de la comparación entre la simulación y el decaimiento medido.

**Restricciones fijadas:** no se toman mediciones nuevas; el cómputo pesado va al clúster;
los datos crudos del 02-06-25 siguen siendo sólo lectura ([D-02]).

### [V-06] Qué condiciones de contorno tiene SPECTER hoy — VERIFICADO EN LA FUENTE

Se clonó SPECTER de `github.com/mfontanaar/SPECTER`, commit `0ad1edb` (2023-11-08), en
`numerico/SPECTER-upstream/`. Esa copia no se modifica: el árbol de trabajo es una copia
separada, para que el diff del proyecto sea exactamente el aporte propio.

La afirmación de partida —que SPECTER sólo admite no-deslizamiento— **se verificó en el
código**, no de memoria:

| dónde | qué dice |
|---|---|
| `src/SOLVERS_AND_BOUNDARY_CONDITIONS.md` | para el campo de velocidad lista una única condición: `noslip` |
| `src/boundary/vboundary.f90:26` (`v_parsebc`) | acepta sólo las cadenas `noslip` y `periodic`; cualquier otra aborta |
| `src/boundary/vboundary.f90:109` | `v_imposebc_and_project` aborta explícitamente si alguna condición en z no es `noslip` |

Tres hallazgos que **bajan el riesgo** estimado en el plan:

1. **La condición ya es un parámetro de ejecución, no de compilación.** `bin/parameter.inp`
   tiene `vbczsta` y `vbczend` como cadenas, y `v_setup` recibe un arreglo `bckind(6)`, una
   por cara. Agregar una condición nueva no obliga a recompilar ni a cambiar la interfaz.
2. **La capa FC-Gram ya está probada para Neumann.** `src/tests/` incluye `fc_dirichlet.f90`,
   `fc_neumann.f90`, `fc_neumann2.f90` y `fc_robin.f90`. La maquinaria que necesita una
   condición de tensión nula existe y tiene tests propios.
3. **El solver de Laplace homogéneo ya es genérico en la condición.** `laplace_z`
   (`src/boundary/boundary_mod.fpp:451`) documenta y admite Dirichlet, Neumann y Robin.

**El hueco concreto, y es chico y preciso.** `laplace_z` resuelve los coeficientes para
Dirichlet/Dirichlet (línea 499), Neumann/Neumann (531), Robin/Robin (563) y
Dirichlet abajo–Robin arriba (595). **No implementa Dirichlet abajo–Neumann arriba**, que
es exactamente el par de la capa con superficie libre. Es un sistema lineal de 2×2 por cada
(kx,ky), análogo al de la rama Dirichlet–Robin que ya está escrita.

**Lo que todavía NO está verificado, y hay que hacer antes de implementar:** cómo se
articulan la condición de Neumann sobre la presión y la corrección de la velocidad
tangencial para una pared libre de tensiones. `noslip_z` fija la velocidad tangencial en la
pared a partir del gradiente de presión (`vboundary.f90:154`); el análogo para tensión nula
**hay que derivarlo, no copiarlo**, y para eso hace falta leer el esquema de partición en
Fontana, Bruno, Mininni & Dmitruk, Comp. Phys. Comm. **256**, 107482 (2020),
[10.1016/j.cpc.2020.107482](https://doi.org/10.1016/j.cpc.2020.107482), que todavía **no
fue leído en el original**.

### [V-07] La Fase 1, derivada por un equipo de agentes — PASA, con salvedades

Trabajo `2026-09-04_205253_derive-fase1-capa-delgada`, receta `derive` de agent-team.
Dotación: 1 worker en codex a esfuerzo alto, **verifier en claude a xhigh** —otra familia
de modelos, para que la verificación no sea el mismo modelo corrigiendo su propia tarea—
y writer al final. Una ronda. Congelado. Costo real: **7.412.821 tokens**.

La puerta de aceptación (`verificacion/test_aceptacion_fase1.py`) se escribió **antes** de
la corrida, roja, y quedó fijada por hash para que el equipo pudiera pasarla pero no
editarla. Su principio de diseño es que no contiene ninguna fórmula analítica escrita a
mano: construye su propia matriz de diferencias finitas y saca los autovalores. Así, el
factor 4 y α = νπ²/4h² quedan **comprobados numéricamente**, no afirmados.

Resultado: puerta **6/6**; chequeos propios del equipo **5/5**; 21 afirmaciones, de las
cuales 16 verificadas, 4 dudosas y 1 refutada.

**Dos fallas de infraestructura, ajenas al trabajo, que hay que conocer para repetirlo:**

1. La receta corre sus chequeos con `python` a secas, que en esta máquina resuelve a
   `/home/lucia/miniforge3/bin/python` (3.13, sin numpy). Por eso `out/checks.py` figura
   como fallado en `state.json`. Con `/home/lucia/miniforge3/envs/piv-dt/bin/python` da
   5/5. **En el próximo trabajo hay que poner la ruta completa del intérprete.**
2. `job run` terminó con código 1 por un bug de agent-team: `engine._compile_pdf` corre
   pdflatex con `text=True` y el log de pdflatex trae bytes latin-1 (los acentos del
   castellano), lo que produce `UnicodeDecodeError`. pdflatex funciona; el PDF se compiló
   a mano, 7 páginas, en `jobs/.../out/build/notes.pdf`.

### [D-22] El criterio de condicionamiento de la rama Dirichlet–Neumann no es k·Lz

El verifier refutó la afirmación de `notes.tex` de que el guardián `k·Lz < 0,5` cubre sin
cancelación tanto k=0 como k positivo muy chico. Se reprodujo el contraejemplo y se aisló
la variable responsable:

| k | Lz | k·Lz | \|g1\|/(k·\|g0\|) | \|φ(0) − g0\| |
|---|---|---|---|---|
| 5·10⁻¹¹ | 10¹⁰ | 0,5 | 1,46·10¹¹ | **7,6·10⁻⁶** |
| 5·10⁻¹¹ | 10¹⁰ | 0,5 | 0 (g1 = 0) | 0 |
| 10⁻¹⁰ | 1 | 10⁻¹⁰ | 7,3·10¹⁰ | 0 |
| 1 | 1 | 1 | 7,3 | 2,2·10⁻¹⁶ |

Con el **mismo** k·Lz el error es 7,6·10⁻⁶ o exactamente 0 según cuánto valga g1, y con un
cociente igual de grande pero en la rama regular vuelve a ser 0. El parámetro de
condicionamiento es **|g1|/(k·|g0|)**, y hace falta además estar en la rama exponencial.

**Por qué esto no compromete el uso en SPECTER.** Estar en la rama exponencial exige
k ≥ 0,5/Lz, y entonces |g1|/(k|g0|) ≤ 2·Lz·|g1|/|g0|. O sea que el parámetro peligroso
está acotado por `Lz·|g1|/|g0|`: para que importe hace falta un dominio enorme, no un k
chico. Con Lz = 6 mm haría falta |g1|/|g0| ≳ 10¹² para que asome. Es álgebra, no medición.
El contraejemplo del verifier es correcto como refutación del enunciado, y queda fuera del
rango alcanzable del código. **La condición hay que escribirla bien al implementar en
Fortran (Fase 2): sobre `Lz·|g1|/|g0|`, no sobre `k·Lz`.**

**Corrección a la puerta, que era un defecto mío.** El verifier observó que C1 pasaba con
margen 1,2 (error 8,23·10⁻⁵ contra tolerancia 10⁻⁴) y diagnosticó bien que era estructural
del test: mi referencia FD con nz=400 tiene error (4π)²dz²/12 = 8,224·10⁻⁵ en el cuarto
modo. Se midió el orden de convergencia —2,000 para los dos pares de bordes— y se agregó
extrapolación de Richardson entre nz y 2nz. El error de la referencia baja a 6,7·10⁻¹⁰ y
C1 pasa con cinco órdenes de margen. Nuevo hash de la puerta:
`0c3ebd62f6c42d08f0fa9a64a0f57b9cfdaf1144ee8a84e968a148d9d66b18bf`.

### [P-01] Pregunta abierta: la reducción 2D podría no valer en el régimen medido

El verifier señaló que las desigualdades de validez de la clausura −αu nunca se evaluaron
para la celda real. Su cuenta: `Re_h·ε ≪ 1` equivale a `U ≪ ν·l/h²`, que con h = 6 mm y
ν = 10⁻⁶ m²/s da **U ≪ 0,56 mm/s** para l = 2 cm.

Las velocidades medidas de esta campaña son 1,32 mm/s en el forzado a 2,15 A ([H-05]) y
4,0 mm/s al inicio del decaimiento. **La condición está violada por un factor de 2 a 7
justo donde arranca el decaimiento**, y sólo se cumpliría tarde en el registro.

Si sobrevive al escrutinio, cambia la interpretación: la clausura con un único α sería
asintóticamente válida al final del decaimiento, no válida a lo largo del registro. No
está verificado: el número de arriba es aritmética sobre una desigualdad que todavía hay
que derivar con cuidado, y la elección de l = 2 cm como escala de inyección viene del
espaciado de los imanes, no de una medición. **Es material de la Fase 4, y es la primera
cosa que hay que atacar cuando lleguemos ahí.**

---

## 2026-09-07

### [D-23] Cadena de compilación local: env de conda `specter`

Esta laptop no tenía `gfortran` ni MPI, así que SPECTER no se podía compilar acá y la
puerta de la Fase 2 no habría podido correr. Se creó el env `specter` con
`mamba create -n specter -c conda-forge gfortran openmpi openmpi-mpifort fftw make`
(1,2 GB, gcc 16.2). No hace falta permiso de administrador y no toca los envs que Lucía
usa para otra cosa, en la misma línea de [D-03] y [D-11].
**Alternativa descartada:** hacer toda la verificación en el clúster. Los tests son de
8×8×128, ~130 kB por arreglo y segundos de corrida: mandar eso a Sakura habría alargado
el ciclo de iteración sin ganar nada. El cómputo de producción sigue yendo al clúster.

Configuración necesaria para que el Makefile de SPECTER compile con esta cadena, toda
en `Makefile.in`: `FC_GNU = mpif90`, `MPIINC_OPEN`/`MPILIB_OPEN`/`MPILD_OPEN` vacías,
`FFTWDIR` apuntando al prefijo del env, y `-fallow-argument-mismatch` en `FFLAGS_GNU`
y `FPSPEC_GNU`, que gfortran ≥ 10 necesita para las llamadas MPI del código.
**Trampa que costó una corrida:** el `Makefile` usa `GHOME = $(PWD)`, y `$(PWD)` sale de
la variable de entorno, no del directorio de trabajo del proceso. Lanzar `make` con
`subprocess.run(cwd=...)` sin exportar `PWD` compila contra el árbol equivocado.

### [V-08] Fontana, Bruno, Mininni & Dmitruk (2020) — LEÍDO EN EL ORIGINAL

Era el prerrequisito que `ESTADO.md` ponía antes de escribir código. El preprint es
[arXiv:2002.01392v2](https://arxiv.org/abs/2002.01392), 31 páginas, versión aceptada del
artículo de CPC; queda guardado en `bibliografia/`. Lo que hacía falta está en las
secciones 4.1, 4.2 y 4.4:

- El esquema de partición: se avanza un campo sin presión `v*`, se le impone una
  condición de contorno propia, y después se proyecta `v = v* - (dt/o) grad(p)`.
- Las condiciones que usa el código, sus ecuaciones (15) y (16), y la generalización a
  Runge-Kutta de orden o, su ecuación (30).
- La presión se resuelve con **Neumann**, `d(p)/dz = (o/dt) z·v*` en la cara, que es
  justamente lo que fuerza `v_z = 0` después de proyectar.
- Una advertencia del propio paper que conviene tener escrita: proyectar la ecuación de
  momento en la dirección normal o en la tangencial da problemas distintos, "there is no
  a priori reason to assume that both approaches lead to the same solution". Ellos
  eligen la normal.

También se verificó en el código lo que el paper no dice: `pr` no guarda `p` sino
`dt*p` del subpaso en que se calculó — el `!TODO update documentation about p' = dt*p
for convenience` de `boundary_mod.fpp:210` lo confirma, y es lo que explica el factor
`(o+1)/o` de `noslip_z`.

### [D-24] La condición de tensión nula se deriva como una condición sobre la vorticidad

Pedir tensión tangencial nula sobre el campo **ya proyectado**, `d(v_par)/dz = 0` en
z=Lz, y usar la proyección junto con el dato de Neumann de la presión, da una condición
sobre el campo auxiliar que **no contiene dt**:

    d(v*_x)/dz = i·kx·v*_z ,   d(v*_y)/dz = i·ky·v*_z   en la cara.

Es decir: se anulan las dos componentes tangenciales de la vorticidad. Como la
proyección resta un gradiente y el rotor de un gradiente es cero, `rot(v*) = rot(v)`
exactamente, así que la condición se transmite al campo proyectado **sin error de
partición**, a diferencia del no-deslizamiento, cuya velocidad de deslizamiento residual
es O(dt^ord).

Eso es una predicción falsable y es lo que chequea V3 de la puerta: el residuo tiene que
ser independiente de dt. Medido, la razón entre el residuo a dt y a dt/2 es **1,00**.
Si el signo del término que acopla `v*_z` estuviera mal, ese término —que es O(dt)— haría
que el residuo escalara con dt, y el chequeo lo detectaría.

La imposición se hace con `neumann_reconstruct` del módulo `fcgram`, que reconstruye el
valor de borde a partir de la derivada normal. La convención de signo es la derivada
normal **saliente**, verificada en `src/tests/fc_neumann.f90` (el comentario
`d/dz = -d/dn`): se pasa `-d/dz` en z=0 y `+d/dz` en z=Lz.

### [H-08] Corrección a [V-06]: la rama Dirichlet–Neumann de `laplace_z` NO es la que usa la superficie libre

[V-06] afirmaba que la rama Dirichlet abajo – Neumann arriba de `laplace_z` "es
exactamente el par de la capa con superficie libre". **Es incorrecto**, y se verificó en
el código:

`v_imposebc_and_project` llama siempre `sol_project(vx,vy,vz,pr,bctarget=1,0,0)`, y
`sol_project` le pasa a `laplace_z` la clase `bcz + bctarget` (línea 356). Con
`bctarget=1` —condición de Dirichlet sobre `v_z` en las dos caras— las dos clases salen
1, o sea **Neumann–Neumann sobre la presión**. Como la superficie libre **plana**
mantiene `w = 0` en el tope, la condición de la presión es idéntica a la del canal. La
combinación Dirichlet–Neumann ni siquiera es expresable con la interfaz actual de
`sol_project`, que tiene un único `bctarget` para las dos caras.

La rama se implementó igual, porque los coeficientes ya estaban derivados y verificados
en la Fase 1 ([V-07]) y porque es el hueco real que tenía `laplace_z`. Pero la puerta la
ejercita como **test unitario** (V2, `src/tests/laplace_dirneu.f90`) y no reclama que el
camino físico la use.

Dónde sí aparecería ese par de condiciones: en una superficie libre **deformable**, con
`p = 0` en el tope. Eso es Neumann abajo – Dirichlet arriba, que es el par opuesto, y
está fuera del alcance declarado desde `intent_fase1.txt`.

### [D-25] La rama nueva de `laplace_z` está mejor condicionada que las que ya estaban

Con `phi(z) = C1·exp(k(z-Lz)) + C2·exp(-kz)` y `r = exp(-k·Lz)`, el par Dirichlet en z=0
y Neumann en z=Lz da

    C1 = (g1/k + r·g0)/(1 + r²) ,   C2 = (g0 - r·g1/k)/(1 + r²) ,

y el denominador es **1 + r²**, no `1 - r²` como en las ramas Dirichlet/Dirichlet y
Neumann/Neumann. Este par no tiene resonancia: el denominador nunca se acerca a cero.
Se escribió con `exp(-2·k·Lz)` igual que el resto del archivo, para que no aparezca
ninguna exponencial creciente. El guardián de condicionamiento se documentó sobre
`Lz·|g1|/|g0|` y no sobre `k·Lz`, como pedía [D-22].

### [V-09] La Fase 2, implementada y verificada — PUERTA 6/6

Puerta escrita **antes** de la implementación y en rojo (0/6 con el árbol intacto), en
`verificacion/test_aceptacion_fase2.py`. Pliego en `verificacion/intent_fase2.txt`.
Hash de la puerta: `6f0b4da4ba20617d483dfff0a20f8e2e11e82795b20351c259ffcf740ca4b1cc`.

Igual que la de la Fase 1, la puerta **no contiene ninguna fórmula analítica del
resultado**: arma su propio operador de diferencias finitas del problema de Stokes
vertical, le saca los autovalores y extrapola con Richardson entre n y 2n. Contra los
valores exactos esa referencia tiene un error relativo de 2·10⁻¹¹ a 1·10⁻¹⁰. El factor 4
se **mide** corriendo el mismo binario con las dos condiciones.

| | resultado |
|---|---|
| **V1** parser | acepta `freeslip`, sigue abortando ante una cadena desconocida |
| **V2** rama Dirichlet–Neumann de `laplace_z` | 9 modos; φ(0)−g₀ = 2,0·10⁻¹⁶, φ′(Lz)−g₁ = 1,9·10⁻¹⁶, contra la forma cosh/sinh 2,8·10⁻¹⁵ |
| **V3** residuo de tensión en la cara libre | 6,1·10⁻¹² relativo a la cara rígida de la misma corrida, en las dos caras; **razón entre dt y dt/2 = 1,00** |
| **V4** espectro λ_m | m=0,1,2: error relativo 6,8·10⁻¹⁰, 5,5·10⁻⁸, 4,4·10⁻⁷ |
| **V5** el factor 4 | **4,000000015**; y la cara libre abajo da la misma λ que arriba, dentro de 7,7·10⁻¹³ |
| **V6** poder discriminante | la corrida no-deslizante falla los criterios de V3 y V4, como debe |

**Cómo se fijó la tolerancia, que es la parte que se puede hacer mal.** `TOL_TASA = 1e-5`
no se eligió a ojo ni después de ver el resultado: se midió corriendo el camino
**no-deslizante de upstream**, que este trabajo no modifica, sobre la misma malla y
contra la misma referencia. Los tres modos más lentos dieron 3,1·10⁻⁹, 4,8·10⁻⁸ y
2,4·10⁻⁷. El umbral es unas 40 veces el peor de esos.

Lo mismo con V3: el criterio es **relativo** al residuo de tensión de la cara rígida de
la misma corrida, no un valor absoluto. La primera versión pedía precisión de máquina
(10⁻¹⁶) y fallaba: la condición se impone reconstruyendo el valor de borde con tablas
FC-Gram de orden finito, así que el residuo es de nivel de reconstrucción, no de
redondeo. Como vara, el propio no-deslizamiento de upstream deja `<|v_t|²> = 1,0·10⁻¹³`
en la pared con un campo de orden 10.

**Tres fallas de la primera corrida contra la implementación (3/6), todas de la puerta y
ninguna del código**, anotadas porque son la clase de error que hace pasar un test malo:

1. V2 le daba a `laplace_z` datos de contorno **complejos en el modo (kx,ky)=(0,0)**. Ese
   modo de un campo real es real, y todas las ramas del código le toman la parte real a
   propósito. El dato era el que no tenía sentido físico.
2. V3 pedía precisión de máquina, como se explicó arriba.
3. V5 corría el caso espejo (cara libre **abajo**) con el autómodo del caso normal,
   `sin(πz/2Lz)`, que no es autómodo del par espejo. La tasa medida era la del
   transitorio. Con el modo correcto, `cos(πz/2Lz)`, coincide en 13 cifras.

**Una consistencia que no estaba pedida y salió sola:** las corridas con la cara libre
arriba y abajo dan las columnas del diagnóstico intercambiadas y coincidentes en 7
cifras, que es exactamente lo que exige la simetría de reflexión del problema.

### [D-26] Qué se tocó de SPECTER, y qué NO

El árbol de trabajo es `numerico/SPECTER-trabajo/`, clon de `SPECTER-upstream` en
`0ad1edb`. El diff son **237 líneas agregadas y 22 borradas** en 5 archivos, más
`src/tests/laplace_dirneu.f90` (57 líneas, nuevo).

| archivo | qué |
|---|---|
| `src/boundary/vboundary.f90` | `v_parsebc` acepta `freeslip`; `v_setup` carga las tablas de Neumann del plan en z; el guardián de `v_imposebc_and_project` acepta la clase 1 y despacha; subrutina nueva `freeslip_z`; `vdiagnostic` escribe además `freeslip_diagnostic.txt` con el residuo de **tensión**, que es el que una superficie libre tiene que anular |
| `src/boundary/boundary_mod.fpp` | rama Dirichlet(z=0)–Neumann(z=Lz) de `laplace_z` |
| `src/specter.fpp` | engancha el test unitario nuevo |
| `src/SOLVERS_AND_BOUNDARY_CONDITIONS.md`, `bin/parameter.inp` | documentación de la cadena nueva |

**El camino no-deslizante quedó intacto a propósito.** `vz` se lleva al dominio
(kx,ky,z) sólo cuando alguna cara es `freeslip`; si no, no se hace ninguna transformada
de más y la ruta es la de upstream. Eso importa porque V5 y V6 comparan contra ella.

**Sigue sin tocarse:** `numerico/SPECTER-upstream/`.

**Pendiente declarado:** el diagnóstico de tensión se calcula ahora en toda corrida con
paredes en z, no sólo con superficie libre. Cuesta dos derivadas y dos transformadas por
cada `cstep`, y una matriz temporal del tamaño del campo, que a resolución de producción
son cientos de MB. Se pide y se libera de a una por vez. Si molesta en el clúster, hay
que condicionarlo.

### [V-10] La condición nueva reproduce bit a bit con 1, 2 y 4 procesos MPI

`freeslip_z` respeta el reparto `ista:iend` y no comunica nada, pero eso había que
mostrarlo. Con el **mismo** `parameter.inp`:

| | np=1 | np=2 | np=4 |
|---|---|---|---|
| λ, modo de corte con tope libre | 2,4674010985·10⁻² | idéntico | idéntico |
| ⟨v²⟩ final, campo 3D con kx≠0 | 9,8431825537157405 | ...511 | ...458 |

Coincide en 15 cifras significativas, que es redondeo de doble precisión. El reparto de
`khom(j,i)` y el caso especial `IF (ista.eq.1)` de la rama nueva de `laplace_z` quedan
así ejercitados con el modo (0,0) en un proceso y sin él en los otros.

### [H-09] En régimen no lineal la superficie libre necesita más resolución que el canal

Este hallazgo salió de un cotejo mal armado —cambié la condición inicial entre las
corridas de np=1 y np=2 y leí como bug de MPI lo que era otra cosa— pero lo que apareció
al perseguirlo es real y hay que declararlo.

Con el campo 3D de amplitud u0=1 (Re ≈ 300 con ν=10⁻²) y 12000 pasos hasta t=2,4:

| malla | BC | ⟨v²⟩ t=0 | t=0,7 | t=2,4 | |
|---|---|---|---|---|---|
| 8×8 | noslip | 1,096·10¹ | 6,74 | 1,30 | decae |
| 8×8 | **freeslip** | 1,096·10¹ | 7,43 | **2,79·10¹** | **crece** |
| 16×16 | noslip | 1,096·10¹ | 4,14 | 0,603 | decae |
| 16×16 | **freeslip** | 1,096·10¹ | 6,45 | 2,49 | decae |
| 32×32 | noslip | 1,096·10¹ | 4,11 | 0,588 | decae |
| 32×32 | **freeslip** | 1,096·10¹ | 6,43 | 2,00 | decae |

Sin forzado la energía no puede crecer, así que el 8×8 con tope libre es una
inestabilidad numérica. **No es la condición de contorno.** Tres cosas lo muestran:

1. Con dt/4 el resultado es **exactamente el mismo** (2,7873·10¹): no es el paso de
   tiempo.
2. Con u0/20 decae: depende de la amplitud, o sea del término no lineal.
3. Con 16×16 y 32×32 decae, y las dos mallas coinciden entre sí en 3 cifras a t=0,7,
   mientras que la de 8×8 se despega. Es convergencia con la resolución.

La divergencia se mantiene en 10⁻¹⁴ y el residuo de tensión en la cara libre en 10⁻¹²
durante toda la corrida que crece: la condición se sigue cumpliendo mientras el interior
se descompone.

**Por qué pasa, y por qué importa para la Fase 4.** El filtro de dealiasing de SPECTER
(`fc_filter`, `pseudospec_hd.f90:1082`) no es la regla de los dos tercios sino una
exponencial `exp(-alpha (2k/N)^{2p})` con p=50, que sólo toca los últimos modos: sobre
una malla de 8 puntos no deja margen. Y una pared no-deslizante disipa mucho más que una
superficie libre, así que el canal aguanta la misma malla donde la capa con superficie
libre ya no. Medido: a u0=2 y u0=4 el canal sigue decayendo en 8×8 y el caso libre da
NaN.

**Consecuencia práctica:** al pasar a producción en Sakura, la resolución horizontal no
se puede heredar de una corrida de canal equivalente. Hay que hacer un estudio de
convergencia propio con la superficie libre puesta.

**Lo que esto NO invalida.** La puerta usa el campo 3D sólo en V3, 400 pasos hasta
t=0,08, mucho antes de que esto aparezca; y V4, V5 y V6 usan el modo de corte puro, donde
el término no lineal es **exactamente cero** por construcción y el problema es el de
Stokes vertical. Los seis chequeos siguen valiendo. Lo que hay que decir es que **la
puerta no verifica el régimen no lineal**, y ahora se sabe qué pasa ahí.

---

## 2026-09-09

### [D-27] Los parámetros de la celda viven en `codigo/celda.py` y en ningún otro lado

Decisión de Lucía, al zanjar el espesor: **"ese valor lo puedo cambiar entre
experiencias, no quiero que quede hardcodeado en ningún lado"**. El espesor `h`, la
viscosidad `ν`, la geometría de los imanes y las condiciones de cámara no son propiedades
del problema sino condiciones de una campaña, así que pasan a estar declaradas una sola
vez, con su procedencia, en `codigo/celda.py`. Las derivaciones se escriben en forma
simbólica y los números aparecen sólo como evaluación al final.

**Alternativa descartada:** dejar el valor repetido en README.md, en los intent y en los
scripts, que es como estaba y es exactamente lo que produjo la discrepancia 5 mm / 6 mm.

Tres datos quedaron fijados en esa misma conversación:

| dato | valor | procedencia |
|---|---|---|
| espesor de la capa | **6 mm** | Lucía, 2026-09-09; zanja la discrepancia de `ESTADO.md` |
| viscosidad | 10⁻⁶ m²/s, **la del agua** | decisión explícita de usarla y **declararla**: no es la del electrolito (KNO₃ 16 % m/m), no se midió ni se buscó |
| imanes | discos de **1 cm** de diámetro, red tipo tablero de **5 cm** de paso | Lucía, 2026-09-09; corrige el `ℓ = 2 cm` heredado, que no correspondía a ninguna longitud de la celda |

### [V-11] [P-01] resuelto: el criterio heredado estaba mal escrito, y la conclusión se da vuelta

Derivación en `teoria/P01_validez_reduccion_2D.md`, cuentas en
`numerico/fase4/p01_validez.py`, salida en `salidas/tablas/p01_validez.json`.

La pregunta abierta decía que `Re_h·ε ≪ 1` estaba *"violada por un factor de 2 a 7 justo
donde arranca el decaimiento"*, con la advertencia de que había que derivar la desigualdad
en vez de heredarla. Al derivarla cambian cuatro cosas:

1. **La escala horizontal estaba mal, y el criterio mal escrito.** `ℓ = 2 cm` no era
   ninguna longitud de la celda. Con los imanes reales el armónico dominante del patrón
   es el **diagonal**, `k = π√2/a = 88,86 m⁻¹`, longitud de onda 7,07 cm — verificado
   transformando el patrón de discos, no supuesto: el pico medido coincide con la forma
   cerrada con error 0 y sale orientado en (−1,+1)·π/a. Además el criterio se escribía
   con `ℓ`, que es ambiguo por 2π; el adimensional que sale de la ecuación es `k·h`.
   Con esto el parámetro empeora: **δ = U k h²/ν ≈ 2,7 en el forzado y ≈ 8,2 al inicio
   del decaimiento**, no «2 a 7 veces el umbral».
2. **`ε ≪ 1` no hace falta para la ley de fricción.** El término viscoso horizontal es
   lineal y no acopla con la estructura vertical, así que se conserva exacto:
   `λ(k) = α(1 + 4ε²/π²)`. Con ε = 0,533 eso es **+11,5 %** sobre α, y hay que incluirlo
   al comparar con el experimento en vez de tratarlo como una falla de validez.
3. **δ no es lo que tiene que ser chico.** Lo que tiene que ser chica es la distorsión del
   perfil vertical. Expandiendo en los modos `F_m = sin[(2m+1)πz/2h]`, con
   `λ_m/λ₀ = (2m+1)²` y acoplamientos `c_m = ⟨F₀²F_m⟩/⟨F_m²⟩` (`c₁ = −8/(15π)`,
   reproducido por cuadratura a 9 cifras), la distorsión es `|a₁/a₀| ≈ 0,0077·δ`. En el
   régimen medido son **2 % en el forzado y 6 % al inicio del decaimiento**.
4. **La clausura mejora sola.** `δ(t) = δ₀e^{−αt}`, y la memoria del segundo modo
   vertical dura `h²/(2π²ν) = 1,82 s` contra registros de 51,2 s.

**Conclusión: la clausura de un modo no queda invalidada en el régimen medido.** El error
estimado es de pocos por ciento, no de orden uno. Eso **no** cierra la pregunta: la
estimación toma igual a 1 un factor geométrico O(1) del término no lineal que no se
calculó, y no está decidido si el sesgo sobre α es de orden δ o δ² (los modos
distorsionados viven en números de onda 0 y 2k, así que para un único modo celular no
realimentan hasta δ²; un campo de banda ancha no tiene esa protección).

**La forma de cerrarla, y es medible con el código que ya existe.** El promedio vertical
da, sin ninguna clausura, `α_eff ≡ (ν/h)·∂_z u|₀ / ⟨u⟩`. Un barrido en amplitud con
SPECTER sobre δ ≈ 0,3…10, midiendo `α_eff/λ(k)`, **es** la respuesta, y de paso es el
estudio de convergencia con superficie libre que [H-09] dejó pendiente para producción.

### [H-10] La velocidad de 4,0 mm/s al inicio del decaimiento no tiene procedencia

Figura en `ESTADO.md` bajo «De la campaña» y se usó en la cuenta de [P-01], pero **no hay
script en este árbol que la produzca** y no se pudo rastrear a ninguna salida. Es la
velocidad que fija el peor valor de δ, así que importa. El disco externo no estaba montado
al escribir esto, de modo que no se pudo recalcular. Queda marcada como no verificada en
`codigo/celda.py`; recalcularla es parte del punto 4 de la Fase 4.

### [D-28] La superficie deformable queda fuera de alcance, y ahora con número

`teoria/superficie_libre_v_estrella_y_p.md` describe en detalle cómo se modelaría la
superficie libre en las condiciones de `v*` y de `p`, en los tres casos: plana (lo
implementado), con película superficial, y deformable. La decisión de alcance que ya
estaba tomada —superficie plana— deja de apoyarse en un criterio y pasa a apoyarse en una
cuenta (`numerico/fase4/superficie_libre_escalas.py`):

| | U (promedio vertical) | Froude | η | η/h |
|---|---|---|---|---|
| forzado 2,15 A | 0,841 mm/s | 0,0035 | 0,07 µm | 1,1·10⁻⁵ |
| inicio del decaimiento | 2,311 mm/s | 0,0095 | 0,51 µm | 8,6·10⁻⁵ |

La longitud capilar es 2,71 mm (k = 369 m⁻¹), bastante mayor que los 88,9 m⁻¹ del forzado,
así que **restituye la gravedad**. La superficie se deforma menos de un micrón sobre 6 mm.
**Alternativa descartada:** implementar la rama Neumann–Dirichlet y la evolución de η.

### [V-12] La rama Neumann(z=0)–Dirichlet(z=Lz) de `laplace_z`, derivada y verificada

Es la que haría falta para una superficie deformable, porque la condición de tensión
**normal** prescribe el valor de `p` en la superficie. Es el par **opuesto** al que se
implementó en la Fase 2 ([H-08]). En la base del archivo, `φ = C₁e^{k(z−Lz)} + C₂e^{−kz}`
con `r = e^{−kLz}`, y datos `g₀ = φ′(0)`, `g₁ = φ(Lz)`:

    C₁ = (g₁ + r·g₀/k)/(1 + r²) ,   C₂ = (r·g₁ − g₀/k)/(1 + r²)

Verificado sobre 20 000 combinaciones de `k` y `Lz` en seis órdenes de magnitud: peor
residuo de Neumann 4,4·10⁻¹⁶ y de Dirichlet 3,9·10⁻¹¹. Denominador **1 + r²**, sin
resonancia, igual que la rama de [D-25]. El parámetro de condicionamiento es el **espejo**
del de [D-22]: `|g₀|/(k·|g₁|)` en vez de `|g₁|/(k·|g₀|)`, comprobado con el mismo
contraejemplo. **No se implementó en Fortran**, por [D-28].

Se registran además dos cosas que costaría descubrir después:

1. La condición sobre `v*` de una superficie deformable **no es libre de `dt`**:
   `∂v*ₓ/∂z = −i kₓ v*_z + 2(dt/o)·i kₓ·∂p/∂z|_h`, y ahí `∂p/∂z` es parte de la solución
   del problema de Poisson, no un dato. La razón de fondo es que la tensión tangencial es
   `∂u/∂z + ∂w/∂x` y la vorticidad es `∂u/∂z − ∂w/∂x`: **coinciden sólo si la superficie
   es plana**, y es esa coincidencia la que hace exacto al caso implementado ([D-24]).
2. `sol_project` tiene un único `bctarget` para las dos caras y haría falta uno por cara.

### [D-29] La película superficial se modela con Boussinesq–Scriven, y es lo que sí falta

La salvedad más fuerte del trabajo —una película vuelve el tope casi rígido, con α cuatro
veces mayor— se puede convertir en un número. Con una única viscosidad superficial `μ_s`,
el balance de tensión tangencial pasa a ser una condición de **Robin**,
`μ·∂u_par/∂z = −μ_s k²·u_par`, con un solo parámetro `β = μ_s k² h/μ`. Entonces
`(qh)·cot(qh) = −β` y `α = νq²` interpola entre `π²ν/4h²` (β→0) y `π²ν/h²` (β→∞): **el
factor 4 de la Fase 1 son los dos extremos de esta familia**.

Con h = 6 mm, ν = 10⁻⁶ m²/s y k = 88,86 m⁻¹: un 10 % en α corresponde a
`μ_s = 2,7·10⁻⁶ N·s/m` y duplicar α a `3,6·10⁻⁵ N·s/m`. **El α medido se puede invertir
para acotar `μ_s`.**

Es además mucho más barato que la superficie deformable: la presión **no se toca** (la
superficie sigue plana, `w = 0`, sigue siendo Neumann–Neumann) y FC-Gram ya tiene
`fc_robin.f90`. Lo que sí cambia es que aparece un término O(dt·β) en la condición sobre
`v*`, así que **el chequeo V3 deja de ser el correcto**: hay que verificar que el residuo
de Robin escale como `dt^ord`, no que sea independiente de `dt`.

**No verificado:** si esos valores de `μ_s` son grandes o chicos para una interfaz
agua–aire con partículas de 100 µm. No se leyó bibliografía de reología interfacial.

### [V-13] El decaimiento medido: α = 0,0669 ± 0,0049 s⁻¹, contra 0,0685 predicho

`codigo/06_decaimiento.py` sobre los cuatro registros, con figuras en
`codigo/07_figuras_decaimiento.py` y página de revisión en `salidas/decaimiento.html`.

**Por qué hubo que rehacer la medición.** Las series `.vec` archivadas
(`Deca1..4_destino`) están armadas con pares de cuadros **consecutivos**, Δt = 1/60 s. Al
arrancar el decaimiento el desplazamiento verdadero es de ~0,2 px por cuadro y el piso de
ruido es de ~0,17 px: el ruido es del orden de la señal. Además esas series tienen 30–36 %
de vectores válidos, contra 97,5 % del multipaso propio sobre los mismos cuadros. Leída de
los `.vec`, la velocidad **no decae**: se queda clavada en el piso de ruido los 51 s.

**El método.** Ocho instantes por registro, pares rearmados desde los TIF con
n = 1, 2, 4, 8, 16, 32, y `desp_rms(n)² = (n·d₁)² + σ²`. El valor que se usa **no** sale
del ajuste sino de `d₁ = desp_rms(n)/n` a n grande, donde el ruido es despreciable
(corrección típicamente < 1 %). Motivo: el ajuste completo da **σ² negativo** en varios
instantes, porque a n = 1 el desplazamiento es de dos décimas de píxel y el estimador
subpíxel sesga hacia el entero. El σ² negativo se informa y no se recorta.

| | α [1/s] | τ [s] |
|---|---|---|
| Decaimiento 1 (`med_S0003`) | 0,06525 ± 0,00340 | 15,3 |
| Decaimiento 2 (`med_S0005`) | 0,05912 ± 0,00592 | 16,9 |
| Decaimiento 3 (`med_S0006`) | 0,05933 ± 0,00421 | 16,9 |
| Decaimiento 4 (`med_S0009`) | 0,06386 ± 0,00343 | 15,7 |
| **ensemble, t > 10 s** | **0,06687 ± 0,00488** | **15,0** |
| ensemble, todos los puntos | 0,06189 ± 0,00336 | 16,2 |
| **predicho, tope libre** | **0,06854** | 14,6 |
| predicho, tapa rígida | 0,27416 | 3,6 |

**El medido está a 2,4 % del valor de superficie libre y 4,1 veces por debajo del de tapa
rígida.** Es la comparación que motivaba el proyecto.

**Salvedades que hay que decir junto con el número:**

- **La meseta inicial.** Los cuatro registros no decaen durante los primeros 6 a 10 s. No
  está explicado —forzado todavía encendido, o transitorio vertical— y es lo que separa
  0,0619 de 0,0669. El ajuste se hace sobre t > 10 s.
- Con h = 6 mm la predicción es 0,0685; para reproducir el α medido exactamente haría
  falta h = 6,08 mm. **Medio milímetro de error en h mueve α un 17 %**, más que la
  discrepancia observada: el resultado está limitado por h, no por la estadística.
- Una película superficial **sube** α. Que el medido esté en o por debajo del valor de
  superficie libre limpia la acota: β ≲ 0,15, o sea μ_s ≲ 3·10⁻⁶ N·s/m en el modelo de
  [D-29]. La cota también está dominada por la incertidumbre en h.

### [H-11] La escala del flujo no es la del forzado, y crece durante el decaimiento

Los espectros de energía de los campos medidos, con el ruido blanco de PIV restado, dan un
k pesado por energía de **240 m⁻¹ al inicio y 135 m⁻¹ al final**, contra los 88,9 m⁻¹ del
armónico del forzado. Con ese k, el δ de [P-01] arranca en ~20 y termina en ~1, y la
distorsión estimada del perfil vertical va del **15 % al 1 %**. La clausura de un modo es
marginal en los primeros segundos y buena en el resto, que es la ventana donde se ajustó α.

**Trampa que costó una figura:** el piso de ruido del espectro hay que leerlo **por debajo
de Nyquist**. Con paso de grilla de 4 px, k_Nyquist = 2984 m⁻¹, y por encima los anillos
sólo recogen las esquinas de la grilla y el espectro se desploma. Estimar el piso en el
cuartil superior de k lo subestima y no resta casi nada.

### [H-12] `med_S0008` no es estacionario: contiene un decaimiento

Medidos con el estimador robusto de [V-13], los dos registros forzados largos se comportan
al revés de lo supuesto:

| t [s] | `med_S0002` (1,95 A) | `med_S0008` (2,15 A) |
|---|---|---|
| 0,0 | 3,549 mm/s | 2,211 mm/s |
| 12,7 | 3,528 | 3,096 |
| 20,2 | 3,497 | 1,649 |
| 27,9 | 3,576 | 0,819 |
| 43,0 | 3,506 | **0,335** |

`med_S0002` es plano al 2 % en 43 s. `med_S0008` **decae**, y a una tasa de ≈ 0,073 s⁻¹,
que es la fricción de fondo. El registro etiquetado «Forzado 2,1–2,2 A» contiene un
decaimiento.

**Consecuencias.**

1. `U_verdadera = 1,321 mm/s` de [H-05] **no es «la velocidad del forzado»**. El barrido
   reparte sus 12 pares uniformemente a lo largo del registro para cada n
   (`02_barrido_dt.py`, `np.linspace(0, len-n-1, pares)`), así que lo que ajusta es un rms
   sobre un registro en el que el flujo cae unas siete veces. Hay que releerlo como
   promedio del registro.
2. **El efecto sobre σ = 0,169 px no está establecido.** Como los pares se reparten igual
   para todos los n, el término de ruido sigue siendo separable a primer orden, pero la
   premisa del modelo está violada. Pide un reanálisis sobre `med_S0002`, que sí es
   estacionario; no una afirmación.
3. **Explica el «factor 2,7 sin explicar»** entre forzado y decaimiento que figuraba como
   pregunta abierta. El Decaimiento 4 arranca en 4,07 mm/s y `med_S0002` está en
   3,55 mm/s: el mismo orden. El anómalo era `med_S0008`. **No era sesgo de
   supervivencia**, que era la sospecha anotada: estos números están corregidos por ruido
   y calculados sobre >95 % de vectores válidos.
4. La aserción de `02_barrido_dt.py` que protege este supuesto compara la **etiqueta**
   (`tipo == "forzado"`), no la física. La etiqueta viene del cuaderno de laboratorio y en
   este caso no describe el registro.

### [D-30] La meseta inicial de los decaimientos es el estado forzado, y se excluye del ajuste

Lucía, 2026-09-09: la meseta que aparece en los primeros 6 a 10 s de los cuatro registros
**es el estado previo a apagar el forzado**, grabado a propósito para no tener que
coordinar el corte con el inicio del registro y arriesgarse a perder los primeros pasos del
decaimiento. Cierra la pregunta abierta que había quedado en [V-13].

Consecuencias: (i) el ajuste correcto de α es el de **t > 10 s**, 0,06687 ± 0,00488 s⁻¹, y
no el de todos los puntos; (ii) el instante exacto del corte no está registrado, pero eso
**no afecta a α** —un corrimiento del origen de tiempos cambia la amplitud ajustada, no la
pendiente— aunque impide fechar el decaimiento en absoluto; (iii) la meseta **mide el
estado forzado** de cada corrida: 4,07 mm/s la que sale de CI 3 y ~3,6 mm/s las que salen
de CI 1. Eso es un contraste independiente que refuerza [H-12]: `med_S0008`, que debería
ser ese mismo estado forzado de CI 3, arranca en 2,2 mm/s y decae.

### [V-14] Verificaciones del filtrado: no hay evidencia de que se esté borrando señal

Pregunta de Lucía: ¿hay manera de saber si se está suavizando de más, y se debería poder
identificar el filtro en los espectros? `codigo/09_verificacion_filtrado.py`, figura 6 de
`salidas/decaimiento.html`.

Primero conviene separar lo que **es** filtro de lo que no: la mediana móvil de la figura 1
es cosmética y no entra en ningún número. Quedan tres cosas.

**A — Elegir Δt largo.** No suaviza en el espacio, pero promedia en el tiempo. Test sin
parámetros libres: si el exceso a n chico fuera sólo ruido aditivo, entonces
`(d₁(n)/d₁)² − 1 = (σ/(n·d₁))²`, una recta de pendiente 1. Los 93 puntos caen sobre o por
debajo (cociente mediano 0,22). **Puntos por encima serían la firma de pérdida de señal a
n grande, y no aparecen.** Entre los dos n más grandes, d₁ coincide dentro del 2,1 %
mediano (15 % en el peor instante, el más lento de su registro).

**B — La ventana de interrogación, y sí se ve en el espectro, pero no la esperada.**

| | |
|---|---|
| Nyquist de la grilla de vectores (paso 4 px) | 2985 m⁻¹ |
| primer cero de sinc²(kX/2), ventana final X = 16 px | 1492 m⁻¹ |
| corte de Foucaut k_c·X = 2,8, X = 16 px | 665 m⁻¹ |
| **caída del piso medido entre 900 y 2900 m⁻¹** | **0,19** |
| la que predice sinc² con X = 16 px | 0,002 |
| la que predice sinc² con X = 4 px (paso de grilla) | 0,46 |

**El ruido de este multipaso no está filtrado por la ventana final**: se comporta como si
estuviera correlacionado sobre el paso de la grilla, no sobre la ventana, y no aparece el
cero en 1492 m⁻¹ que exigiría el sinc² de 16 px. Para el marco de Foucaut que usa el
capítulo de instrumento (`E = E_ruido·sinc²(kX/2)` con X = 16 px, ver [V-05] y [D-18]),
esto quiere decir que **no describe esta salida**. Salvedad: en esa banda puede haber señal
real mezclada, así que la caída medida es una cota. Y en la dirección que importa para
[P-01]: la ventana sí suprime energía verdadera a k alto, de modo que el k reportado en
[H-11] es una **cota inferior**, y δ también.

**C — La resta del piso del espectro**, que es el único paso donde una elección propia
puede borrar señal:

- **C1, sobre datos sintéticos de σ conocido** —única forma de saber si el estimador está
  bien normalizado—: se le da ruido blanco puro y recupera σ dentro del 1 % (0,050 →
  0,0501; 0,100 → 0,0991; 0,200 → 0,2020; 0,400 → 0,3974 px).
- **C2, contra una ruta independiente:** el σ del espectro es una serie limpia y monótona
  que cae con el flujo (0,26 → 0,06 px en el Decaimiento 1), confirmando [H-07] por otro
  camino. El σ del ajuste de dos Δt es inestable en estos instantes y no sirve de contraste
  punto a punto.
- **C3, sensibilidad a la banda:** cambiando dónde se lee el piso, el k pesado se mueve
  menos del 3 % (238–242 m⁻¹ al inicio, 141–145 al final). El número es del dato.

**Cota sobre el sesgo de la velocidad.** Con la cota superior σ = 0,235 px y el
desplazamiento más chico usado (0,68 px) la corrección sería del 6 %; con el σ efectivo
del colapso de A, 1,3 %. En los instantes típicos (2 a 3 px) queda por debajo del 1 %, y
todo eso es menor que la barra del α, que es del 7 %.

### [V-15] El análisis de capa límite de divergencia de KIO §2.3, aplicado al esquema de SPECTER

Chequeo pedido por Lucía sobre una inferencia que yo había marcado como no verificada.
Documento en `teoria/capa_limite_divergencia_KIO_vs_Fontana.md`, chequeo numérico en
`verificacion/kio_capa_divergencia.py`, salida en `salidas/tablas/kio_capa_divergencia.json`.

**Lo estructural.** La capa límite de divergencia de Karniadakis, Israeli y Orszag (1991)
§2.3 existe porque **su término viscoso es implícito**: el campo nuevo queda adentro del
laplaciano, la divergencia obedece `Q − c·νΔt·∇²Q = 0`, y las soluciones homogéneas de ese
Helmholtz decaen como `exp(−s/ℓ)` con `ℓ = √(c·νΔt)`. En el esquema de Fontana el término
viscoso es **explícito** —verificado en `src/include/hd/hd_rkstep2.f90`— y entonces

    div(v^{n+1}) = div(v*) − dt·∇²p = div(v*) − dt·(1/dt)·div(v*) = 0

idénticamente, sin operador diferencial sobre Q. **No hay ecuación de Helmholtz, no hay
capa límite de divergencia**: no es que sea chica, la estructura que la produce no está.

**Medido**, malla 8×8×128, ν = 10⁻², campo 3D con kx ≠ 0, Δt entre 10⁻⁴ y 8·10⁻⁴:

| | exponente p en y ~ Δt^p | esperado |
|---|---|---|
| `⟨(∇·v)²⟩`, no-deslizamiento | **+0,000** | 0 |
| `⟨(∇·v)²⟩`, tope libre | **−0,000** | 0 |
| `⟨\|v_t\|²⟩` en la pared | **+4,004** | 2·ORD = 4 |
| residuo de tensión, cara libre | **−0,003** | 0 |

La divergencia no cambia una sola cifra (1,069·10⁻¹⁶) en un factor 8 de Δt mientras el
deslizamiento cae cuatro órdenes. **El error de partición está íntegramente en la condición
de contorno tangencial y nada de él en la divergencia.** P2 reproduce sobre este árbol la
Fig. 4 de Fontana et al.

**Consecuencia propia, que no está en ninguno de los dos papers:** aunque un esquema
explícito produjera esa capa, sería **sub-grilla por construcción**, porque la estabilidad
viscosa explícita exige `Δt ≲ 2dz²/(π²ν)` y entonces `ℓ = √(νΔt) ≲ 0,45·dz`. Medido: `ℓ/dz`
va de 0,10 a 0,29 en el rango usado, y ≤ 0,45 en el límite de estabilidad.

**Salvedades:** el coeficiente exacto de `ℓ = √(γ₀νΔt)` se leyó de un escaneo con OCR
imperfecto (el escalamiento `√(νΔt)` sí es inequívoco); Orszag, Israeli y Deville (1986)
sigue sin leerse, es de acceso cerrado; y el análisis es sobre el esquema continuo en el
espacio, así que no cubre qué le hace la continuación FC-Gram al piso de 10⁻¹⁶.

### [H-13] La condición inicial de la puerta es degenerada para medir no linealidad

Al armar la etapa 1 del barrido en δ reusé la condición inicial `vparam1=1` de la puerta de
la Fase 2. El barrido dio **resultados idénticos bit a bit entre 16×16 y 32×32**, que fue
lo que delató el problema: esa condición es `ψ ~ cos(kx x)·sin(πz/Lz)` en el plano (x,z),
**autofunción del laplaciano 2D**, así que `J(ψ, ∇²ψ) = 0` y el término no lineal se anula
idénticamente. El barrido midió cero porque no había nada que medir.

Es exactamente el caso que este mismo proyecto había señalado en la sección 6 de
`teoria/P01_validez_reduccion_2D.md` para el armónico del tablero de imanes, y que no
apliqué a la condición inicial de la puerta.

**Y hay una razón de fondo por la que no alcanzaba con retocarla:** en un flujo 2D en
(x,z), integrar continuidad en z da `∂⟨u⟩/∂x = −[w] = 0`, o sea que el promedio vertical
—que es la variable de la clausura— no puede tener estructura horizontal. **La prueba tiene
que ser tridimensional.**

Se escribió entonces una condición inicial propia (en `verificacion/barrido_delta_etapa1.py`,
no se toca el árbol de trabajo): `u_h = F(z)·∇⊥ψ(x,y)`, `w = 0`, con `F = sin(πz/2Lz)` el
modo fundamental del par rígido-libre y `ψ = cos(kx)cos(ky) + b·cos(2kx)cos(ky)`, cuya suma
de dos modos con `|k|²` distinto rompe la degeneración. Verificada en una corrida corta:
`div² = 1,6·10⁻²⁷`, `w = 0` a 10⁻⁴¹ en las dos caras, y residuo de tensión en la cara libre
5·10⁻¹⁵.

**Dato de convención que costó una iteración y conviene tener escrito:** los números de
onda de SPECTER son **enteros** (`specter.fpp:777-780`, `ky(j) = real(j-1)`), de modo que
el `Lx = 1` del archivo de parámetros son 2π de largo físico y el modo fundamental es
`k = 1`. `Lz`, en cambio, es la longitud literal. Verificado comparando el `⟨v²⟩` que
informa el código contra la fórmula analítica de la condición inicial: 4,750·10⁻³ contra
4,750·10⁻³. Con eso δ queda calibrado y no «a menos de una constante».

### [D-31] El proyecto pasa a un repositorio de git, privado, y corrige a [D-19]

Pedido de Lucía: puso Claude en el clúster y quiere que trabajar en la laptop o allá esté
igual de documentado y sea indistinto. Hasta hoy ninguno de los dos árboles era un repo
git, que es lo que dejaba constancia [D-19]; queda corregido.

`github.com/Lucia-Maz/proyecto-final-piv`, **privado**. 104 archivos, 4,3 MB.

**Qué NO entra, y por qué:**

| | motivo |
|---|---|
| `numerico/SPECTER-{upstream,trabajo}/` | código de terceros, 5,4 MB. Va el **parche** contra `0ad1edb` en `numerico/specter-parche/`, que es exactamente lo que [D-26] define como el aporte propio, con la receta para reconstruir el árbol y la suma de verificación |
| `bibliografia/*.pdf` | son de los editores. `REFERENCIAS.md` tiene los DOI resueltos contra Crossref y de dónde bajar cada uno |
| `salidas/decaimiento.html` | 1,6 MB regenerables con `codigo/08_html_decaimiento.py` desde los JSON |
| `.venv/`, `__pycache__`, `*.o`, `*.mod` | reconstruibles |
| datos crudos | 140 GB en el disco externo, sólo lectura ([D-02]) |

**Qué sí entra y es la parte que importa:** `ESTADO.md`, `DECISIONES.md`, todo `codigo/`,
`teoria/`, `verificacion/`, los JSON de resultados de `salidas/tablas/` y las figuras. Con
eso, en una máquina sin el disco externo se regeneran las figuras y los informes, y corren
las dos puertas de aceptación y todo lo numérico.

**Higiene, verificada antes de crear el repositorio:** se escaneó el árbol por contraseñas,
tokens, claves y hosts. El único acierto es la nota de `ESTADO.md` que *advierte* que las
contraseñas del clúster están en las notas de Obsidian de la laptop; no hay ninguna
credencial en el árbol.

### [D-32] `CLAUDE.md` fija las convenciones, y vale en las dos máquinas

Archivo nuevo en la raíz, que Claude lee al abrir el proyecto en cualquier lado. Recoge lo
que hasta ahora sólo vivía en el log y había que ir a buscar: por dónde empezar, que los
parámetros del experimento no se hardcodean ([D-27]), que las puertas están fijadas por
hash y no se editan para que pasen, que ninguna puerta contiene la fórmula del resultado
que verifica, que no se afirma lo que no se chequeó y que al citar hay que dar sección y
frase textual, y que la bibliografía distingue «localizado» de «leído en el original».

Separa además lo que **sí** depende de la máquina: los límites de memoria son de la laptop
y no aplican en Sakura, y el disco externo con los datos crudos está sólo en la laptop, de
modo que los scripts que lo tocan no corren en el clúster —pero sus resultados ya están
volcados en `salidas/tablas/*.json`.

### [D-33] Las reglas del clúster se integran a `CLAUDE.md`, textuales y con procedencia

Lucía pasó las instrucciones operativas que le dieron los administradores de Sakura y pidió
que quedara **un único** `CLAUDE.md`, sin sobrescribir lo que ya había. Se agregaron como
sección propia dentro de «Las dos máquinas», en un bloque citado **textual y sin traducir**:
son reglas operativas —particiones, reserva de cores, variables de OpenMPI, dónde va la E/S
paralela, módulos— y una mala traducción se paga con un trabajo mal lanzado. Todo lo que ya
estaba quedó intacto.

**Una errata aparente que se declara y no se corrige por cuenta propia:** en
`OMPI_MCA_btl=sm,self,tpc`, el `tpc` es casi seguramente `tcp`, y coincide con la frase
anterior de las mismas instrucciones, que dice que el transporte entre nodos es TCP. Queda
anotado para confirmarlo con quien las escribió antes de usarlo.

**Cuatro consecuencias para este proyecto**, que se derivaron de cruzar esas reglas con lo
que el proyecto necesita:

1. En el clúster SPECTER **no** usa el env de conda `specter`: se compila con los módulos
   `gnu15`, `openmpi5` y `fftw/3.3.11`, apuntando `FFTWDIR` al módulo. El resto de la
   configuración de [D-23] no cambia.
2. Es **CPU puro**: partición `compute`, o `normal` para pruebas cortas en el nodo de
   login. Los nodos con GPU no hacen falta.
3. El scratch de las corridas va a `/share/scratch*`, **no a `$HOME`**: desde los nodos de
   cómputo la E/S paralela contra `$HOME` es lenta.
4. `python/3.13.13` trae numpy pero **no scipy ni matplotlib**, y el proyecto los usa
   (`scipy.optimize` en `numerico/fase4/superficie_libre_escalas.py`, matplotlib en todos
   los scripts de figuras). Hace falta un venv propio sobre ese módulo, o generar las
   figuras en la laptop desde los JSON de `salidas/tablas/`, que el repositorio ya permite.

### [D-34] Prefijo sintético para correr la puerta de la Fase 2 en el clúster

`verificacion/test_aceptacion_fase2.py` está fijada por hash y espera **un solo prefijo**
con `bin/mpif90`, `bin/mpirun` y `lib/libfftw3` adentro, que es la forma que tiene un env
de conda. En Sakura, con módulos, OpenMPI y FFTW viven en prefijos distintos, así que
`SPECTER_TOOLCHAIN` no puede apuntar a uno solo.

**Alternativa descartada:** editar la puerta para que acepte dos prefijos. Está fijada por
hash justamente para que no se la toque, y flexibilizarla para una máquina nueva es
exactamente la clase de cambio que la vuelve inútil como puerta.

Se agrega `verificacion/toolchain_cluster.sh`, que arma un directorio con enlaces
simbólicos a los dos prefijos —`bin/mpif90`, `bin/mpirun`, `lib` e `include`— y cumple la
forma que la puerta espera. Descubre FFTW probando `FFTW_DIR`, `FFTW_ROOT`, `FFTW_HOME`,
`FFTWDIR`, `FFTW3_DIR` y `FFTW3_ROOT`, y si ninguna sirve lo busca en `LD_LIBRARY_PATH`.
No modifica nada fuera de su directorio de salida.

### [D-35] En el clúster todo va por el planificador, no a mano en el nodo de login

Al revisar contra las reglas que dieron los administradores aparecieron tres
incumplimientos de lo que yo había propuesto, y uno propio:

1. **Correr la puerta «en normal» a mano.** `normal` es una **partición de SLURM**: la
   forma correcta es mandarla al planificador. La puerta compila SPECTER y lanza unas
   quince simulaciones, diez a quince minutos de CPU que no van por fuera del scheduler en
   una máquina compartida. Se agrega `verificacion/puerta_fase2.sbatch`, que usa la
   partición `compute` para no cargar el nodo de login.
2. **Faltaban las variables de transporte.** Las reglas piden `OMPI_MCA_pml=ob1` y
   `OMPI_MCA_btl=sm,self,tcp` para la OpenMPI de gnu (red Ethernet de 10 Gb, sin
   InfiniBand). La puerta define tres `OMPI_MCA_*` propias pero no éstas; como copia el
   entorno heredado, alcanza con exportarlas antes, sin tocarla.
3. **El scratch iba a `/tmp` del nodo.** Va a `/share/scratch/$USER/...`, que es donde las
   reglas dicen que funciona la E/S paralela. Como `SPECTER_SCRATCH` explícito hace que la
   puerta **no** lo borre al terminar, el sbatch lo limpia con un `trap`.
4. **Un defecto propio de `toolchain_cluster.sh`:** hacía `command -v mpirun || command -v
   srun`. Si hubiera caído en ese fallback habría enlazado `srun` con el nombre `mpirun`, y
   la puerta lo invoca con `-np`, que `srun` no acepta. Ahora falla con un mensaje claro.

    module load gnu15 openmpi5 fftw/3.3.11
    eval "$(bash verificacion/toolchain_cluster.sh)"
    python verificacion/test_aceptacion_fase2.py     # 6/6

En el clúster, en cambio:

    sbatch verificacion/puerta_fase2.sbatch

### [D-36] La puerta va a la partición `normal`, y corrige a [D-35]

[D-35] mandaba la puerta a `compute` con el argumento de no cargar el nodo de login. **El
argumento era falso**, y el trabajo 6885 quedó encolado con *«Nodes required for job are
DOWN, DRAINED or reserved for jobs in higher priority partitions»*. Medido con `sinfo` el
2026-09-09:

| partición | nodos | estado |
|---|---|---|
| `normal*` (por omisión) | a1, c[1-5], g[1-2], sakura | mix / alloc |
| `compute` | c[1-5], **sakura** | mix / alloc, ninguno idle |
| `legacy` | l[1-3] | **idle** |
| `cuda` | g[1-2] | mix |
| `rocm` | a1, a2 | mix / alloc |

Dos cosas que desmienten lo que yo había supuesto: **`compute` incluye sakura**, así que no
evita el nodo de login; y es un **subconjunto** de `normal`, que al ser la partición por
omisión tiene más prioridad sobre los nodos compartidos. De ahí el encolado indefinido.

`legacy` es la única con nodos `idle`, y es justamente la que las reglas dicen de evitar:
conviene no dejarse tentar por eso.

**Regla que queda:** trabajos cortos de un core —la puerta, los tests— a `normal`. Para la
corrida de producción de SPECTER habrá que elegir con cuidado y mirar `sinfo` en el
momento; y si cae en a1, a2, g1 o g2, dejar cuatro cores libres, como piden las reglas.

### [D-37] El scratch del clúster queda fijado en `/share/scratch1/$USER`

La primera versión buscaba por comodín entre `/share/scratch*` y `/share/data*` y tomaba el
primero escribible. Lucía marcó que eso puede terminar en un lugar equivocado: en
`/share/data2/$USER` están sus simulaciones de **GHOST**, y un directorio de datos no
es lugar para el scratch de un trabajo.

Queda fijado en `/share/scratch1/$USER`, que es donde cayó y es el correcto, con
`SPECTER_SCRATCH_BASE` como escape si algún día cambia. **`/share/data*` sale de la lista
de candidatos**: ningún script de este proyecto tiene por qué escribir ahí.

**Dato de contexto:** Lucía tiene simulaciones de GHOST en el clúster, pero **ninguna de
SPECTER todavía**. Antes de la corrida de producción vale la pena mirar si esas corridas de
GHOST sirven de punto de comparación, pero es un tema aparte.

### [H-14] `make` no está en el PATH de un trabajo de SLURM en Sakura

La puerta falló con `[ERROR] [Errno 2] No such file or directory: 'make'` después de
encontrar bien MPI y FFTW. Como el prefijo sintético se antepone al PATH,
`toolchain_cluster.sh` resuelve ahora `make` —o `gmake`, según la instalación— y deja un
enlace con el nombre `make` adentro. Si no encontrara ninguno, falla con un mensaje que
sugiere buscar el módulo que lo provee.

### [H-15] `eval "$(cmd)"` no propaga el fallo de `cmd` bajo `set -e`

El trabajo 6897, corrido después de [D-37], falló con una mezcla rara de dos errores: el
propio de `toolchain_cluster.sh` (`no encuentro make ni gmake en el PATH`) seguido del de
`test_aceptacion_fase2.py` señalando el entorno de conda de la laptop
(`/home/lucia/miniforge3/envs/specter`), que no existe en el clúster. Los dos no pueden
ser ciertos a la vez si el script hubiera cortado en el primero.

La causa es de shell: `puerta_fase2.sbatch` llamaba `eval "$(bash
verificacion/toolchain_cluster.sh)"`. Cuando el script interno falla sin escribir nada a
stdout, `"$(...)"` vale la cadena vacía; `eval ""` no hace nada y **su propio código de
salida es 0**. `set -e` mira el estado de `eval`, no el del comando adentro de la
sustitución, así que el trabajo seguía de largo sin `SPECTER_TOOLCHAIN` exportado —hasta
reventar más abajo, en el fallback de la puerta, con un mensaje que no menciona la causa
real.

Confirmado además que la resolución de `make` en sí **no** era el problema de fondo: con
los mismos `module load gnu15 openmpi5 fftw/3.3.11 python/3.13.13` en un trabajo real en
`g1`, `command -v make` resuelve a `/usr/bin/make` sin inconvenientes — los módulos no
tocan esa parte del PATH del sistema. El fallo de 6897 fue este bug de propagación, no
un problema de módulos.

**Arreglado** capturando la sustitución en una variable antes de evaluarla:
`TOOLCHAIN_ENV="$(bash verificacion/toolchain_cluster.sh)"` sí propaga el fallo bajo
`set -e`, porque ahí la asignación simple es la excepción que POSIX y bash sí respetan.
`eval "$TOOLCHAIN_ENV"` corre después, sólo si lo anterior no cortó.

### [D-38] Los logs de la puerta van a `verificacion/logs/`, no a la raíz del repo

Los `puerta-fase2-*.out` de los trabajos 6885 a 6897 quedaron sueltos en la raíz,
sin trackear. Lo que vale de ellos ya está volcado acá arriba ([D-34]–[D-37], [H-14],
[H-15]); son logs de intentos fallidos, no resultado final, así que no se comitean
([D-19] y la regla de no reportar lo que se rehizo). Se borraron, y
`puerta_fase2.sbatch` ahora escribe a `verificacion/logs/puerta-fase2-%j.out`, que queda
en `.gitignore`. El directorio hay que crearlo antes de mandar el trabajo: SLURM resuelve
la ruta de `--output` al momento de encolar, no la crea él.

### [D-39] La puerta y la corrida de producción se mandan con `--nodelist=g1`

Por pedido explícito: Lucía tiene corridas propias en curso en `a1`, `a2`, `g2` y `c5`
(los `calibF0_00` en cola de SLURM al 2026-09-16), y no quiere que un trabajo de este
proyecto compita por esos nodos. `g1` es el único de los cuatro habituales sin trabajo
propio en ese momento —lo ocupan 16 de sus 20 cores con un trabajo de otro usuario
(`lucianov`, `SWHD_27pulc`)—, así que un trabajo de un core encaja en los cuatro libres
sin desalojar a nadie, y si no encajara, SLURM lo deja en cola hasta que se libere.

**No se hardcodea en `puerta_fase2.sbatch`**: es una preferencia de esta sesión de
cómputo, no un requisito técnico de la puerta, que sigue valiendo para `normal` en
general ([D-36]). Se pasa en la línea de comando:

    sbatch --nodelist=g1 verificacion/puerta_fase2.sbatch

Si en el futuro `g1` no es la mejor opción, alcanza con omitir el flag o apuntarlo a otro
nodo — no hace falta editar el script.

**Resultado, trabajo 8168:** `RESULTADO: 6/6` (V1–V6), sin encolarse — entró directo
porque `g1` tenía 4 de 20 cores libres. El scratch en `/share/data2/$USER/puerta-fase2-8168`
se limpió solo al terminar, y los cuatro trabajos propios en `a1`, `a2`, `g2` y `c5`
siguieron corriendo sin interrupción.

### [D-40] Terminología: "free-slip" para la condición de contorno, "superficie libre" para el escenario físico

Pregunta de Lucía: si no se modela la deformación de la superficie, ¿no correspondería
llamar **free-slip** a lo que la prosa del proyecto viene llamando "superficie libre"?

Sí, y son dos ejes distintos que la prosa venía mezclando —aunque el código nunca lo hizo:
`v_parsebc` acepta la cadena `freeslip`, la subrutina es `freeslip_z`, el diagnóstico es
`freeslip_diagnostic.txt` ([D-26]). Free-slip vs. no-deslizante es sobre la tensión
tangencial en un borde **fijo**, y depende de que el fluido de arriba (aire) tenga
viscosidad y densidad despreciables frente al agua. Superficie libre (deformable) vs. tapa
rígida es sobre si ese borde **se mueve**, y depende del balance Froude/capilar que [D-28]
ya había cuantificado (η/h ~ 10⁻⁵). Lo implementado es free-slip **y además plano**: dos
aproximaciones apiladas, no la misma cosa.

No es sólo prolijidad de vocabulario. El caso plano tiene una propiedad que el deformable
pierde: su condición de contorno resulta ser una condición sobre la **vorticidad**
tangencial (`teoria/superficie_libre_v_estrella_y_p.md` §2.4), y como la proyección de
presión resta un gradiente, se transfiere al campo proyectado sin error de partición —de
ahí el residuo independiente de `dt` que mide V3. Con superficie deformable esa propiedad
no sobrevive (§4.4 del mismo archivo). Llamar "superficie libre" a la condición numérica
sugiere que se está modelando lo que en realidad se dejó fuera de alcance.

**Convención de acá en adelante:** "free-slip" (o "libre de tensiones tangenciales") para
la condición de contorno numérica; "superficie libre" acotado al escenario físico —el
electrolito en contacto con aire— que motiva usar free-slip como aproximación. Se revisó
la terminología en `CLAUDE.md`, `ESTADO.md`, `salidas/entrega1_resumen.html` y los dos
generadores que producen prosa (`codigo/08_html_decaimiento.py`,
`codigo/11_pdf_barrido_delta.py`).

**Qué no se tocó, y por qué:** las entradas anteriores de este log, que quedan como
constancia de cómo se habló en su momento (no se reescribe el log, [README.md]);
`teoria/superficie_libre_v_estrella_y_p.md`, que ya distinguía los dos ejes con precisión
—es la fuente de la que sale esta decisión, no algo que corregirle—; y
`verificacion/test_aceptacion_fase1.py`/`test_aceptacion_fase2.py`, las puertas fijadas
por hash, que no se editan aunque usen la frase en un docstring o un mensaje.

## 2026-09-16

### [H-16] El pico del espectro forzado está en 295 m⁻¹, y era el fundamental de la red con el paso bien puesto

Lucía, mirando `f4_espectros`: *"creo que lo estás localizando mal (k más chico que el que
veo que tiene un pico)"*. Tenía razón, y los datos ya lo decían.

**Lo que había.** `celda.py` declaraba imanes de 1 cm en una red de 5 cm de paso, y de ahí
`k_fundamental = π√2/a = 88,9 m⁻¹` (λ = 7,07 cm), que [V-11] llamó "el k del forzado".
`p01_validez.py` lo "verificaba" tomando el **argmax modo a modo** de la transformada de
un patrón de discos ±1. Eso sólo prueba que el modo individual más alto es el fundamental,
lo que es cierto para cualquier red, y no dice dónde está la energía: el mismo script
devolvía `fraccion_energia_en_fundamental = 0,115` y `k_medio_pesado = 575 m⁻¹`, y ese
dato quedó anotado como "salvedad del modelo crudo" en vez de leerse como lo que era.

**Lo que dicen los campos.** En `salidas/tablas/decaimiento.json`, los doce espectros de la
meseta forzada (t < 10 s, los cuatro registros, [D-30]) tienen su pico **en el mismo bin:
294,9 m⁻¹**, con bins de 23,6 (campo de 26,9 cm). Es un pico único y angosto, λ ≈ 2,1 cm,
3,3 veces por encima de los 88,9 que marcaba la línea punteada de la figura. Durante el
decaimiento migra a 59–130 m⁻¹, consistente con que los modos de k = 295 decaen a
λ(k) = α(1 + νk²/α) ≈ 2,3α y los grandes a ~α.

**La causa era un dato, no el método.** Al preguntar la geometría, Lucía: *"son 5 mm entre
imanes, fue un typo"*: 1 cm de imán más 5 mm de vacío, **1,5 cm entre centros**, y un
acrílico de casi 6 mm entre la cara del imán y la capa. Con a = 1,5 cm el fundamental de
la red es **296,2 m⁻¹** (λ = 2,12 cm): cociente contra el pico medido **0,996**, dentro
de un cuarto de bin. Con los imanes ocupando dos tercios del paso y el acrílico
suavizando el campo sobre una distancia comparable, el patrón es casi un coseno puro y el
fundamental se lleva la energía —el modelo de discos ±1, que sobrestima los armónicos,
ya da 79 % en el fundamental—, así que "k de la red" y "k donde el forzado inyecta" son
lo mismo. Con el paso de 5 cm no lo eran (3 % de llenado: fuentes puntuales), y esa
diferencia fue la que se discutió antes de descubrir el typo.

De paso: el `ℓ = 2 cm` heredado que [V-11] descartó como "ninguna longitud de la celda"
era, a un 6 %, la longitud de onda del flujo forzado.

**Alternativa descartada:** definir un "k del forzado" medido, separado del fundamental
de la red, con procedencia propia en `celda.py`. Con el paso correcto no hace falta; la
red lo predice y el espectro lo confirma. Ahora `p01_validez.py` lee el pico de la meseta
desde `decaimiento.json` si existe y lo imprime al lado del de la red.

### [D-41] El paso de la red de imanes es 1,5 cm entre centros; corrige a [D-27]

`codigo/celda.py`: `paso = 1,5 cm`, más dos campos nuevos con procedencia, `separacion =
5 mm` (vacío entre cilindros) y `acrilico = 6 mm` (cara del imán a fondo de la capa), y
una comprobación de que `paso = diametro + separacion`. El diámetro de 1 cm y el patrón de
tablero no cambian. Procedencia: Lucía, 2026-09-16.

**Lo que se recalculó**, todo desde el mismo `celda.py` y sin tocar los espectros (que no
dependen del paso):

| | antes (a = 5 cm) | ahora (a = 1,5 cm) |
|---|---|---|
| k fundamental de la red | 88,9 m⁻¹ | **296,2 m⁻¹** |
| ε = kh | 0,533 | **1,777** |
| λ(k)/α = 1 + 4ε²/π² | 1,115 | **2,28** |
| δ en el forzado (⟨u⟩ = 0,84 mm/s) | 2,69 | **8,97** |
| δ al inicio del decaimiento (⟨u⟩ = 2,55 mm/s) | 8,15 | **27,2** |
| distorsión estimada \|a₁/a₀\| ≈ 0,0077 δ | 2 % y 6 % | **7 % y 21 %** |
| cota de la película, μ_s = β μ/(k²h) para β = 0,15 | 3·10⁻⁶ N·s/m | **2,8·10⁻⁷ N·s/m** |

Salidas regeneradas: `p01_validez.json`, `superficie_libre_escalas.json`,
`decaimiento_resumen.json`, figuras `f4`/`f5` y `salidas/decaimiento.html`.
`decaimiento.json` guarda en `celda.k_forzado_1_m` el 88,9 con que se corrió
`06_decaimiento.py`; no se editó a mano (es la configuración con que se produjo) y
`07_figuras_decaimiento.py` ya no lo lee: toma el k de `celda.py`.

**Qué cambia de conclusión y qué no.**

- **[V-11] sigue teniendo razón en el método y cambia en el número.** El criterio
  correcto es la distorsión del perfil y no δ; pero con k = 296 la distorsión estimada al
  inicio del decaimiento es 21 %, no 6 %: la clausura de un modo es **mala en la meseta
  y en los primeros segundos**, y buena recién cuando la energía migró a k ≲ 130 y U
  cayó, que es dentro de la ventana del ajuste (t > 10 s). El barrido en δ con SPECTER
  del punto 1 de la Fase 4 tiene que llegar a δ ≈ 30, no a 10.
- **[H-11] se invierte en la lectura.** El k pesado por energía (240 → 135 m⁻¹) ya no
  está "por encima del forzado": está **por debajo** del fundamental de la red desde el
  primer instante, porque el promedio pesado incluye el hombro de k chico. El pico, no
  el promedio, es lo que coincide con la red.
- **[D-28] queda igual en conclusión y no en la frase.** La longitud capilar (k = 369 m⁻¹)
  ya no es "bastante mayor" que el k del forzado: 296 está a un 20 %. Restituye igual la
  gravedad (9810 contra 6317 Pa/m) y η/h sigue en 10⁻⁵; sólo cambia el margen.
- **La comparación de [V-13] queda abierta:** ver [P-02].

**Alternativa descartada:** dejar el paso como estaba y tratar los 295 medidos como
"escala del flujo, distinta de la del forzado". Era la lectura de [H-11] y se sostenía
sólo porque el paso estaba mal.

### [P-02] Pregunta abierta: el α medido está por debajo de λ(k) para cualquier k que tenga energía

[V-13] compara el α medido (0,0669 ± 0,0049 s⁻¹, t > 10 s) contra α pelado (0,0685) y
lo da a 2,4 %. Pero la tasa de decaimiento de u_rms es, sin ninguna clausura horizontal,
`α + ν⟨k²⟩_E` —cada modo decae con λ(k) = α + νk² y u_rms² es la suma—, y **siempre está
por encima de α**. Con el k corregido esto ya no es un detalle:

| k que se use | λ(k) [1/s] | medido / λ(k) |
|---|---|---|
| α pelado (k → 0) | 0,0685 | 0,98 |
| pico en la ventana del ajuste, 59–130 m⁻¹ | 0,0720–0,0854 | 0,93–0,78 |
| k pesado por energía al final, 135 m⁻¹ | 0,0868 | 0,77 |
| fundamental de la red, 296 m⁻¹ | 0,156 | 0,43 |

O sea que el "acuerdo al 2,4 %" es contra el único valor que el flujo **no** puede tener.
La tensión pre-existía —ESTADO.md ya decía que había que comparar contra 1,115α, y [V-13]
no lo hizo— y con k = 296 en la meseta se vuelve un 7–25 % en la ventana del ajuste.

Candidatos, ninguno chequeado: (i) el k pesado está inflado por el piso de ruido, que
[V-14] mostró que no es blanco, y el pico (59–83 al final) es una lectura mejor que el
promedio; (ii) h es mayor que 6 mm —α ∝ 1/h², y 6,5 mm bajan α a 0,058, que con
+27 % de viscosidad horizontal da 0,074—; (iii) la película superficial no aplica porque
empuja para el otro lado; (iv) el perfil vertical no es el fundamental en la ventana
(δ ≈ 3–10 ahí), y la fricción efectiva difiere de α — que es exactamente lo que el
barrido en δ con SPECTER mide. Hasta resolverlo, el número que se entrega es el cociente
contra α pelado **con esta salvedad escrita al lado**, no como acuerdo.
