# Marco teórico

**Estado: borrador del 2026-09-02.** Todas las fórmulas y valores que se atribuyen a un
paper fueron verificados contra el PDF que está en `bibliografia/`, no citados de
memoria. Donde el original no está disponible o no fue leído, se dice explícitamente.

---

## 1. El régimen del problema

PIV mide un **desplazamiento** entre dos exposiciones y lo divide por el intervalo entre
ellas. Todo lo que sigue depende de una sola relación: cuánto se mueve una partícula
entre cuadros, comparado con la precisión con la que el algoritmo puede localizar el
pico de correlación.

En la campaña del 02-06-25 esa relación es adversa. El desplazamiento entre cuadros
consecutivos es de **0,09 a 0,28 px**, y el piso de ruido típico de PIV es de orden
0,1 px. La relación señal/ruido *en desplazamiento* es del orden de 1, y en buena parte
de un decaimiento es menor que 1.

Ese es el régimen que este trabajo estudia. No es el régimen habitual de la literatura
de PIV, que trabaja con desplazamientos de varios píxeles y donde el ruido es una
corrección; acá es comparable a la señal. Dos consecuencias inmediatas:

- La velocidad medida está **sobrestimada**, porque un error de media cero en el
  desplazamiento no es de media cero en el módulo de la velocidad.
- La sobrestimación **crece cuando el flujo se frena**, porque el error es constante en
  píxeles y la señal no. En un decaimiento eso produce un piso aparente en la curva de
  velocidad que no es físico.

La palanca experimental es que cada corrida tiene 3072 cuadros consecutivos a 60 fps.
Los pares se rearman con cualquier separación Δt = n/60 s sin volver a medir, lo que
convierte a Δt en una variable de análisis en vez de una constante del experimento.

## 2. Qué mide PIV y de dónde sale el error

El algoritmo correlaciona dos ventanas de interrogación y ubica el máximo de la
correlación con precisión subpíxel, típicamente ajustando una gaussiana de tres puntos.
El desplazamiento medido en una ventana es

    D_med = D_verdadero + ε

con ε el error de desplazamiento. La velocidad se obtiene dividiendo por Δt y por la
escala s en px/m:

    U_med = D_med / (Δt · s) = U_verdadera + ε / (Δt · s)

**La observación que organiza todo el trabajo** es que ε está en píxeles y no depende de
Δt, mientras que la conversión a m/s sí. Un error de desplazamiento fijo se traduce en
un error de velocidad que es inversamente proporcional a Δt. Por eso el error en
velocidad se puede reducir arbitrariamente aumentando Δt — hasta que la decorrelación lo
impide.

## 3. El modelo espectral de Foucaut, Carlier & Stanislas

Referencia central del trabajo:
[10.1088/0957-0233/15/6/003](https://doi.org/10.1088/0957-0233/15/6/003).
La versión de MST 2004 es de acceso pago; lo que hay en `bibliografia/` y lo que se cita
acá es la versión de ISPIV 2003, **leída en el original**.

### 3.1 El experimento de desplazamiento nulo

Su método de partida es medir la respuesta del sistema a una entrada nula: un flujo sin
movimiento, obtenido en un acuario de vidrio, en cuarto oscuro para que el nivel de
fondo sea cercano a cero, analizado por correlación cruzada sin desplazamiento y con
ajuste gaussiano de tres puntos. Si el desplazamiento verdadero es cero, todo lo que
mide PIV es su propio ruido.

Su tabla 1, sobre 125 campos:

| ventana | desp. medio en x | desp. medio en y | σ_u | σ_v |
|---|---|---|---|---|
| 16 × 16 | −0,016 px | 0,000 px | **0,087 px** | 0,050 px |
| 32 × 32 | −0,023 px | 0,000 px | 0,055 px | 0,029 px |
| 64 × 64 | −0,026 px | 0,000 px | 0,033 px | 0,015 px |

Dos detalles del texto que importan acá y que no suelen citarse:

- Hay un **error sistemático** de orden 0,02 px en x, que crece levemente con el tamaño
  de ventana. No es sólo dispersión.
- Sobre el efecto que este trabajo también encontró: *"The displacement being zero, the
  results were weakly affected by the peak-locking effect which tends to reduce the
  noise intensity."* Es decir, ellos mismos señalan que a desplazamiento nulo el peak
  locking **deprime** el ruido medido. Su experimento de cero es, en su propia lectura,
  una cota inferior — débilmente sesgada en su montaje, pero sesgada.
- σ_u es sistemáticamente unas dos veces σ_v. Lo atribuyen a coincidencia con Raffel
  et al. para ventana de 16×16.

### 3.2 El modelo espectral

El espectro de potencia del desplazamiento medido sobre esos campos sin movimiento es un
ruido blanco multiplicado por un sinc al cuadrado (su ecuación 1):

    E_ii(k) = E_noise · sinc²(k X / 2)

El sinc es la transformada de Fourier de la función puerta que representa la ventana de
interrogación: correlacionar sobre una ventana de lado X equivale a promediar sobre ella,
y promediar es multiplicar por un sinc en el espacio de Fourier. Miden
E_noise = 17,5·10⁻³ px³ para E₁₁ y 4,4·10⁻³ px³ para E₂₂.

De ahí sale el resultado que usamos como criterio de resolución: el número de onda de
corte es **k_c · X = 2,8**, y es **universal** — no depende del tamaño de ventana, que es
lo que permite usarlo como criterio general.

Variando el otro lado de la ventana encuentran que la densidad espectral escalada es
constante (su ecuación 2):

    Y · E_ii = ζ,   con   ζ = E_noise · Y

y como σ_u es la raíz de la integral del espectro, su ecuación 3:

    σ_u = √( 4 ζ I / (X Y) ),   con   I = ∫ sinc²(u) du = 1,492

(con 50 % de solape la integral se limita a π y da I = 1,418). De acá sale la predicción
estructural: **σ_u es inversamente proporcional a la raíz del área de la ventana de
interrogación.**

### 3.3 Lo que este paper deja pendiente, y que es justo nuestro término

En la sección de caracterización del ruido, el texto dice: *"This should be sensitive at
least to the effect of the background, speckle and camera noise together with the effect
of the two pulse shape difference. **The noise due to particle motion will be studied
afterwards.**"*

En la versión de ISPIV 2003 que tenemos, ese "afterwards" no llega: el trabajo tiene seis
secciones y ninguna desarrolla el ruido debido al movimiento de las partículas.
Presumiblemente está en la versión larga de MST 2004, que es de pago y no fue leída.

**Esto ubica exactamente nuestra contribución.** El término que Foucaut difiere es el que
en nuestros datos resulta dominante (§ 7).

## 4. La descomposición de George & Stanislas

[arXiv 2010.10768](https://arxiv.org/abs/2010.10768), v2 del 23-12-2020. PDF en
`bibliografia/`. **Leídos el resumen, la lista de símbolos y las secciones 3.2–3.3**; el
desarrollo completo, no.

Su tesis es que cuando PIV se usa para medir turbulencia hay que tratarlo como una señal
dependiente del tiempo, y que la velocidad de salida tiene **tres contribuciones**:

1. la velocidad dependiente del tiempo;
2. un ruido de **cuantización o pixelización**;
3. un ruido que viene de que **la velocidad no es uniforme dentro del volumen de
   interrogación**.

Las varianzas de (2) y (3) dependen inversamente del número medio de partículas en el
volumen de interrogación, y las tres contribuciones están filtradas espacialmente por la
extensión finita de la ventana.

La frase que sostiene nuestro resultado principal está en su sección 3.2, sobre las
desviaciones respecto de la velocidad verdadera promediada en el volumen:

> *"any departures from the true volume-averaged instantaneous velocity must be
> interpreted as either a bias (if the mean difference is non-zero) or noise (on the
> statistics of the fluctuating quantities) on the data, **albeit noise resulting from
> the flow itself**."*

Es decir: hay una parte del "ruido" de PIV que **es producida por el flujo**, no por el
instrumento. Con pocas partículas por ventana, la salida es una aproximación cruda al
promedio volumétrico euleriano instantáneo, y la diferencia se comporta como ruido.

**Su advertencia metodológica, que hay que declarar.** En la sección 3.3 señalan que en
PIV la pixelización ocurre *antes* de la integración sobre el volumen, de modo que tanto
la pixelización como el efecto del número finito de partículas quedan filtrados
espacialmente, *"and hence difficult to distinguish"*. Es una advertencia dirigida a los
métodos espectrales. El método que usamos acá los separa de otra manera — quitando
físicamente el flujo (§ 7) — que sortea esa dificultad pero introduce las suyas.

**Un supuesto suyo que en esta celda no es gratuito.** Su desarrollo asume que las
partículas siguen exactamente al flujo, y advierten que si no lo hacen, lo único que se
puede decir es que la estadística es la de un híbrido euleriano-partícula. Con partículas
de 100 µm en una capa de 5 mm, eso es una hipótesis que habría que justificar y que hoy
no está justificada.

## 5. El modelo del barrido en Δt

Este es el modelo propio, y no está tomado de ninguno de los papers anteriores.

Para una corrida **forzada**, que es estacionaria, la velocidad verdadera no depende de
la separación entre cuadros. Tomando el valor cuadrático medio sobre el campo, y
suponiendo que ε no está correlacionado con U_verdadera:

    ⟨U_med²⟩ = ⟨U_verdadera²⟩ + σ² / (Δt · s)²

y con Δt = n / fps:

    U_med(n)² = U_verdadera² + ( σ · fps / (n · s) )²

Es una **recta** al graficar U_med² contra n⁻². La ordenada al origen da la velocidad
verdadera y la pendiente da σ. Es un método de dos o más Δt: no aparece en el catálogo de
los cuatro métodos que compara Sciacchitano et al. 2015, lo que es a la vez su novedad y
su pasivo, porque no está caracterizado por terceros.

### 5.1 Los supuestos, explícitos

El modelo se apoya en cinco hipótesis. Conviene tenerlas separadas porque hoy sabemos que
una es falsa:

| | supuesto | estado |
|---|---|---|
| A1 | σ es constante: no depende de n ni del flujo | **falso** — ver § 7 |
| A2 | ε no está correlacionado con U_verdadera | no verificado |
| A3 | el flujo está congelado en el mayor Δt usado, de modo que D_verdadero ∝ n | no verificado; a n = 128, Δt = 2,13 s |
| A4 | el flujo es estadísticamente estacionario | razonable en una corrida forzada |
| A5 | descartar vectores inválidos no sesga el estadístico | dudoso — es el sesgo de supervivencia, pregunta abierta |

## 6. Los dos extremos del rango útil de Δt

El modelo anterior da el extremo **inferior**: por debajo de cierto Δt, el ruido domina.
El extremo **superior** lo dan dos criterios de la literatura:

- **La regla del cuarto**, de Keane & Adrian 1990
  ([10.1088/0957-0233/1/11/013](https://doi.org/10.1088/0957-0233/1/11/013)): el
  desplazamiento no debe superar 1/4 del lado de la ventana de interrogación, o se pierde
  la correlación por partículas que salen de la ventana. **No leído en el original**;
  citado de literatura secundaria, y hay que verificarlo antes de usarlo como criterio.
- **La caída de la fracción de vectores válidos**, que es la manifestación medible de lo
  anterior y no requiere creerle a ninguna regla.

Ambos extremos salen de la misma figura del barrido, que es lo que lo hace un método de
selección de Δt y no sólo una medición de σ.

Del lado de la resolución espacial, el límite es el filtrado pasabajos de la ventana:
**Willert & Gharib 1991** ([10.1007/BF00190388](https://doi.org/10.1007/BF00190388)) dan
una longitud de onda de corte igual al doble del tamaño de ventana. Es la fuente que
Foucaut cita para la respuesta tipo sinc. **No leído en el original.**

El multipaso 64 → 32 → 16 px con deformación simétrica de imagen, que es lo que corren
nuestros scripts y lo que corría la GUI, se apoya en **Scarano 2002**
([10.1088/0957-0233/13/1/201](https://doi.org/10.1088/0957-0233/13/1/201)). Hace falta
para justificar que iterar con deformación no invalida el modelo de ventana única sobre
el que se apoyan Foucaut y el barrido. **No leído en el original.**

## 7. Lo que la medición dice sobre estos modelos

Tres mediciones de σ sobre la misma cámara y la misma configuración de PIV, que se
diferencian sólo en qué degradación física incluyen (`DECISIONES.md`, [V-04] y [H-07]):

| medición | σ por componente | qué agrega respecto de la anterior |
|---|---|---|
| par sintético, corrimiento conocido | ≤ 0,011 px | correlación y estimador subpíxel solos |
| fluido en reposo, par real | 0,022 px (σ_u), 0,021 px (σ_v) | ruido de lectura, iluminación, movimiento propio de partículas |
| forzado 2,15 A, ajuste del barrido | ≈ 0,12 px | no uniformidad de la velocidad dentro de la ventana, movimiento fuera del plano |

El término que agrega el flujo es √(0,120² − 0,022²) ≈ **0,118 px: el 97 % de la
varianza.** El instrumento aporta el 3 %.

**Eso refuta A1**, y con ella la lectura habitual del modelo. σ no es una constante del
instrumento: su término dominante es el ruido *producido por el flujo* de George &
Stanislas, y es justamente el término que Foucaut difiere (§ 3.3). En un decaimiento,
donde U cae un factor ~10, σ no es constante a lo largo del registro, y la corrección a
la curva no puede usar un solo número.

### 7.1 Una advertencia de convención, que corrige a [D-18]

El barrido mide `desp_rms = √⟨u² + v²⟩`, que es un σ **de módulo**. Foucaut informa σ
**por componente**. Comparar los dos directamente es un error de √2.

`DECISIONES.md` [D-18] compara los 0,169 px del ajuste contra los 0,087 px de Foucaut y
concluye "un factor 2". Convertido a la misma convención, el σ por componente del ajuste
es 0,169/√2 ≈ **0,120 px**, y la comparación correcta contra 0,087 px es un **factor
1,4**, no 2. La tabla de arriba ya está en convención por componente. La conclusión
cualitativa de [D-18] no cambia — el orden de magnitud coincide y la diferencia se
explica por ζ, que es propiedad de la cadena de registro — pero el número sí.

### 7.2 Dos contrastes contra Foucaut que salen de nuestros datos

- **La asimetría x/y no se reproduce.** Foucaut mide σ_u ≈ 2 σ_v. Nuestro registro en
  reposo da σ_u = 0,0216 px y σ_v = 0,0213 px: iguales dentro del error. Es consistente
  con que la asimetría sea propiedad de su cadena óptica —YAG pulsado, imágenes de
  partícula de 2,5 px— y no un rasgo general del método.
- **El error sistemático coincide en orden de magnitud.** Foucaut mide un sesgo de
  ~0,02 px, independiente de la ventana. Nuestra calibración con corrimientos conocidos
  ([V-04]) mide un sesgo máximo de 0,018 px en el rango 0–1,5 px. Es el contraste externo
  más limpio que tenemos hasta ahora.

## 8. Ubicación frente a la literatura

**El antecedente más cercano** es el PIV espacialmente adaptativo de **Theunissen,
Scarano & Riethmuller** ([10.1007/s00348-009-0782-7](https://doi.org/10.1007/s00348-009-0782-7)):
ventana variable en tamaño, forma y orientación según gradientes y densidad de sembrado.
La diferencia que delimita qué sería nuevo acá: **ellos optimizan la calidad local de la
correlación, no un objetivo global sobre un estadístico de turbulencia, y no tocan Δt.**
**No leído en el original.**

**El marco de incertidumbre** lo dan **Sciacchitano et al. 2015**
([10.1088/0957-0233/26/7/074004](https://doi.org/10.1088/0957-0233/26/7/074004)), que
comparan cuatro métodos de cuantificación; la revisión de **Sciacchitano 2019**
([10.1088/1361-6501/ab1db8](https://doi.org/10.1088/1361-6501/ab1db8)); y el *image
matching* de **Sciacchitano, Wieneke & Scarano 2013**
([10.1088/0957-0233/24/4/045302](https://doi.org/10.1088/0957-0233/24/4/045302)), que es
la vía independiente para contrastar σ a posteriori sobre las mismas imágenes.
**Ninguno leído en el original.**

**Dos referencias de higiene.** **Garcia 2010**
([10.1016/j.csda.2009.09.020](https://doi.org/10.1016/j.csda.2009.09.020)) es el
`smoothn` que OpenPIV aplica cuando `smoothn_each_pass=True`, activo en la familia de
sesiones `180226`, y cuya función de transferencia falta caracterizar ([D-07]).
**Mendes et al. 2020** ([10.1016/j.softx.2020.100537](https://doi.org/10.1016/j.softx.2020.100537))
es un generador de imágenes sintéticas, **revisado y descartado**: es MATLAB y sólo
admite flujos analíticos, no un campo arbitrario.

**Pendiente de revisar:** **Manovski et al. 2025**
([10.1007/s00348-025-04059-0](https://doi.org/10.1007/s00348-025-04059-0)), sobre la
respuesta espectral de PIV resuelto en tiempo. Es el trabajo más reciente de la lista y
hace falta leerlo para saber si el hueco que identificamos sigue abierto.

## 9. La herramienta que el resultado vuelve necesaria

Si σ depende del flujo (§ 7), el barrido en Δt deja de ser suficiente: para caracterizar
σ a lo largo de un decaimiento hace falta un método que dé σ **a Δt fijo**, porque en un
decaimiento no se puede barrer Δt sin que el flujo cambie mientras se barre.

El método espectral de Foucaut (§ 3.2) es exactamente eso: separa ruido de señal dentro
del espectro de un solo campo, a un Δt. Aplicado a nuestros campos daría σ en función del
tiempo a lo largo del registro, que es lo que la corrección necesita y el barrido no
puede dar. Es la continuación natural del marco teórico, y hoy no está implementado.

## 10. Estado de la bibliografía

El proyecto distingue entre **fuente localizada** y **contenido inspeccionado**, y no
apoya ninguna afirmación en un paper que no haya sido leído. El balance actual:

| trabajo | estado |
|---|---|
| Foucaut et al. 2003 (ISPIV) | **leído en el original**; PDF local. La versión MST 2004 es de pago y no se leyó |
| George & Stanislas 2020 | resumen, símbolos y secciones 3.2–3.3 leídos; el desarrollo completo, no |
| Keane & Adrian 1990 | **no leído** |
| Willert & Gharib 1991 | **no leído** |
| Scarano 2002 | **no leído** |
| Theunissen et al. 2009 | **no leído** |
| Sciacchitano et al. 2013, 2015, 2019 | **no leídos** |
| Garcia 2010 | **no leído** |
| Mendes et al. 2020 | revisado y descartado |
| Manovski et al. 2025 | **no leído — pendiente** |

Todos los DOI fueron resueltos contra la API de Crossref el 2026-09-02; la salida cruda
está en `bibliografia/_crossref_raw.json`.

**El riesgo concreto que esto deja abierto** es que el resultado principal (§ 7) se apoya
en la descomposición de George & Stanislas, que es el trabajo del que menos leímos. La
frase que lo sostiene está verificada en su texto, pero no el desarrollo que la produce,
ni el régimen de sembrado y resolución para el que fue formulada. Antes de que § 7 entre
en el informe final, ese paper hay que leerlo entero.
