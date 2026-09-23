# Estado — dónde estamos y cómo retomar

**Instantánea del 2026-09-23.** Este archivo existe para que se pueda retomar el trabajo
sin la conversación que lo produjo: si la sesión se cae, o pasa una semana, alcanza con
leer esto, `DECISIONES.md` y `bibliografia/REFERENCIAS.md`.

Si retomás con un agente, la forma más corta de ponerlo al día es:

```
lee ESTADO.md, DECISIONES.md y bibliografia/REFERENCIAS.md de este directorio,
y seguí desde "Lo que sigue"
```

---

## En una línea

El decaimiento de una capa delgada está dominado por la fricción con el fondo, que nace
de la asimetría entre el fondo no-deslizante y la superficie libre. En las simulaciones
del grupo esa fricción entra como un parámetro ajustado, porque los códigos usados o no
tienen bordes o los tienen no-deslizantes de los dos lados. Este proyecto implementa en
SPECTER la condición **free-slip** —plana y sin tensión tangencial, la aproximación de esa
superficie libre válida cuando su deformación es despreciable ([D-28])— para **calcular**
α en vez de ajustarlo, y contrasta el resultado contra el decaimiento medido por PIV.

**Terminología, para no confundir dos ejes distintos.** Free-slip vs. no-deslizante es
sobre la tensión tangencial en un borde fijo; superficie libre (deformable) vs. tapa
rígida es sobre si ese borde se mueve. Lo implementado es free-slip **y además plano**, dos
aproximaciones apiladas y no la misma cosa. De acá en más: "free-slip" para la condición de
contorno (como ya la llama el código: `freeslip`, `freeslip_z`), "superficie libre" para el
escenario físico. Ver [D-40].

## El proyecto cambió de objeto el 2026-09-04

Hasta el 02-09 el objeto era el **Δt óptimo de PIV**. Se cambió por [D-21]: esa pregunta
se puede resolver tomando mejores mediciones, así que no deja instalada una herramienta
para el doctorado, y el plazo se extendió a un mes.

**Nada de aquel trabajo se tira.** El piso de ruido de PIV medido pasa a ser el capítulo
de instrumento y la barra de error de la comparación con el experimento. Las entradas
[V-04], [H-05], [H-06] y [H-07] de `DECISIONES.md` siguen valiendo tal como están.

## Restricciones fijadas, que no se renegocian sin decisión explícita

- **No se toman mediciones nuevas.** Todo corre sobre la campaña del 02-06-25 o es puro
  cálculo.
- **El cómputo pesado va al clúster** (Sakura). Esta laptop tiene 5,7 GB: sirve para
  desarrollo, tests 1D y verificación a resolución chica.
- Los datos crudos del disco USB son **sólo lectura** ([D-02]).
- **Los parámetros de la celda no se hardcodean.** Espesor, viscosidad y geometría
  de los imanes cambian entre experiencias y viven sólo en `codigo/celda.py`, con su
  procedencia ([D-27]).

---

## Estado de las fases

| fase | qué es | estado |
|---|---|---|
| **0** | Conseguir SPECTER y leer cómo impone las condiciones de contorno | **hecha** — [V-06], corregida en un punto por [H-08] |
| **1** | La derivación analítica y un solver de referencia no espectral | **hecha** — [V-07], puerta 6/6 |
| **2** | Implementar el tope free-slip (libre de tensiones tangenciales) en Fortran | **hecha** — [V-09], [D-24], [D-26] |
| **3** | Verificación V1–V6, con un caso que falla a propósito | **hecha** — puerta 6/6, [V-09] |
| **4** | La física, la comparación con el experimento y la entrega | **en curso** — decaimiento medido ([V-13]) y puente local 32² ([V-19]); falta convergencia 64², producción y entrega |

## Números ya establecidos

### De la Fase 1, verificados numéricamente

| | |
|---|---|
| Espectro con no-deslizamiento en las dos paredes | λ_n = ν(nπ/h)², n = 1, 2, … |
| Espectro con fondo no-deslizante y tope free-slip | λ_m = ν((m+½)π/h)², m = 0, 1, … |
| Fricción de fondo, tope free-slip | **α = νπ²/(4h²)** |
| Fricción de fondo, canal | α = νπ²/h² |
| Cociente entre ambas | **exactamente 4** — por tres rutas independientes |
| Sesgo del PIV de superficie | u(h)/⟨u⟩ = **π/2 = 1,5708** (exceso del 57 %) |
| Coeficientes de Laplace Dirichlet–Neumann | C₁ = (g₁/k + r·g₀)/(1+r²), C₂ = (g₀ − r·g₁/k)/(1+r²), con r = exp(−k·Lz) |
| Convención de base de SPECTER | φ = C₁·exp(k(z−Lz)) + C₂·exp(−kz) |

Con los parámetros de `codigo/celda.py` (h = 6 mm, ν = 10⁻⁶ m²/s): **α = 0,06854 s⁻¹**,
τ = 1/α = **14,59 s**, contra los 51,2 s que dura cada registro. Un modo de número de
onda k decae con **λ(k) = α(1 + 4ε²/π²)**, ε = kh, no con α pelado; para el fundamental
de la red de imanes (k = 296 m⁻¹, [D-41]) eso es **2,28 α**, y ahí está la pregunta
abierta [P-02].

### De la Fase 2, medidos con SPECTER ya modificado

Malla 8×8×128 (Cz=25, 103 puntos físicos en z), ν = 10⁻², Lz = 1, Runge-Kutta de orden 2.

| | |
|---|---|
| λ₀ con fondo rígido y tope free-slip | 2,467401099·10⁻², contra 2,467401100·10⁻² de la referencia — **error 6,8·10⁻¹⁰** |
| λ₁ (m=1), primer armónico vertical | 2,220660867·10⁻¹, contra 2,220660990·10⁻¹ — error 5,5·10⁻⁸ |
| λ₂ (m=2), segundo armónico vertical | 6,168500038·10⁻¹, contra 6,168502750·10⁻¹ — error 4,4·10⁻⁷ |
| λ(canal) / λ(free-slip) | **4,000000015** |
| Cara free-slip abajo en vez de arriba | misma λ dentro de 7,7·10⁻¹³ |
| Residuo de tensión en la cara free-slip | 6,1·10⁻¹² relativo a la cara rígida, **independiente de dt** |

Ese último número es la parte falsable de la derivación: la condición se impone sobre la
vorticidad tangencial, y como la proyección resta un gradiente, se conserva sin error de
partición. Ver [D-24].

> **Zanjado el 2026-09-09:** el espesor de esta campaña es **6 mm** (Lucía). Junto con la
> viscosidad y la geometría de los imanes quedó declarado en `codigo/celda.py` con su
> procedencia, y no vuelve a escribirse en ningún otro lado ([D-27]).

### Del capítulo de instrumento — el piso de ruido de PIV

Es la barra de error de la Fase 4. Convención: σ por componente.

| medición | σ | qué agrega |
|---|---|---|
| par sintético, corrimiento conocido [V-04] | ≤ **0,011 px** | correlación y estimador subpíxel solos |
| fluido en reposo, par real (campaña 05-09-24) | **0,030 px** | ruido de lectura, iluminación, movimiento propio |
| forzado 2,15 A, ajuste del barrido [H-05] | **0,169 px** | gradiente dentro de la ventana, movimiento fuera del plano |

El término que agrega el flujo es 0,166 px: el **98 % de la varianza** ([H-07]).

### De la Fase 4 — [P-01], la validez de la reducción 2D ([V-11])

Derivación en `teoria/P01_validez_reduccion_2D.md`, cuentas en
`numerico/fase4/p01_validez.py`.

| | |
|---|---|
| Parámetro que gobierna la clausura | **δ = U k h²/ν = Re_h·ε**, con ε = k h |
| ε para esta celda (k = 296 m⁻¹, [D-41]) | 1,78 → la viscosidad horizontal agrega **128 %** a α para ese modo |
| δ en el forzado (⟨u⟩ = 0,84 mm/s) | **8,97** |
| δ al inicio del decaimiento (⟨u⟩ = 2,55 mm/s) | **27,2** |
| Distorsión del perfil vertical, \|a₁/a₀\| ≈ 0,0077 δ | **7 %** y **21 %** |
| Memoria del segundo modo vertical, h²/(2π²ν) | 1,82 s, contra registros de 51,2 s |

**Lo que hay que retener:** el criterio heredado estaba mal escrito (ℓ en vez de k, ambiguo
por 2π), y el que lo reemplazó tenía el paso de los imanes mal por un factor 3,3 ([D-41]).
Con los números de hoy, δ no es lo que tiene que ser chico, pero la distorsión que sí lo
es vale **21 % al inicio del decaimiento**: la clausura de un modo es **mala en la meseta
y en los primeros segundos**, y buena recién dentro de la ventana del ajuste (t > 10 s),
cuando la energía migró a k ≲ 130 m⁻¹ y U cayó. El barrido en δ con SPECTER tiene que
llegar a δ ≈ 30.

### De la Fase 4 — el decaimiento medido ([V-13])

Página de revisión con figuras: `salidas/decaimiento.html`
(https://claude.ai/code/artifact/6363ddc3-404f-4aa0-993d-26bb1d67c17b).

| | α [1/s] | τ [s] |
|---|---|---|
| **ensemble de los 4, t > 10 s** | **0,06687 ± 0,00488** | **15,0** |
| ensemble, todos los puntos | 0,06189 ± 0,00336 | 16,2 |
| **predicho, free-slip** | **0,06854** | 14,6 |
| predicho, tapa rígida | 0,27416 | 3,6 |

**A 2,4 % del valor free-slip y 4,1 veces por debajo del de tapa rígida.** Es la
comparación que motivaba el proyecto. Limitado por `h`, no por la estadística: medio
milímetro de error en el espesor mueve α un 17 %.

> **Salvedad abierta desde el 2026-09-16 ([P-02]):** α pelado es la tasa de un modo con
> k → 0, y la de u_rms es α + ν⟨k²⟩, siempre mayor. Con el k que sí tiene el flujo en la
> ventana del ajuste (pico en 59–130 m⁻¹) la predicción es 0,072–0,085 s⁻¹ y el medido
> queda un 7–22 % **por debajo**. El "2,4 %" es contra el único valor que el flujo no
> puede tener. No está resuelto.

Dos cosas más que salieron de ahí:

- **El pico del espectro forzado es el fundamental de la red de imanes** ([H-16]): los
  doce espectros de la meseta pican en 295 m⁻¹, y la red con el paso correcto da 296
  ([D-41]; el 88,9 anterior venía de un paso de 5 cm que era un typo por 5 mm). El k
  pesado por energía va de 240 a 135 m⁻¹ durante el decaimiento, y con él δ arranca en
  ~20 y termina en ~1.
- **`med_S0008` contiene un decaimiento** aunque está etiquetado como forzado: cae de
  3,10 a 0,335 mm/s. `med_S0002` sí es estacionario. Ver [H-12].

### De la Fase 4 — el barrido en δ con SPECTER ([V-16])

`verificacion/barrido_delta_etapa1.py`, dos casos ([D-42]), dos mallas cada uno, corrido
en la laptop el 2026-09-16. Informes en `informe/barrido_delta_etapa1_{ventana,forzado}.pdf`.
Dos medidas independientes: λ/λ_lin (pendiente de log⟨v²⟩ contra la corrida lineal) y
α_eff/α (tensión en el fondo sobre promedio vertical, medida en el campo, contra νπ²/4h²).

| δ en la ventana | ε = 0,71 (ventana del ajuste): λ/λ_lin · α_eff/α | ε = 1,78 (meseta forzada): λ/λ_lin · α_eff/α |
|---|---|---|
| 0,3 | 1,000 · 1,001 | 1,000 · 1,001 |
| 1 | 1,005 · 1,004 | 1,002 · 1,004 |
| 3 | 1,044 · 1,025 | 1,014 · 1,028 |
| 5–7 | 1,143 · 1,083 | 1,007 · 1,155 |
| 8–10 | 1,27–1,31 · 1,18–1,28 | no convergido |

Convergido entre 16² y 32² dentro de 1,2 % en ε = 0,71 hasta δ = 10, y hasta δ = 7 en
ε = 1,78 (los dos puntos de más amplitud no convergen: [H-09]).

**Lo que hay que retener:** el término no lineal **sube** la fricción de fondo — 2,5 % en
δ = 3, 8 % en δ = 5, 15–18 % en δ ≈ 7–8, 26–28 % en δ = 10 — y la estimación de orden de
[V-11] la subestimaba por un factor ~3. En ε = 1,78 la tasa total no lo muestra porque
la energía migra a escalas grandes que disipan menos horizontalmente; en ε = 0,71 sí. En
la ventana del ajuste experimental (δ de ~10 a ~1,3) la corrección es del 10–30 % en la
primera mitad y despreciable en la segunda. **[P-02] no se explica por esto: empeora**;
el candidato nuevo es que el PIV ajusta u(h) y no ⟨u⟩, y con perfil distorsionado
u(h)/⟨u⟩ ≠ π/2 y cambia con el tiempo. Medido en [V-17], abajo.

### De la Fase 4 — la tasa de la velocidad de superficie ([V-17])

`barrido_delta_etapa1.py --superficie`, las mismas corridas de [V-16] escribiendo cinco
campos en la ventana; 16² solamente. λ_sup es la pendiente de log u(h) —lo que ajusta un
PIV de superficie—, λ_prom la de log⟨u⟩; cada una dividida por la misma tasa en la
referencia lineal. JSON en `salidas/tablas/barrido_delta_etapa1_*_superficie.json`.

| δ | ε = 0,71: λ_sup/λ_lin · λ_prom/λ_lin · u(h)/⟨u⟩ | ε = 1,78: λ_sup/λ_lin · λ_prom/λ_lin · u(h)/⟨u⟩ |
|---|---|---|
| 1 | 1,004 · 1,005 · 1,57 | 1,000 · 1,001 · 1,57 |
| 3 | 1,033 · 1,040 · 1,56 | 0,996 · 1,011 · 1,55 |
| 5 | 1,092 · 1,128 · 1,52 | — |
| 7–8 | 1,118 · 1,242 · 1,45 | **0,891** · 1,013 · 1,44 |
| 10 | 1,043 · 1,305 · 1,37 | no convergido |

**Lo que hay que retener:** el perfil distorsionado tiene menos exceso de superficie
(u(h)/⟨u⟩ baja de π/2 = 1,571 a 1,37 en δ ≈ 10) y, como δ cae en el tiempo, el cociente
se recupera dentro de la ventana; por eso u(h) decae más lento que ⟨u⟩ (hasta 20 % en
δ ≈ 10). En ε = 0,71 eso casi cancela la subida no lineal de la fricción: la tasa que ve
un PIV queda entre +3 % y +12 % de la lineal hasta δ = 10, **nunca por debajo**. En
ε = 1,78 no hay qué cancelar (la tasa total no sube, [V-16]) y la tasa de u(h) cae **11 %
por debajo de la lineal** en δ = 7. Para [P-02]: el candidato (v) es real y cambia el
modelo de comparación (la predicción para el PIV es λ_sup, no λ_prom), pero con el ε de la
ventana del ajuste no baja la predicción por debajo de λ(k); no alcanza.

### De la Fase 4 — prueba espectral del decaimiento ([V-18])

Los 32 espectros ya calculados permiten una prueba por banda sin volver a los TIF. Después
de restar el piso de PIV y corregir las distintas separaciones entre cuadros, la amplitud
espectral total decae a **0,0711 ± 0,0024 s⁻¹**, consistente con la ruta independiente de
dos Δt ([V-13]). Entre 11 y 40 s la fracción de energía en 59–106 m⁻¹ sube de 20 % a
58 %, mientras la de 224–342 m⁻¹ baja de 44 % a 17 %: la escala crece durante el ajuste.

En ocho bins entre 83 y 248 m⁻¹ las tasas aparentes siguen
**λ(k) = (0,0481 ± 0,0038) s⁻¹ + (1,026 ± 0,095)·10⁻⁶ m²/s·k²** (incerteza entre los
cuatro registros). La pendiente recupera ν sin imponerla; el déficit queda casi entero en
una ordenada menor que α = 0,0685 s⁻¹. Esto descarta que "poco turbulento" equivalga sin
más al límite lineal, pero no identifica la ordenada con fricción: cada banda intercambia
energía, el PIV mide la superficie y el piso es correlacionado. Script y salida:
`codigo/13_decaimiento_espectral.py`, `salidas/tablas/decaimiento_espectral.json`, figura
`f9`.

### De la Fase 4 — puente numérico de banda ancha ([V-19])

La prueba SPECTER local usa la forma de los diez anillos medidos a `t ≈ 11,1 s` con
`delta_0 = 8,887`, y la compara contra la misma condición a amplitud `10^-4`. En 32², el
control lineal recupera **`alpha = 0,06777 s^-1`, `nu = 0,987·10^-6 m²/s`** (`R² =
0,99992`). La banda ancha da **`alpha = 0,06544 s^-1`, `nu = 1,171·10^-6 m²/s`**: baja
la ordenada sólo 0,00233 s⁻¹, contra el déficit PIV de 0,02048 ± 0,00378 s⁻¹. La tasa
global de superficie baja sólo 1,5 %.

El efecto no lineal existe y redistribuye tasas entre anillos, pero a esta resolución es
unas nueve veces menor que lo que habría que explicar. La 32² no es el cierre: el anillo
10 está cerca del corte de dealiasing. El siguiente paso inmediato es el par 64² por SLURM;
el runner, el `sbatch` y el handoff autocontenido están listos. JSON en
`salidas/tablas/decaimiento_banda_ancha.json`, figura `f10`.

### De la campaña

| | |
|---|---|
| Cámara | Photron FASTCAM-1024PCI, 1024×1024, 10 bits, 60 fps, shutter 1/60 |
| Escala | 3800 px/m → 0,263 mm/px, campo de 26,9 cm |
| Corridas | 3 forzados (0,95 / 1,95 / 2,15 A), **4 decaimientos**, 3 quenchings (pulso breve de forzado, máximo en t = 12–14 s; [H-18]) |
| Cada corrida | 3072 cuadros = 51,2 s |
| Capa de electrolito | 6 mm, KNO₃ al 16 % m/m, partículas de 100 µm |
| Imanes | discos de 1 cm, 5 mm de vacío entre ellos → red tipo tablero de **1,5 cm** entre centros; acrílico de 6 mm hasta la capa ([D-41]) |
| Armónico dominante del forzado | **diagonal**, λ = 2,12 cm, k = 296,2 m⁻¹ — medido en la meseta forzada: 295 ([H-16]) |
| Velocidad del forzado a 2,15 A | 1,321 mm/s, IC68 [1,306–1,346] |
| Velocidad al inicio de un decaimiento | **4,07 mm/s**, medida y corregida por ruido ([V-13]); cierra [H-10] |

Los cuatro decaimientos (`med_S0003`, `S0005`, `S0006`, `S0009`) son un **ensemble**: las
barras de error de la comparación con la simulación salen de ahí, no de una sola curva.

---

## Lo que sigue: la Fase 4, la física y la comparación con el experimento

La herramienta ya existe: SPECTER puede calcular α en vez de ajustarlo. Lo que queda es
usarla y contrastarla.

El árbol de trabajo es `numerico/SPECTER-trabajo/`, clon de `SPECTER-upstream` en
`0ad1edb`. El diff contra upstream son **237 líneas agregadas y 22 borradas** en 5
archivos, más `src/tests/laplace_dirneu.f90`. **Todavía no está commiteado**: hacer
`git -C numerico/SPECTER-trabajo diff` muestra exactamente el aporte propio.

Para volver a correr la verificación:

```bash
/home/lucia/miniforge3/envs/piv-dt/bin/python verificacion/test_aceptacion_fase2.py
```

Compila SPECTER en un scratch temporal con el env de conda `specter` ([D-23]) y corre
unas quince simulaciones chicas. Tarda del orden de diez minutos y no pasa de unos pocos
cientos de MB. `SPECTER_KEEP=1` deja el scratch para inspeccionar.

**Hechos el 2026-09-09:** [P-01] queda resuelto ([V-11]) y el espesor, la viscosidad y la
geometría de los imanes quedan declarados en `codigo/celda.py` ([D-27]).

**Hechos el 2026-09-16:** el paso de los imanes corregido a 1,5 cm ([D-41]) a partir de
que Lucía notó que el pico del espectro no estaba donde la figura marcaba el forzado
([H-16]); ε, δ y λ(k) recalculados, y la comparación de [V-13] abierta como [P-02]. El
barrido en δ rediseñado con esa geometría ([D-42]) y lanzado en la laptop, secuencial,
con `ulimit -v`; logs en `verificacion/logs/barrido_delta_etapa1_*.log`.

**Hechos el 2026-09-17:** la etapa 1b del barrido cerró y está en el log como [V-17]. El
repositorio se reencuadró para hacerse público centrado en lo numérico ([D-43]): README
nuevo, receta del clúster en `numerico/CLUSTER.md`, `notas` fuera del árbol, `jobs/README.md`.

**Hecho el 2026-09-21:** la prueba espectral rápida de [V-18] confirma migración hacia k
bajos y recupera el término νk², pero deja una ordenada común de 0,048 s⁻¹. [P-02] no se
cierra: ahora el déficit está localizado en la parte casi independiente de k.

**Hecho el 2026-09-22:** la prueba puente [V-19] se corrió localmente en 16² y 32². El
control lineal recupera `alpha + nu k²`; la banda ancha con `delta_0 = 8,887` reduce la
ordenada sólo 0,00233 s⁻¹, unas nueve veces menos que el déficit PIV. Falta convergencia
64² en Sakura; el handoff y el `sbatch` están escritos, pero no se envió ningún trabajo.

**Hecho el 2026-09-23:** se trajo al log lo que el análisis de los quenching para el póster
de tesis (fuera del repo, `~/Desktop/PhD/Docs/poster_powerpoint_first_layout/experimental/`)
tenía de útil acá: un fondo de movimiento en la banda de los imanes, anterior al forzado y
que decae a 0,5 α ([H-18]). Restado como piso constante sube el α medido 4–5 %: candidato
parcial para [P-02], no cierre.

**Lo que sigue, en este orden:**

1. ~~Medir α_eff con SPECTER en un barrido en δ~~ — **hecho** ([D-42], [V-16], y la
   etapa 1b de superficie en [V-17]). La prueba puente de banda ancha de [V-18] ya está
   hecha localmente hasta 32² ([V-19]); **lo inmediato es convergerla en 64² en Sakura**.
   Después quedan una malla más fina para los dos puntos no convergidos del caso forzado
   (δ ≥ 11), estudio de [H-09] para producción, y u(h) comparado entre mallas, que sólo se
   midió en 16² en aquel barrido.
2. **Una corrida de producción en Sakura** con la geometría real de la celda, y de ahí el
   α calculado. La resolución sale de 1: 32² alcanza hasta δ ≈ 7 en ε = 1,78 y hasta
   δ ≈ 10 en ε = 0,71.
3. ~~El decaimiento medido~~ — **hecho** ([V-13], [H-11], [H-12]).
4. ~~Reanalizar el capítulo de instrumento a la luz de [H-12]~~ — **se cae del plan**
   ([D-43]): es instrumento, no objeto. Queda anotado como pregunta abierta de la zona
   exploratoria; si sobra tiempo, es rehacer el barrido en Δt sobre `med_S0002`.
5. **La entrega**: el repositorio se hace **público** con el README nuevo ([D-43], [D-44];
   la visibilidad la cambia Lucía en GitHub). `salidas/entrega1_resumen.html` ya está
   reordenada en las tres zonas, con la sección de superficie libre (no deformable: por
   qué vale; deformable: el siguiente paso), [V-17] y la figura `f8`; publicada en
   https://claude.ai/artifact/21j1baf1T54uYqdkF99Xwi. El PDF explicado está en
   `informe/entrega.pdf` (`codigo/12_pdf_entrega.py`), 12 páginas. **El repositorio es
   público desde el 2026-09-17 a las 18:19.** Lo que queda de la entrega es integrar la
   corrida de producción cuando salga. Plazo: completo unos días antes de las
   presentaciones, que son en dos semanas (~2026-10-01).

Los puntos 2 y 5 son independientes. La convergencia 64² del punto 1 y la corrida del
punto 2 son los cómputos pesados que quedan.

**Lo que la Fase 2 dejó pendiente, y hay que resolver antes de la corrida de producción:**

- El diagnóstico de tensión se calcula ahora en toda corrida con paredes en z, y pide una
  matriz temporal del tamaño del campo. A resolución de producción son cientos de MB por
  cada `cstep`. Si molesta, hay que condicionarlo ([D-26]).
- **La resolución horizontal no se hereda de una corrida de canal.** Con el tope free-slip
  hay menos disipación, y sobre la malla chica de los tests el término no lineal se
  desestabiliza donde el canal aguanta ([H-09]). Antes de producción, estudio de
  convergencia propio con el tope free-slip puesto.
- **La puerta no verifica el régimen no lineal.** V4, V5 y V6 usan el modo de corte puro,
  donde el término no lineal es exactamente cero; V3 usa un campo 3D pero sólo 400 pasos.
- Sólo se probó en doble precisión y con `ORD=2`. Con 1, 2 y 4 procesos MPI sí está
  verificado, y coincide en 15 cifras ([V-10]).
- El borde superior es **free-slip y plano**: `w = 0` en el tope, sin tensión tangencial.
  La deformación de la superficie libre real está fuera de alcance y cuantificada como
  despreciable ([D-28], η/h ~ 10⁻⁵); si alguna vez hiciera falta modelarla, ese sería el
  caso que sí necesitaría la rama Neumann–Dirichlet de `laplace_z`, que no está.

## Preguntas abiertas

- **[P-02] El α medido está por debajo de λ(k) para cualquier k con energía.** Es la
  comparación principal del proyecto y hoy se lee como acuerdo sólo contra α pelado.
  [V-16] descartó la no linealidad como explicación (empuja para el otro lado) y [V-17]
  midió el candidato (v): el PIV ajusta u(h), no ⟨u⟩, y con perfil distorsionado u(h)
  decae más lento; es real y hay que incluirlo en la comparación, pero con el ε de la
  ventana del ajuste no baja la predicción por debajo de λ(k). [V-18] localizó el déficit
  en una ordenada casi independiente de k y [V-19] muestra que, a 32², la transferencia de
  una banda ancha con la intensidad medida explica sólo alrededor de una décima; falta su
  convergencia 64². Quedan además h > 6 mm, los sesgos del PIV (piso correlacionado y
  calibración temporal, [V-14]) y el fondo en la banda de los imanes de [H-18], que
  explicaría 4–5 % de la tasa si se confirma con el estimador de [V-13].
- **[P-01] cerrado en lo cuantitativo ([V-16]).** La clausura de un modo tiene un error
  de pocos por ciento para δ ≲ 3 y de 15–30 % en δ ≈ 7–10, medido con el código y no
  estimado; el sesgo sobre α es hacia arriba. [V-19] ya probó localmente un campo de banda
  ancha; falta su convergencia 64². La meseta como estado forzado no se probó.
- ~~La menor escala con energía del flujo medido~~ — **medida** ([H-11], [H-16]): pico
  en 295 m⁻¹ en la meseta, que es el fundamental de la red; migra a 59–130 durante el
  decaimiento.
- **La viscosidad es la del agua, no la del electrolito.** Decisión explícita de usarla y
  declararla ([D-27]); entra a primer orden en α.
- ~~La superficie libre real podría no ser free-slip ideal por una película superficial~~
  — **no aplica** ([D-46]): no ocurre en este montaje ni es relevante para el caso de
  estudio (Lucía). Queda modelada en [D-29] y se menciona como salvedad breve en README
  y entrega, sin presentarla como la salvedad principal.
- ~~El factor 2,7 entre forzado y decaimiento~~ — **explicado** ([H-12]): el registro
  anómalo era `med_S0008`, que contiene un decaimiento pese a estar etiquetado como
  forzado. El forzado que sí es estacionario (`med_S0002`, con menos corriente) da
  3,55 mm/s, del mismo orden que los 4,07 mm/s con que arranca el decaimiento. **No era
  sesgo de supervivencia.**
- ~~La meseta inicial de los decaimientos~~ — **explicada** ([D-30]): es el estado
  forzado previo al corte, grabado a propósito para no tener que coordinar el corte con
  el inicio del registro. Se excluye del ajuste, y el α que vale es el de t > 10 s.
- **σ = 0,169 px se midió sobre un registro no estacionario** ([H-12]). Cuánto lo
  contamina no está establecido. Es la barra de error de la comparación principal.
- **El marco de Foucaut no describe el ruido de este multipaso** ([V-14]): el piso del
  espectro no tiene el cero en 1492 m⁻¹ que exigiría el sinc² de la ventana de 16 px, y
  se comporta como si el ruido estuviera correlacionado sobre el paso de la grilla. Toca
  al capítulo de instrumento ([V-05], [D-18]).
- **La rama Dirichlet–Neumann de `laplace_z` quedó implementada pero sin uso en el camino
  físico.** [H-08] corrige a [V-06] en este punto. Está verificada como test unitario. Si
  alguna vez se ataca la superficie libre deformable, el par que hace falta es el
  opuesto, Neumann abajo – Dirichlet arriba, y no está.
- El método de dos Δt sobre `Decaimiento_4` daba 0,083 px sobre `.vec` con reemplazo de
  vectores; [H-06] muestra que ese reemplazo invierte el signo de la dependencia con n,
  así que no es comparable. Contrastable rehaciéndolo sobre campos sin filtrar.

## Lo que no se toca

- Los datos crudos del 02-06-25, en el disco externo: **sólo lectura** ([D-02]). Es la
  única copia de una campaña de un día que no se puede repetir.
- `numerico/SPECTER-upstream/`: copia intacta de upstream. Todo el trabajo va en
  `numerico/SPECTER-trabajo/`.
- `verificacion/test_aceptacion_fase1.py` y `verificacion/test_aceptacion_fase2.py`: las
  puertas están fijadas por hash. Se pueden leer, no se editan para que pasen.
  Fase 1: `0c3ebd62f6c42d08f0fa9a64a0f57b9cfdaf1144ee8a84e968a148d9d66b18bf`.
  Fase 2: `6f0b4da4ba20617d483dfff0a20f8e2e11e82795b20351c259ffcf740ca4b1cc`.
- `jobs/2026-09-04_205253_derive-fase1-capa-delgada/out/`: es el registro de lo que
  produjo el equipo de agentes. Las correcciones van en `DECISIONES.md`, no ahí adentro.
- El otro árbol, `~/GW-AI-course/piv-s0008-forzado`, tampoco se toca ([D-19]).

> **Higiene para la entrega:** hay contraseñas del clúster en texto plano en las notas de
> Obsidian (`~/Desktop/PhD/Obsidian/`). El repositorio es público. No deben entrar.
