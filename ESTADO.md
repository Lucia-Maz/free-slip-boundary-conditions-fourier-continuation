# Estado — dónde estamos y cómo retomar

**Instantánea del 2026-09-09 (segunda del día).** Este archivo existe para que se pueda retomar el trabajo
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
| **4** | La física, la comparación con el experimento y la entrega | **en curso** — [P-01] resuelto ([V-11]) y **decaimiento medido ([V-13])**; falta la corrida de producción y la entrega |

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
τ = 1/α = **14,59 s**, contra los 51,2 s que dura cada registro. Para el armónico
fundamental del forzado hay que comparar contra **λ(k) = α(1 + 4ε²/π²) = 1,115 α**, no
contra α pelado ([V-11]).

### De la Fase 2, medidos con SPECTER ya modificado

Malla 8×8×128 (Cz=25, 103 puntos físicos en z), ν = 10⁻², Lz = 1, Runge-Kutta de orden 2.

| | |
|---|---|
| λ₀ con fondo rígido y tope free-slip | 2,467401099·10⁻², contra 2,467401100·10⁻² de la referencia — **error 6,8·10⁻¹⁰** |
| λ₁, λ₂ | error relativo 5,5·10⁻⁸ y 4,4·10⁻⁷ |
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
| ε para esta celda | 0,533 → la viscosidad horizontal agrega **11,5 %** a α |
| δ en el forzado (⟨u⟩ = 0,84 mm/s) | **2,69** |
| δ al inicio del decaimiento (⟨u⟩ = 2,55 mm/s) | **8,15** |
| Distorsión del perfil vertical, \|a₁/a₀\| ≈ 0,0077 δ | **2 %** y **6 %** |
| Memoria del segundo modo vertical, h²/(2π²ν) | 1,82 s, contra registros de 51,2 s |

**Lo que hay que retener:** el criterio heredado estaba mal escrito (usaba ℓ = 2 cm, que
no es ninguna longitud de la celda, y ℓ en vez de k, que es ambiguo por 2π). Bien escrito,
δ es **peor** que lo que decía la pregunta abierta; pero δ no es lo que tiene que ser
chico, y la distorsión que sí lo es queda en pocos por ciento. **La clausura no queda
invalidada en el régimen medido.**

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

Dos cosas más que salieron de ahí:

- **El k del flujo no es el del forzado** ([H-11]): pesado por energía va de 240 a
  135 m⁻¹ contra los 88,9 del forzado, y con él δ arranca en ~20 y termina en ~1.
- **`med_S0008` contiene un decaimiento** aunque está etiquetado como forzado: cae de
  3,10 a 0,335 mm/s. `med_S0002` sí es estacionario. Ver [H-12].

### De la campaña

| | |
|---|---|
| Cámara | Photron FASTCAM-1024PCI, 1024×1024, 10 bits, 60 fps, shutter 1/60 |
| Escala | 3800 px/m → 0,263 mm/px, campo de 26,9 cm |
| Corridas | 3 forzados (0,95 / 1,95 / 2,15 A), **4 decaimientos**, 3 quenchings |
| Cada corrida | 3072 cuadros = 51,2 s |
| Capa de electrolito | 6 mm, KNO₃ al 16 % m/m, partículas de 100 µm |
| Imanes | discos de 1 cm, red tipo tablero de 5 cm de paso |
| Armónico dominante del forzado | **diagonal**, λ = 7,07 cm, k = 88,86 m⁻¹ ([V-11]) |
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

**Lo que sigue, en este orden:**

1. **Medir α_eff con SPECTER en un barrido en δ.** Es lo que cierra el fleco que dejó
   [P-01]: el promedio vertical da, sin ninguna clausura,
   `α_eff = (ν/h)·∂_z u|₀ / ⟨u⟩`, y la curva de `α_eff/λ(k)` contra δ ≈ 0,3…10 **es** la
   respuesta. Con dos resoluciones horizontales, porque además es el estudio de
   convergencia con el tope free-slip puesto que [H-09] dejó pendiente para producción. El diseño
   está escrito en la sección 7 de `teoria/P01_validez_reduccion_2D.md`.
2. **Una corrida de producción en Sakura** con la geometría real de la celda, y de ahí el
   α calculado. Depende de 1 para la resolución.
3. ~~El decaimiento medido~~ — **hecho** ([V-13], [H-11], [H-12]).
4. **Reanalizar el capítulo de instrumento a la luz de [H-12].** `med_S0008` no es
   estacionario, y el barrido en Δt supone que lo es. Hay que rehacerlo sobre
   `med_S0002`, que sí lo es, y ver qué le pasa a σ = 0,169 px, que es la barra de error
   de todo el trabajo.
5. **La entrega**: repositorio de GitHub, página HTML y PDF ([D-19] deja constancia de
   que hoy ninguno de los dos árboles es un repo git). Ya existe
   `salidas/decaimiento.html` como página de revisión del decaimiento.

Los puntos 1, 2 y 4 son independientes entre sí; el 1 y el 4 son cómputo local.

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

- **[P-01] quedó resuelto ([V-11]), pero con dos flecos.** La clausura no queda
  invalidada —la distorsión del perfil es de pocos por ciento—, pero la estimación toma
  igual a 1 un factor geométrico O(1) del término no lineal que no se calculó, y no está
  decidido si el sesgo sobre α es de orden δ o δ². Los dos se cierran midiendo α_eff.
- **La menor escala con energía del flujo medido.** δ ∝ k, así que el criterio lo fija esa
  escala y no el fundamental del forzado. Hay que sacarla del espectro o de la
  autocorrelación de los campos PIV. Necesita el disco montado.
- **La viscosidad es la del agua, no la del electrolito.** Decisión explícita de usarla y
  declararla ([D-27]); entra a primer orden en α.
- **La superficie libre real no es free-slip ideal.** Un electrolito con
  partículas flotando puede desarrollar una película superficial que la vuelve casi
  rígida — y ese es justo el caso no-deslizante, con α cuatro veces mayor. Es la salvedad
  más fuerte del trabajo. Hay que declararla, y si se puede, acotarla con los datos.
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
