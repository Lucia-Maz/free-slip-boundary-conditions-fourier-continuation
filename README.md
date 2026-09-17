# Fricción de fondo de una capa delgada: la condición free-slip en SPECTER

Proyecto final del curso *Ondas Gravitacionales e Investigación Asistida por IA* (2026).
Lucía Mazaira, INFINA (CONICET-UBA).

## La pregunta

Una capa delgada de fluido que fluye sobre un fondo sólido se frena por la fricción con ese
fondo. En las simulaciones de turbulencia cuasi-2D esa fricción entra como un término
lineal −α**u** con un α que **se ajusta**, porque los códigos usados o no tienen bordes
(periódicos) o los tienen no-deslizantes de los dos lados. Pero el α nace de una asimetría
concreta: fondo no-deslizante abajo, superficie libre arriba. La pregunta del proyecto es
si α se puede **calcular** en vez de ajustar, y qué le hace la no linealidad.

Para eso se implementó en [SPECTER](https://github.com/mfontanaar/SPECTER) —pseudoespectral
con continuación de Fourier, no periódico en z— la condición **free-slip** en la cara
superior: plana y sin tensión tangencial, que es la aproximación de una superficie libre
válida cuando su deformación es despreciable (acá, η/h ~ 10⁻⁵). Después se verificó, se
usó para medir α en régimen no lineal, y se contrastó contra un decaimiento medido.

**Terminología.** "Free-slip" es la condición de contorno (así la llama el código:
`freeslip`, `freeslip_z`); "superficie libre" es el escenario físico que la motiva. No son
sinónimos: lo implementado es free-slip *y además plano*.

## Lo que se estableció

Cada número apunta a lo que lo produce: un chequeo ejecutable, un script y su JSON, y la
entrada del log donde está discutido.

| | | procedencia |
|---|---|---|
| Fricción de fondo con tope free-slip | **α = νπ²/(4h²)**; con tapa rígida, νπ²/h²: **exactamente 4 veces más** | `verificacion/test_aceptacion_fase1.py` (referencia por diferencias finitas y Richardson); derivación en `jobs/…/out/build/notes.pdf` con `out/provenance.json`; [V-07] |
| SPECTER con el parche reproduce la tasa del modo lento | error **6,8·10⁻¹⁰** contra la referencia; residuo de tensión en la cara free-slip 6,1·10⁻¹², independiente de dt | `verificacion/test_aceptacion_fase2.py`, casos V4 y V3 (compila y corre SPECTER); [V-09], [D-24] |
| La no linealidad **sube** la fricción de fondo | α_eff/α = 1,025 en δ = 3, 1,08 en δ = 5, 1,15–1,18 en δ ≈ 7–8, 1,26–1,28 en δ = 10 (medido, convergido entre dos mallas) | `verificacion/barrido_delta_etapa1.py --caso {ventana,forzado}` → `salidas/tablas/barrido_delta_etapa1_{ventana,forzado}.json`; informes `informe/barrido_delta_etapa1_*.pdf`; [V-16] |
| La velocidad de superficie decae más lento que el promedio vertical | u(h)/⟨u⟩ baja de π/2 a 1,37 en δ ≈ 10 y se recupera en el tiempo; para un PIV de superficie en ε = 0,71 la predicción lineal queda buena al 12 % aunque cada ingrediente no lineal valga 20–30 %; en ε = 1,78 la tasa de u(h) cae 11 % por debajo de la lineal | el mismo script con `--superficie` → `salidas/tablas/barrido_delta_etapa1_*_superficie.json`; [V-17] |

Los parámetros que gobiernan: ε = kh (esbeltez del modo) y δ = U k h²/ν = Re_h·ε. Las
entradas `[V-16]` y `[V-17]` de `DECISIONES.md` tienen las tablas completas.

## Dónde termina el resultado

### El siguiente paso: la superficie deformable

La propuesta habla de *superficie libre*, y lo implementado es free-slip **plano**: es el
orden cero de esa superficie libre en el número de Froude, y en esta celda el orden
siguiente es despreciable con número, no con criterio — la presión dinámica ρU² contra la
restitución ρg + σk² deforma la superficie **0,04–0,33 µm sobre 6 mm, η/h ≲ 6·10⁻⁵**
(`numerico/fase4/superficie_libre_escalas.py` → `salidas/tablas/superficie_libre_escalas.json`; [D-28]). Por eso el α de este proyecto no cambia con
la superficie deformable, y por eso no se implementó acá.

Es, en cambio, **la continuación del trabajo**, y ya está preparada:

- Las condiciones linealizadas están derivadas: cinemática `w(h) = ∂η/∂t`, tensión
  tangencial completa, y tensión normal que prescribe `p` en la superficie
  (`teoria/superficie_libre_v_estrella_y_p.md` §4).
- La rama de presión que hace falta, Neumann abajo – Dirichlet arriba, está derivada y
  verificada en Python sobre 20 000 combinaciones de k y Lz, con su criterio de
  condicionamiento ([V-12]). No está en Fortran.
- Lo que falta, dicho sin adornos: esa rama en `laplace_z`; un `bctarget` por cara en
  `sol_project`; la condición sobre `v*` queda **acoplada implícitamente** a `∂p/∂z` en la
  superficie (iterar dentro del subpaso, o aceptar un error de partición O(dt)); el campo
  η(x, y, t) avanzado con el mismo Runge-Kutta; y una puerta nueva, porque la propiedad
  que hace exacto al caso plano —residuo independiente de dt— no sobrevive. La puerta
  natural es la relación de dispersión de ondas gravito-capilares amortiguadas.

### Fuera de alcance de este proyecto, a propósito

- **La película superficial.** Un electrolito con partículas flotando puede desarrollar
  una película que vuelve el tope casi rígido, y ese es justo el caso con α cuatro veces
  mayor: es **la salvedad más fuerte del trabajo**. Está modelada (Boussinesq–Scriven:
  condición de Robin con un solo parámetro β, y el factor 4 de la Fase 1 son los dos
  extremos de esa familia, [D-29]) y no implementada. El α medido acota β por arriba.
- **El régimen forzado como estado estacionario.** El barrido en δ es de decaimiento
  libre; la meseta forzada del experimento no se simuló como tal.
- **La puerta no verifica el régimen no lineal.** V4–V6 usan el modo de corte puro, donde
  el término no lineal es cero; la convergencia en régimen no lineal se estudió aparte
  ([V-16], [H-09]) y sólo hasta δ ≈ 7–10.
- **Sólo doble precisión y `ORD=2`.** Con 1, 2 y 4 procesos MPI sí está verificado ([V-10]).
- **La viscosidad es la del agua**, no la del electrolito al 16 % ([D-27]); entra a primer
  orden en α.

## Qué hay acá

El repositorio tiene tres zonas, por madurez:

1. **Lo numérico** — el objeto del proyecto, cerrado y verificado.
2. **El contraste con el experimento** — un número medido, con su salvedad.
3. **Explorando** — lo que salió de preguntas disparadas por los resultados; en proceso.

### 1. Lo numérico

| fase | qué | dónde |
|---|---|---|
| 0 | Qué condiciones de contorno tiene SPECTER, leído en la fuente | `DECISIONES.md` [V-06], [H-08] |
| 1 | La derivación analítica (espectros, α, sesgo π/2 del PIV de superficie) y un solver de referencia **no espectral** | `jobs/2026-09-04_…/out/build/notes.pdf` (la derivación), `out/capa.py` (la referencia, lo que carga la puerta); `teoria/` |
| 2 | La implementación en Fortran: 237 líneas agregadas y 22 borradas en 5 archivos, más un test unitario | `numerico/specter-parche/` (el parche contra upstream `0ad1edb` y la receta) |
| 3 | La verificación, seis casos con uno que falla a propósito | `verificacion/test_aceptacion_fase{1,2}.py` |
| 4 | La física: el barrido en δ (fricción efectiva, tasa del promedio vertical y de la superficie), y la validez de reducir a 2D | `verificacion/barrido_delta_etapa1.py`, `numerico/fase4/`, `informe/barrido_delta_etapa1_*.pdf`, `teoria/P01_validez_reduccion_2D.md` |

Las **puertas de aceptación** (`verificacion/test_aceptacion_fase*.py`) están fijadas por
hash: se escribieron antes de la implementación, arrancaron en rojo, y **no contienen la
fórmula analítica** de lo que verifican — construyen su propia referencia por diferencias
finitas y extrapolación de Richardson, para que el resultado quede comprobado y no
afirmado.

### 2. El contraste con el experimento

Con los parámetros de la celda (`codigo/celda.py`: h = 6 mm, ν = 10⁻⁶ m²/s) el α
calculado es **0,0685 s⁻¹** (τ = 14,6 s). El decaimiento medido por PIV en la campaña del
02-06-25, ensemble de cuatro registros y t > 10 s, da **0,0669 ± 0,0049 s⁻¹** — a 2,4 %
del valor free-slip y **4,1 veces por debajo** del de tapa rígida (0,274 s⁻¹).

**Salvedad, que va siempre al lado del número ([P-02]):** α pelado es la tasa de un modo con
k → 0. Para el k que el flujo sí tiene en la ventana del ajuste (59–130 m⁻¹) la tasa
lineal es 0,072–0,085 s⁻¹, y el medido queda un 7–22 % **por debajo**. El barrido en δ
descartó la no linealidad como explicación (empuja para el otro lado) y mostró que la
predicción correcta para lo que ve un PIV de superficie es la tasa de u(h), no la de ⟨u⟩;
con eso la comparación sigue abierta. El detalle está en `DECISIONES.md` [V-13], [P-02],
[V-16], [V-17].

Código: `codigo/06_decaimiento.py` (necesita el disco externo) y `07`/`08` (figuras y
página, desde `salidas/tablas/decaimiento*.json`, sin el disco). Página de revisión:
`salidas/entrega1_resumen.html`.

### 3. Explorando: preguntas disparadas por los resultados

Esta zona **está en proceso** y no es parte del objeto del proyecto. Se conserva porque el
contraste de arriba depende de ella y porque documenta lo que se intentó.

- **El capítulo de instrumento: parámetros óptimos de PIV.** Fue el objeto original del
  proyecto hasta el 04-09 ([D-21]). Dejó medido el piso de ruido del PIV — ≤ 0,011 px
  (par sintético), 0,030 px (fluido en reposo), 0,169 px (forzado) — y el hallazgo de que
  el 98 % de la varianza la agrega el flujo, no el instrumento ([H-07]). Abierto: el
  0,169 px se midió sobre un registro que resultó no estacionario ([H-12]), y el marco de
  Foucaut no describe el ruido de este multipaso ([V-14]). Código: `codigo/01`–`05`, `09`;
  teoría: `informe/01_marco_teorico.md`; resultados: `salidas/tablas/*med_S0008*`,
  `corrimiento_sintetico`, `desplazamiento_nulo`, `verificacion_filtrado`.

### Los archivos pensados para humanos y los pensados para agentes

| | |
|---|---|
| `ESTADO.md` | dónde estamos y qué sigue: **empezar por acá** |
| `DECISIONES.md` | el log completo, una entrada por decisión: qué se decidió, por qué, qué alternativa se descartó, qué salió mal. Se agrega al final, no se reescribe. Las etiquetas `[D-nn]` decisión, `[H-nn]` hallazgo, `[V-nn]` verificación, `[P-nn]` pregunta abierta resuelven ahí |
| `CLAUDE.md` | las convenciones del proyecto, que lee el agente al abrirlo en cualquier máquina |
| `bibliografia/REFERENCIAS.md` | los papers con DOI resueltos contra Crossref, y cuáles se leyeron en el original |
| `salidas/tablas/*.json` | todo resultado numérico, con la configuración que lo produjo adentro del mismo archivo |
| `salidas/figuras/` | las figuras, en versión clara y oscura |
| `jobs/` | el registro del equipo de agentes que hizo la Fase 1 |

## Reproducir

Todo lo numérico corre en cualquier máquina sin datos externos. Tres cosas quedan fuera
del repositorio a propósito: el árbol de SPECTER (código de terceros, se reconstruye desde
upstream con el parche), los PDF de la bibliografía (son de los editores) y los datos
crudos del PIV (140 GB en un disco externo, sólo lectura).

```bash
git clone https://github.com/Lucia-Maz/proyecto-final-piv.git
cd proyecto-final-piv

# 1. entornos
mamba create -n piv-dt -c conda-forge python=3.11 numpy scipy matplotlib \
      tifffile imageio scikit-image pip && pip install openpiv
mamba create -n specter -c conda-forge gfortran openmpi openmpi-mpifort fftw make

# 2. SPECTER, desde upstream más el parche propio
#    (la receta completa está en numerico/specter-parche/README.md)
git clone https://github.com/mfontanaar/SPECTER.git numerico/SPECTER-upstream
cd numerico/SPECTER-upstream && git checkout 0ad1edb && cd ../..
cp -r numerico/SPECTER-upstream numerico/SPECTER-trabajo
cd numerico/SPECTER-trabajo
git apply ../specter-parche/superficie-libre.patch
cp ../specter-parche/laplace_dirneu.f90 src/tests/ && cd ../..

# 3. que la verificación pase, que es la prueba de que quedó bien
python verificacion/test_aceptacion_fase1.py   # 6/6, segundos
python verificacion/test_aceptacion_fase2.py   # 6/6, ~10 min: compila SPECTER en un
                                               # scratch temporal y corre 15 simulaciones

# 4. la física (opcional: ~13 min por malla de 16², ~1 h por malla de 32²)
python verificacion/barrido_delta_etapa1.py --caso ventana --resolucion 16
python verificacion/barrido_delta_etapa1.py --caso ventana --resolucion 16 --superficie
```

`python` es el del entorno `piv-dt`. Con 5,7 GB de RAM alcanza: la puerta de la Fase 2
no pasa de unos pocos cientos de MB y el barrido de 130 MB de pico; correrlos de a uno y
con tope duro (`ulimit -v 2500000`).

**Qué corre sin el disco externo:** las dos puertas, el barrido en δ, el chequeo de la
capa límite de divergencia (`verificacion/kio_capa_divergencia.py`), y todo lo que lee de
`salidas/tablas/`, o sea todas las figuras y los informes. Los scripts `codigo/0[1-9]_*.py`
que rearman PIV desde los TIF crudos necesitan el disco montado (`piv_comun.RAIZ`); sin él
fallan con `FileNotFoundError`, a propósito, en vez de seguir con datos parciales.

## Las dos máquinas

El desarrollo, los tests 1D y la verificación a resolución chica son en la laptop; el
cómputo pesado va al clúster (Sakura, SLURM). Las convenciones que valen en las dos están
en `CLAUDE.md`; la primera vez en el clúster —clave de acceso, módulos, scratch, cómo
mandar la puerta por SLURM— está en `numerico/CLUSTER.md`. El repositorio es la única
fuente de verdad entre ambas: no se copian archivos por `scp`.

## Licencia

MIT (`LICENSE`) para todo lo propio: código, derivaciones, texto y figuras. El parche sobre
SPECTER es obra derivada de un código sin licencia publicada; ver la nota en
`numerico/specter-parche/README.md`.

## SPECTER

El código base es SPECTER, de M. Fontana, O. P. Bruno, P. D. Mininni y P. Dmitruk. Sus
autores piden citar: Fontana, Bruno, Mininni & Dmitruk, *Fourier continuation method for
incompressible fluids with boundaries*, Comp. Phys. Comm. 256, 107482 (2020),
[10.1016/j.cpc.2020.107482](https://doi.org/10.1016/j.cpc.2020.107482). Lo que este
repositorio aporta sobre SPECTER es exactamente el parche de `numerico/specter-parche/`.
