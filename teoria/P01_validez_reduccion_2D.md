# [P-01] Cuándo vale la reducción bidimensional, derivado y no heredado

**2026-09-09.** Primera cosa de la Fase 4. La pregunta abierta decía: la condición
`Re_h·ε ≪ 1` da `U ≪ 0,56 mm/s`, y las velocidades medidas son de 1,32 a 4,0 mm/s, así
que estaría *"violada por un factor de 2 a 7 justo donde arranca el decaimiento"*.

Esa cuenta se rehizo de cero. **Dos cosas cambian, y en direcciones opuestas.**

- La escala horizontal estaba mal. `ℓ = 2 cm` no correspondía a ninguna longitud de la
  celda; el espaciado real de los imanes es 5 cm y el armónico que dominan es el
  **diagonal**, de 7,07 cm. Además el criterio se escribió con `ℓ` en vez de con el
  número de onda `k`, y esa ambigüedad vale un factor 2π. Con `k` medido sobre el patrón
  real, el parámetro es **peor** de lo que decía la cuenta heredada: δ ≈ 2,7 en el
  forzado y ≈ 8 al inicio del decaimiento, no 2 a 7 veces uno.
- Pero δ no es lo que tiene que ser chico. Lo que tiene que ser chica es la **distorsión
  del perfil vertical**, y esa lleva un prefactor calculable de ≈ 0,008. En el régimen
  medido vale **2 % en el forzado y 6 % al inicio del decaimiento**.

La conclusión, entonces, es distinta de la que sugería la pregunta abierta: la clausura
de un modo **no queda invalidada** en el régimen medido. Lo que queda es un error
estimado de pocos por ciento, un factor geométrico O(1) que no se calculó, y una manera
concreta de medirlo con el código que ya existe.

Todos los números de acá salen de `numerico/fase4/p01_validez.py`; los parámetros de la
celda, de `codigo/celda.py`. La salida cruda queda en `salidas/tablas/p01_validez.json`.

---

## 1. La parte que es exacta, y conviene separarla

Sea la capa `0 < z < h`, no-deslizante en `z = 0` y libre de tensiones en `z = h`
(superficie plana: `w = 0` arriba). El promedio vertical `⟨·⟩ = h⁻¹∫₀ʰ · dz` aplicado a
la ecuación horizontal de Navier–Stokes da, **sin ninguna hipótesis sobre el perfil**,

    ∂_t⟨u⟩ + ⟨(u·∇)u⟩ = −ρ⁻¹∇_h⟨p⟩ + ν∇_h²⟨u⟩ + (ν/h)[∂_z u]₀ʰ .

El último término es la fricción, y con tope libre `∂_z u|_h = 0` exactamente, así que
se reduce a `−(ν/h)·∂_z u|₀`. De ahí sale una definición que **no depende de ninguna
clausura** y que es directamente medible en una simulación:

> **α_eff ≡ (ν/h)·∂_z u|₀ / ⟨u⟩** ,   modo horizontal por modo horizontal.

Con el perfil del modo fundamental `F₀(z) = sin(πz/2h)` se obtiene
`α_eff = π²ν/(4h²) = α`, que es el resultado de la Fase 1. **Toda la pregunta de [P-01]
es cuánto se aparta α_eff de α**, y la formulación de arriba la convierte en algo que se
mide, no que se argumenta.

Para esta celda (h = 6 mm, ν del agua): α = 0,06854 s⁻¹, τ = 14,59 s, contra los 51,2 s
que dura cada registro.

## 2. El parámetro pequeño, escrito sin ambigüedad

Adimensionalizando con `h` en la vertical, `k⁻¹` en la horizontal, `U` para la velocidad
horizontal y **el tiempo de difusión vertical `τ_ν = h²/ν`** como unidad de tiempo —esta
última elección es la que hay que justificar, y se justifica a posteriori porque deja el
operador vertical a orden cero—, la ecuación horizontal queda

    ∂_T û + δ·[(û·∇)û + ŵ ∂_Z û] = −∇p̂ + ε²∇²_h û + ∂²_Z û + f̂ ,

con exactamente dos parámetros:

| | definición | qué mide |
|---|---|---|
| **ε** | `k·h` | espesor relativo a la escala horizontal |
| **δ** | `U·k·h²/ν = Re_h·ε` | inercia frente a difusión vertical |

De continuidad, `w ~ ε·U`, y entonces `w ∂_z u` y `(u·∇)u` entran **al mismo orden**:
el ordenamiento es consistente y el término de advección vertical no se puede tirar.

Dos observaciones que la versión heredada no hacía:

**(a) El criterio se escribe con `k`, no con `ℓ`.** «La menor escala horizontal» es
ambigua por un factor 2π según se entienda longitud de onda o inverso del número de
onda. El adimensional que sale de la ecuación es `k·h`, sin ambigüedad. La cuenta
heredada usaba `ℓ = 2 cm`, es decir `k = 50 m⁻¹`; el valor correcto para esta celda es
`k = 88,86 m⁻¹` (sección 4).

**(b) `ε ≪ 1` no hace falta para la ley de fricción.** El término viscoso horizontal es
lineal y diagonal en el modo horizontal: no acopla con la estructura vertical. Se puede
conservar **exacto**, y entonces la tasa de decaimiento lineal del modo `k` es

    λ(k) = ν[k² + (π/2h)²] = α·(1 + 4ε²/π²) .

Para esta celda ε = 0,533 y ε² = 0,284, o sea que la viscosidad horizontal agrega un
**11,5 %** a la tasa de decaimiento del armónico fundamental del forzado. Eso no es una
falla de validez: es una corrección que hay que **incluir** al comparar con el
experimento, y que el modelo 2D ya incluye porque conserva `ν∇²U`. El único parámetro
genuinamente pequeño que hace falta es δ.

## 3. Qué es lo que realmente tiene que ser chico

Expandir en los modos verticales del problema con estas condiciones de contorno,

    F_m(z) = sin[(2m+1)πz/2h] ,    λ_m = ν[(2m+1)π/2h]² ,    λ_m/λ₀ = (2m+1)² ,

y escribir `u = Σ_m a_m(x,y,t)·F_m(z)`. A orden cero sólo vive `a₀`, y como
`∇_h·a₀ = 0` resulta `w = 0` **exactamente**: el modo fundamental solo no genera
velocidad vertical.

El término no lineal de ese campo es `F₀²(z)·(a₀·∇)a₀`. Su proyección sobre `F₀` es la
no linealidad 2D que el modelo ya tiene, con el factor `β = ⟨F²⟩ = π²/8` de la Fase 1.
Lo que el modelo **no** tiene es el resto, que cae sobre los modos `m ≥ 1` con los
acoplamientos

    c_m = ⟨F₀²F_m⟩/⟨F_m²⟩ ,    c₁ = −8/(15π) = −0,169765 ,

y `c_m` cae rápido: −0,0243, −0,0081, −0,0037… para m = 2, 3, 4. (El valor cerrado de
`c₁` y los de β = π²/8 y u(h)/⟨u⟩ = π/2 se reproducen por cuadratura a nueve y once
cifras respectivamente; es el control de que la base y las normalizaciones están bien.)

En balance cuasi-estacionario —lícito porque `λ_m ≫ λ₀` ya para m = 1, donde el cociente
es 9— la amplitud del modo distorsionado es

    |a_m/a₀| ≈ |c_m|·(4/π²)·δ / (2m+1)² ,

o sea **|a₁/a₀| ≈ 0,0077·δ**, dominante sobre todos los demás. Ese, y no δ, es el número
que tiene que ser chico: es la fracción del perfil vertical que se aparta del modo
fundamental supuesto.

> **Acá está el factor que no se calculó.** El paso `|(a₀·∇)a₀| ~ k·U²` lleva un factor
> geométrico O(1) que depende de la estructura horizontal del campo y de la proyección
> de la presión, y se tomó igual a 1. Es una estimación de orden, no un cálculo. Ver la
> sección 6.

## 4. El forzado: dónde está `k`, medido y no elegido

Los imanes son discos de 1 cm de diámetro en una red cuadrada de 5 cm de paso, con la
polaridad alternada como un tablero. Para el signo `(−1)^(i+j)` sobre una red de paso
`a`, el primer armónico está en `k = (π/a)(±1,±1)`: es **diagonal**, y su longitud de
onda es `a√2`, no `a`.

`p01_validez.py` no da eso por sentado: arma el patrón con discos de diámetro finito, lo
transforma y busca el pico. Resultado:

| | |
|---|---|
| k dominante | **88,858 m⁻¹**, contra π√2/a = 88,858 (error 0) |
| orientación | (−1, +1)·π/a — **diagonal**, como decía Lucía |
| longitud de onda | **7,071 cm** = 5 cm × √2 |
| ε = k·h | 0,5331 |

**Salvedad sobre este modelo, que es deliberadamente crudo.** El campo se modeló como
±1 dentro de cada disco y 0 afuera. Eso fija la **posición** de los picos exactamente,
porque es una propiedad de la red y no del imán, pero **sobrestima el contenido en
armónicos**: un imán real tiene un campo suave a la altura del fluido, y el patrón
verdadero se parece mucho más a una sinusoide que este tablero de discos. Por eso el
«11,5 % de la energía en el fundamental» que devuelve el script no debe usarse como si
fuera del experimento. El contenido espectral real del **flujo** —que es lo que importa,
porque δ ∝ k— hay que sacarlo de los campos de velocidad medidos, y eso queda pendiente
(sección 6).

## 5. Evaluación sobre la celda

El PIV mide la superficie, porque las partículas flotan, y la clausura vive sobre el
promedio vertical: hay que dividir por `u(h)/⟨u⟩ = π/2` antes de evaluar δ. Con
h = 6 mm, ν = 10⁻⁶ m²/s y k = 88,86 m⁻¹, resulta `δ = 3199 · U[m/s]`:

| u_rms superficie | ⟨u⟩ | Re_h | **δ** | **distorsión \|a₁/a₀\|** |
|---|---|---|---|---|
| 0,25 mm/s | 0,159 | 0,95 | 0,51 | 0,4 % |
| 0,50 | 0,318 | 1,91 | 1,02 | 0,8 % |
| 1,00 | 0,637 | 3,82 | 2,04 | 1,6 % |
| **1,32** (forzado 2,15 A) | 0,840 | 5,04 | **2,69** | **2,1 %** |
| 2,00 | 1,273 | 7,64 | 4,07 | 3,1 % |
| **4,00** (inicio de decaimiento) | 2,546 | 15,28 | **8,15** | **6,2 %** |
| 6,00 | 3,820 | 22,92 | 12,22 | 9,3 % |

Dos condiciones más, y las dos se cumplen con holgura:

- **Memoria modal.** Para que el perfil *sea* el fundamental hay que esperar
  `1/(λ₁−λ₀) = h²/(2π²ν) = 1,82 s`. Los registros duran 51,2 s: después de unos 5 s el
  segundo modo vertical está 15 veces atenuado. La clausura no vale en el primer
  segundo del decaimiento, y sí en todo el resto.
- **Decaimiento de δ.** δ ∝ U, así que `δ(t) = δ₀·e^{−αt}`. Aun arrancando de δ₀ = 8, a
  los 30 s vale 1,0 y la distorsión baja al 0,8 %. La clausura **mejora sola** a lo
  largo del registro; el peor momento es el primero.

**La velocidad de 4,0 mm/s no está verificada.** Figura en `ESTADO.md` sin script que la
produzca y no se pudo rastrear a ninguna salida de este árbol. El disco externo no está
montado, así que no se pudo recalcular. Queda anotada como tal en `codigo/celda.py`.

## 6. Lo que no quedó resuelto

1. **El factor geométrico O(1) de la sección 3.** La estimación de la distorsión toma
   `|(a₀·∇)a₀| ~ kU²` con prefactor 1. Puede ser bastante menor: para el armónico
   fundamental del tablero, `ψ ∝ cos(k_x x)cos(k_y y)` es autofunción del laplaciano,
   de modo que `J(ψ, ∇²ψ) = 0` y la parte solenoidal del término no lineal **se anula
   idénticamente**. Es decir que el patrón ideal de los imanes maneja un flujo que, en
   su fundamental, no se distorsiona a este orden. Un campo de banda ancha —que es lo
   que hay en un decaimiento real— no tiene esa protección. No sé cuál de los dos
   describe mejor los registros, y no lo afirmo.
2. **A qué orden en δ se sesga α_eff.** Los modos distorsionados viven en números de
   onda horizontales 0 y 2k, no en k, así que para un único modo celular no realimentan
   la tasa del modo k hasta orden δ². Para un campo de banda ancha hay tríadas que sí lo
   hacen a orden δ. No está decidido analíticamente.
3. **La menor escala presente en el flujo.** δ ∝ k, así que el criterio lo fija la
   escala más chica con energía apreciable, no el fundamental del forzado. Hay que
   sacarla del espectro o de la autocorrelación de los campos PIV medidos. Requiere el
   disco externo montado.
4. **La viscosidad.** Se usa la del agua y no la del electrolito (KNO₃ al 16 % m/m), por
   decisión explícita de Lucía del 2026-09-09, declarada en `codigo/celda.py`. Entra a
   primer orden: `δα/α = δν/ν − 2δh/h`.
5. **La superficie libre real.** Sigue en pie la salvedad más fuerte del trabajo: una
   película superficial vuelve el tope casi rígido, y ese es el caso con α cuatro veces
   mayor. Es independiente de todo lo de acá.

## 7. Cómo se cierra: medir α_eff, que ahora se puede

La sección 1 deja la pregunta en una forma medible, y SPECTER con la condición de la
Fase 2 puede producirla. El diseño mínimo:

- Caja delgada, periódica en x e y, con `freeslip` arriba y `noslip` abajo, sin forzado.
- Condición inicial `u = F₀(z)·∇^⊥ψ(x,y)` con `ψ` el patrón celular del forzado, que ya
  cumple `w = 0` y divergencia horizontal nula.
- Barrido en la amplitud inicial para recorrer δ ≈ 0,3 … 10, que es el rango de la
  tabla de la sección 5.
- Se mide `α_eff = (ν/h)·∂_z u|₀ / ⟨u⟩` y se lo compara con `λ(k) = α(1+4ε²/π²)`. La
  curva `α_eff/λ(k)` contra δ **es** la respuesta a [P-01], y de paso da el factor
  geométrico del punto 1 y el orden del punto 2.
- Hay que hacerlo con dos resoluciones horizontales: [H-09] mostró que con superficie
  libre el término no lineal se desestabiliza donde el canal aguanta, así que la malla
  no se hereda. Este barrido es también el estudio de convergencia que la corrida de
  producción necesita.

## 8. Consecuencia para el resto de la Fase 4

- Al comparar con el decaimiento medido hay que contrastar contra **λ(k) = α(1+4ε²/π²)**
  y no contra α pelado: son 11,5 % en el fundamental, muy por encima de la barra de
  error del capítulo de instrumento.
- El primer segundo de cada registro está fuera de la clausura por memoria modal, y el
  arranque es además donde δ es mayor. Ajustar sobre el registro completo y mirar si el
  residuo tiene tendencia sistemática en los primeros ~10 s.
- Al pasar de lo medido al modelo, dividir por π/2. El sesgo de superficie es del 57 % y
  entra en δ y en cualquier comparación de amplitudes.
