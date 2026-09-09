# La superficie libre en las condiciones de v\* y de p

**2026-09-09.** Descripción detallada de cómo se modela una superficie libre dentro del
esquema de partición de SPECTER: qué condición cumple el campo auxiliar `v*`, qué
condición cumple la presión, y por qué las dos juntas producen la condición física sobre
el campo proyectado.

Tres casos, en orden de lo que está hecho a lo que no:

| caso | estado | condición sobre p | condición sobre v\* |
|---|---|---|---|
| **plana, libre de tensiones** | **implementado y verificado** ([V-09]) | Neumann, igual que el canal | vorticidad tangencial nula, **sin dt** |
| **con película superficial** | no implementado; **es el que le falta a este experimento** | Neumann, sin cambios | Robin en las componentes tangenciales |
| **deformable** | no implementado; **no hace falta acá**, ver §4.1 | **Dirichlet arriba** — rama que no existe | acoplada a ∂p/∂z, **con dt** |

Todo lo que se afirma del código está verificado leyéndolo, y las cuentas están en
`numerico/fase4/superficie_libre_escalas.py`.

---

## 1. Por qué la pregunta se plantea sobre v\* y sobre p, y no sobre v

SPECTER no avanza la velocidad directamente. En cada subpaso de Runge-Kutta de orden `o`
([V-08], secciones 4.1 y 4.4 de Fontana, Bruno, Mininni y Dmitruk 2020):

1. se avanza un campo **sin presión** `v*`, con todos los demás términos;
2. se le impone a `v*` una condición de contorno **propia**;
3. se resuelve la presión, y
4. se proyecta

       v = v* − (dt/o)·∇p ,

   donde además —y esto no está en el paper, se verificó en el código— la variable
   `pr` no guarda `p` sino `dt·p` del subpaso en que se calculó
   (`boundary_mod.fpp:210`, el `!TODO update documentation about p' = dt*p`).

De modo que «modelar la superficie libre» son tres decisiones acopladas:

- **(a)** qué condición se le pide a la presión en la cara;
- **(b)** qué condición se le pide a `v*` en la cara;
- **(c)** verificar que (a) y (b) juntas den, **sobre el campo ya proyectado**, la
  condición física que se quería.

El paso (c) es el que se puede hacer mal en silencio: una condición correcta sobre `v*`
con la presión equivocada da un campo proyectado que no cumple nada. Y hay una advertencia
del propio paper que se hereda: proyectar la ecuación de momento en la dirección normal o
en la tangencial da problemas de contorno distintos para la presión, y *«there is no a
priori reason to assume that both approaches lead to the same solution»*. Ellos eligen la
normal; este trabajo también.

## 2. El caso implementado: superficie plana y libre de tensiones

### 2.1 Las condiciones físicas

Superficie **plana** en `z = h`, sin deformación:

    w = 0 ,        ∂u/∂z = ∂v/∂z = 0        en z = h .

La primera es impermeabilidad; las otras dos son tensión tangencial nula. Nótese que como
`w ≡ 0` **sobre todo el plano**, sus derivadas tangenciales `∂w/∂x` y `∂w/∂y` se anulan
idénticamente ahí, y eso es lo que hace especial a este caso (§2.4).

### 2.2 La presión: Neumann, exactamente la misma que el canal

La condición de impermeabilidad se impone pidiendo `v_z = 0` **después** de proyectar. Con
`v_z = v*_z − (dt/o)·∂p/∂z`, eso equivale al dato de Neumann

    ∂p/∂z |_cara = (o/dt)·v*_z |_cara .

Verificado en el código, no supuesto: `v_imposebc_and_project` llama siempre
`sol_project(vx,vy,vz,pr,bctarget=1,0,0)` (`vboundary.f90:193`), y `sol_project` le pasa a
`laplace_z` la clase `bcz + bctarget` (`boundary_mod.fpp:357`). Con `bctarget=1` —«la
condición la tiene que cumplir la componente normal después de proyectar»— las dos clases
salen **1**, es decir **Neumann–Neumann sobre la presión**, idéntico al canal.

Esto es [H-08], que corrigió a [V-06]: la rama Dirichlet–Neumann de `laplace_z` **no** es
la que usa la superficie libre plana. La superficie plana no toca el camino de la presión
en absoluto.

### 2.3 El campo auxiliar: vorticidad tangencial nula, y sin dt

Pedir la condición física sobre el campo **proyectado** y sustituir la proyección:

    ∂v_x/∂z = ∂v*_x/∂z − (dt/o)·i·k_x·∂p/∂z ,

y usar el dato de Neumann de arriba, `∂p/∂z = (o/dt)·v*_z`, hace que el `dt` se cancele
**exactamente**:

    ∂v*_x/∂z = i·k_x·v*_z ,     ∂v*_y/∂z = i·k_y·v*_z       en la cara.

Es decir: se anulan las dos componentes tangenciales de la **vorticidad** de `v*`. Así
está escrito en `freeslip_z`, con el signo de la cara (`sgn = −1` abajo, `+1` arriba)
porque `neumann_reconstruct` trabaja con la derivada normal **saliente** —convención
verificada en `src/tests/fc_neumann.f90`, el comentario `d/dz = -d/dn`.

### 2.4 Por qué esta condición no tiene error de partición, y qué la hace especial

La componente `y` del tensor de tensiones tangenciales es `∂u/∂z + ∂w/∂x`, mientras que la
componente `y` de la vorticidad es `ω_y = ∂u/∂z − ∂w/∂x`. **Coinciden sólo cuando
`∂w/∂x = 0`**, que es justamente lo que garantiza una superficie plana.

Y como la proyección resta un gradiente, y el rotor de un gradiente es cero,

    rot(v*) = rot(v)     exactamente,

de modo que una condición escrita sobre la vorticidad se transfiere al campo proyectado
**sin error de partición**, a diferencia del no-deslizamiento, cuya velocidad de
deslizamiento residual es O(dt^ord). Esa es la predicción falsable de [D-24], y es lo que
mide V3 de la puerta: el residuo de tensión es 6,1·10⁻¹² relativo a la cara rígida y la
razón entre `dt` y `dt/2` es **1,00**.

**Retener esto**, porque es lo que se pierde en el caso deformable: la superficie plana es
exacta en el esquema de partición *porque su condición de contorno resulta ser una
condición sobre la vorticidad*.

## 3. La película superficial: el caso que este experimento sí necesita

Es la salvedad más fuerte del proyecto: un electrolito con partículas flotando puede
desarrollar una película que vuelve la superficie casi rígida, y ese es el caso
no-deslizante, con α **cuatro veces mayor**.

### 3.1 La condición física

El balance de tensión tangencial en una interfaz con una película deja de ser «cero» y
pasa a igualar la divergencia superficial de las tensiones de la película. En su forma más
simple —modelo de Boussinesq-Scriven con una única viscosidad superficial `μ_s`, sin
gradientes de tensión superficial (sin Marangoni)— y para un modo horizontal `k`:

    μ·∂u_par/∂z |_h = −μ_s·k²·u_par |_h .

Es una condición de **Robin** sobre las componentes tangenciales, con un único parámetro
adimensional

    β = μ_s·k²·h / μ ,

que interpola entre los dos casos ya conocidos: `β → 0` es el tope libre y `β → ∞` es el
canal.

### 3.2 Qué le pasa a α, y por qué esto convierte la salvedad en una medición

Con `u = sin(qz)`, que ya cumple `u(0) = 0`, la condición de arriba da

    (q·h)·cot(q·h) = −β ,        α = ν·q² ,

y `q·h` va de `π/2` a `π`, o sea α de `π²ν/4h²` a `π²ν/h²`. El factor 4 de la Fase 1 sale
como los dos extremos de esta familia. Con h = 6 mm, ν = 10⁻⁶ m²/s y k = 88,86 m⁻¹:

| β | q·h/π | α [1/s] | α/α_libre | μ_s [N·s/m] |
|---|---|---|---|---|
| 0 | 0,5000 | 0,06854 | 1,000 | 0 |
| 0,03 | 0,5060 | 0,07020 | 1,024 | 6,3·10⁻⁷ |
| **0,127** | — | 0,07539 | **1,10** | **2,7·10⁻⁶** |
| 0,3 | 0,5543 | 0,08424 | 1,229 | 6,3·10⁻⁶ |
| **1,69** | — | 0,13708 | **2,00** | **3,6·10⁻⁵** |
| 10 | 0,9112 | 0,22765 | 3,321 | 2,1·10⁻⁴ |
| 1000 | 0,9990 | 0,27361 | 3,992 | 2,1·10⁻² |

Leído al revés: **el α medido en el decaimiento se puede invertir para acotar `μ_s`**. Eso
transforma la salvedad —que hoy se declara y no se acota— en un resultado con número.

> **Lo que no se puede decir todavía.** Si `μ_s = 2,7·10⁻⁶ N·s/m` es un valor grande o
> chico para una interfaz agua–aire con partículas de 100 µm no lo sé: no se leyó
> bibliografía de reología interfacial y no se midió. La tabla dice qué `μ_s` implica cada
> α; para decir si ese `μ_s` es plausible hace falta una fuente y todavía no la hay.

### 3.3 Qué cambia en las condiciones de v\* y de p — muy poco

**La presión no cambia en absoluto.** La superficie sigue siendo plana, `w = 0` sigue
valiendo, y por lo tanto el camino Neumann–Neumann de §2.2 queda intacto.

**Sobre `v*`**: repitiendo la cuenta de §2.3 con el miembro derecho no nulo,

    ∂v*_x/∂z + (μ_s k²/μ)·v*_x = i·k_x·v*_z + (μ_s k²/μ)·(dt/o)·i·k_x·p |_h .

Aparece `p` en la cara —no su derivada—, así que **el `dt` ya no se cancela solo**: queda
un término O(dt) proporcional a `β`. En el límite `β → 0` se recupera exactamente §2.3.
Dos consecuencias honestas: (i) el caso con película **sí** tiene error de partición, de
orden `dt·β`, y (ii) por eso mismo el diagnóstico de V3 —residuo independiente de `dt`—
deja de ser el chequeo correcto y hay que reemplazarlo por uno que mida el residuo de
Robin y verifique que escala como `dt^ord`.

La maquinaria de reconstrucción ya existe: `src/tests/` tiene `fc_robin.f90`, así que
FC-Gram sabe imponer Robin.

## 4. La superficie deformable

### 4.1 Primero: en esta celda no hace falta, y es un número

La presión dinámica de un flujo de velocidad `U` es del orden de `ρU²`, y la superficie se
restituye con `ρg·η` (gravedad) más `σk²·η` (capilaridad), de donde
`η ~ ρU²/(ρg + σk²)`.

| | U (promedio vertical) | Froude | η | η/h | restituye |
|---|---|---|---|---|---|
| forzado 2,15 A | 0,841 mm/s | 0,0035 | **0,07 µm** | 1,1·10⁻⁵ | gravedad |
| inicio del decaimiento | 2,311 mm/s | 0,0095 | **0,51 µm** | 8,6·10⁻⁵ | gravedad |

La longitud capilar es 2,71 mm, o sea `k_capilar = 369 m⁻¹`, bastante mayor que los
88,9 m⁻¹ del forzado: **manda la gravedad**. Y la superficie se deforma menos de un micrón
sobre una capa de 6 mm. La aproximación de superficie plana es buena a una parte en 10⁵, y
eso justifica con número —no con criterio— la decisión de alcance que ya estaba tomada.

*(Estos números usan ρ y σ del agua limpia. El electrolito es más denso y una interfaz con
partículas tiene menos tensión; ninguna de las dos correcciones puede mover el resultado
cuatro órdenes de magnitud.)*

### 4.2 Las condiciones físicas, linealizadas

Superficie en `z = h + η(x,y,t)`, con `η ≪ h` y pendiente chica. A primer orden, evaluando
en `z = h`:

- **Cinemática:** `w(h) = ∂η/∂t` — ya no es `w = 0`. Este es el cambio de fondo.
- **Tensión tangencial:** `∂u/∂z + ∂w/∂x = 0` y `∂v/∂z + ∂w/∂y = 0`. La *forma* es la misma
  que en el caso plano, pero ahora el segundo término **no** se anula.
- **Tensión normal:** con `p` la desviación de la hidrostática,

      p(h) = ρ·g·η − σ·∇²_h η + 2μ·∂w/∂z |_h .

  El signo del término capilar depende de la convención de curvatura y **hay que fijarlo
  con un test**, no de memoria: ver §5.

### 4.3 La presión pasa a ser Dirichlet arriba, y esa rama no existe

La condición de tensión normal **prescribe el valor de p en la superficie**. El fondo
rígido sigue dando Neumann. De modo que el problema de la presión pasa a ser

    Neumann en z = 0 ,   Dirichlet en z = Lz ,

que es exactamente el par que `laplace_z` **no** tiene. La rama que se implementó en la
Fase 2 es la **opuesta**, Dirichlet abajo – Neumann arriba ([H-08]).

Los coeficientes, en la misma base que usa el archivo, `φ(z) = C₁e^{k(z−Lz)} + C₂e^{−kz}`
con `r = e^{−k·Lz}`, y datos `g₀ = φ′(0)`, `g₁ = φ(Lz)`:

    C₁ = (g₁ + r·g₀/k)/(1 + r²) ,     C₂ = (r·g₁ − g₀/k)/(1 + r²) .

Verificado sobre 20 000 combinaciones aleatorias de `k` y `Lz` en seis órdenes de magnitud:
el residuo de Neumann en `z=0` es 4,4·10⁻¹⁶ y el de Dirichlet en `z=Lz` es 3,9·10⁻¹¹ en el
peor caso. El denominador es **1 + r²**, igual que la rama ya implementada: no tiene
resonancia, a diferencia de las ramas Dirichlet/Dirichlet y Neumann/Neumann, que llevan
`1 − r²` ([D-25]).

El parámetro de condicionamiento es el **espejo** del de [D-22]: allá era `|g₁|/(k·|g₀|)`,
acá es `|g₀|/(k·|g₁|)`. Verificado con el mismo contraejemplo: con el mismo `k·Lz = 0,5` el
error de Dirichlet es 1,9·10⁻⁶ o exactamente 0 según cuánto valga `g₀`, y con un cociente
igual de grande pero fuera de la rama exponencial vuelve a ser 0. El guardián hay que
escribirlo sobre `Lz·|g₀|/|g₁|`.

**Y hay un problema de interfaz, no sólo de fórmula.** `sol_project` tiene un **único**
`bctarget` para las dos caras, y acá hacen falta distintos: abajo la condición la cumple
`v_z` después de proyectar (`bctarget = 1`, que sumado a `bczsta = 0` da clase 1 =
Neumann), y arriba la condición es sobre la presión misma (`bctarget = 0`, clase 0 =
Dirichlet). Hay que abrir la firma a un `bctarget` por cara.

Dos cosas mejoran, de paso: el problema Neumann–Neumann actual determina `p` sólo a menos
de una constante y exige una condición de compatibilidad entre el flujo por las dos caras;
el par Neumann–Dirichlet es **unívocamente soluble** y no la exige.

### 4.4 La condición sobre v\* se acopla a la presión, y pierde el dt

Ahora hay que imponer la condición completa `∂v_x/∂z + ∂v_z/∂x = 0` sobre el campo
proyectado. Sustituyendo `v = v* − (dt/o)∇p` en los dos términos:

    ∂v*_x/∂z + i·k_x·v*_z − 2·(dt/o)·i·k_x·∂p/∂z |_h = 0 ,

o sea

    ∂v*_x/∂z = − i·k_x·v*_z + 2·(dt/o)·i·k_x·∂p/∂z |_h .

Se reduce correctamente al caso plano: si además se impone `v_z = 0` en la cara, entonces
`i k_x v*_z − (dt/o) i k_x ∂p/∂z = 0`, y lo que queda es `∂v*_x/∂z = i k_x v*_z`, que es
§2.3.

**Las dos dificultades, dichas sin adornos:**

1. **`∂p/∂z|_h` ya no es un dato.** Con Dirichlet arriba, la derivada normal de la presión
   en la cara es parte de la *solución* del problema de Poisson, no del dato de contorno.
   La condición sobre `v*` y el cálculo de `p` quedan **implícitamente acoplados**. Las dos
   salidas son iterar dentro de cada subpaso (resolver `p`, reimponer la condición sobre
   `v*`, reproyectar, hasta converger) o retrasarla un subpaso, que mete un error de
   partición O(dt).
2. **Se pierde la exactitud de §2.4.** Con superficie deformable `∂w/∂x ≠ 0`, así que la
   tensión tangencial deja de coincidir con la vorticidad tangencial, y el argumento
   `rot(v*) = rot(v)` ya no aplica. La propiedad que hace elegante al caso plano —residuo
   independiente de `dt`— **no sobrevive**, ni siquiera con la iteración.

### 4.5 El campo η y el dominio

`η(x,y,t)` es un campo 2D, periódico en x e y, así que vive en Fourier puro y no necesita
FC-Gram. Se avanza con el mismo Runge-Kutta, `∂η/∂t = w(h) − u_h·∇_h η`.

El dominio **no se mueve**: las condiciones se aplican en `z = h` con `η` como campo
aparte. Es la formulación linealizada, y es la coherente con `η/h ~ 10⁻⁵`. Una coordenada
mapeada tipo σ, que seguiría la superficie, cambia la estructura de todos los operadores en
z y es un proyecto de otra escala.

## 5. Cómo se verificaría cada uno

Siguiendo el criterio de las Fases 1 y 2: la puerta se escribe **antes**, en rojo, y no
contiene la fórmula del resultado.

- **Película:** que α(β) reproduzca los dos límites conocidos (`β→0` da `π²ν/4h²`, `β→∞` da
  `π²ν/h²`, cociente 4) corriendo el mismo binario, y que el residuo de la condición de
  Robin escale como `dt^ord` — que es lo contrario de lo que se verifica hoy en V3, y por
  eso hay que cambiar el chequeo, no reusarlo.
- **Deformable:** la relación de dispersión de ondas de gravedad-capilaridad,
  `ω² = (g·k + σ·k³/ρ)·tanh(k·h)`, es **exacta, lineal y sin parámetros libres**. Es el
  test correcto y además fija el signo del término capilar de §4.2: con el signo cambiado
  las ondas cortas crecen en vez de oscilar, y el test lo detecta enseguida.
- **Rama Neumann–Dirichlet de `laplace_z`:** test unitario igual que `laplace_dirneu.f90`,
  contra la forma cerrada `cosh/sinh` y contra sus propiedades definitorias, con cuidado de
  darle datos **reales** al modo (0,0) — que fue una de las tres fallas de la primera
  corrida de la puerta de la Fase 2.

## 6. Salvedades

- El signo del término de curvatura de §4.2 no está fijado; se fija con el test de §5.
- El término viscoso normal `2μ·∂w/∂z` de la condición de tensión normal se suele
  despreciar; acá no se evaluó si se puede.
- El costo de la iteración de §4.4 no se estimó.
- Se hereda la advertencia del paper sobre proyectar en la normal o en la tangencial.
- ρ y σ de §4.1 son de agua limpia, no del electrolito ni de una interfaz con partículas.
- Nada de §3 y §4 está implementado. Esto es una descripción de diseño.
