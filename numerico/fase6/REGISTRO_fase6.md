# Fase 6 — superficie deformable de orden N en SPECTER: registro de reproducción

**Fecha:** 2026-09-28. **Pedido de Lucía, textual:** *«Try to design a further implementation
of the deformable surface on SPECTER, allowing for higher order deformations, to the best of
your abilities. Please remember to verify key steps, and keep a separated record of your
decisions, problems, bibliography consulted, etc.; any info you deem important for the
reproduction of this task.»*

Este archivo es ese registro separado. Las entradas de `DECISIONES.md` ([D-51], [D-52],
[H-23]–[H-26], [V-21], [P-04]) lo resumen; acá está el detalle que hace falta para rehacer
el trabajo sin la conversación. Convención: **(§)** marca una derivación, estimación o
hipótesis de esta sesión, no validada contra bibliografía.

---

## 1. Resumen

**Qué se construyó.** Un camino general en SPECTER (`src/boundary/fsorder.f90`, nuevo) que
impone en la frontera fija z = Lz las condiciones exactas de una superficie gráfica
z = h + η(x,y,t) **desarrolladas en serie de Taylor en η y truncadas a orden N en la
amplitud** (N = 1, 2, 3; los términos de velocidad llevan hasta η^(N−1)), con el volumen de
Navier–Stokes completo. N = 1 reproduce la superficie lineal de la Fase 5
([D-49]) a 2,5·10⁻¹⁵. Se activa con `fsorder = 2` o `3` en `&freesurface`.

**Qué funciona.**
- **N = 2.** Estable con η > −0,667 Δz en todas partes, al ν dt/Δz² = 0,104 de la puerta; el
  límite depende de ese número (P3). Para cualquier η > 0, a cualquier η/Δz probado (hasta 9),
  es estable y la iteración converge a 10⁻¹⁰ en ≤ 4–7 pasadas.
- **N = 3.** Estable para todo η probado (−4 Δz a +9 Δz) con ν dt/Δz² ≤ 0,19 (P14). Pasadas
  ≤ 7–17 según el caso; residuo 10⁻¹⁰–10⁻⁸.
- Las funcionales de orden N tienen el orden correcto contra la geometría exacta (H3):
  1, 2 y 3, y 2, 2 y 4 para la curvatura.

**Qué se encontró que no funciona, y por qué.** Cuatro hallazgos, todos con evidencia
reproducible (§4):
1. **[H-23]** El cierre que congelaba el intent es **inestable en el tiempo**. Extrapolaba
   las derivadas normales de orden ≥ 2 desde nodos interiores, y es inestable para
   |η| > 0,5 Δz: radio espectral del paso 2,8 en η = Δz y 39 en η = 4 Δz. Se reemplazó por
   trazas espectrales resueltas de forma implícita, con un precondicionador exacto del
   auto-acople.
2. **[H-24]** **El modelo de orden 2 está mal planteado donde η < 0**, ya en el continuo
   **(§)**. La condición transferida ∂z u + η ∂zz u = 0 admite una capa espuria
   exp((z − Lz)/|η|) que crece como ν/η²; en el problema 2D linealizado completo la raíz es
   exactamente λ = ν(1/η² − k²) (verificador).
   - En la malla, el paso explícito se vuelve inestable donde 1 + σ baja a ~0,7 ν dt/Δz². Eso
     es η = −0,667 Δz con los parámetros de la puerta; el cero de 1 + σ, −0,716 Δz, es sólo el
     límite dt → 0.
   - A N = 3 los modos espurios son neutros para η uniforme, pero con pendiente pueden crecer
     débilmente, Re λ ≈ ν(2s²/η² − k²). Esto lo encontró el verificador y no se reprodujo aparte.
3. **[H-25]** **Para una elevación media uniforme, N = 2 ya tiene error O(η³)**. Los términos
   O(η²) se cancelan exactamente **(§)**. La escalera de profundidad congelada (H5) esperaba
   orden 2 en N = 2: el criterio estaba mal a priori.
4. **[H-26]** **La traza espectral de ∂z³u limita a N = 3.**
   - Su error FC-Gram tiene un mínimo cerca de nz = 128 y después **crece** como Δz⁻³: un piso
     determinista de la continuación, ~2·10⁻¹¹·|f|/Δz³. Lo encontró el verificador y se
     reprodujo (P14).
   - Con estas tablas, N = 3 no mejora refinando en z.
   - Además, la iteración tiene a N = 3 un piso de redondeo de unos 10⁻⁹ con η = 4 Δz, de
     origen cuantitativo no establecido.
   - A nz = 128, N = 3 no mejora a N = 2 en la escalera de profundidad. En la onda no lineal,
     donde no hay degeneración, sí tiene orden 3.

**Puerta congelada: 4/10** (H1, H3, H4, H6). Los seis fallos se explican en §6: dos por
criterios fijados mal a priori, tres porque piden a N = 2 el régimen donde está mal planteado,
y uno por el piso de redondeo de N = 3. Regresión: Fase 2 6/6 y Fase 5 (§6).

**Para la aplicación física no aplica ninguno de los límites**: en la celda η/h ~ 10⁻⁵ y, a
nz = 128, η/Δz ~ 10⁻³. Ahí 1 + σ ≈ 1, la iteración converge en 2–3 pasadas, y los errores
de traza entran multiplicados por η ~ 10⁻⁵.

---

## 2. El modelo — [D-51]

Con Tr_M[f] = Σ_{m=0}^{M} η^m/m! ∂z^m f (evaluado en z = Lz; 0 si M < 0), s = ∇_h η,
S = ∇v + ∇vᵀ y p la presión cinemática menos la hidrostática del reposo:

```
cinemática   ∂η/∂t = Tr_{N-1}[w] − s·Tr_{N-2}[u_h]
tangencial   Tr_{N-1}[S_az] − s_b Tr_{N-2}[S_ab] + s_a Tr_{N-2}[S_zz] − s_a s_b Tr_{N-3}[S_zb] = 0
normal       Tr_{N-1}[p] = g η − γ κ_N + ν B_N
             κ_N = ∇·( s Σ_{2j+1≤N} binom(−1/2, j) |s|^{2j} )
             B_N = Σ_{2j≤N−1} (−|s|²)^j ( Tr_{N−1−2j}[S_zz] − 2 s_a Tr_{N−2−2j}[S_az]
                                            + s_a s_b Tr_{N−3−2j}[S_ab] )
```

**De dónde sale.** Son los polinomios de Taylor de las condiciones exactas de Wang, Tice &
Kim (2014), DOI 10.1007/s00205-013-0700-2, §1.1, ecs. (1.4)–(1.8), leídas en el original
(arXiv:1109.1798v2) en la Fase 5. El desarrollo y el truncamiento son de esta sesión **(§)**.

**Verificado a redondeo** con `numerico/fase6/chequeo_expansion.py`, sin SPECTER. Compara
los coeficientes de Taylor de las condiciones exactas, extraídos por integrales de Cauchy
numéricas, contra las fórmulas truncadas:
- para N = 1..5, peor error relativo entre 3,4·10⁻¹⁵ y 1,05·10⁻¹⁴;
- control negativo con el signo de la pendiente cambiado: 1,8.

Salida en `salidas/tablas/fase6_chequeo_expansion.json`.

**Alcance.** Es un desarrollo en la **deformación**, con el campo de velocidad completo. No
es un desarrollo en la amplitud del flujo. Exige η chico frente a las escalas verticales del
flujo, incluida la capa de vorticidad de la superficie. [H-24] agrega una condición que no
se sabía al congelar: a N = 2, η > −0,67 Δz en todas partes (con ν dt/Δz² = 0,1; P3).

---

## 3. El esquema — [D-52] (versión final)

Por subpaso `o` del Runge–Kutta de SPECTER (`fs_general_imposebc`, `fsorder.f90:568`):

1. **Explícito, del subpaso anterior:** la cinemática (η^(o) = η⁰ + (dt/o)·K) y la parte de la
   tensión normal sin presión, E = gη − γκ_N + νB_N. Esta da el dato Dirichlet del potencial,
   d_top = (dt/o)E − P. Es igual que en la Fase 5.
2. **Implícito, iterado:** el estado nuevo tiene que cumplir la condición tangencial
   completa de orden N y la transferencia de presión P = Σ_{m≥1} η^m/m! ∂z^m d. Incógnitas
   de superficie: x = (W, Dx, Dy, P).
   - W = w(Lz) es exacto con B_k, como en la Fase 5.
   - D = Tr_{N−1}[S_az] + … − S_az(m = 0) son las correcciones de orden ≥ 2 de la tensión.
   - El dato de Neumann de u* es ikₓ(w* − 2W) − Dx.
3. **Trazas** (`fs_traces_general`, `:146`):
   - ∂z^m u, ∂z^m v, m = 0..N, por **sumas espectrales** sobre la columna reconstruida;
   - ∂z^m w por continuidad;
   - las del potencial, exactas: la parte armónica analítica (d'' = k²d) y la no homogénea
     por sumas espectrales.
   - El término de Nyquist en kz se omite en las derivadas impares (`:114`), por la
     convención del interpolante simétrico (§4, P8).
4. **Precondicionador del auto-acople** (`:735-750`). La traza espectral de ∂z^m u
   responde al dato de Neumann como β_m por unidad de dato. β_m se mide en `fs_gen_setup`
   (`:133`) con una sonda: columna nula más dato unitario. En unidades de Δz^(m−1) los
   valores son β = (0,4800; 1,0052; 1,3967; 1,1515), iguales a nz = 64, 128 y 256.
   - Entonces Dx responde a sí mismo con el factor puntual σ = β₂η + β₃η²/2, este último
     término sólo si N ≥ 3. Medido en SPECTER: ganancia −15,28 = −σ en cada modo, con
     acoples cruzados ≤ 0,018.
   - La iteración usa G = (F − 2σ∂ₓ(W_nuevo − W) + σ·Dx)/(1 + σ), en el espacio físico
     acolchado. G tiene los mismos puntos fijos que F mientras 1 + σ > 0.
   - Se acelera con Anderson tipo II, memoria 5 y **coeficientes reales** (§4, P6).
5. **Criterio** (`:694-729`). Cada bloque se mide contra el dato que corrige:
   - W, contra sí mismo;
   - Dx y Dy, contra el dato de Neumann completo;
   - P, contra el dato Dirichlet sin la media, porque la media es la estática g⟨η⟩.

   Se para en `fstol`, o por estancamiento. Estancamiento es: desde la sexta pasada, el
   residuo actual ≤ 10²·fstol y el mejor residuo de las últimas cuatro pasadas no llega a la
   mitad del mejor anterior (§4, P4 y P13).
6. **Guardias** (`:642-670`). La corrida aborta con un mensaje explícito en tres casos
   ([H-24], P14):
   - a N = 2, si min(1 + σ) ≤ 0,05 + ν dt/Δz²; con los parámetros de la puerta eso es
     η ≤ −0,607 Δz, antes del umbral −0,667 Δz;
   - a N = 3, si min(1 + σ) ≤ 0,05, que no se alcanza nunca porque min(1 + σ) = 0,153;
   - a N = 3, si ν dt/Δz² > 0,19. El límite explícito es 0,205, y a partir de 0,99 de él los
     modos cerca de η = −1,2 Δz son inestables.

Namelist, diagnóstico `freesurface_order_diagnostic.txt` (13 columnas) y tests: los que
congela `verificacion/intent_fase6.txt`, sin cambios. **Lo que sí cambió respecto del
intent** es el cierre de las trazas y el criterio de parada. El intent es un archivo
congelado y no se editó; la desviación está registrada en [D-52] con su motivo ([H-23]).

---

## 4. Problemas encontrados, en orden, y qué se hizo

**P1 — NaN en la primera corrida N = 2 de la escalera de profundidad** (η0 = 0,04 ≈ 4,1 Δz).
- Una copia instrumentada, sólo en el scratch, mostró que dentro de cada subpaso la
  iteración converge geométricamente.
- Pero de subpaso a subpaso las trazas crecían ×3–10: Dx0 = 3,4·10⁻⁹ → 5,4·10⁻⁷ → … →
  7,4·10⁻⁵.
- Hipótesis: el dato ∂z u = −η·∂zz u (extrapolada) convierte al nodo de borde en una
  extrapolación con pesos ~(η/Δz)·324/2 ≈ 640 sobre el interior.

**P2 — confirmación y rediseño.** `numerico/fase6/estabilidad_cierre.py`: réplica 1D del
paso de SPECTER (difusión FC-Gram, neumann_reconstruct, RK2, ν = 0,01, dt = 10⁻³, nz = 128).
Radio espectral del paso completo; el piso de estabilidad es 0,999975:

| cierre | N = 2 inestable en η/Δz ∈ | N = 3 inestable en η/Δz ∈ |
|---|---|---|
| extrapolación interior, q = 8 (el del intent) | fuera de (−0,5, +0,5) | fuera de (−1,9, +0,4) |
| extrapolación interior, q = 5 | ≤ −1,1 y ≥ 3,3 | ≥ 1,8 |
| extrapolación interior, q = 3 | −4 (y es de primer orden) | −4 |
| **traza espectral implícita (final)** | **≤ −0,667 con D = 0,104** (depende de D, P3) | **ninguno** con D ≤ 0,97 D_max (P14) |

- Figura: `salidas/figuras/f_fs6_estabilidad.png`; tabla:
  `salidas/tablas/fase6_estabilidad_cierre.json`.
- En SPECTER, con la versión final, η0 = +0,04 corre hasta t = 4 sin crecer, en 5–7
  pasadas por subpaso. Con la versión anterior daba NaN en cinco subpasos.

**P3 — el régimen η < 0 de N = 2 ([H-24]).** Análisis local de modos **(§)**:
- u ~ e^{λt} e^{κ(z−Lz)} con ∂z u + η ∂zz u = 0 da κ(1 + ηκ) = 0. La raíz κ = −1/η decae
  hacia el interior cuando η < 0, con λ = νκ² = ν/η² > 0.
- A N = 3: 1 + ηκ + (ηκ)²/2 = 0 da ηκ = −1 ± i y λ = ∓2iν/η², neutro.
- Es el polinomio de Taylor truncado de e^{κη}, que no tiene ceros. El modo espurio vive en
  la escala |η|, donde el desarrollo ya no vale.
- La réplica tiende a la tasa del continuo a medida que la capa se resuelve: en η/Δz = −4, −3 y
  −2, ρ − 1 vale 6,68·10⁻³, 1,22·10⁻² y 3,01·10⁻², contra exp(ν dt/η²) − 1 = 6,52·10⁻³,
  1,16·10⁻² y 2,64·10⁻². Son diferencias del 2, 5 y 14 %.
- En SPECTER, con la guardia desactivada (sólo en el scratch), N = 2 con η0 = −0,01 da NaN
  antes de t = 0,3. N = 2 con −0,005 y N = 3 con −0,01 y −0,04 son estables.
- El umbral real depende de D = ν dt/Δz². Lo da `estabilidad_cierre.py` por bisección
  (entrada `umbral_vs_D` del JSON), y el verificador lo reprodujo aparte con su propio código
  (en el
  verificador). 1 + σ en el umbral ≈ 0,7 D, y N = 3 es estable en todos los casos:

  | D | umbral η/Δz | 1 + σ ahí |
  |---|---|---|
  | 0,021 (nz = 128, dt = 2·10⁻⁴) | −0,707 | 0,013 |
  | 0,104 (la puerta: nz = 128, dt = 10⁻³) | −0,667 | 0,069 |
  | 0,132 (nz = 256, dt = 2,5·10⁻⁴) | −0,651 | 0,091 |
  | 0,177 (nz = 128, dt = 1,7·10⁻³) | −0,622 | 0,131 |

- Decisión: la guardia de §3.6, con umbral 0,05 + D a N = 2. **No** se cambió el modelo de
  N = 2. Una variante bien planteada queda como [P-04].

**P4 — piso del residuo a N = 3.** Con la iteración corregida, Dx se estancaba en unos
10⁻⁹ relativo; W y P bajaban a 10⁻¹⁵ y 10⁻²¹. La cadena de controles, en una copia
instrumentada:
- F es exactamente determinista: evaluada dos veces, diferencia 0.
- La ganancia medida es la modelada: −σ.
- Picard puro sobre G tiene el mismo piso, así que no es Anderson.
- Las segundas diferencias de F valen ~8·10⁻¹⁶ para pasos ≥ 10⁻⁶·esc y ~2·10⁻²⁰ por debajo.
  Es redondeo que se descorrelaciona, no curvatura.
- El piso vive sólo en el modo de la onda, y en SPECTER se ve ya en la traza m = 3: su
  segunda diferencia vale 1,1·10⁻¹² con h = 10⁻⁵·esc, contra 4·10⁻¹⁷ por debajo del ulp. Por
  η²/2 = 8·10⁻⁴ eso da los ~8·10⁻¹⁶ de Dx.
- **La atribución cuantitativa no está establecida.** Estimé ε·max|û|·kz³ ≈ 10⁻¹², pero la
  estimación omitía el 1/nz de `fs_ph`, como encontró el verificador (C9). El ruido blanco de
  los coeficientes da ~10⁻¹⁴ en la traza.
- En su réplica de la cadena en z, el ruido medido de la traza m = 3 es 2,5–9·10⁻¹⁴. Lo domina
  la cancelación de la continuación FC, con sumas de filas |dir| ≤ 9,5·10³. Eso queda 10–40
  veces por debajo de lo medido en SPECTER. El resto vendría de pasos que la réplica no tiene:
  la proyección y las transformadas horizontales. No se aisló.
- Lo establecido: F es determinista y el piso está en la traza m = 3; es ∝ η² y no es un
  error de la iteración. Por eso el criterio de estancamiento de §3.5.

**P5 — el acople W → Dx.** El dato de u contiene −2ikₓW, así que Dx depende de W con −2σikₓ.
La primera pasada de Anderson saltaba a un residuo 11,6. Se incorporó al precondicionador:
el residuo de la pasada 2 bajó a 6·10⁻⁵.

**P6 — Anderson con coeficientes complejos.** Las incógnitas son coeficientes de Fourier de
campos reales. F es lineal sobre ℝ pero no sobre ℂ, porque tiene productos en el espacio
físico. Un γ complejo predice residuos que F no produce. Se pasó a γ real, con la parte real
de la matriz de Gram. No era la causa de P4, pero es lo correcto.

**P7 — la escala del bloque P.** Estaba dominada por el modo (0,0), la presión media g·η0,
unas 10⁴ veces mayor que la de la onda. Eso hacía el criterio laxo en los modos que importan.
Se excluye la media.

**P8 — parte imaginaria espuria a N = 3.** Im η(1,1) ~ 3·10⁻¹³ contra ~10⁻¹⁶ a N = 2.
- Causa: con kz(Nyquist) = −nz/2, las derivadas impares de una columna real reciben un
  término imaginario, amplificado por kz³ en m = 3.
- El interpolante trigonométrico real tiene en Nyquist un coseno, cuyas derivadas impares
  se anulan en los nodos, y z = Lz es un nodo.
- Se omite ese término para m impar. Resultado: Im η(1,1) ≤ 4·10⁻¹⁸.

**P9 — el binario TESTS no carga las tablas de Neumann.** `setup_bc` no corre con
`-DTESTS_`, y la sonda de `fs_gen_setup` las necesita. Ahora el test las carga con
`load_neumann_tables`, la misma llamada de `v_setup`, y `fs_gen_setup` aborta con un mensaje
si faltan.

**P10 — mensajes de aborto sin `FLUSH`.** Se agregó `FLUSH(6)` antes de cada `MPI_ABORT`
nuevo, como en la Fase 5.

**P11 — el Makefile de `src/boundary` no sigue los `INCLUDE`.** Una compilación
incremental usa objetos viejos si sólo cambió `fsorder.f90` o `vboundary.f90`. No se
modificó el Makefile de upstream. Hay que borrar `boundary/*.o` o compilar en un árbol
limpio, que es lo que hacen las puertas.

**P13 — el primer criterio de estancamiento cortaba convergencias lentas.** La corrida de
H7 hecha fuera de la puerta mostró residuos aceptados de 10⁻⁷. La historia por pasada
(onda 2D, N = 3, ε = 0,04):
- Anderson no es monótono. Una secuencia 7,8·10⁻⁴ → 1,2·10⁻⁶ → 6,8·10⁻⁸ → 1,0·10⁻⁷ disparaba
  el criterio en la pasada 4, cuando ese mismo subpaso llega a 5·10⁻¹¹ en la pasada 10.
  Lo mismo pasaba en la pasada 7 con un solo rebote después de un mínimo nuevo.
- La versión final compara el **mejor** residuo de las últimas cuatro pasadas con el mejor
  anterior, desde la sexta pasada, con umbral 10²·fstol.
- Resultado en esa onda: los subpasos que terminan por encima de 10⁻⁹ pasan de 1798 a 17
  (de 8002), y las muestras del diagnóstico de ≤ 10⁻⁷ a ≤ 1,3·10⁻¹⁰.
- En la escalera de profundidad (el piso real de P4) terminan por encima de 10⁻⁹ 481 de 8002
  (≤ 9,2·10⁻⁹, 8–16 pasadas); las muestras del diagnóstico, ≤ 1,8·10⁻⁹.
- Que el piso escala como η² es una **inferencia** del factor η²/2 de la traza m = 3; no se
  midió con un segundo η, porque los residuos aceptados están acotados por el umbral 10⁻⁸. En la
  celda (η/Δz ~ 10⁻³) quedaría unas 10⁷ veces por debajo. La
  estimación ε·kz³·(η²/2)·(nz/k)·C tiene una constante libre C y no es una comprobación (C9).

**P14 — la revisión independiente** (rol «verifier»; informe y código en
`numerico/fase6/verificador/`).

Verificó por su cuenta:
- el modelo, exacto hasta grado N contra las condiciones exactas, para N = 1–4, en 2D y 3D;
- la réplica, contra una transliteración literal del Fortran;
- β_m, la degeneración de la profundidad y el precondicionador;
- Anderson real y la convención de Nyquist;
- la aritmética de la guardia y la lectura del código (estado entre pasadas, estado aceptado,
  reducciones MPI, índices del criterio nuevo).

Encontró seis cosas: cinco se corrigieron y una queda abierta.
- **B1.** La guardia a 1 + σ ≤ 0,05 dejaba una ventana inestable sin mensaje, entre −0,680 y
  −0,667 Δz con el dt de la puerta, y más ancha con dt mayores. Se reprodujo con bisección en
  la réplica (tabla en P3). **Se corrigió:** umbral 0,05 + ν dt/Δz² a N = 2.
- **B2.** El residuo de [H-25] no es −(ε³/2)/(…) sino +(ε³/3)/(…): el numerador exacto es
  ε cosh ε − sinh ε. **Se corrigió** en todos los documentos; la conclusión O(η³) no cambia.
- **B3.** La cabecera de `fsorder.f90` escribía G sin el término de W, y el umbral como
  −0,72 Δz. **Se corrigió.**
- **B4.** Llegar a `fsmaxit` sin converger era silencioso. **Se agregó** un aviso por salida
  estándar en cada salida del diagnóstico, con cuántos subpasos no convergieron y el peor
  residuo. El archivo congelado no cambia.
- **C2 a N = 3.** El orden 3 es neutro sólo con η uniforme. Con pendiente congelada s, el par
  espurio crece con Re λ ≈ ν(2s²/η² − k²); por ejemplo +0,48 en η = −0,01, s = 0,05. Es un
  hallazgo del verificador, **no reproducido aparte**, de relevancia práctica incierta: los k
  inestables son los de la escala en que varía η. No se vio en las corridas (T = 4,
  ε ≤ 0,04). También afecta a la opción de [P-04]. Queda registrado en [H-24].
- **C9**, el origen del piso: corregido en P4.

**Segunda ronda, sobre las correcciones** (`informe_verificador_fase6_ronda2.md` en la misma
carpeta).

Verificó:
- que la guardia nueva de N = 2 cubre toda la ventana inestable, con cualquier D estable;
- las cantidades que usa, que no dispara a N = 3 ni en la celda, y la cabecera nueva;
- la lógica del contador de `fsmaxit`.

Encontró, y se hizo:
- **N1.** A N = 3 el paso es inestable cerca de η = −1,2 Δz, pero sólo a ≥ 0,97–0,99 del límite
  explícito (D_max = 0,2048 a nz = 128). Reproducido: ρ = 1,0097 a 0,99 D_max y estable a 0,97.
  **Se agregó** una guardia: N = 3 exige ν dt/Δz² ≤ 0,19.
- **N2.** El error de la traza m = 3 no es O(Δz²) más allá de nz ≈ 128. Reproducido con el
  operador de la réplica aplicado a sin(1,2z): 1,8·10⁻⁴ (nz = 128), 2,4·10⁻⁴ (256),
  2,0·10⁻³ (512), 1,7·10⁻² (1024). Con 1 + 0,3z², cuya ∂z³ exacta es 0, crece ×10 por
  duplicación desde nz = 64. Es un piso determinista de las tablas FC. **Se corrigió** [H-26]:
  con estas tablas N = 3 no mejora refinando en z, y el ×3,4 de 128 a 256 en la profundidad no
  lo explica la traza.
- **Documentación.** Todavía se afirmaba el umbral viejo (−0,70 / −0,68 Δz), el conteo de
  pasadas no coincidía, se afirmaba una corrida de la puerta antes de terminarla y había
  referencias de línea viejas. **Corregido.**

**P15 — la condición inicial de banda ancha de [V-19] sólo es correcta con un proceso MPI.**
- Lo encontró la primera comparación np = 1 contra np = 2 de los campos de superficie: η
  difería un 300 %.
- `condicion_inicial` de `verificacion/decaimiento_banda_ancha.py` escribe todos los modos
  dentro de `IF (ista .eq. 1)`, con índices globales de kx. Con más de un proceso, los modos
  de las columnas de otros procesos se pierden o pisan memoria.
- [V-19] corrió siempre con un proceso, así que sus números no se ven afectados, y su
  `sbatch` pide `--ntasks=1`.
- El camino general, en cambio, está bien en paralelo. Con la condición inicial de la puerta,
  que sí respeta la partición, np = 1 y np = 2 dan diferencia 0 a N = 2 y 7·10⁻¹⁸ a N = 3.
- `numerico/fase6/campos_superficie.py` usa una versión que escribe cada modo en el proceso
  que tiene su columna. Con ella np = 1 y np = 2 dan η idéntico (`--comparar-np`).

**P16 — las transformadas de superficie del orden N son O(nxy³) y no escalan con MPI.**
- `fs_to_phys` y `fs_to_spec` son DFT explícitas sobre la grilla acolchada (2nx × 2ny), y cada
  proceso hace la transformada completa después de juntar los modos.
- Medido: 16² tarda 33 ms por paso con un proceso, y 64² ~1,2 s por paso con tres. A 64²
  dominan.
- Antes de una producción a 128² hay que pasarlas a FFT (FFTW 2D sobre la grilla acolchada) y
  repartirlas. No se hizo: queda anotado en ESTADO.

**P12 — antes de congelar la puerta** (sesión anterior, resumido en su registro):
- el Newton de la referencia no convergía con tol = 10⁻¹³, por debajo del piso de redondeo;
  se pasó a 10⁻¹¹ con detección de estancamiento;
- Nz = 24 no resolvía el arranque impulsivo; se usa Nz = 40;
- la guardia de Nyquist en y rechazaba ky = 0 con ny = 1; corregido en `vboundary.f90`;
- errores del arnés de la puerta: normalización, piso, 101 muestras y escala de la deriva.

---

## 5. Verificación en SPECTER fuera de la puerta

**Escalera de profundidad** (η0 uniforme, modo (1,1), γ = 0,2, ν = 0,01; error contra la
referencia MAC a profundidad 1 + η0, con el piso de malla restado). El cálculo es el de H5,
con dt y nz libres: `numerico/fase6/escaleras.py profundidad`, que escribe
`salidas/tablas/fase6_escalera_profundidad.json`. Tardó 27 min y el pico fue 281 MB. Las filas
a nz = 128 son de corridas exploratorias de la sesión.

| | η0 = 0,01 | 0,02 | 0,04 | órdenes |
|---|---|---|---|---|
| truncamiento no viscoso esperado, N = 2, T = 2 **(§)** | 2,0·10⁻⁷ | 1,5·10⁻⁶ | 1,2·10⁻⁵ | 2,96 · 2,92 |
| SPECTER N = 2, nz = 256, dt = 2,5·10⁻⁴, T = 2 | 5,2·10⁻⁷ | 2,25·10⁻⁶ | 1,48·10⁻⁵ | 2,11 · 2,72 |
| truncamiento no viscoso esperado, N = 3, T = 2 **(§)** | 9,9·10⁻⁸ | 7,7·10⁻⁷ | 5,8·10⁻⁶ | 2,96 · 2,92 |
| SPECTER N = 3, nz = 256 | 8,1·10⁻⁷ | 2,76·10⁻⁶ | 1,20·10⁻⁵ | 1,77 · 2,12 |
| SPECTER N = 3, nz = 128, dt = 2,5·10⁻⁴ | — | — | 4,1·10⁻⁵ | — |

- El truncamiento esperado sale de ω_N² = (gk + γk³)·S_N/C_N, con S_N, C_N las trazas de
  orden N − 1 de sinh y cosh (`truncacion_h5.py` en el scratch).
- Para N = 1 predice 3,05·10⁻² a T = 4, y SPECTER da 3,2·10⁻².
- **dt no cambia nada:** los errores de N = 2 y N = 3 a nz = 128 son idénticos con dt = 10⁻³
  y dt = 2,5·10⁻⁴.
- **nz sí:** el piso de malla, con T = 2, cae de 1,2·10⁻⁶ a 3,6·10⁻⁸, un factor 34 (quinto
  orden FC-Gram). El 2,2·10⁻⁶ de H5 es el mismo piso con T = 4.
- N = 3 cae ×3,4, pero eso **no** lo explica la traza m = 3, cuyo error para perfiles suaves
  es parecido a nz = 128 y 256 (P14). El origen del exceso de N = 3 a nz = 128 no está
  establecido.

**Onda estacionaria no lineal 2D** (la de H7, 16×1×128, contra la geometría exacta con la
métrica de H7). La puerta no llega a medirla porque aborta en su primer caso N = 2. La mide
`numerico/fase6/escaleras.py no_lineal`, que escribe
`salidas/tablas/fase6_escalera_no_lineal.json` (3 min, 313 MB):

| | ε = 0,01 | 0,02 | 0,04 | órdenes | pasadas ≤ | residuo ≤ |
|---|---|---|---|---|---|---|
| N = 1 | 4,96·10⁻³ | 9,92·10⁻³ | 1,98·10⁻² | 1,00 · 1,00 | 2 | 2,4·10⁻¹⁵ |
| N = 3 | 6,87·10⁻⁷ | 3,45·10⁻⁶ | 2,81·10⁻⁵ | 2,33 · 3,02 | 13 | 1,3·10⁻¹⁰ |

La figura de las dos escaleras es `salidas/figuras/f_fs6_escaleras.png`.

**η0 negativo** (nz = 128, T = 1):

| caso | resultado |
|---|---|
| N = 2, η0 = −0,005 (−0,51 Δz) | estable, 4 pasadas, residuo 5·10⁻¹¹ |
| N = 2, η0 = −0,01 (−1,02 Δz) | la guardia aborta: «1 + sigma = −4,248E−01 <= 0,05 …» |
| N = 2, η0 = −0,01, sin guardia (scratch) | NaN antes de t = 0,3 |
| N = 3, η0 = −0,01 y −0,04 | estables, ≤ 6 pasadas, residuo ≤ 6,5·10⁻¹⁰ |

---

## 6. La puerta congelada — [V-21]

**Cómo se escribió.** `verificacion/test_aceptacion_fase6.py`, `referencia_fase6.py` e
`intent_fase6.txt` se escribieron antes del Fortran y se corrieron en rojo contra el árbol
sin la implementación: **1/10**, sólo H6, que verifica la referencia
(`verificacion/logs/fase6_rojo_2026-09-28.log`). Después se congelaron por hash (§9).

**Corridas.** Cuatro corridas completas, con el código en cuatro estados, todas con **4/10** y
los mismos chequeos y números:
1. el primer criterio de estancamiento: `logs/fase6_primera_4de10_2026-09-28.log`;
2. el criterio definitivo de P13: `logs/fase6_segunda_4de10_2026-09-28.log`, 9 min 21 s;
3. la guardia de N = 2 que depende de ν dt/Δz² y el aviso de `fsmaxit`:
   `logs/fase6_tercera_4de10_2026-09-28.log`, 9 min 18 s;
4. **final**, además con la guardia de N = 3 sobre ν dt/Δz² (segunda ronda del verificador):
   `logs/fase6_2026-09-28.log`, 9 min 16 s, pico 294 MB.

Después de la cuarta sólo cambió un comentario de la cabecera de `fsorder.f90`, no código.

Cada log tiene dos bytes NUL justo después del cartel de `prterun` de las dos corridas
abortadas, cuya salida la puerta copia. Vienen del lanzador de MPI y están en todos los logs en
el mismo lugar relativo: no es corrupción del disco.

| | resultado | número | por qué |
|---|---|---|---|
| H1 parser | **pasa** | acepta 2, rechaza 0 y 4 | — |
| H2 trazas | falla | m = 0: 2,5·10⁻¹⁶; m = 1: 3,5·10⁻⁹ (umbral 10⁻⁹); m = 2: 5,6·10⁻⁷ (10⁻⁶); m = 3: 2,9·10⁻⁴ (10⁻⁴) | umbrales fijados para el cierre descartado; son el piso FC-Gram de trazas espectrales ([H-26]) |
| H3 escalera algebraica | **pasa** | N = 1: 0,97–1,16 (κ 1,99); N = 2: 1,99–2,23 (κ 1,99); N = 3: 2,98–3,16 (κ 3,98) | — |
| H4 regresión N = 1 | **pasa** | Fase 5 contra MAC 1,4·10⁻⁶; general N = 1 contra Fase 5 2,5·10⁻¹⁵ | — |
| H5 profundidad | falla | N = 1: 0,98, 0,96; N = 2: 0,81, 1,80; N = 3: 1,49, 1,58 | criterio mal a priori ([H-25]), más el piso de malla a nz = 128 y el error de ∂z³u ([H-26]) |
| H6 referencia | **pasa** | lineal 6,6·10⁻⁹; balance 2,8·10⁻⁷; Nz 32 contra 40 5,1·10⁻⁸ | — |
| H7 escalera no lineal | falla | aborta en `h7_N2_0.01`: 1 + σ = −0,42 | N = 2 con η hasta −1,02 Δz ([H-24]). Fuera de la puerta, N = 1 y N = 3 en §5 |
| H8 orden temporal N = 2 | falla | aborta en `h8_0.001`: 1 + σ = −1,85 | ídem, η hasta −2,04 Δz |
| H9 iteración | falla | pasadas ≤ 12 (tope 30), **residuo 1,8·10⁻⁹** (≤ 10⁻⁹), deriva 9,5·10⁻¹⁷, max\|η\|/Δz 4,08 | sólo el residuo: el piso de redondeo de N = 3 con η = 4 Δz ([H-26], P4) |
| H10 discriminante | falla | «H7 no produjo la corrida N=2» | consecuencia de H7 |

**Lectura.** Ningún chequeo falla por un error de programación encontrado. Hay cuatro
aprobados que prueban el modelo, la regresión y la referencia. De los que fallan:
- **H2 y H5** tenían criterios fijados antes de saber cosas que no se sabían: la
  degeneración de [H-25] y el cambio de cierre de [H-23];
- **H7, H8 y H10** piden a N = 2 un régimen donde el modelo está mal planteado ([H-24]);
- **H9** pide un residuo por debajo del piso de precisión de ∂z³u a η = 4 Δz.

La puerta no se editó.

---

## 7. Bibliografía consultada

Los DOI se resolvieron contra la API de Crossref el 2026-09-28.

| trabajo | DOI | estado | para qué |
|---|---|---|---|
| Wang, Tice & Kim, *The viscous surface-internal wave problem: global well-posedness and decay*, ARMA 212, 1–92 (2014) | 10.1007/s00205-013-0700-2 | **leído en el original**, §1.1, ecs. (1.3)–(1.12) (arXiv:1109.1798v2, en la Fase 5) | las condiciones exactas de la superficie que se desarrollan en Taylor ([D-51]) |
| Kamrin, Bazant & Stone, *Effective slip boundary conditions for arbitrary periodic surfaces: the surface mobility tensor*, J. Fluid Mech. 658, 409–437 (2010) | 10.1017/s0022112010001801 | **leído en el original, sólo §3.1** (arXiv:0911.1328v2) | precedente del método: *«Using the method of domain perturbations, the expression of the boundary condition along the periodic surface becomes much easier to deal with by expanding in a Taylor series about z = 0»* (§3.1). Ahí se desarrolla **la solución** en ε (u = u₀ + εu₁ + …) para Stokes con no deslizamiento; acá se transfieren las condiciones con el campo completo. No se usa ninguno de sus resultados. |
| Walker & Ni, *Anderson acceleration for fixed-point iterations*, SIAM J. Numer. Anal. 49, 1715–1735 (2011) | 10.1137/10078356x | **localizado, no leído**: el reporte abierto de OSTI (10.2172/1240289) no se pudo descargar desde la laptop | nombre y variante de la aceleración. La implementación (`fs_anderson`) no se apoya en ese texto: está escrita completa en el código y se verificó por su convergencia medida (§4) |
| Fornberg, *Generation of finite difference formulas on arbitrarily spaced grids*, Math. Comp. 51, 699–706 (1988) | 10.1090/s0025-5718-1988-0935077-0 | **localizado, no leído** | pesos del cierre interior descartado y de la réplica. Verificados numéricamente: reproducen monomios hasta z⁷ en las derivadas 0–3 a 1,2·10⁻⁷, que es el redondeo con h = 0,01 |
| Trefethen, *Spectral Methods in MATLAB*, SIAM (2000) | 10.1137/1.9780898719598 | **localizado, no leído en esta sesión** | la matriz de Chebyshev de la referencia congelada; la referencia se verifica numéricamente en H6 |
| Hairer & Wanner, *Solving Ordinary Differential Equations II*, Springer (1996) | 10.1007/978-3-642-05221-7 | **localizado, no leído** | Radau IIA de 3 etapas en la referencia. Coeficientes verificados numéricamente: c = (0,155051; 0,644949; 1), condiciones B(5) y C(3) a 5,6·10⁻¹⁷ |

**Búsqueda que no se hizo, y por eso no se afirma nada.** No se buscó si la mala
postura de la transferencia de Taylor de orden 2 de la tensión tangencial viscosa ([H-24])
está descrita en la literatura. [H-24] queda como derivación de esta sesión, con su
verificación numérica. No es una afirmación de novedad.

---

## 8. Cómo reproducir

Todo en la laptop, en secuencia, con `ulimit -v 2500000`. Picos medidos: ≤ 280 MB al
compilar y ≤ 150 MB al correr.

```bash
PY=/home/lucia/miniforge3/envs/piv-dt/bin/python
# 0. árbol de SPECTER desde upstream + parche (numerico/specter-parche/README.md)
# 1. el modelo, sin SPECTER (segundos)
$PY numerico/fase6/chequeo_expansion.py          # -> salidas/tablas/fase6_chequeo_expansion.json
# 2. estabilidad de los cierres, réplica 1D (unos 20 s, 150 MB)
$PY numerico/fase6/estabilidad_cierre.py         # -> tabla y f_fs6_estabilidad.png
# 3. las dos escaleras fuera de la puerta (compilan en un temporal)
$PY numerico/fase6/escaleras.py profundidad --nz 256 --dt 2.5e-4 --T 2   # 27 min, 281 MB
$PY numerico/fase6/escaleras.py no_lineal                                 # 3 min, 313 MB
$PY numerico/fase6/escaleras.py figura
# 4. la puerta congelada (compila en un scratch; unos 9 min, 294 MB)
ulimit -v 2500000
SPECTER_SCRATCH=/ruta/scratch $PY verificacion/test_aceptacion_fase6.py
# 5. regresión
$PY verificacion/test_aceptacion_fase2.py        # 6/6 (6 min, 281 MB)
$PY verificacion/test_aceptacion_fase5.py        # 8/8 (10 min, 282 MB)
```

Para correr el camino de orden N en una simulación propia:

```
&freesurface
fsgrav = …  fstens = …  fsamp = …  fskx = …  fsky = …  fscoup = 2
fsorder = 2          ! 1, 2 o 3
fsmean  = 0.0        ! elevación media inicial
fstol   = 1e-10      ! tolerancia relativa por subpaso
fsmaxit = 30
fsgen   = 0
/
```

`freesurface_order_diagnostic.txt` informa, en cada `cstep`:
- la tensión tangencial completa de orden N (columna 10);
- las pasadas y el residuo final (columnas 11 y 12);
- max|η|/Δz (columna 13).

**Revisar la 12 y la 13.** La corrida aborta sola en tres casos:
- a N = 2, si 1 + σ ≤ 0,05 + ν dt/Δz² en algún punto (η ≤ −0,607 Δz con ν dt/Δz² = 0,104);
- a N = 3, si ν dt/Δz² > 0,19;
- en cualquier orden, avisa por salida estándar si algún subpaso agotó `fsmaxit`.

---

## 9. Hashes

| archivo | SHA-256 |
|---|---|
| `verificacion/test_aceptacion_fase6.py` (congelado) | `e9724eb1b21f92a7064d224a597234df3fa9721c14c16c34b2903148a10017b2` |
| `verificacion/referencia_fase6.py` (congelado) | `9b3574099287957e2cdfc1576c3eb38f6f195216ed1d442689e4ea5c5123cc23` |
| `verificacion/intent_fase6.txt` (congelado) | `f34e0ae028b0407d4d1a9b0d41a20db8f22ffecbe93e2e980bcd4cc5c6746022` |
| Fases 1, 2 y 5: sin cambios | ver `ESTADO.md` |

El hash del parche regenerado está en `numerico/specter-parche/README.md`.

---

## 10. Lo que queda abierto

- **[P-04]** Un modelo de orden 2 bien planteado para η < 0. Por ejemplo, llevar sólo la
  transferencia tangencial a orden 3 (consistente a orden 2, sin el cero de 1 + σ) o una
  forma de Padé. Cualquiera de las dos cambia el modelo congelado y necesita decisión.
- Una ∂z³u más precisa en el borde, para que N = 3 mejore a N = 2 a nz = 128. Por ejemplo,
  usar la ecuación de momento en la frontera para reemplazar derivadas normales altas.
- La aplicación física (el efecto O(Fr²) de la deformación sobre α, N = 2 contra `freeslip`
  a g y γ físicos) **no se corrió** en esta sesión: con 4–17 pasadas por subpaso, un 16² como
  el de [V-19] es del orden de una hora en la laptop, y corresponde al clúster.
