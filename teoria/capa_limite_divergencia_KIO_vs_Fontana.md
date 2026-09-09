# La capa límite de divergencia de KIO §2.3, aplicada al esquema de SPECTER

**2026-09-09.** Chequeo pedido por Lucía. En la respuesta anterior yo había dicho que el
error de partición de SPECTER es "el mismo error contabilizado en otra variable" que el de
Karniadakis, Israeli y Orszag, y había marcado que eso lo estaba **infiriendo**. Acá está
hecho el análisis, y la inferencia resulta correcta, pero por una razón estructural más
nítida de lo que yo había supuesto.

Convención de marcas: **[K]** = Karniadakis, Israeli & Orszag, *J. Comput. Phys.* **97**,
414–443 (1991), [10.1016/0021-9991(91)90007-8](https://doi.org/10.1016/0021-9991(91)90007-8),
PDF en `bibliografia/karniadakis1991.pdf`. **[F]** = Fontana, Bruno, Mininni & Dmitruk,
*Comput. Phys. Commun.* **256**, 107482 (2020), preprint en `bibliografia/`. **[C]** =
verificado en el código de `numerico/SPECTER-trabajo/`. **[D]** = derivado acá.

---

## 1. De dónde sale la capa límite en [K]

**El supuesto que la produce es que la viscosidad se integra de manera implícita.** [K] lo
declara en §3, al enumerar los parámetros del esquema:

> *"the set of parameters (J_e, J_p, J_i) denoting the individual accuracy of the explicit
> integration (nonlinear terms), the pressure boundary condition, and the implicit
> integration (viscous terms), respectively."*

Con el término viscoso implícito, el campo nuevo `v^{n+1}` aparece **dentro** del
laplaciano. Al tomar la divergencia de la ecuación semidiscreta y anular el miembro
derecho —que es lo que impone la ecuación de la presión—, [K] §2.3 llega a su ecuación
(10b),

    Q − γ₀·ν·Δt·∇²Q = 0 ,      Q ≡ ∇·v^{n+1} ,

y de ahí:

> *"It is clear, therefore, that there exists a (numerical) boundary layer of thickness
> ℓ = √(γ₀νΔt), so that Q = Q_w e^{−s/ℓ}, and thus the boundary divergence is
> Q_w = −ℓ(∂Q/∂s)_w."*

*(El coeficiente γ₀ se lee de un escaneo con OCR imperfecto; lo que importa para todo lo
que sigue es el escalamiento ℓ ∝ √(νΔt), que sí es inequívoco.)*

La cadena causal que cierra el argumento, también textual:

> *"This relation demonstrates that the time-differencing error of the velocity field is
> one order smaller in Δt than the corresponding error in the boundary divergence."*

> *"The above argument clearly demonstrates that the time accuracy of the global solution
> is directly dependent on the boundary values of the divergence (i.e., ∂Q/∂n) and,
> therefore (see Eq. (6a)), on the treatment of pressure boundary conditions."*

En una frase: **el operador de Helmholtz `1 − c·νΔt·∇²` es la firma de la viscosidad
implícita, y sus soluciones homogéneas —que decaen como e^{−s/ℓ}— son la capa límite.**

## 2. El mismo paso, sobre el esquema de [F]

**El término viscoso de [F] es explícito.** Está en su ecuación (12),
`v* = v^t + Δt[(v·∇)v + ν∇²v^t + f^t]`, y se verificó en el código [C],
`src/include/hd/hd_rkstep2.f90`:

```fortran
vx(k,j,i) = C1(k,j,i) + dt*(nu*vx(k,j,i)-C4(k,j,i)+fx(k,j,i))*rmp
```

donde `vx` ya fue reemplazado por su laplaciano en la llamada anterior a `laplak`. No hay
ninguna resolución implícita.

Ahora se repite el paso de [K] §2.3. Las ecuaciones (13) y (14) de [F] son

    ∇²p = (1/Δt)·∇·v* ,        v^{n+1} = v* − Δt·∇p .

Tomando la divergencia de la segunda y sustituyendo la primera [D]:

    ∇·v^{n+1} = ∇·v* − Δt·∇²p = ∇·v* − Δt·(1/Δt)·∇·v* = 0 .

**Idénticamente, punto a punto, y sin ningún operador diferencial actuando sobre Q.** El
paso que en [K] lleva de su (10a) a su (10b) no tiene análogo: allá sobrevive el término
`−γ₀νΔt∇²Q` precisamente porque `v^{n+1}` está adentro del término viscoso; acá no está.

**Conclusión 1.** En el esquema de [F] no hay ecuación de Helmholtz para la divergencia,
no hay solución homogénea, y por lo tanto **no hay capa límite de divergencia**. No es que
sea chica: la estructura que la produce no está.

Dos condiciones que hay que declarar para que esto valga, y las dos se cumplen:

- **La ecuación de Poisson tiene que resolverse exactamente.** [F] la resuelven como
  `p = p^I + p^H`, con `p^H` una solución **homogénea** construida analíticamente con
  exponenciales (sus ecuaciones 19–25). Como `∇²p^H = 0` exactamente, sumarla no perturba
  `∇²p`, y la identidad de arriba sobrevive. El residuo queda fijado por la representación
  FC-Gram y el redondeo, **no por Δt**.
- **La imposición sobre `v*` es por inyección** —[F] la llama *"strong imposition (also
  known as 'injection')"*—, o sea que los valores de borde de `v*` se reemplazan. Eso
  cambia `∇·v*` en la pared, pero (13) se resuelve para el `v*` **ya reemplazado**, así
  que la cancelación sigue siendo exacta.

## 3. Adónde va, entonces, el error de partición

A la condición de contorno tangencial. Tomando la parte tangencial de (14) en la pared e
insertando la condición (15) de [F], `v*_∥ = Δt·∇_∥p^t` [D]:

    v_∥^{n+1}|pared = Δt·∇_∥(p^t − p^{n+1}) = O(Δt²)

para Euler hacia adelante, y O(Δt^o) para Runge-Kutta de orden o. El origen es que (15)
necesita la presión **al nivel nuevo**, que todavía no se calculó, y usa la del nivel
anterior. [F] lo dice y lo mide:

> *"these boundary conditions result in a O(∆t²) slip velocity at the walls, a fact that
> will be numerically verified."*

**Conclusión 2.** Los dos esquemas pagan el mismo error de partición, pero en variables
distintas: **[K] satisface la condición de contorno y paga con una capa de divergencia de
espesor √(νΔt); [F] satisface la divergencia exactamente y paga con un defecto O(Δt^o) en
la condición de contorno.** La inferencia que yo había marcado como no verificada era
correcta, y ahora está derivada.

Es coherente con lo que [F] miden en su Fig. 4 y su texto: `⟨w²⟩|pared` ≈ 10⁻³⁵ e
independiente de Δt (la condición normal es exacta, porque la ecuación (16) se deriva
pidiendo `v_z^{n+1} = 0` **al nivel nuevo**, sin extrapolar), `⟨(∇·v)²⟩` ≈ 10⁻¹⁶ e
independiente de Δt, y `⟨u²⟩|pared` ∝ Δt⁴.

## 4. Una consecuencia que no está en ninguno de los dos papers

**[D], y la marco como propia:** aunque un esquema explícito llegara a producir una capa
de divergencia, sería **sub-grilla por construcción**. La estabilidad del término viscoso
explícito exige `Δt ≲ 2/(ν·k_max²)` con `k_max ~ π/dz`, es decir `Δt ≲ 2dz²/(π²ν)`. Entonces

    ℓ = √(νΔt) ≲ √(2/π²)·dz ≈ 0,45·dz ,

o sea que el espesor de la capa nunca supera medio paso de grilla. La capa límite de [K]
es un fenómeno de esquemas con viscosidad implícita, que son justamente los que pueden
tomar Δt grande.

## 5. El caso de la superficie libre

Con tope libre la condición se impone sobre la **vorticidad** tangencial, y como la
proyección resta un gradiente, `rot(v*) = rot(v)` exactamente ([D-24]). Entonces el error
de partición **no aparece ni en la divergencia ni en la condición de contorno**: las dos
quedan a nivel de representación. Es una predicción más fuerte que la del canal, y es la
que V3 de la puerta ya midió (razón entre el residuo a Δt y a Δt/2 igual a 1,00).

## 6. El chequeo numérico

`verificacion/kio_capa_divergencia.py`, sobre el árbol de trabajo, malla 8×8×128, ν = 10⁻²,
campo 3D con kx ≠ 0, cuatro pasos de tiempo entre 10⁻⁴ y 8·10⁻⁴ hasta t = 0,08. Tres
predicciones falsables:

| | predicción | exponente medido, y ~ Δt^p | veredicto |
|---|---|---|---|
| **P1** | `⟨(∇·v)²⟩` independiente de Δt, no-deslizamiento | **p = +0,000** | **pasa** |
| **P1** | `⟨(∇·v)²⟩` independiente de Δt, tope libre | **p = −0,000** | **pasa** |
| **P2** | `⟨\|v_t\|²⟩` en la pared ∝ Δt^(2·ORD) = Δt⁴ | **p = +4,004** | **pasa** |
| **P3** | residuo de tensión en la cara libre, independiente de Δt | **p = −0,003** | **pasa** |

Los números crudos, con Δt entre 10⁻⁴ y 8·10⁻⁴:

| Δt | `⟨(∇·v)²⟩` noslip | `⟨\|v_t\|²⟩` pared | `⟨(∇·v)²⟩` freeslip | residuo de tensión, cara libre |
|---|---|---|---|---|
| 8·10⁻⁴ | 1,069·10⁻¹⁶ | 3,208·10⁻¹¹ | 3,340·10⁻¹⁷ | 6,763·10⁻⁸ |
| 4·10⁻⁴ | 1,069·10⁻¹⁶ | 1,994·10⁻¹² | 3,340·10⁻¹⁷ | 6,775·10⁻⁸ |
| 2·10⁻⁴ | 1,069·10⁻¹⁶ | 1,243·10⁻¹³ | 3,340·10⁻¹⁷ | 6,785·10⁻⁸ |
| 1·10⁻⁴ | 1,069·10⁻¹⁶ | 7,759·10⁻¹⁵ | 3,340·10⁻¹⁷ | 6,802·10⁻⁸ |

**La divergencia no cambia una sola cifra** en un factor 8 de Δt, mientras el deslizamiento
cae cuatro órdenes de magnitud. Eso es exactamente la conclusión 2: el error de partición
está íntegramente en la condición de contorno y nada de él está en la divergencia. P2
reproduce sobre este árbol la Fig. 4 de [F].

Dos lecturas honestas de la tabla:

- El residuo de tensión de la cara libre **sube** un 0,6 % al bajar Δt ocho veces, de donde
  sale el exponente −0,003. No es un escalamiento con Δt: es deriva de nivel de
  representación. Lo que P3 afirma —que no escala como Δt^o— se cumple con margen de
  cuatro órdenes.
- En la corrida con tope libre, la columna de z = 0 vale 1,196·10⁴ e independiente de Δt.
  No es un error: es el **corte físico** en la pared rígida de abajo, que es contra lo que
  se normaliza el residuo de la cara libre. El cociente da 5,7·10⁻¹², coherente con el
  6,1·10⁻¹² que había medido V3 de la puerta.

Y el espesor que tendría la capa, si existiera: con dz = 0,0098, `ℓ/dz` va de 0,10 a 0,29
en el rango de Δt usado, y el límite de estabilidad viscosa explícita (Δt ≲ 1,9·10⁻³) da
`ℓ/dz ≲ 0,45`. Es la comprobación numérica del argumento de la sección 4.

## 7. Salvedades

- El coeficiente exacto de `ℓ = √(γ₀νΔt)` se leyó de un escaneo con OCR imperfecto.
- No se leyó Orszag, Israeli & Deville (1986), que es de acceso cerrado; lo que [K] §2.1
  dice de él es que *"it was shown in [2] that such a boundary condition leads to
  instabilities"* para la versión ingenua de la condición de presión.
- El análisis de §2 es sobre el esquema continuo en el espacio. No cubre lo que la
  continuación FC-Gram le hace al residuo, que es lo que fija el piso de 10⁻¹⁶.
- No se estudió el caso de [K] §2.1 —escribir el término viscoso como `−ν∇×(∇×v)` en la
  condición de presión— porque en un esquema explícito ese término no aparece en la
  condición de contorno de la presión.
