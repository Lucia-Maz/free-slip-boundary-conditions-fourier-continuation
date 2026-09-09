# Bibliografía

Todos los DOI de esta lista fueron **resueltos contra la API de Crossref** el 2026-09-02,
no citados de memoria. La columna de metadatos es la que devolvió Crossref, así que
cualquiera puede repetir la verificación con:

```bash
curl -s https://api.crossref.org/works/<DOI> | python3 -m json.tool | head -40
```

La salida cruda de esa verificación está en `_crossref_raw.json`.

**Leído en el original** significa que bajé el texto completo y lo cité de ahí.
El resto está listado porque hace falta leerlo, y está marcado como tal: no se debe
apoyar ninguna afirmación en ellos hasta haberlos leído.

---

## PIV optimization for the study of turbulent flow using spectral analysis

- **DOI:** [`10.1088/0957-0233/15/6/003`](https://doi.org/10.1088/0957-0233/15/6/003) — verificado: **OK**
- **Autores:** Foucaut, J M; Carlier, J; Stanislas, M
- **Publicación:** Measurement Science and Technology **15**, 1046-1058 (2004)
- **Para qué se usa acá:** **Referencia central.** De acá salen: E_ii = E_noise·sinc²(kX/2); corte universal k_c·X = 2.8; ζ = E_noise·Y constante; σ_u = sqrt(4ζI/XY) con I=1.492 (1.418 con 50% de solape); σ medido de 0.087 px para ventana 16×16, 0.055 para 32×32, 0.033 para 64×64 (su tabla 1). El valor de 0.087 px es contra el que se contrasta el σ que despejamos de los datos de Lucía. PDF de la versión ISPIV 2003 guardado acá y **leído en el original** (verificado formula por formula, ver [V-05]); la versión larga de MST 2004 es de acceso pago y **no fue leída**. Atención: esta versión anuncia que el ruido debido al movimiento de las partículas "will be studied afterwards" y nunca lo desarrolla — y ese es justamente el termino que en nuestros datos resulta dominante.

## Optimization of particle image velocimeters. I. Double pulsed systems

- **DOI:** [`10.1088/0957-0233/1/11/013`](https://doi.org/10.1088/0957-0233/1/11/013) — verificado: **OK**
- **Autores:** Keane, R D; Adrian, R J
- **Publicación:** Measurement Science and Technology **1**, 1202-1215 (1990)
- **Para qué se usa acá:** Origen de la 'regla del cuarto': el desplazamiento no debe superar 1/4 del lado de la ventana de interrogación. Fija el extremo superior del rango útil de Δt. **No verificado en el texto** — citado de la literatura secundaria, hay que leer el original antes de usarlo como criterio.

## Digital particle image velocimetry

- **DOI:** [`10.1007/BF00190388`](https://doi.org/10.1007/BF00190388) — verificado: **OK**
- **Autores:** Willert, C. E.; Gharib, M.
- **Publicación:** Experiments in Fluids **10**, 181-193 (1991)
- **Para qué se usa acá:** Efecto de filtrado pasabajos de la ventana de interrogación; longitud de onda de corte igual al doble del tamaño de ventana. Es la fuente que Foucaut cita para la respuesta tipo sinc. **No leído en original.**

## Iterative image deformation methods in PIV

- **DOI:** [`10.1088/0957-0233/13/1/201`](https://doi.org/10.1088/0957-0233/13/1/201) — verificado: **OK**
- **Autores:** Scarano, F
- **Publicación:** Measurement Science and Technology **13**, R1-R19 (2001)
- **Para qué se usa acá:** Métodos iterativos con deformación de ventana, que es lo que usa el multipaso 64→32→16 de las sesiones de Lucía (`deformation_method = symmetric`). Necesario para justificar que el multipaso no invalida el modelo de ventana única. **No leído en original.** Crossref fecha la publicación en 2001 aunque el volumen 13 de MST es de 2002.

## Spatially adaptive PIV interrogation based on data ensemble

- **DOI:** [`10.1007/s00348-009-0782-7`](https://doi.org/10.1007/s00348-009-0782-7) — verificado: **OK**
- **Autores:** Theunissen, Raf; Scarano, Fulvio; Riethmuller, Michel L.
- **Publicación:** Experiments in Fluids **48**, 875-887 (2009)
- **Para qué se usa acá:** PIV adaptativo: ventana variable en tamaño, forma y orientación según gradientes y sembrado. Es el antecedente más cercano a un 'recomendador de configuración', y sirve para delimitar qué es nuevo acá: ellos optimizan calidad local de correlación, no un objetivo global sobre un estadístico de turbulencia, y no tocan Δt. **No leído en original.** Crossref lo fecha en 2009 (online first); el número impreso es de 2010.

## Robust smoothing of gridded data in one and higher dimensions with missing values

- **DOI:** [`10.1016/j.csda.2009.09.020`](https://doi.org/10.1016/j.csda.2009.09.020) — verificado: **OK**
- **Autores:** Garcia, Damien
- **Publicación:** Computational Statistics &amp; Data Analysis **54**, 1167-1178 (2010)
- **Para qué se usa acá:** El algoritmo `smoothn` que OpenPIV aplica cuando `smoothn_each_pass=True`, activo en la familia de sesiones `180226`. Hace falta para caracterizar la función de transferencia extra que ese suavizado introduce. Ver decisión [D-07].

## Collaborative framework for PIV uncertainty quantification: comparative assessment of methods

- **DOI:** [`10.1088/0957-0233/26/7/074004`](https://doi.org/10.1088/0957-0233/26/7/074004) — verificado: **OK**
- **Autores:** Sciacchitano, Andrea; Neal, Douglas R; Smith, Barton L; Warner, Scott O; Vlachos, Pavlos P; Wieneke, Bernhard
- **Publicación:** Measurement Science and Technology **26**, 074004 (2015)
- **Para qué se usa acá:** Comparación de los cuatro métodos de cuantificación de incertidumbre en PIV. Marco de referencia para ubicar el método de dos Δt que usamos acá, que no es ninguno de esos cuatro. **No leído en original.**

## Uncertainty quantification in particle image velocimetry

- **DOI:** [`10.1088/1361-6501/ab1db8`](https://doi.org/10.1088/1361-6501/ab1db8) — verificado: **OK**
- **Autores:** Sciacchitano, A
- **Publicación:** Measurement Science and Technology **30**, 092001 (2019)
- **Para qué se usa acá:** Revisión general de incertidumbre en PIV. Punto de entrada bibliográfico. **No leído en original.**

## PIV uncertainty quantification by image matching

- **DOI:** [`10.1088/0957-0233/24/4/045302`](https://doi.org/10.1088/0957-0233/24/4/045302) — verificado: **OK**
- **Autores:** Sciacchitano, Andrea; Wieneke, Bernhard; Scarano, Fulvio
- **Publicación:** Measurement Science and Technology **24**, 045302 (2013)
- **Para qué se usa acá:** Método de disparidad de partículas / image matching: estima el error a posteriori usando las imágenes y el campo medido. Alternativa independiente para contrastar el σ que obtenemos. **No leído en original.**

## piv-image-generator: An image generating software package for planar PIV and Optical Flow benchmarking

- **DOI:** [`10.1016/j.softx.2020.100537`](https://doi.org/10.1016/j.softx.2020.100537) — verificado: **OK**
- **Autores:** Mendes, Luís; Bernardino, Alexandre; Ferreira, Rui M.L.
- **Publicación:** SoftwareX **12**, 100537 (2020)
- **Para qué se usa acá:** Generador de imágenes sintéticas de PIV. Revisado y **descartado por ahora**: es MATLAB y sólo admite flujos analíticos (uniforme, corte, punto de estancamiento, vórtice de Rankine), no un campo arbitrario de DNS. Se registra para no volver a evaluarlo.

## The spectral response of time-resolved PIV in a turbulent boundary layer

- **DOI:** [`10.1007/s00348-025-04059-0`](https://doi.org/10.1007/s00348-025-04059-0) — verificado: **OK**
- **Autores:** Manovski, Peter; Abu Rowin, Wagih; Ng, Henry; Gulotta, Paul; Giacobello, Matteo; de Silva, Charitha
- **Publicación:** Experiments in Fluids **66**, — (2025)
- **Para qué se usa acá:** Trabajo reciente sobre respuesta espectral de PIV resuelto en tiempo. Relevante para ver si el hueco que identificamos sigue abierto. **No leído — pendiente.**

## Fourier continuation method for incompressible fluids with boundaries

- **DOI:** [`10.1016/j.cpc.2020.107482`](https://doi.org/10.1016/j.cpc.2020.107482) — verificado contra Crossref el 2026-09-07: **OK**
- **Autores:** Fontana, Mauro; Bruno, Oscar P.; Mininni, Pablo D.; Dmitruk, Pablo
- **Publicación:** Computer Physics Communications **256**, 107482 (2020)
- **Para qué se usa acá:** **Referencia central de la parte numérica.** Es el paper de SPECTER, el código en el que se implementa la condición de superficie libre. De acá salen: el esquema de partición temporal con campo auxiliar sin presión y proyección posterior (sección 4.1, ecuaciones 12–14); las condiciones de contorno que el código impone sobre `v*` y sobre la presión (ecuaciones 15 y 16); la solución analítica de la ecuación de Laplace homogénea en la caja y sus coeficientes (ecuaciones 19–25), que es la que extiende `laplace_z`; y la generalización a Runge-Kutta de orden o (sección 4.4, ecuaciones 27–31). Preprint [arXiv:2002.01392v2](https://arxiv.org/abs/2002.01392) (13 jul 2020, versión aceptada, 31 páginas) guardado acá y **leído en el original** en las secciones 2, 4.1, 4.2 y 4.4; el resto —FC-Gram en detalle, paralelización, y los casos de Poiseuille y Rayleigh-Bénard— **no fue leído**. Atención a una advertencia del propio paper que hay que declarar al usar el método: proyectar la ecuación de momento en la dirección normal o en la tangencial da problemas de contorno distintos para la presión, y "there is no a priori reason to assume that both approaches lead to the same solution"; ellos eligen la normal, y este trabajo hereda esa elección. Ver [V-08].

---

## Sin DOI

### On the noise in statistics of PIV measurements

- **arXiv:** [`2010.10768`](https://arxiv.org/abs/2010.10768) (v2, 23 dic 2020)
- **Autores:** George, William K.; Stanislas, Michel
- **Para qué se usa acá:** teoría moderna del ruido de PIV. Descompone la salida en
  velocidad, ruido de pixelización y ruido por no uniformidad de la velocidad dentro
  del volumen de interrogación; las dos últimas con varianza inversamente proporcional
  al número de partículas, y las tres filtradas espacialmente por la ventana.
  PDF guardado acá. **Leídos el resumen, la lista de símbolos y las secciones 3.2–3.3**; el desarrollo completo, no. De la 3.2 sale la frase que sostiene [H-07]: las desviaciones respecto del promedio volumétrico son "noise resulting from the flow itself". De la 3.3, la advertencia de que pixelización y número finito de partículas son "difficult to distinguish" espectralmente. **Es el trabajo del que menos leímos y sobre el que más se apoya el resultado principal: hay que leerlo entero antes del informe final.**

---

## La cadena de la condición de contorno de SPECTER

Estas cuatro entradas se agregaron el 2026-09-09 al rastrear de dónde salen las
condiciones sobre `v*` y sobre `p`. Los DOI se resolvieron contra la API de Crossref ese
mismo día. Los números entre corchetes son los de la lista de referencias del preprint de
Fontana et al., para poder seguir el hilo en el original.

### Boundary conditions for incompressible flows

- **DOI:** [`10.1007/BF01061454`](https://doi.org/10.1007/BF01061454) — verificado: **OK**
- **Autores:** Orszag, Steven A.; Israeli, Moshe; Deville, Michel O.
- **Publicación:** Journal of Scientific Computing **1**, 75-111 (1986)
- **Para qué se usa acá:** **Es el origen de las condiciones de contorno que usa SPECTER.**
  Fontana et al. lo citan como su ref. [51] y dicen textualmente: *"We also employ the
  stability enhancing modification to the wall-normal velocity presented in (51), where the
  following boundary conditions were introduced"*, y a continuación van sus ecuaciones (15)
  y (16), que son `v*_∥ = Δt ∇_∥ p` y `∂p/∂z = (1/Δt) ẑ·v*` en las caras. O sea que la
  condición de Neumann sobre la presión y la condición tangencial sobre el campo auxiliar
  **no son de SPECTER ni de su paper: son de acá**. Es la referencia que hay que leer para
  entender por qué la presión lleva Neumann y de dónde sale la velocidad de deslizamiento
  residual O(Δt^ord). **Localizado y verificado en Crossref; NO leído en el original.**

### Application of a fractional-step method to incompressible Navier-Stokes equations

- **DOI:** [`10.1016/0021-9991(85)90148-2`](https://doi.org/10.1016/0021-9991(85)90148-2) — verificado: **OK**
- **Autores:** Kim, J; Moin, P
- **Publicación:** Journal of Computational Physics **59**, 308-323 (1985)
- **Para qué se usa acá:** ref. [58] de Fontana et al. Es el esquema global O(Δt²) que ellos
  usan y que satisface las condiciones de contorno con error O(Δt). **Localizado; NO leído.**

### Numerical solution of the Navier-Stokes equations

- **DOI:** [`10.1090/S0025-5718-1968-0242392-2`](https://doi.org/10.1090/S0025-5718-1968-0242392-2) — verificado: **OK**
- **Autores:** Chorin, Alexandre Joel
- **Publicación:** Mathematics of Computation **22**, 745-762 (1968)
- **Para qué se usa acá:** ref. [57] de Fontana et al. El método de proyección original: el
  campo auxiliar sin presión y su proyección posterior. Es el eslabón más viejo de la
  cadena. **Localizado; NO leído.**

### High-order splitting methods for the incompressible Navier-Stokes equations

- **DOI:** [`10.1016/0021-9991(91)90007-8`](https://doi.org/10.1016/0021-9991(91)90007-8) — verificado: **OK**
- **Autores:** Karniadakis, George Em; Israeli, Moshe; Orszag, Steven A
- **Publicación:** Journal of Computational Physics **97**, 414-443 (1991)
- **Para qué se usa acá:** **No** lo cita Fontana et al. Se agrega porque es la continuación
  directa de Orszag, Israeli y Deville: introduce condiciones de contorno de presión de
  orden alto en el tiempo, pensadas para suprimir la capa límite numérica espuria que
  inducen los esquemas de partición. Es el lugar donde mirar si alguna vez molesta el error
  O(Δt^ord) del no-deslizamiento. **Localizado; NO leído.**

### Spectral Methods in Fluid Dynamics

- **Libro:** Canuto, C.; Hussaini, M. Y.; Quarteroni, A.; Zang, T. A., Springer, 1988.
  [10.1007/978-3-642-84108-8](https://doi.org/10.1007/978-3-642-84108-8)
- **Para qué se usa acá:** ref. [5] de Fontana et al., citada junto con Chorin como
  referencia general del esquema de partición. **Localizado; NO leído.**

> **Lo que NO se encontró.** Se buscó una formulación publicada de la condición sobre `v*`
> para una pared **libre de tensiones** dentro de este esquema de partición, y no se
> encontró. Eso no significa que no exista: la búsqueda fue con buscador general, no
> sistemática, y hay literatura de elementos espectrales sobre condiciones stress-free que
> no se leyó. La derivación propia está en [D-24] y en
> `teoria/superficie_libre_v_estrella_y_p.md`; lo que sí está verificado es que
> **SPECTER upstream no la tiene**: la documentación del repositorio, consultada el
> 2026-09-09, sigue listando `noslip` como única condición para la velocidad.
