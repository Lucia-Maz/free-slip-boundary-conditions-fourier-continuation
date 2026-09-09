# Ronda 01 — verifier (independiente, adversarial)

Leídos `spec.json`, `state.json`, `reports/r00-pi-staffing.md`, `r01-lead.md`,
`r01-worker-1.md`, la puerta completa, `out/{notes.tex,capa.py,checks.py,provenance.json}`
y `SPECTER-upstream/src/boundary/boundary_mod.fpp`. `inbox.jsonl` no existe: cero
direcciones humanas. No había claims previos `refuted` ni `unclear`.

Nada se dio por bueno leyendo el argumento del worker: cada resultado se rehízo con un
método distinto. Scripts propios en `work/verif/` (`indep.py`, `indep2.py`, `indep3.py`,
`indep4.py`, `prov.py`), ejecutados **en serie**; RSS máximo observado 68 MB.

## Rutas independientes usadas

| resultado | ruta del worker | mi ruta independiente |
|---|---|---|
| espectros | fórmula cerrada + FD | colocación de Chebyshev (N=64) **y** raíces de la secular por `brentq` |
| alfa | autovalor lento | `-(nu/h)[dF/dz]_0^h` **y** decaimiento asintótico de `<u>` bajo Crank–Nicolson desde dato **genérico** |
| pi/2, beta | `trapezoid` en 2e5 puntos | `scipy.integrate.quad` adaptativa |
| Laplace DN | residuos de borde | resolución del BVP tridiagonal `solve_banded` punto a punto + barrido de 2404 casos |
| solver_1d | `expm_multiply` | Crank–Nicolson propio sobre la misma malla |
| convenciones SPECTER | lectura | relectura y reconstrucción algebraica de las ramas DD y NN existentes |

## Reproducido (números míos)

- **Espectro DD** `nu(n pi/h)^2`, n>=1: Chebyshev err rel 2.66e-14; secular por `brentq` err 0.0.
- **Espectro DN** `nu((m+1/2)pi/h)^2`, m>=0: Chebyshev 2.45e-12; `brentq` 0.0.
  Escalado en `nu` y `h` exacto a 3.4e-16 con `nu=3.7, h=0.013`.
- **alfa**: `alfa_DD=pi^2 nu/h^2`, `alfa_DN=pi^2 nu/(4h^2)`, cociente 4. Tercera ruta:
  Crank–Nicolson (nz=2000, T=2) desde un dato inicial que **no** es autovector; la tasa
  asintótica de `<u>` da 2.467400671 (DN, err rel 1.7e-07) y 9.870051915 (DD, 4.5e-05)
  contra `capa.alfa`. Esto verifica el término `-alfa*U` en sí, no sólo el autovalor.
- **Sesgo PIV** `u(h)/<u>=pi/2=1.5707963267948966` por cuadratura adaptativa; |dif| 2.2e-16.
  Exceso 57.0796326795 %.
- **beta**: `<F^2>=pi^2/8` para **ambos** perfiles normalizados; `beta_s=beta*(2/pi)=pi/4`
  reproducido por el cambio de variable `U_s=(pi/2)U`.
- **Brechas**: `lam1-lam0=2pi^2` (DN), `lam2-lam1=3pi^2` (DD); identidades exactas.
- **Escalas**: `tau_nu/tau_H = Re_l eps^2 = Re_h eps` idénticas a precisión de máquina.
- **Coeficientes de expansión** `a_n=b_m=(2/h)\int u_0 Z`: matriz de Gram `= (h/2)I` a
  1.11e-16 para las dos familias; reconstrucción con 200 modos converge.
  *Esto no lo cubre ni `checks.py` ni la puerta; lo verifiqué yo.*
- **Laplace DN**: rederivé el sistema `r C1 + C2 = g0`, `k(C1 - r C2) = g1` y obtuve los
  mismos `C1=(g1/k + r g0)/(1+r^2)`, `C2=(g0 - r g1/k)/(1+r^2)`. Coincide punto a punto con
  un BVP tridiagonal independiente. Barrido de 2404 casos `(k,Lz)` cruzando el umbral
  `kLz=0.5`: residuo máximo **3.886e-15**. Sin desborde hasta `kLz=1e8` (`C1,C2` finitos,
  residuo 0). Rama `k=0` exacta. La forma regular
  `g0[cosh kz - tanh(kLz) sinh kz] + g1 sinh(kz)/(k cosh kLz)` satisface ambos bordes
  **exactamente** (verificación algebraica) y `np.sinc(1j kz/pi).real = sinh(kz)/(kz)`
  sin warnings; datos complejos correctos; continua a través del umbral.
- **Convenciones SPECTER** (todas las de `notes.tex` §Laplace son correctas):
  `a = coef(1)e^{k(z-Lz)} + coef(2)e^{-kz}`, `b = k[coef(1)e^{k(z-Lz)} - coef(2)e^{-kz}]`
  (derivada hacia `+z`); fila `k=0`: `a=coef(1)z+coef(2)`, `b=coef(1)` -> pendiente/ordenada;
  `bc(1)`=dato en `z=0`, `bc(2)`=dato en `z=Lz` (deducido de la rama DD con `k=0`);
  `khom=sqrt(kx^2+ky^2)` (`specter.fpp:798`). Denominador DD `1-e^{-2kLz}` vs DN `1+e^{-2kLz}`.
- **`solver_1d` integra de verdad**: coincide con mi Crank–Nicolson en la misma malla a
  2.8e-10 (DD) y 4.1e-11 (DN). Su error de tasa es exactamente el error espacial FD
  `lam*kappa^2 dz^2/12`: predicho 2.2020e-04 a nz=48 (DN), observado 2.202e-04. Orden 2.000.
- **Puerta**, corrida por mí desde la raíz del proyecto: **6/6**, exit 0, módulo cargado
  `jobs/.../out/capa.py`. C1 8.23e-05, C2 6.30e-05, C3 4.000000 vs 3.999985,
  C4 1.570796 vs 1.570798, C5 orden 2.00, C6 5 casos. `sha256 =
  2db98ae5a286c64ff42e6bcbdf9e9d69ccbcae60cd631c16d28a93e9855c5527` (coincide). RSS 68 MB.
- **`out/checks.py`**: **5/5**, exit 0, RSS 67 MB; todos los números idénticos a los del worker.
- **Provenance**: 25 ocurrencias `\src`, 27 claves citadas, 26 únicas; 26 entradas en
  `provenance.json`; **0 faltantes, 0 huérfanas**. Los 19 literales de >=6 decimales de
  `notes.tex` resuelven todos en `checks.py --values` y/o su forma cerrada.
- **LaTeX**: compila en 2 pasadas, **6 páginas**, 0 errores, 0 referencias indefinidas,
  2 `Overfull` tipográficos.
- **Higiene**: ningún archivo modificado fuera del subárbol del trabajo desde las 20:52;
  SPECTER intacto.

## Refutado

**La frase de `notes.tex` §Laplace: "Así se cubren sin cancelación tanto k=0 como valores
positivos muy pequeños" es falsa tal como está escrita.** El guardián `k*Lz < 0.5` no es el
criterio correcto. Contraejemplo hallado por barrido dirigido:

    k = 5e-11, Lz = 1e10  (kLz = 0.5 -> rama exponencial), g0 = 1.0, g1 = 7.3
      |phi(0) - g0| = 7.629e-06        <-- 4 ordenes por encima de la tol 1e-10 de la puerta

Mecanismo: en la rama exponencial `C1 ~ g1/k` y `C2 ~ -r g1/k` se cancelan al formar
`phi(0) = C1 r + C2`, con pérdida relativa `~ eps * (e^{-kLz}/(kLz)) * Lz|g1|/|g0|`. El grupo
peligroso es `Lz|g1|/|g0|`, no `kLz`. Modelo y medición coinciden en orden (1.9e-05 vs 7.6e-06).

Alcance honesto: **esto no toca nada de lo contratado**. El pedido era estabilidad para
`k*Lz` grande, y eso se cumple (verificado hasta `kLz=1e8`, residuo 0). En SPECTER `k` está
cuantizado como múltiplo de `2 pi/Lx`, así que `Lz|g1|/|g0| ~ 1e10` no ocurre. Es la
**redacción** la que promete más de lo que el código entrega. Arreglo mínimo: reescribir la
frase acotándola a `kLz` moderado/grande, o cambiar el umbral por uno sobre `Lz|g1|/|g0|`.

## No resuelto (unclear)

1. **Las desigualdades de validez nunca se evalúan para la celda de ~6 mm que motiva el
   trabajo.** `Re_h eps << 1` equivale a `U << nu*l/h^2`; con `h=6 mm` y `nu=1e-6 m^2/s`
   eso es `U << 0.028*l` (SI), o sea `U << 0.56 mm/s` para `l = 2 cm`. En el trabajo no hay
   ningún `U` ni `l` experimental, así que **no puedo decidir** si la clausura aplica al
   experimento — pero la pregunta central del `intent` ("¿se puede calcular alfa?") depende
   de esa aritmética y la nota no la hace. Hace falta un número humano (`U` típica y la menor
   escala horizontal que el modelo quiera conservar).
2. **Zona de continuación FC-Gram, para la fase Fortran.** `laplace_z` recorre `DO k=1,nz`
   con `z(k)=rmp*(k-1)`, `rmp=Lz/(nz-Cz-1)` (`specter.fpp:744-749`): el dominio físico es
   `k=1..nz-Cz` y los `Cz` puntos restantes tienen `z > Lz`. `capa.py::_validar_z` **rechaza**
   `z > Lz`, y la afirmación de estabilidad de la nota está acotada a `[0,Lz]`. Allí
   `C1 e^{k(z-Lz)}` crece como `e^{kLz*Cz/(nz-Cz-1)}`. Es la misma estructura que las ramas
   DD/NN ya existentes, así que no es un defecto nuevo, pero ni la nota ni el módulo lo
   mencionan y **no lo verifiqué**.
3. **Punteros de provenance más anchos que el check citado.** `espectro_dd` y `espectro_dn`
   apuntan a `checks.py::test_spectra_against_fd`, que sólo prueba las cuatro tasas más
   lentas: no ejercita las autofunciones ni los coeficientes `a_n`, `b_m` que el mismo
   `statement` cubre. Yo los verifiqué aparte (Gram + reconstrucción), así que la física
   está bien; lo que sobra es la promesa del puntero.
4. **Margen de C1 en la puerta: 8.23e-05 contra tol 1e-04 (factor 1.2).** Es estructural de
   la puerta, no del módulo: la referencia FD con `nz=400` tiene error relativo
   `(4 pi)^2 dz^2/12 = 8.224e-05` en el cuarto modo DD, que reproduje analíticamente. Es
   determinista y siempre pasará, pero no hay margen si alguien sube `n`.
5. **Fuera de alcance por contrato, no verificado por mí**: nada experimental, ningún Fortran,
   ningún MPI, ninguna DNS 3D, ninguna deformación de superficie. La sección de salvedades de
   `notes.tex` es sustantiva y correcta; no la encontré inflada ni vacía.

## Veredicto

La entrega es correcta en todo lo contratado: los dos espectros completos, `alfa` y el
factor cuatro, el sesgo `pi/2`, la rama Laplace DN estable y las convenciones de SPECTER
resisten reproducción por métodos independientes. Puerta 6/6 y checks 5/5 confirmados por
mí, no heredados. El único enunciado que falla es una frase sobreextendida sobre la
ausencia de cancelación; y la brecha más importante para el proyecto no es un error sino
una omisión: la nota no evalúa sus propias desigualdades de validez con los parámetros de
la celda de 6 mm.
