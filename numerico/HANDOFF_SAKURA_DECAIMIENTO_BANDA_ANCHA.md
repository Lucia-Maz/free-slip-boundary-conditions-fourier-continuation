# Instrucciones para retomar en Sakura: convergencia del decaimiento de banda ancha

Este archivo está escrito para pegarlo completo como primer mensaje de una sesión nueva de
un agente abierta **dentro del clon de este repositorio en Sakura**. Fecha del estado local:
2026-09-22.

## Tu misión

Continuá la prueba puente [V-19] entre el decaimiento PIV y SPECTER. La laptop ya completó
16² y 32²; en Sakura tenés que medir primero el costo y después correr la convergencia 64².
No corras 128² salvo que 64² no sostenga la conclusión y hayas presentado el costo y
recibido aprobación explícita de Lucía.

No abras conexiones hacia otra máquina ni mandes mensajes externos. Todo cómputo pesado se
hace mediante SLURM. No ejecutes SPECTER en el nodo de login.

## Antes de tocar nada

1. Leé completos `CLAUDE.md`, `ESTADO.md`, `DECISIONES.md` desde [V-18] en adelante,
   este archivo, `numerico/CLUSTER.md`, `verificacion/decaimiento_banda_ancha.py` y
   `verificacion/decaimiento_banda_ancha.sbatch`.
2. Corré `git status --short --branch`, `git remote -v` y `git log -3 --oneline`. No
   descartes ni sobrescribas cambios locales. La historia se reescribió el 2026-09-17
   ([D-47]); si el clon de Sakura todavía tiene la historia vieja, `git pull` no puede
   reconciliarla. Hacé `git fetch origin`, mostrale a Lucía la divergencia y pedí permiso
   antes de cualquier `git reset --hard`. Si sólo falta tracking, fijalo con
   `git branch --set-upstream-to=origin/main main` después de verificar que ésa es la rama
   correcta.
3. La revisión que contiene esta campaña debe incluir, como mínimo:

   - `verificacion/decaimiento_banda_ancha.py`
   - `verificacion/decaimiento_banda_ancha.sbatch`
   - `salidas/tablas/decaimiento_banda_ancha.json`
   - `codigo/14_figura_decaimiento_banda_ancha.py`
   - `numerico/HANDOFF_SAKURA_DECAIMIENTO_BANDA_ANCHA.md`
   - la entrada [V-19] de `DECISIONES.md`

   Si no están, no los reconstruyas de memoria: el trabajo local todavía no llegó al
   remoto. Informalo y esperá a que Lucía sincronice el repositorio.
4. Las puertas `verificacion/test_aceptacion_fase1.py` y
   `verificacion/test_aceptacion_fase2.py` están congeladas por hash: se leen e importan,
   pero **no se editan**. Tampoco se toca `numerico/SPECTER-upstream/` ni se regeneran
   datos PIV crudos.
5. Verificá antes de enviar un trabajo:

   ```bash
   python -m py_compile verificacion/decaimiento_banda_ancha.py
   bash -n verificacion/decaimiento_banda_ancha.sbatch
   python verificacion/higiene.py --autotest
   mkdir -p verificacion/logs
   squeue -u "$USER"
   ```

   No compartas un nodo con otro trabajo propio. Evitá `l1-l3`; si SLURM fuera a elegir
   un nodo inconveniente, proponé una selección compatible con las reglas de `CLAUDE.md`.

## Qué se está probando

La condición inicial es un campo horizontal solenoidal, de fases deterministas, con los
pesos de los diez anillos medidos por PIV a `t ≈ 11,13 s`:

`u_h(x,y,z) = sin(pi z / 2h) grad_perp psi(x,y)`, `w = 0`.

Se corre exactamente la misma forma espectral en dos amplitudes:

- `lineal`: amplitud `10^-4` de la experimental, control del propio discretizado;
- `experimental`: `delta_0 = 8,887`, para permitir transferencia no lineal.

Las condiciones de borde están fijadas y no son una decisión pendiente:

- `x,y`: periódicas;
- `z=0`: no deslizamiento, `u=v=w=0`;
- `z=h`: free-slip plano, `w=0`, `∂z u=∂z v=0`;
- presión en `z`: Neumann–Neumann, como requiere el paso de proyección implementado.

La malla vertical queda en `NZ=64`, `CZ=25`, es decir 39 puntos físicos. El control lineal
ya muestra que alcanza para esta pregunta; la convergencia pendiente es horizontal. A 32²
los anillos 1–10 entran, pero el anillo 10 está cerca del corte de dealiasing. A 64² los
mismos anillos quedan holgados dentro del rango retenido.

## Resultado local que debés poder reproducir leyendo el JSON

El par 32² tardó 58 min 28 s en total en la laptop, secuencial, con pico de 152 MB. No lo
repitas en Sakura salvo que el JSON falte o esté corrupto.

| magnitud | control lineal 32² | banda ancha 32² |
|---|---:|---:|
| `lambda_sup` [s⁻¹] | 0,085282 | 0,083970 |
| intercepto del ajuste modal [s⁻¹] | 0,067767 | 0,065438 |
| pendiente `nu` [m²/s] | 0,9868·10⁻⁶ | 1,1715·10⁻⁶ |
| R² de `lambda(k)` | 0,99992 | 0,97842 |

El control recupera `alpha = 0,068539 s⁻¹` y `nu = 1,0·10⁻⁶ m²/s`. La no linealidad baja
el intercepto sólo `0,00233 s⁻¹` y la tasa global de superficie sólo 1,5 %. En cambio, el
PIV da un intercepto `0,04805 ± 0,00378 s⁻¹`, o sea un déficit de
`0,02048 ± 0,00378 s⁻¹` respecto del nominal. En 32², por lo tanto, la transferencia de
una banda ancha con la intensidad medida es real pero demasiado pequeña para explicar
[P-02]. Falta demostrar que esa conclusión no depende del corte de 32².

Controles ya satisfechos en 32²: energía de superficie monótona; pesos iniciales de los
diez anillos iguales a los prescritos a precisión de máquina; residuo free-slip menor que
`4·10^-37`; `delta` del caso experimental baja de 8,887 a 0,447 durante la ventana; la
energía vertical permanece menor que `10^-3` al primer campo ajustado y luego decrece.

## Ejecución autorizada en Sakura

El `sbatch` ya aplica los módulos, arma el prefijo sintético del toolchain, usa un scratch
propio bajo `/share`, limpia sólo ese subdirectorio y hace que el runner lance SPECTER con
`srun --mpi=pmix`. No copies resultados por `scp` ni escribas la corrida contra `$HOME`.

### 1. Prueba corta de 64²

Enviá:

```bash
PRUEBA=1 RESOLUCION=64 sbatch verificacion/decaimiento_banda_ancha.sbatch
```

Registrá el job ID. Al terminar, revisá el log y `sacct` (estado, elapsed, MaxRSS). La
prueba corta sólo valida compilación, lanzamiento, E/S, seis campos y costo; sus tasas no
son resultado físico. Comprobá que apareció `n64_prueba` en
`salidas/tablas/decaimiento_banda_ancha.json` y que no alteró `n16`, `n16_prueba` ni `n32`.

### 2. Par completo de 64²

Si la prueba corta termina limpia y el walltime de 8 h es suficiente con margen, enviá:

```bash
RESOLUCION=64 sbatch verificacion/decaimiento_banda_ancha.sbatch
```

Los casos son secuenciales y el JSON se guarda después de cada uno. Si el trabajo termina
por walltime después del control lineal, no uses `--repetir`: reenviar el mismo comando
reanuda sólo el caso que falta. No lances dos trabajos que puedan escribir el JSON a la
vez.

La extrapolación cruda desde 32² es unas 4,5–5 h para el par 64², pero el dato válido será
el de la prueba corta en el nodo asignado. Informá si la predicción cambia de forma
material antes de gastar más.

### 3. Auditoría de la salida

No te limites a los dos escalares impresos. Comprobá por código y reportá:

- que los dos casos están completos, sin NaN/Inf;
- energía de superficie estrictamente decreciente en los seis campos;
- pesos iniciales por anillo contra `config.pesos_usados`;
- `residuo_freeslip_max`, `w2_sobre_v2`, `delta`, `k_rms_fisico_1_m` y
  `alfa_eff_sobre_alfa` en el tiempo;
- tasas de cada anillo y dispersión del ajuste temporal;
- ajuste `lambda(k)=intercepto+pendiente*k²` de los anillos 3–10;
- `delta_alpha = alpha_banda-alpha_lineal` y
  `lambda_sup_banda/lambda_sup_lineal`;
- diferencias 64²−32² de esas dos magnitudes y de cada tasa por anillo;
- control lineal contra `alpha=0,06853891945 s⁻¹` y `nu=10^-6 m²/s`.

Como regla para decidir si hace falta 128²: 64² sostiene la conclusión si el control
lineal conserva R² > 0,995, su pendiente queda a menos de 3 % de `nu`, y el cambio de
`delta_alpha` entre 32² y 64² es menor que `0,002 s⁻¹` (aproximadamente media SEM del
intercepto PIV). Mirá también los anillos 9–10: un cambio localizado allí puede ser sólo la
liberación del corte de dealiasing de 32² y no invalida por sí solo el intercepto.

Si esa regla no se cumple, **no envíes 128² todavía**. Estimá memoria y walltime desde el
64² real, explicá qué observable no convergió y pedí aprobación. La extrapolación local
sugiere del orden de 20–24 h por el par, pero no es una reserva ni una autorización.

## Cierre y entrega de la sesión

Si 64² converge:

1. agregá una entrada [V-20] a `DECISIONES.md` con procedencia, costo, controles y cifras;
2. actualizá `ESTADO.md` para cerrar la prueba puente y decir con precisión qué queda
   abierto de [P-02];
3. no intentes generar `f10` en Sakura si el módulo Python no trae matplotlib: el script
   `codigo/14_figura_decaimiento_banda_ancha.py` selecciona automáticamente la resolución
   completa más fina y se regenerará en la laptop;
4. corré la higiene del árbol y revisá `git diff --check`, `git diff --stat` y
   `git status --short`;
5. mostrá el resultado y los cambios a Lucía antes de cualquier 128². Sólo hacé commit y
   push si ella te lo pide en esa sesión.

Tu informe final debe distinguir tres afirmaciones: (a) qué demostró el control lineal,
(b) cuánto cambia la transferencia no lineal el observable a resolución convergida y
(c) si esa magnitud alcanza o no para explicar el déficit experimental. No cierres [P-02]
por descarte: aunque esta hipótesis falle, el espesor real de la capa y los sesgos del PIV
siguen siendo candidatos independientes.
