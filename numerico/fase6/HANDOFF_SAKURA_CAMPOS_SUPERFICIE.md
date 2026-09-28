# Instrucciones para Sakura: campos de superficie con la superficie deformable de orden 2

Este archivo está escrito para pegarlo completo como primer mensaje de una sesión nueva de
un agente abierta **dentro del clon de este repositorio en Sakura**. Estado local: 2026-09-28.

## Tu misión, y sólo ésa

Correr **una** simulación con la superficie deformable de orden 2 a la escala de la celda, y
subir sus resultados chicos al repositorio:

- `salidas/tablas/fase6_campos_superficie_n64_m4.json`
- `salidas/campos/fase6_campos_superficie_n64_m4.npz`

La corrida de la laptop, de 1,5 memorias (`…_n64_m1.5.*`), ya está en el repositorio: no la
toques.

La corrida la define `numerico/fase6/campos_superficie.py`. Es la condición inicial de banda
ancha de [V-19], con el tope `freesurface`, `fsorder = 2`, y g y σ/ρ físicos. Ya se corrió
en la laptop con 1,5 memorias. Acá se pide la versión más larga: 4 memorias, 12 campos. No
cambies la física ni el código de SPECTER, y no corras nada más.

**Reglas:**
- No abras conexiones hacia otras máquinas ni mandes mensajes externos.
- Todo el cómputo va por SLURM. No ejecutes SPECTER en el nodo de login.
- No edites las puertas `verificacion/test_aceptacion_fase*.py`, que están congeladas por hash.

## Pasos

1. **Leé** `CLAUDE.md` (sobre todo «Reglas del clúster»), este archivo,
   `numerico/fase6/campos_superficie.py`, `numerico/fase6/campos_superficie.sbatch` y
   `numerico/specter-parche/README.md`.
2. **Traé el estado.** Corré `git status --short --branch` y `git fetch origin`, y hacé
   `git pull --ff-only`. Si no avanza en fast-forward, **pará** y mostrale la divergencia a
   Lucía: no hagas `reset --hard`. El commit tiene que contener
   `numerico/fase6/campos_superficie.py` y este archivo.
3. **Reconstruí el árbol de SPECTER desde upstream más el parche actual.** El de Sakura puede
   estar viejo:
   ```bash
   cd numerico
   [ -d SPECTER-upstream ] || git clone https://github.com/mfontanaar/SPECTER.git SPECTER-upstream
   (cd SPECTER-upstream && git checkout 0ad1edb && git status --short | head)   # tiene que estar limpio
   [ -d SPECTER-trabajo ] && mv SPECTER-trabajo SPECTER-trabajo.previo-$(date +%Y%m%d-%H%M)   # no se borra
   cp -r SPECTER-upstream SPECTER-trabajo
   cd SPECTER-trabajo && git apply ../specter-parche/superficie-libre.patch
   cp ../specter-parche/laplace_dirneu.f90 ../specter-parche/laplace_neudir.f90 \
      ../specter-parche/fs_orden.f90 src/tests/
   cp ../specter-parche/fsorder.f90 src/boundary/
   cd ../..
   sha256sum numerico/specter-parche/superficie-libre.patch   # a1d31136b1543cec… (README del parche)
   ```
4. **Chequeos rápidos,** en el login, que no corren SPECTER:
   ```bash
   module load python/3.13.13
   python -m py_compile numerico/fase6/campos_superficie.py
   bash -n numerico/fase6/campos_superficie.sbatch
   python verificacion/higiene.py --autotest
   mkdir -p verificacion/logs && squeue -u "$USER"
   ```
5. **Una prueba corta primero,** que compila y corre un minuto con el mismo `sbatch`:
   ```bash
   PRUEBA=1 sbatch --time=00:30:00 numerico/fase6/campos_superficie.sbatch
   ```
   Tiene que terminar con una línea `corrida: … avisos: 0`. Si falla, mostrale el log a
   Lucía y pará. No corrijas el código.
6. **La corrida:**
   ```bash
   NXY=64 MEMORIAS=4 CAMPOS=12 sbatch numerico/fase6/campos_superficie.sbatch
   ```
   - Estimación: unas 10 000 pasos a ~1 s por paso, del orden de 3 h. Las transformadas de
     superficie no se aceleran con más procesos: el sbatch pide 2, no pidas más.
   - Seguila con `squeue` y con `tail` del log. Es normal que el log quede casi vacío hasta
     el final.
   - Si aparece `[WARNING] freesurface` o `[ERROR]`, anotalo tal cual.
7. **Verificá el resultado** antes de subir nada:
   - el JSON existe, `avisos` está vacío o documentado, y `instantes` tiene 13 entradas;
   - `eta_rms_m` es del orden de 10⁻⁷–10⁻⁶ m;
   - `max_eta_sobre_dz` es menor que 0,05;
   - `serie_diagnostico.pasadas` es menor que 30 y `residuo` menor que 10⁻⁸ en toda la serie.
8. **Subí** sólo esos dos archivos, más una línea en `ESTADO.md` bajo «Hecho el 2026-09-28,
   Fase 6» con el trabajo de SLURM, el tiempo y el pico de memoria:
   ```bash
   git add salidas/tablas/fase6_campos_superficie_n64_m4.json salidas/campos/fase6_campos_superficie_n64_m4.npz ESTADO.md
   git commit -m "Campos de superficie de orden 2, 64², 4 memorias, en Sakura (trabajo <id>)"
   git push      # corre el chequeo de higiene antes
   ```
   Si el push falla por el hook de higiene, no lo saltees: mostrale a Lucía qué encontró.

## Qué no hace falta

La figura (`numerico/fase6/figura_campos_superficie.py`) se hace en la laptop, porque necesita
matplotlib. Tampoco hace falta la puerta de la Fase 6: da 4/10 por diseño, ver [V-21]. Para
confirmar que el árbol quedó bien alcanza la prueba del paso 5.
