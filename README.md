# Δt óptimo en PIV de turbulencia cuasi-2D

Proyecto final del curso *Ondas Gravitacionales e Investigación Asistida por IA*,
sobre datos propios: la campaña de PIV del **2 de junio de 2025** en la celda de
forzado electromagnético.

## La pregunta

En esos registros el desplazamiento de las partículas entre cuadros consecutivos es de
**0,09 a 0,28 píxeles**, y el piso de ruido típico de PIV es de orden 0,1 px. La
relación señal/ruido en desplazamiento es del orden de 1, y en buena parte de un
decaimiento es menor que 1.

La pregunta es entonces cuál es el Δt correcto, cómo se lo elige a partir de las
imágenes, y cuánto cambia lo que ya está medido. El punto de partida es que **no hace
falta medir de nuevo**: cada corrida tiene 3072 cuadros consecutivos a 60 fps, así que
los pares se pueden rearmar con cualquier separación Δt = n/60 s.

## Qué hay acá

| | |
|---|---|
| `DECISIONES.md` | **el log de decisiones.** Qué se decidió, por qué, y qué alternativa se descartó. Se agrega al final, no se reescribe |
| `bibliografia/` | los papers y sus DOI, **resueltos contra Crossref**, no citados de memoria. `REFERENCIAS.md` dice para qué se usa cada uno y cuáles todavía no fueron leídos en el original |
| `codigo/` | los scripts. `piv_comun.py` lee la campaña; `02_barrido_dt.py` hace el barrido |
| `salidas/tablas/` | resultados numéricos en JSON, con la configuración que los produjo adentro del mismo archivo |
| `salidas/figuras/` | las figuras |
| `datos-derivados/` | caches intermedios, reconstruibles |

## Los datos

No están en este repositorio: son ~140 GB en un disco externo, y se tratan como sólo
lectura ([D-02]).

```
/mnt/usb-.../Lucía Mazaira/02-06-25/
```

Once mediciones crudas de 3072 cuadros (salvo `med_S0001`, la grilla de calibración,
y `med_S0010`, guardada desde el cuadro 500), más seis carpetas con salidas previas de
OpenPIV. El mapa carpeta → medición está en `piv_comun.MEDICIONES` y sale de la entrada
del 2 de junio del cuaderno de laboratorio; la verificación es que los timestamps de
los archivos siguen exactamente ese orden.

**Condiciones**: Photron FASTCAM-1024PCI, 1024×1024, 10 bits, 60 fps, shutter 1/60 s.
Escala 3800 px/m (0,263 mm/px, campo de 26,9 cm). Capa de 6 mm de KNO₃ al 16 % m/m con
partículas de 100 µm, sobre imanes de 1 cm en red tipo tablero de 1,5 cm entre centros
(5 mm de vacío entre imanes, acrílico de 6 mm hasta la capa; [D-41]).
Corrientes de forzado: 0,9–1,0 A (16 V), 1,9–2,0 A y 2,1–2,2 A (32 V).

Estos valores **cambian entre experiencias** y por eso están declarados una sola vez, con
su procedencia, en `codigo/celda.py` ([D-27]). Lo de acá es una copia de lectura: el valor
que vale es el de ese archivo.

## Reproducir

```bash
mamba create -n piv-dt -c conda-forge python=3.11 numpy scipy matplotlib \
      tifffile imageio scikit-image pip
conda activate piv-dt && pip install openpiv

python codigo/02_barrido_dt.py --carpeta med_S0008 --pares 12
```

Con el disco desmontado nada de esto corre: `piv_comun.RAIZ` apunta al punto de
montaje y los scripts fallan con `FileNotFoundError`, a propósito, en vez de seguir
con datos parciales.

**Memoria**: la máquina tiene 5,7 GB. Los scripts procesan un par de imágenes por vez
y no retienen ningún campo completo, sólo escalares. Conviene igual correrlos de a uno
y con tope duro:

```bash
systemd-run --user --scope --quiet -p MemoryMax=3G -p MemorySwapMax=0 \
  python codigo/02_barrido_dt.py
```

## El modelo que se está poniendo a prueba

Para una corrida **forzada**, que es estacionaria, la velocidad verdadera no depende de
la separación entre cuadros. Si el error de PIV es un error de desplazamiento de
desviación σ en píxeles, independiente de Δt porque es propiedad de la correlación y no
del flujo, entonces

$$U_{med}(n)^2 = U_{verdadera}^2 + \left(\frac{\sigma}{\Delta t(n)\, s}\right)^2,
\qquad \Delta t(n) = n/\text{fps}$$

o sea una **recta** al graficar $U_{med}^2$ contra $n^{-2}$: la ordenada al origen es la
velocidad verdadera y la pendiente da σ. A n grande la recta se rompe por pérdida de
correlación, y ahí crecen los vectores inválidos. Los dos extremos del rango útil de Δt
salen de la misma figura.

El σ obtenido así se contrasta contra **0,087 px**, el valor que Foucaut, Carlier &
Stanislas miden para ventana de 16×16 px
([10.1088/0957-0233/15/6/003](https://doi.org/10.1088/0957-0233/15/6/003)).
Un cálculo preliminar sobre `med_S0009`, usando los dos procesamientos que Lucía ya
tenía con distinto Δt, dio **0,083 px**.

---

## Corrección al modelo (2026-09-02)

La sección anterior dice que σ es "independiente de Δt porque es propiedad de la
correlación y no del flujo". **Eso resultó falso**, y es el hallazgo principal hasta
ahora. Tres mediciones de σ sobre la misma cámara y la misma configuración de PIV,
que se diferencian sólo en qué degradación física incluyen:

| medición | σ | qué agrega respecto de la anterior |
|---|---|---|
| par sintético, corrimiento conocido | ≤ 0,011 px | correlación y estimador subpíxel solos |
| fluido en reposo, par real | 0,030 px | ruido de lectura, iluminación, movimiento de partículas |
| forzado 2,15 A | 0,169 px | gradiente de velocidad dentro de la ventana, movimiento fuera del plano |

El término que agrega el flujo es √(0,169² − 0,030²) = 0,166 px: el **98 % de la
varianza**. Es la descomposición de George & Stanislas (arXiv 2010.10768) — ruido de
pixelización más ruido por no uniformidad de la velocidad dentro del volumen de
interrogación — y el registro en reposo mide sólo el primer término.

Consecuencia: en un decaimiento, donde U cae un factor ~10, σ tampoco es constante a lo
largo del registro, y la corrección a la curva no puede usar un solo número. El ajuste
del barrido sigue siendo válido como descripción de *esta* corrida, pero su σ es un
valor efectivo, no una constante del instrumento.

Ver `DECISIONES.md`, entradas [V-04], [H-06] y [H-07], con las salvedades — que son
fuertes: el registro en reposo es de otra campaña, y el 0,169 px sale de un ajuste cuyo
modelo esta misma corrección declara inválido.

---

## Reproducir en otra máquina

El repositorio es autocontenido salvo por tres cosas que, a propósito, no están adentro:
los datos crudos (140 GB, disco externo), el árbol de SPECTER (código de terceros) y los
PDF de la bibliografía (son de los editores).

```bash
git clone git@github.com:Lucia-Maz/proyecto-final-piv.git
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
python verificacion/test_aceptacion_fase1.py   # 6/6
python verificacion/test_aceptacion_fase2.py   # 6/6
```

**Qué corre sin el disco externo:** todo lo numérico y lo teórico —las dos puertas, el
barrido en δ, el chequeo de la capa límite de divergencia, y todo lo que lee de
`salidas/tablas/`—. Los scripts `codigo/0[1-9]_*.py` que rearman PIV desde los TIF crudos
necesitan el disco montado; sus resultados ya están volcados en `salidas/tablas/*.json`,
así que las figuras y los informes se regeneran sin él.

**Convenciones del proyecto:** están en `CLAUDE.md`, y valen igual en la laptop y en el
clúster.

### La primera vez en el clúster

El repositorio es **privado**, así que hay que autenticarse. Se usa una **deploy key**:
una clave SSH cuyo alcance es este repositorio y nada más. Sakura es una máquina
compartida, y una clave de cuenta sin passphrase ahí daría acceso a todo el GitHub; la
deploy key, si se filtra, sólo alcanza a este proyecto.

```bash
# 1. ¿Sale tráfico hacia GitHub? Muchos clústeres bloquean el 22 de salida.
ssh -T git@github.com          # "Permission denied (publickey)" YA ES BUENA SEÑAL:
                               # significa que llegó. Si queda colgado o da timeout,
                               # ver el punto 5.

# 2. Generar la clave, sin passphrase para no tipearla en cada push
ssh-keygen -t ed25519 -C "sakura-proyecto-final-piv" -f ~/.ssh/id_ed25519 -N ""
cat ~/.ssh/id_ed25519.pub

# 3. Cargarla en el REPOSITORIO (no en la cuenta):
#    github.com/Lucia-Maz/proyecto-final-piv -> Settings -> Deploy keys
#    -> Add deploy key -> pegar -> MARCAR "Allow write access" -> Add
#    Sin esa marca queda de sólo lectura y no se puede pushear.
#    Desde una máquina con gh autenticado es equivalente:
#      gh repo deploy-key add ~/.ssh/id_ed25519.pub \
#         --repo Lucia-Maz/proyecto-final-piv --title sakura --allow-write

# 4. Probar. Con deploy key el mensaje nombra al REPO, no al usuario, y eso es correcto:
ssh -T git@github.com
#    Hi Lucia-Maz/proyecto-final-piv! You've successfully authenticated, but GitHub
#    does not provide shell access.

# 5. Si el puerto 22 está bloqueado, GitHub escucha SSH también en el 443:
cat >> ~/.ssh/config <<'CFG'
Host github.com
  Hostname ssh.github.com
  Port 443
  User git
CFG
```

**No editar `~/.ssh/config` en sakura.** Ese archivo lo genera Warewulf, el sistema de
aprovisionamiento del clúster, y trae

```
Host *
   IdentityFile ~/.ssh/cluster
   StrictHostKeyChecking=no
```

Dos motivos para no tocarlo: en cuanto hay un `IdentityFile` explícito, ssh **deja de
probar las claves por omisión**, que es exactamente por qué la clave nueva no se ofrecía;
y si reimaginan el nodo, cualquier edición se pierde. La configuración se pone **por
repositorio**, que vive en `.git/config` y sobrevive.

Clonar **por SSH** —las deploy keys no funcionan por HTTPS— indicando la clave sólo para
ese comando, y después fijarla en el repo:

```bash
GIT_SSH_COMMAND="ssh -i ~/.ssh/id_ed25519 -o IdentitiesOnly=yes" \
  git clone git@github.com:Lucia-Maz/proyecto-final-piv.git

cd proyecto-final-piv
git config core.sshCommand "ssh -i ~/.ssh/id_ed25519 -o IdentitiesOnly=yes"
git config user.name  "Lucia-Maz"       # local al repo: la deploy key autentica la
git config user.email "<tu correo>"   # máquina, no firma la autoría

module load gnu15 openmpi5 fftw/3.3.11        # ver CLAUDE.md
# y después la reconstrucción de SPECTER de la sección anterior
```

Desde ahí `git pull` y `git push` andan solos dentro de ese directorio, sin afectar las
conexiones a los nodos de cómputo ni a ningún otro repositorio.

Si alguna vez querés clonar **otro** repo tuyo en Sakura, hace falta otra clave distinta
—GitHub rechaza la misma deploy key en dos repositorios— más un alias por repo en
`~/.ssh/config`.

### Trabajar en las dos máquinas

El repositorio es la única fuente de verdad; no se copian archivos por `scp`.

```bash
git pull                       # SIEMPRE antes de empezar
# ... trabajar ...
git add -A && git commit -m "..." && git push
```

Lo que **no** viaja por git y hay que tener en cuenta: los datos crudos (disco externo,
sólo en la laptop), el árbol de SPECTER (se reconstruye con el parche) y los PDF de la
bibliografía. Los resultados sí viajan, en `salidas/tablas/*.json`, así que las figuras y
los informes se regeneran de un lado o del otro indistintamente.
