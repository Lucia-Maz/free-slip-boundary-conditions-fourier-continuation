# El clúster: la primera vez, y cómo trabajar entre las dos máquinas

Esto es operativo, no es para el lector del proyecto. Las reglas del clúster tal como las
dieron los administradores, y qué implican para SPECTER, están en `CLAUDE.md`; la puerta
se manda con `verificacion/puerta_fase2.sbatch` y la cadena de compilación con módulos
está en `verificacion/toolchain_cluster.sh` ([D-33] a [D-39] en `DECISIONES.md`).

## Acceso al repositorio desde Sakura

Mientras el repositorio fue privado hizo falta autenticarse. Se usó una **deploy key**: una
clave SSH cuyo alcance es este repositorio y nada más. Sakura es una máquina compartida, y
una clave de cuenta sin passphrase ahí daría acceso a todo el GitHub; la deploy key, si se
filtra, sólo alcanza a este proyecto. Con el repositorio público se puede clonar por HTTPS
sin nada de esto, pero para **pushear** desde el clúster la deploy key sigue siendo la
forma prudente.

```bash
# 1. ¿Sale tráfico hacia GitHub? Muchos clústeres bloquean el 22 de salida.
ssh -T git@github.com          # "Permission denied (publickey)" YA ES BUENA SEÑAL:
                               # significa que llegó. Si queda colgado o da timeout,
                               # ver el punto 5.

# 2. Generar la clave, sin passphrase para no tipearla en cada push
ssh-keygen -t ed25519 -C "sakura-free-slip-boundary-conditions-fourier-continuation" -f ~/.ssh/id_ed25519 -N ""
cat ~/.ssh/id_ed25519.pub

# 3. Cargarla en el REPOSITORIO (no en la cuenta):
#    github.com/Lucia-Maz/free-slip-boundary-conditions-fourier-continuation -> Settings -> Deploy keys
#    -> Add deploy key -> pegar -> MARCAR "Allow write access" -> Add
#    Sin esa marca queda de sólo lectura y no se puede pushear.
#    Desde una máquina con gh autenticado es equivalente:
#      gh repo deploy-key add ~/.ssh/id_ed25519.pub \
#         --repo Lucia-Maz/free-slip-boundary-conditions-fourier-continuation --title sakura --allow-write

# 4. Probar. Con deploy key el mensaje nombra al REPO, no al usuario, y eso es correcto:
ssh -T git@github.com
#    Hi Lucia-Maz/free-slip-boundary-conditions-fourier-continuation! You've successfully authenticated, but GitHub
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
  git clone git@github.com:Lucia-Maz/free-slip-boundary-conditions-fourier-continuation.git

cd free-slip-boundary-conditions-fourier-continuation
git config core.sshCommand "ssh -i ~/.ssh/id_ed25519 -o IdentitiesOnly=yes"
git config user.name  "Lucia-Maz"     # local al repo: la deploy key autentica la
git config user.email "<tu correo>"   # máquina, no firma la autoría

module load gnu15 openmpi5 fftw/3.3.11        # ver CLAUDE.md
# y después la reconstrucción de SPECTER del README (numerico/specter-parche/README.md)
```

Desde ahí `git pull` y `git push` andan solos dentro de ese directorio, sin afectar las
conexiones a los nodos de cómputo ni a ningún otro repositorio.

Si alguna vez hay que clonar **otro** repositorio en Sakura, hace falta otra clave
distinta —GitHub rechaza la misma deploy key en dos repositorios— más un alias por repo
en `~/.ssh/config`.

## Trabajar en las dos máquinas

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

## Para retomar con un agente en el clúster (2026-09-17; archivado)

> **No pegar este prompt para la campaña actual.** Conserva el plan anterior como registro.
> El handoff vigente, autocontenido y listo para copiar, es
> `numerico/HANDOFF_SAKURA_DECAIMIENTO_BANDA_ANCHA.md`.


`CLAUDE.md` se lee solo al abrir el proyecto. Lo que sigue es el pedido de la sesión de
producción; se pega tal cual como primer mensaje.

> El repositorio se llama ahora `free-slip-boundary-conditions-fourier-continuation`
> ([D-45]) y **su historia se reescribió el 2026-09-17** ([D-47]), así que el clon de acá
> no se actualiza con `git pull`: primero
> `git remote set-url origin git@github.com:Lucia-Maz/free-slip-boundary-conditions-fourier-continuation.git`,
> después `git fetch origin && git reset --hard origin/main` (no hay cambios locales
> que conservar allá), y activá el chequeo de higiene como dice `CLAUDE.md`
> (`git config core.hooksPath verificacion/hooks` y `verificacion/higiene.local` con tu
> usuario del clúster). En todo lo que escribas, el usuario del clúster es `$USER`. Leé `ESTADO.md` (sección "Lo que sigue") y, de
> `DECISIONES.md`, las entradas [D-42], [V-16], [V-17], [H-09] y [D-34] a [D-39]. Estás en
> Sakura: aplican las reglas del clúster de `CLAUDE.md` — todo por SLURM, módulos
> `gnu15 openmpi5 fftw/3.3.11 python/3.13.13`, `OMPI_MCA_pml=ob1 OMPI_MCA_btl=sm,self,tcp`,
> scratch en un subdirectorio propio bajo `/share/data2/$USER` que se limpia al terminar.
> Antes de elegir nodo mirá `squeue -u $USER`: mis trabajos propios no se comparten
> nodo; en la sesión anterior eso dio `--nodelist=g1` ([D-39]), hoy puede ser otro. La
> puerta de la Fase 2 ya dio 6/6 acá (trabajo 8168); no hace falta repetirla salvo que el
> árbol de SPECTER cambie.
>
> Objetivo: el punto 2 de "Lo que sigue", la corrida de producción, en tres pasos. No
> pasás al siguiente sin mi aprobación explícita.
>
> 1. **Convergencia del caso forzado a δ ≥ 11** ([H-09], [V-16]): correr
>    `verificacion/barrido_delta_etapa1.py --caso forzado --resolucion 64` (y 128 si 64 no
>    cierra) como trabajo de SLURM, con un `sbatch` nuevo calcado de
>    `verificacion/puerta_fase2.sbatch` — mismo toolchain sintético, mismo scratch, misma
>    limpieza con `trap`. El script importa los helpers de la puerta, así que hereda
>    `SPECTER_TOOLCHAIN` y `SPECTER_SCRATCH` y corre con un solo proceso MPI; si hace falta
>    más de uno, proponelo con el costo, no lo cambies solo. Sólo necesita numpy. Antes de
>    lanzar en serio, medí el costo con `--prueba` y decime cuánto va a tardar. El JSON va
>    a `salidas/tablas/` y se commitea; las figuras se hacen en la laptop, acá no hay
>    matplotlib.
> 2. **Diseño de la corrida de producción**, con el costo cuantificado y esperando mi OK
>    porque es una decisión científica: condición inicial (dos modos como el barrido, o
>    banda ancha), ε y δ de la celda real (ε = 1,78 en el fundamental de la red; δ ≈ 27
>    al inicio del decaimiento; la ventana del ajuste experimental va de δ ≈ 10 a ≈ 1,3
>    con ε bajando a ≈ 0,7), la resolución que la convergencia del paso 1 justifique, la
>    duración en memorias modales, y qué se mide: α_eff, λ_prom y λ_sup como en [V-16] y
>    [V-17]. El diagnóstico de tensión de `vboundary` pide una matriz del tamaño del campo
>    por `cstep` ([D-26]); si a esa resolución molesta, proponé condicionarlo antes.
> 3. **Lanzarla y registrarla**: entrada [V-18] en `DECISIONES.md` con procedencia, JSON en
>    `salidas/tablas/`, `ESTADO.md` al día, commit y push.
>
> Restricciones que no se renegocian: las puertas no se editan; los parámetros de la celda
> viven sólo en `codigo/celda.py`; no se escribe en `/share/data2/$USER` fuera del
> subdirectorio del trabajo; no se afirma nada que no se haya chequeado.
