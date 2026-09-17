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
