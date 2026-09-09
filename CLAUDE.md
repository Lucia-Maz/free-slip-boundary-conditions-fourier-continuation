# Convenciones de este proyecto

Este archivo lo lee Claude al abrir el proyecto, en cualquier máquina. La idea es que
trabajar en la laptop o en el clúster sea indistinto.

## Por dónde empezar

Leer, en este orden: `ESTADO.md` (dónde estamos y qué sigue), `DECISIONES.md` (el log
completo, una entrada por decisión) y `bibliografia/REFERENCIAS.md`. Con eso alcanza para
retomar sin la conversación que lo produjo.

Las etiquetas `[D-nn]`, `[H-nn]`, `[V-nn]`, `[P-nn]` que aparecen por todos lados resuelven
en `DECISIONES.md`: decisión, hallazgo, verificación, pregunta abierta.

## En una línea

El decaimiento de una capa delgada está dominado por la fricción con el fondo, que nace de
la asimetría entre el fondo no-deslizante y la superficie libre. Este proyecto implementa
la condición de superficie libre en SPECTER para **calcular** α en vez de ajustarlo, y lo
contrasta contra el decaimiento medido por PIV.

## Reglas que no se renegocian sin decisión explícita

- **Los parámetros del experimento no se hardcodean.** El espesor de la capa, la
  viscosidad y la geometría de los imanes cambian entre experiencias y viven **sólo** en
  `codigo/celda.py`, con su procedencia. Las derivaciones se escriben en forma simbólica y
  los números aparecen sólo como evaluación al final ([D-27]).
- **Las puertas de aceptación están fijadas por hash.** `verificacion/test_aceptacion_fase1.py`
  y `..._fase2.py` se pueden leer e importar, **no se editan para que pasen**. Una puerta se
  escribe antes de la implementación y arranca en rojo.
- **Ninguna puerta contiene la fórmula analítica del resultado que verifica.** Construye su
  propia referencia numérica —diferencias finitas, extrapolación de Richardson— para que el
  resultado quede comprobado y no afirmado.
- **No afirmar lo que no se chequeó.** Nada de «nadie lo hizo», «no existe» o «es nuevo» sin
  haber hecho la búsqueda y poder mostrarla. Al apoyarse en bibliografía hay que dar la
  **sección** y, si es el nudo del argumento, la **frase textual**. Y marcar en cada paso
  qué viene del paper y qué se derivó acá.
- **La bibliografía distingue «localizado» de «leído en el original».** Ninguna afirmación
  se apoya en un trabajo que no se haya verificado en su texto. Los DOI se resuelven contra
  la API de Crossref, no se citan de memoria.
- **Los datos crudos del 02-06-25 son sólo lectura** ([D-02]). Es la única copia de una
  campaña que no se puede repetir.
- **Todo análisis de los datos experimentales termina en figuras** que se puedan mirar, no
  sólo en escalares. Mostrar el dato crudo junto al corregido.

## Las dos máquinas

**Laptop.** 5,7 GB de RAM y ~4 GB utilizables. Pasarse no significa ir lento: entra en
thrashing y se cuelga. Scripts pesados **secuenciales**, nunca en paralelo; `ulimit -v
2500000` antes de lanzar; medir con `/usr/bin/time -v` y reportar el pico. Sirve para
desarrollo, tests 1D y verificación a resolución chica. Ya mató tres corridas: si hay que
correr algo largo, conviene pedir que se cierre el navegador.

**Clúster (Sakura).** Es donde va el cómputo pesado: corridas de producción de SPECTER y
estudios de convergencia. Los límites de memoria de arriba **no** aplican ahí; sí aplica
todo lo demás de este archivo. Las reglas operativas están abajo.

### Reglas del clúster

Tal como las dieron los administradores, textuales y sin traducir, porque son operativas y
una mala traducción se paga con un trabajo mal lanzado:

> ## Cluster (sakura, SLURM)
> - Nodes: g1,g2: two Tesla P4 (cc61) each, partition `cuda`; a1,a2: two MI210 each,
>   partition `rocm`; c1-c5 and sakura (login) CPU nodes, `compute`/`normal`. Avoid
>   l1-l3 nodes or ask before using. When using CPUs in a1, a2, g1, g2, leave four
>   cores free. When using GPUs request one core per task (other jobs fill the cores,
>   GPUs stay free) and `--gres=gpu:N`; `srun --mpi=pmix`.
> - Network is 10 Gb Ethernet, no InfiniBand: TCP is the inter-node transport.
>   ROCm runs: `OMPI_MCA_pml=ucx UCX_TLS=rocm_copy,rocm_ipc,sm,self,tcp`.
>   On sakura under srun the gnu OpenMPI needs `OMPI_MCA_pml=ob1
>   OMPI_MCA_btl=sm,self,tpc`. Short CPU jobs can run on sakura (partition normal).
> - MPI-IO works on `/share/data*` and `/share/scratch*` (use a directory with
>   the user name there); from compute nodes to $HOME I/O is slow when parallel.
> - Use modules: `gnu15`, `openmpi5`, `fftw/3.3.11`, `cmake/4.1.2`, `python/3.13.13`
>   (has numpy/PyYAML), `rocm/7.2`, `cuda/12.4`, `nvhpc-stack/24.5`.

**Una errata aparente, que no corrijo por mi cuenta:** en `OMPI_MCA_btl=sm,self,tpc`, el
`tpc` es casi seguramente `tcp` — coincide con la frase anterior, que dice que el
transporte entre nodos es TCP. Antes de usarlo conviene confirmarlo con quien las escribió;
si `sm,self,tcp` funciona y `tpc` no, es eso.

### Qué implican para este proyecto

- **SPECTER en el clúster no usa el env de conda `specter`.** Ahí se compila con los
  módulos: `gnu15`, `openmpi5` y `fftw/3.3.11`. `Makefile.in` hay que apuntarlo al
  `FFTWDIR` del módulo en vez de al prefijo de conda ([D-23] describe el resto de la
  configuración, que no cambia).
- **Es CPU puro.** SPECTER en este proyecto se compila y corre en CPU, así que la partición
  es `compute` (o `normal` para pruebas cortas en el nodo de login). Los nodos con GPU no
  hacen falta, y si alguna vez se usan sus CPUs hay que dejar cuatro cores libres.
- **El scratch de las corridas va a `/share/scratch*`, no a `$HOME`.** SPECTER escribe
  salida binaria y desde los nodos de cómputo la E/S paralela contra `$HOME` es lenta. Usar
  un directorio con el nombre de usuario adentro.
- **`python/3.13.13` trae numpy pero no scipy ni matplotlib**, y este proyecto los necesita:
  `scipy.optimize` en `numerico/fase4/superficie_libre_escalas.py`, y matplotlib en todos
  los scripts de figuras. Para eso hace falta un venv propio sobre ese módulo, o correr las
  figuras en la laptop desde los JSON de `salidas/tablas/`, que es lo que el repositorio ya
  permite.
- **Las puertas de aceptación corren igual**, y son la prueba de que el árbol quedó bien
  armado: `verificacion/test_aceptacion_fase2.py` compila SPECTER en un scratch temporal y
  tiene que dar 6/6. Es un trabajo corto de un solo proceso, así que entra en `normal`.
- Al lanzar con `srun`, el binario de SPECTER se invoca con `--mpi=pmix`, y en sakura la
  OpenMPI de gnu necesita además las variables `OMPI_MCA_*` de arriba.

**El disco externo con los datos crudos está sólo en la laptop.** Los scripts que lo tocan
(`codigo/0*_*.py`, vía `piv_comun.RAIZ`) no corren en el clúster. Lo que sí viaja son los
JSON derivados de `salidas/tablas/`, que ya tienen los resultados.

## Entornos

```bash
# análisis de PIV y todos los scripts de codigo/
mamba create -n piv-dt -c conda-forge python=3.11 numpy scipy matplotlib \
      tifffile imageio scikit-image pip && pip install openpiv

# cadena de compilación de SPECTER ([D-23])
mamba create -n specter -c conda-forge gfortran openmpi openmpi-mpifort fftw make
```

En la laptop el intérprete es `/home/lucia/miniforge3/envs/piv-dt/bin/python`. **Usar la
ruta completa**, no `python` a secas: en esa máquina resuelve a un intérprete sin numpy, y
eso ya hizo figurar como fallado un chequeo que estaba bien.

Trampa del Makefile de SPECTER: usa `GHOME = $(PWD)`, y `$(PWD)` sale de la variable de
entorno, no del directorio del proceso. Lanzarlo con `subprocess.run(cwd=...)` sin exportar
`PWD` compila contra el árbol equivocado.

## SPECTER

El árbol no está en el repositorio. Se reconstruye desde upstream con el parche:
ver `numerico/specter-parche/README.md`. Después, `verificacion/test_aceptacion_fase2.py`
tiene que dar 6/6.

Convención que costó una iteración: **los números de onda de SPECTER son enteros**
(`specter.fpp:777-780`), así que el `Lx = 1` del archivo de parámetros son 2π de largo
físico y el modo fundamental es `k = 1`. `Lz`, en cambio, es la longitud literal ([H-13]).

## Higiene

El repositorio es privado hoy, pero la idea es publicarlo con el proyecto. No deben entrar
contraseñas ni credenciales del clúster; las de Sakura están en texto plano en las notas de
Obsidian de la laptop y **no van acá**. Tampoco entran los PDF de la bibliografía, que son
de los editores: `bibliografia/REFERENCIAS.md` dice de dónde bajar cada uno.
