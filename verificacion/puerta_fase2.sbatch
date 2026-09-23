#!/bin/bash
#SBATCH --job-name=puerta-fase2
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:40:00
#SBATCH --output=verificacion/logs/puerta-fase2-%j.out
#
# La puerta de aceptación de la Fase 2 en Sakura. Se manda al planificador y no se
# corre a mano: compila SPECTER y lanza unas quince simulaciones chicas, o sea diez a
# quince minutos de CPU en un solo core.
#
# Partición 'normal', que es la de omisión y la correcta acá. Medido con sinfo el
# 2026-09-09: 'normal' abarca a1, c1-c5, g1-g2 y sakura; 'compute' es un SUBCONJUNTO
# (c1-c5 y sakura, o sea que tampoco evita el nodo de login) y tiene menos prioridad,
# de modo que un trabajo mandado ahí queda esperando huecos que se reservan para
# 'normal'. Ver [D-36].
#
# NO usar 'legacy' (l1-l3) aunque figure 'idle': las reglas del clúster dicen de
# evitar esos nodos o preguntar antes.
#
#   sbatch verificacion/puerta_fase2.sbatch
#
# Tiene que terminar con RESULTADO: 6/6.

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

module load gnu15 openmpi5 fftw/3.3.11 python/3.13.13

# Transporte, según las reglas del clúster: red Ethernet de 10 Gb, sin InfiniBand.
# OJO: las instrucciones originales dicen 'sm,self,tpc'; se usa 'tcp', que es lo que
# tiene sentido y coincide con la frase que dice que el transporte entre nodos es TCP.
# Si esto fallara, probar exactamente como está en las instrucciones. Ver [D-33].
export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=sm,self,tcp

# Prefijo sintético: la puerta está fijada por hash y espera MPI y FFTW bajo un mismo
# prefijo; con módulos están separados. Ver [D-34].
#
# OJO ([H-15]): 'eval "$(cmd)"' NO propaga el fallo de cmd bajo 'set -e'. Si
# toolchain_cluster.sh falla, "$(...)" devuelve stdout vacío, 'eval ""' no hace nada y
# vale 0, y el script seguía de largo sin SPECTER_TOOLCHAIN exportado -- hasta reventar
# más abajo, en test_aceptacion_fase2.py, con un mensaje que apunta al entorno de conda
# de la laptop y no dice nada de lo que realmente falló. Por eso se captura aparte.
TOOLCHAIN_ENV="$(bash verificacion/toolchain_cluster.sh)"
eval "$TOOLCHAIN_ENV"

# La E/S no va a $HOME: desde los nodos de cómputo la E/S paralela contra $HOME es
# lenta.
#
# Destino: /share/data2/$USER, que es donde Lucía trabaja en el clúster (ahí están sus
# simulaciones de GHOST). Se puede pisar con SPECTER_SCRATCH_BASE.
#
# Confirmado desde el clúster el 2026-09-23 ([D-48]): scratch1 está al 100 % y data2
# tiene cientos de GB libres, sin cuota de usuario. El trabajo escribe en un
# subdirectorio propio y lo borra al terminar, así que no se mezcla con lo que ya haya.
elegir_base() {
    local d
    for base in ${SPECTER_SCRATCH_BASE:-} /share/data2 /share/scratch1; do
        [ -d "$base" ] || continue
        d="$base/$USER"
        if mkdir -p "$d" 2>/dev/null && [ -w "$d" ]; then echo "$d"; return 0; fi
    done
    return 1
}

if BASE="$(elegir_base)"; then
    echo "scratch: $BASE"
else
    # Para ESTA puerta no es grave: corre con -np 1 y arreglos de ~130 kB, así que la
    # advertencia sobre MPI-IO paralela no aplica. Para una corrida de producción sí,
    # y ahí hay que resolver el acceso a /share antes de lanzar.
    BASE="${TMPDIR:-/tmp}/$USER"
    mkdir -p "$BASE"
    echo "[aviso] no hay /share escribible; se usa $BASE." >&2
    echo "[aviso] alcanza para la puerta (un proceso, arreglos chicos), NO para producción." >&2
fi

# OJO: <base>/$USER es un directorio que YA EXISTE y donde Lucía guarda simulaciones.
# El trabajo escribe en un subdirectorio propio, con el número de trabajo en el nombre,
# y al terminar borra SÓLO ese. La limpieza verifica la forma del nombre antes de
# borrar: un rm -rf al lado de datos reales no puede depender de que las variables
# tengan el valor que uno espera.
export SPECTER_SCRATCH="$BASE/puerta-fase2-${SLURM_JOB_ID:-manual}"
mkdir -p "$SPECTER_SCRATCH"

limpiar() {
    case "${SPECTER_SCRATCH:-}" in
        */puerta-fase2-?*)
            [ -d "$SPECTER_SCRATCH" ] && rm -rf "$SPECTER_SCRATCH" ;;
        *)
            echo "[aviso] no borro '${SPECTER_SCRATCH:-}': no tiene la forma esperada" >&2 ;;
    esac
}
trap limpiar EXIT   # SPECTER_SCRATCH explícito no se autolimpia

python verificacion/test_aceptacion_fase2.py
