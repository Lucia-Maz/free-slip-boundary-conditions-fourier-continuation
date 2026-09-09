#!/bin/bash
#SBATCH --job-name=puerta-fase2
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:40:00
#SBATCH --output=puerta-fase2-%j.out
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
eval "$(bash verificacion/toolchain_cluster.sh)"

# La E/S no va a $HOME: desde los nodos de cómputo la E/S paralela contra $HOME es
# lenta. Las reglas dicen /share/data* y /share/scratch*, CON ASTERISCO, así que los
# directorios reales tienen algún sufijo y hay que descubrirlos en vez de adivinarlos.
elegir_base() {
    local d
    for base in ${SPECTER_SCRATCH_BASE:-} /share/scratch* /share/data*; do
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

export SPECTER_SCRATCH="$BASE/puerta-fase2-${SLURM_JOB_ID:-manual}"
mkdir -p "$SPECTER_SCRATCH"
trap 'rm -rf "$SPECTER_SCRATCH"' EXIT   # SPECTER_SCRATCH explícito no se autolimpia

python verificacion/test_aceptacion_fase2.py
