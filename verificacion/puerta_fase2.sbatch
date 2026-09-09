#!/bin/bash
#SBATCH --job-name=puerta-fase2
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:40:00
#SBATCH --output=puerta-fase2-%j.out
#
# La puerta de aceptación de la Fase 2 en Sakura. Se manda al planificador y no se
# corre a mano en el nodo de login: compila SPECTER y lanza unas quince simulaciones
# chicas, o sea diez a quince minutos de CPU. Las reglas del clúster permiten trabajos
# CPU cortos en la partición 'normal' (nodo de login); acá se usa 'compute' para no
# cargar el login, y cambiarlo es de una línea.
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

# La E/S va a /share/scratch, no a $HOME: desde los nodos de cómputo la E/S paralela
# contra $HOME es lenta.
export SPECTER_SCRATCH="/share/scratch/$USER/puerta-fase2-${SLURM_JOB_ID:-manual}"
mkdir -p "$SPECTER_SCRATCH"
trap 'rm -rf "$SPECTER_SCRATCH"' EXIT   # SPECTER_SCRATCH explícito no se autolimpia

python verificacion/test_aceptacion_fase2.py
