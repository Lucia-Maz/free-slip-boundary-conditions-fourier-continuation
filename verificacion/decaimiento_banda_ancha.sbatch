#!/bin/bash
#SBATCH --job-name=banda-ancha
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=verificacion/logs/banda-ancha-%j.out
#
# Convergencia 64²/128² de la prueba puente [V-19]. No se ejecuta en el login:
#
#   PRUEBA=1 RESOLUCION=64 sbatch verificacion/decaimiento_banda_ancha.sbatch
#   RESOLUCION=64 sbatch verificacion/decaimiento_banda_ancha.sbatch
#
# El segundo comando corre primero el control lineal y después la amplitud experimental.
# Cada caso se guarda inmediatamente en salidas/tablas/decaimiento_banda_ancha.json, de
# modo que una segunda invocación reanuda lo que falte. Para una eventual 128² hay que
# estimar y pedir el walltime a partir de 64², no lanzarla a ciegas.

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

module load gnu15 openmpi5 fftw/3.3.11 python/3.13.13

export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=sm,self,tcp

TOOLCHAIN_ENV="$(bash verificacion/toolchain_cluster.sh)"
eval "$TOOLCHAIN_ENV"

elegir_base() {
    local d base
    for base in ${SPECTER_SCRATCH_BASE:-} /share/data2 /share/scratch1; do
        [ -d "$base" ] || continue
        d="$base/$USER"
        if mkdir -p "$d" 2>/dev/null && [ -w "$d" ]; then
            echo "$d"
            return 0
        fi
    done
    return 1
}

if ! BASE="$(elegir_base)"; then
    echo "ERROR: no hay scratch escribible en /share; no corro producción contra HOME." >&2
    exit 1
fi

export SPECTER_SCRATCH="$BASE/banda-ancha-${SLURM_JOB_ID:?}"
mkdir -p "$SPECTER_SCRATCH"
limpiar() {
    case "${SPECTER_SCRATCH:-}" in
        */banda-ancha-?*)
            [ -d "$SPECTER_SCRATCH" ] && rm -rf "$SPECTER_SCRATCH" ;;
        *)
            echo "[aviso] no borro '${SPECTER_SCRATCH:-}': forma inesperada" >&2 ;;
    esac
}
trap limpiar EXIT

RESOLUCION="${RESOLUCION:-64}"
SOLO="${SOLO:-ambas}"
PRUEBA="${PRUEBA:-0}"
case "$RESOLUCION" in 64|128) ;; *) echo "RESOLUCION debe ser 64 o 128" >&2; exit 2;; esac
case "$SOLO" in lineal|experimental|ambas) ;; *) echo "SOLO inválido" >&2; exit 2;; esac
case "$PRUEBA" in 0|1) ;; *) echo "PRUEBA debe ser 0 o 1" >&2; exit 2;; esac

ARGS=(--resolucion "$RESOLUCION" --solo "$SOLO")
if [ "$PRUEBA" = 1 ]; then ARGS+=(--prueba); fi

# El timeout interno queda por encima del walltime de omisión. SLURM sigue siendo quien
# limita el trabajo; esta variable evita que Python corte una 128² válida a las cuatro horas.
export SPECTER_RUN_TIMEOUT_S=86400
python verificacion/decaimiento_banda_ancha.py "${ARGS[@]}"
