#!/bin/bash
#SBATCH --job-name=fs6-superficie
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --output=verificacion/logs/fs6-superficie-%j.out
#
# Campos de superficie con la superficie deformable de orden 2 a la escala de la celda
# (numerico/fase6/campos_superficie.py). Instrucciones completas en
# numerico/fase6/HANDOFF_SAKURA_CAMPOS_SUPERFICIE.md. Va a `normal` porque compila adentro y c2–c5 no tienen make ([H-19]). No se ejecuta en el login:
#
#   PRUEBA=1 sbatch --time=00:30:00 numerico/fase6/campos_superficie.sbatch   # cadena, 16²
#   NXY=64 MEMORIAS=4 CAMPOS=12 sbatch numerico/fase6/campos_superficie.sbatch
#
# Costo medido en la laptop: 64², 1,5 memorias, 3 procesos, ~1,2 s por paso (3800 pasos).
# Las transformadas de superficie del orden N son DFT explícitas, O(nxy³), y se repiten en
# cada proceso: más procesos no las aceleran, y 128² cuesta ~8 veces más por paso. No
# pedir --mem (ver decaimiento_banda_ancha.sbatch).

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
    echo "ERROR: no hay scratch escribible en /share; no corro contra HOME." >&2
    exit 1
fi
export SPECTER_SCRATCH="$BASE/fs6-superficie-${SLURM_JOB_ID:?}"
mkdir -p "$SPECTER_SCRATCH"
limpiar() {
    case "${SPECTER_SCRATCH:-}" in
        */fs6-superficie-?*) [ -d "$SPECTER_SCRATCH" ] && rm -rf "$SPECTER_SCRATCH" ;;
        *) echo "[aviso] no borro '${SPECTER_SCRATCH:-}': forma inesperada" >&2 ;;
    esac
}
trap limpiar EXIT

if [ "${PRUEBA:-0}" = 1 ]; then
    python numerico/fase6/campos_superficie.py --nxy 16 --prueba --np 1
    exit 0
fi
NXY="${NXY:-64}"
MEMORIAS="${MEMORIAS:-4}"
CAMPOS="${CAMPOS:-12}"
case "$NXY" in 32|64) ;; *) echo "NXY debe ser 32 o 64 (128 no sin estimar el costo)" >&2; exit 2;; esac
export SPECTER_RUN_TIMEOUT_S=21000
python numerico/fase6/campos_superficie.py --nxy "$NXY" --memorias "$MEMORIAS" \
    --campos "$CAMPOS" --np "${SLURM_NTASKS:-1}"
