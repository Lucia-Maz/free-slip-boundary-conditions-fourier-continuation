#!/usr/bin/env bash
# Arma un prefijo sintético para SPECTER_TOOLCHAIN en el clúster.
#
# POR QUÉ HACE FALTA
# ------------------
# test_aceptacion_fase2.py está fijada por hash y no se edita ([D-26], CLAUDE.md). Esa
# puerta espera UN solo prefijo que tenga adentro bin/mpif90, bin/mpirun y lib/libfftw3,
# que es como queda un env de conda en la laptop. Con módulos, en cambio, OpenMPI y FFTW
# viven en prefijos distintos. Este script arma un directorio con enlaces simbólicos a
# los dos, que cumple esa forma, sin tocar la puerta ni los módulos.
#
# USO
#   module load gnu15 openmpi5 fftw/3.3.11
#   eval "$(bash verificacion/toolchain_cluster.sh)"
#   python verificacion/test_aceptacion_fase2.py     # tiene que dar 6/6
#
# El script NO modifica nada fuera de su directorio de salida, y no borra nada que no
# haya creado él.

set -euo pipefail

DESTINO="${1:-$HOME/.local/share/specter-toolchain}"

err() { echo "ERROR: $*" >&2; exit 1; }

# --- MPI: tiene que estar en el PATH, o sea con el módulo cargado -------------------
MPIF90="$(command -v mpif90 || true)"
MPIRUN="$(command -v mpirun || true)"
[ -n "$MPIF90" ] || err "no encuentro mpif90 en el PATH. ¿Cargaste 'module load openmpi5'?"
# Tiene que ser mpirun y no srun: la puerta lo invoca con '-np', que srun no acepta.
# Antes que enlazar srun con el nombre mpirun y fallar de manera confusa, falla acá.
[ -n "$MPIRUN" ] || err "no encuentro mpirun en el PATH. ¿Cargaste 'module load openmpi5'?"

# --- FFTW: se busca su prefijo, probando lo que suelen definir los módulos ----------
FFTW_PREFIJO=""
for v in "${FFTW_DIR:-}" "${FFTW_ROOT:-}" "${FFTW_HOME:-}" "${FFTWDIR:-}" \
         "${FFTW3_DIR:-}" "${FFTW3_ROOT:-}"; do
    if [ -n "$v" ] && [ -d "$v/lib" ] && ls "$v"/lib/libfftw3* >/dev/null 2>&1; then
        FFTW_PREFIJO="$v"; break
    fi
done

# si ninguna variable sirvió, se busca libfftw3 en el LD_LIBRARY_PATH que dejó el módulo
if [ -z "$FFTW_PREFIJO" ]; then
    IFS=':' read -ra DIRS <<< "${LD_LIBRARY_PATH:-}"
    for d in "${DIRS[@]}"; do
        if [ -n "$d" ] && ls "$d"/libfftw3* >/dev/null 2>&1; then
            FFTW_PREFIJO="$(dirname "$d")"; break
        fi
    done
fi

[ -n "$FFTW_PREFIJO" ] || err "no encuentro el prefijo de FFTW. ¿Cargaste 'module load fftw/3.3.11'?
       Probá 'module show fftw/3.3.11' y pasá el prefijo a mano:
       $0 <destino>  con FFTW_DIR=<prefijo> exportado."

[ -d "$FFTW_PREFIJO/include" ] || err "el prefijo de FFTW ($FFTW_PREFIJO) no tiene include/"

# --- armar el prefijo sintético ------------------------------------------------------
rm -rf "$DESTINO"
mkdir -p "$DESTINO/bin"
ln -sf "$MPIF90" "$DESTINO/bin/mpif90"
ln -sf "$MPIRUN" "$DESTINO/bin/mpirun"
ln -sfn "$FFTW_PREFIJO/lib" "$DESTINO/lib"
ln -sfn "$FFTW_PREFIJO/include" "$DESTINO/include"

{
  echo "# prefijo sintético para SPECTER_TOOLCHAIN"
  echo "#   mpif90 -> $MPIF90"
  echo "#   mpirun -> $MPIRUN"
  echo "#   fftw   -> $FFTW_PREFIJO"
} >&2

echo "export SPECTER_TOOLCHAIN=$DESTINO"
