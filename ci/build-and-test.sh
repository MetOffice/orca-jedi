#!/bin/bash
set -euo pipefail

finally() {
    trap '' ERR
    trap '' EXIT
    if [[ -d "${WORKD:-}" ]]; then
        cd /
        rm -fr "${WORKD}"
    fi
}

HERE="$(cd "$(dirname "$0")" && pwd)"
THIS="$(basename "$0")"
NPROC=${NPROC:-2}
WORKD="$(mktemp -d "${THIS}-XXXXXX" -t)"

trap finally ERR
trap finally EXIT

cd "${WORKD}"

if [[ -f /opt/spack-environment/activate.sh ]]; then
    source /opt/spack-environment/activate.sh
fi

# -- Enable OpenMPI over subscription -----------------------------------------
if command -v ompi_info &>/dev/null; then
    echo "Check support for MPI_THREAD_MULTIPLE"
    ompi_info | grep -i 'thread support'
    ompi_vn=$(ompi_info | awk '/Ident string:/ {print $3}')
    case $ompi_vn in
        4.*) export OMPI_MCA_rmaps_base_oversubscribe=1 ;;
        5.*) export PRTE_MCA_rmaps_default_mapping_policy=:oversubscribe ;;
    esac
fi

echo "
-------------------------------
gcc version $(gcc -dumpversion ||:)
$(ecbuild --version ||:)
atlas version $(atlas --version ||:)
eckit version $(eckit-version ||:)
ectrans version $(ectrans --version ||:)
fckit version $(fckit --version ||:)
fiat version $(fiat --version ||:)
lz4 version $(lz4 --version ||:)
odc version $(odc --version ||:)
-------------------------------
"

rm -f "${HERE}/orca-jedi"
ln -s '..' "${HERE}/orca-jedi"
ecbuild -S "${HERE}"
make -j "${NPROC}"

env OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    ATLAS_TRACE=1 ATLAS_DEBUG=1 \
    LD_LIBRARY_PATH="${HERE}/lib:${LD_LIBRARY_PATH}" \
    ATLAS_DATA_PATH="${HERE}/atlas-data" \
    ctest -j "${NPROC}" -V --output-on-failure --test-dir "./orca-jedi"

exit
