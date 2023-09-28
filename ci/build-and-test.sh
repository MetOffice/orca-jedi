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

source /opt/spack-environment/activate.sh

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

if [[ ! -f share/plugins/atlas-orca.yml ]]; then
 echo "ERROR atlas-orca.yml not found!" | tee >(cat >&2)
 exit 1
fi

env OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    ATLAS_TRACE=1 ATLAS_DEBUG=1 \
    LD_LIBRARY_PATH="${HERE}/lib:${LD_LIBRARY_PATH}" \
    ATLAS_DATA_PATH="${HERE}/atlas-data" \
    ctest -j "${NPROC}" -V --output-on-failure --test-dir "./orca-jedi"

exit
