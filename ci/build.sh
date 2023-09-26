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

source /etc/profile.d/z10_spack_environment.sh

rm -f "${HERE}/orca-jedi"
ln -s '..' "${HERE}/orca-jedi"
ecbuild -S "${HERE}"
make -j "${NPROC}"

echo ls
ls

echo ls share
ls share

echo ls atlas-orca
ls atlas-orca

echo ls atlas
ls atlas

echo ls atlas/share/plugins
ls atlas/share/plugins

echo ls atlas-orca/share/plugins
ls atlas-orca/share/plugins

echo ls share/plugins
ls share/plugins

exit 1

env OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    ATLAS_TRACE=1 ATLAS_DEBUG=1 \
    LD_LIBRARY_PATH="${HERE}/lib:${LD_LIBRARY_PATH}" \
    ATLAS_DATA_PATH="${HERE}/atlas-data" \
    ctest -j "${NPROC}" -V --output-on-failure --test-dir "./orca-jedi"

exit
