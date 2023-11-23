#!/usr/bin/env bash
set -euo pipefail

LINTER="tools/cpplint.py --filter=-runtime/references --quiet"

function cpplint () {
  LOGFILE=$(mktemp)
  find ./ -type f -name '*.h' -exec $LINTER {} >> ${LOGFILE} 2>&1 \;
  find ./ -type f -name '*.cc' -exec $LINTER {}>> ${LOGFILE} 2>&1 \; 
  find ./ -type f -name '*.cpp' -exec $LINTER {} >> ${LOGFILE} 2>&1 \; 

  # if log file not empty
  if [[ -s $LOGFILE ]]; then
    cat $LOGFILE
    rm $LOGFILE
    return 1
  fi

  rm $LOGFILE
  return 0
}

echo "${BASH_SOURCE} ${0}"
if [[ "$0" == "$BASH_SOURCE" ]] ; then
  cpplint
  exit $?
fi

