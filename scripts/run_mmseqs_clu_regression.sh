#!/bin/bash -e
BENCHDIR="${1}"
MMSEQS="${2}"
EVALUATE="${3}"
VERSION="${4}"
RESULTDIR="${5}"
PROFILE="${6:-0}"
NAME="${7:-regression-test}"

BENCHDB="small-benchmark-db"

QUERY="$BENCHDIR/${BENCHDB}/clu.fasta"
QUERYDB="$BENCHDIR/${BENCHDB}/clu"
RESULTS="${RESULTDIR}/mmseqs-${NAME}-${VERSION}"

TIMERS=()
function lap() { TIMERS+=($(date +%s.%N)); }
CLU_PARM="--min-seq-id 0.3 -s 4"

rm -rf "$RESULTS"
mkdir -p "$RESULTS/tmp"
lap
${MMSEQS} cluster "$QUERYDB" "$RESULTS/results_clu" "$RESULTS/tmp" ${CLU_PARM} 1>&2 || exit 125 
lap
${MMSEQS} createtsv "$QUERYDB" "$QUERYDB" "$RESULTS/results_clu" "$RESULTS/results_clu.tsv" ${VERBOSE} 1>&2  || exit 125
lap

${EVALUATE} -d "$QUERY" -c "$RESULTS/results_clu.tsv" > "${EVALPREFIX}.log"

CLUSIZE=$(grep "^# clusters " "${EVALPREFIX}.log" | tail -n 1)
CLUWRONG=$(grep "^#wrong sequences: " "${EVALPREFIX}.log" | awk '{print $3}')
CLUCORRUPTED=$(grep "^#corrupted clusters: " "${EVALPREFIX}.log" | awk '{print $3}')
echo -e "${NAME}\t${VERSION}\t${CLUSIZE}\t${CLUWRONG}\t${CLUCORRUPTED}\t$(printf '%s\t' "${TIMERS[@]}")"
