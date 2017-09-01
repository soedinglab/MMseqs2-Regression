#!/bin/bash -e
QUERY="${1}"
MMSEQS="${2}"
NAME="${3}"
VERSION="${4}"
RESULTDIR="${5}"
LINCLUST="${6}"
PARAMS="${7}"

RESULTS="${RESULTDIR}/mmseqs-${NAME}-${VERSION}"
QUERYDB="$RESULTS/clu"

TIMERS=()
function lap() { TIMERS+=($(date +%s.%N)); }

MODE="cluster"
BENCHID="3"
if [[ "$LINCLUST" == 1 ]]; then
    MODE="linclust"
    BENCHID="6"
fi

rm -rf "$RESULTS"
mkdir -p "$RESULTS/tmp"
lap
${MMSEQS} createdb $QUERY "$QUERYDB" 1>&2 || exit 125
lap
${MMSEQS} $MODE "$QUERYDB" "$RESULTS/results_clu" "$RESULTS/tmp" ${PARAMS} 1>&2 || exit 125
lap
${MMSEQS} createtsv "$QUERYDB" "$QUERYDB" "$RESULTS/results_clu" "$RESULTS/results_clu.tsv" 1>&2  || exit 125
lap

RESULT="$(awk 'BEGIN { l = "" }  l != $1 { l = $1; cnt++; } { t++; } END { print cnt"\t"t"\t"(t/cnt) }' "$RESULTS/results_clu.tsv")"
IFS=$'|' VALUES=(${RESULT//$'\t'/|})
echo -e "${NAME}\t${VERSION}\t$(($BENCHID + 0))\t${VALUES[0]}\t$(printf '%s\t' "${TIMERS[@]}")"
echo -e "${NAME}\t${VERSION}\t$(($BENCHID + 1))\t${VALUES[1]}\t$(printf '%s\t' "${TIMERS[@]}")"
echo -e "${NAME}\t${VERSION}\t$(($BENCHID + 2))\t${VALUES[2]}\t$(printf '%s\t' "${TIMERS[@]}")"

