#!/bin/bash -e
BENCHDIR="${1}"
MMSEQS="${2}"
EVALUATE="${3}"
VERSION="${4}"
RESULTDIR="${5}"
NAME="${6:-regression-nucl}"
THREADS="${7:-4}"

BENCHDB="small-benchmark-db"

QUERY="$BENCHDIR/${BENCHDB}/query.fasta"
QUERYDB="$BENCHDIR/${BENCHDB}/nucl"
TARGETDB="$BENCHDIR/${BENCHDB}/db2"
DBANNOTATION="$BENCHDIR/${BENCHDB}/targetannotation.fasta"

RESULTS="${RESULTDIR}/mmseqs-${NAME}-${VERSION}"

TIMERS=()
function lap() { TIMERS+=($(date +%s.%N)); }
SEARCH_PARM="--min-ungapped-score 15 -e 10000.0 -s 4 --max-seqs 4000"

rm -rf "$RESULTS"
mkdir -p "$RESULTS/tmp"
mkdir -p "$RESULTS/tmp_index"
lap
${MMSEQS} search "$QUERYDB" "$TARGETDB" "$RESULTS/results_aln" "$RESULTS/tmp" --threads $THREADS ${SEARCH_PARM} 1>&2 || exit 125 
lap
${MMSEQS} convertalis "$QUERYDB" "$TARGETDB" "$RESULTS/results_aln" "$RESULTS/results_aln.m8" ${VERBOSE} 1>&2  || exit 125
lap
rm -rf ${TARGETDB}.sk*

EVALPREFIX="${RESULTS}/evaluation"
${EVALUATE} "$QUERY" "$DBANNOTATION" "$RESULTS/results_aln.m8" "${EVALPREFIX}_roc5.dat" 4000 1 > "${EVALPREFIX}.log"
AUC=$(grep "^ROC5 AUC:" "${EVALPREFIX}.log" | cut -d" " -f3)
cat "${EVALPREFIX}.log" >&2
echo -e "${NAME}\t${VERSION}\t11\t${AUC}\t$(printf '%s\t' "${TIMERS[@]}")"

