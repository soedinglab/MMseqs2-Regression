#!/bin/bash -e
BENCHDIR="${1}"
MMSEQS="${2}"
EVALUATE="${3}"
VERSION="${4}"
RESULTDIR="${5}"
PROFILE="${6:-0}"
NAME="${7:-regression-test}"
THREADS="${8:-4}"

BENCHDB="small-benchmark-db"

QUERY="$BENCHDIR/${BENCHDB}/query.fasta"
QUERYDB="$BENCHDIR/${BENCHDB}/query"
TARGETDB="$BENCHDIR/${BENCHDB}/db2.fasta"
DBANNOTATION="$BENCHDIR/${BENCHDB}/targetannotation.fasta"

RESULTS="${RESULTDIR}/mmseqs-${NAME}-${VERSION}"

TIMERS=()
function lap() { TIMERS+=($(date +%s.%N)); }
SEARCH_PARM="--min-ungapped-score 15 -e 10000.0 -s 4 --max-seqs 4000 --split 1"
if [ $PROFILE -ne 0 ]; then
    SEARCH_PARM="--min-ungapped-score 0 -e 10000.0 -s 4 --max-seqs 4000 --num-iterations 2 --split 1 "
fi

rm -rf "$RESULTS"
mkdir -p "$RESULTS/tmp"
lap
${MMSEQS} easy-search "$QUERY" "$TARGETDB" "$RESULTS/results_aln.m8" "$RESULTS/tmp" --threads $THREADS ${SEARCH_PARM} 1>&2 || exit 125 
lap

EVALPREFIX="${RESULTS}/evaluation"
if [ $PROFILE -ne 0 ]; then
    LC_ALL=C sort -k1,1 -k11,11g "$RESULTS/results_aln.m8" > "$RESULTS/results_aln_sorted.m8"
    mv "$RESULTS/results_aln_sorted.m8" "$RESULTS/results_aln.m8"
fi

${EVALUATE} "$QUERY" "$DBANNOTATION" "$RESULTS/results_aln.m8" "${EVALPREFIX}_roc5.dat" 4000 1 > "${EVALPREFIX}.log"
AUC=$(grep "^ROC5 AUC:" "${EVALPREFIX}.log" | cut -d" " -f3)
cat "${EVALPREFIX}.log" >&2
echo -e "${NAME}\t${VERSION}\t${PROFILE}\t${AUC}\t$(printf '%s\t' "${TIMERS[@]}")"
