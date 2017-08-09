#!/bin/bash -e
BENCHDIR="${1}"
MMSEQS="${2}"
EVALUATE="${3}"
VERSION="${4}"
RESULTDIR="${5}"
NAME="${6:-regression-test}"
THREADS="${7:-4}"

BENCHDB="small-benchmark-db"

QUERY="$BENCHDIR/${BENCHDB}/query.fasta"
QUERYDB="$BENCHDIR/${BENCHDB}/query"
TARGETDB="$BENCHDIR/${BENCHDB}/db2"
DBANNOTATION="$BENCHDIR/${BENCHDB}/targetannotation.fasta"

RESULTS="${RESULTDIR}/mmseqs-${NAME}-${VERSION}"

TIMERS=()
function lap() { TIMERS+=($(date +%s.%N)); }
SEARCH_PARM="--min-ungapped-score 15 -e 10000.0 -s 1 --max-seqs 4000 --split 1 --target-profile"

rm -rf "$RESULTS"
mkdir -p "$RESULTS/tmp"
lap
${MMSEQS} mergedbs "$TARGETDB" "${TARGETDB}_fasta" "${TARGETDB}_h" "$TARGETDB" --prefixes ">" 1>&2 || exit 125
lap
${MMSEQS} msa2profile "${TARGETDB}_fasta" "${TARGETDB}_profile" --filter-msa 0 --threads $THREADS 1>&2 || exit 125
lap
${MMSEQS} search "$QUERYDB" "${TARGETDB}_profile" "$RESULTS/results_aln" "$RESULTS/tmp" --threads $THREADS ${SEARCH_PARM} 1>&2 || exit 125 
lap
${MMSEQS} convertalis "$QUERYDB" "$TARGETDB" "$RESULTS/results_aln" "$RESULTS/results_aln.m8" ${VERBOSE} 1>&2  || exit 125
lap

EVALPREFIX="${RESULTS}/evaluation"
${EVALUATE} "$QUERY" "$DBANNOTATION" "$RESULTS/results_aln.m8" "${EVALPREFIX}_roc5.dat" 4000 1 > "${EVALPREFIX}.log"
AUC=$(grep "^ROC5 AUC:" "${EVALPREFIX}.log" | cut -d" " -f3)
cat "${EVALPREFIX}.log" >&2
echo -e "${NAME}\t${VERSION}\t9\t${AUC}\t$(printf '%s\t' "${TIMERS[@]}")"
