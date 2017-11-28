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
TARGETDB="$BENCHDIR/${BENCHDB}/db2"
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
mkdir -p "$RESULTS/tmp_index"
lap
${MMSEQS} createindex "$TARGETDB" "$RESULTS/tmp_index"  --threads $THREADS --split 1 1>&2 || exit 125
lap
${MMSEQS} search "$QUERYDB" "$TARGETDB" "$RESULTS/results_aln" "$RESULTS/tmp" --threads $THREADS ${SEARCH_PARM} 1>&2 || exit 125 
lap
${MMSEQS} convertalis "$QUERYDB" "$TARGETDB" "$RESULTS/results_aln" "$RESULTS/results_aln.m8" ${VERBOSE} 1>&2  || exit 125
lap
rm -rf ${TARGETDB}.sk*

EVALPREFIX="${RESULTS}/evaluation"
if [ $PROFILE -ne 0 ]; then
    LC_ALL=C sort -k1,1 -k11,11g "$RESULTS/results_aln.m8" > "$RESULTS/results_aln_sorted.m8"
    mv "$RESULTS/results_aln_sorted.m8" "$RESULTS/results_aln.m8"

	${MMSEQS} createtsv "$QUERYDB" "$TARGETDB" "$RESULTS/tmp/latest/pref_1" "$RESULTS/results_pref.tsv" ${VERBOSE} 1>&2  || exit 125
	awk '{print $1"\t"$2"\t"0"\t"0"\t"0"\t"0"\t"0"\t"0"\t"0"\t"0"\t"$3"\t"0}' "$RESULTS/results_pref.tsv" > "$RESULTS/results_pref.m8"
    LC_ALL=C sort -k1,1 -k11,11g "$RESULTS/results_pref.m8" > "$RESULTS/results_pref_sorted.m8"
    mv "$RESULTS/results_pref_sorted.m8" "$RESULTS/results_pref.m8"

	${EVALUATE} "$QUERY" "$DBANNOTATION" "$RESULTS/results_pref.m8" "${EVALPREFIX}_pref_roc5.dat" 4000 1 > "${EVALPREFIX}_pref.log"
	AUC=$(grep "^ROC5 AUC:" "${EVALPREFIX}_pref.log" | cut -d" " -f3)
	echo -e "${NAME}_pref\t${VERSION}\t2\t${AUC}"
else
	${MMSEQS} createtsv "$QUERYDB" "$TARGETDB" "$RESULTS/tmp/latest/pref_4.0" "$RESULTS/results_pref.tsv" ${VERBOSE} 1>&2  || exit 125
	awk '{print $1"\t"$2"\t"0"\t"0"\t"0"\t"0"\t"0"\t"0"\t"0"\t"0"\t"$3"\t"0}' "$RESULTS/results_pref.tsv" > "$RESULTS/results_pref.m8"

	${EVALUATE} "$QUERY" "$DBANNOTATION" "$RESULTS/results_pref.m8" "${EVALPREFIX}_pref_roc5.dat" 4000 1 > "${EVALPREFIX}_pref.log"
	AUC=$(grep "^ROC5 AUC:" "${EVALPREFIX}_pref.log" | cut -d" " -f3)
	echo -e "${NAME}_pref\t${VERSION}\t2\t${AUC}"
fi

${EVALUATE} "$QUERY" "$DBANNOTATION" "$RESULTS/results_aln.m8" "${EVALPREFIX}_roc5.dat" 4000 1 > "${EVALPREFIX}.log"
AUC=$(grep "^ROC5 AUC:" "${EVALPREFIX}.log" | cut -d" " -f3)
cat "${EVALPREFIX}.log" >&2
echo -e "${NAME}\t${VERSION}\t${PROFILE}\t${AUC}\t$(printf '%s\t' "${TIMERS[@]}")"
