#!/bin/bash -e
BENCHDIR="${1}"
MMSEQS="${2}"
EVALUATE="${3}"
VERSION="${4}"
RESULTDIR="${5}"
NAME="${6:-cs-regression-test}"
THREADS="${7:-4}"

BENCHDB="small-benchmark-db"

QUERY="$BENCHDIR/${BENCHDB}/query.fasta"
QUERYDB="$BENCHDIR/${BENCHDB}/query"
TARGETDB="$BENCHDIR/${BENCHDB}/db2"
DBANNOTATION="$BENCHDIR/${BENCHDB}/targetannotation.fasta"

RESULTS="${RESULTDIR}/mmseqs-${NAME}-${VERSION}"

TIMERS=()
function lap() { TIMERS+=($(date +%s.%N)); }

rm -rf "$RESULTS"
mkdir -p "$RESULTS/tmp"
lap
"${MMSEQS}" search "${TARGETDB}" "${TARGETDB}" "${RESULTS}/aln_target_profile" "${RESULTS}/tmp" -s 2 --e-profile 0.1 -e 0.1 -a --realign 1>&2 || exit 125
lap
"${MMSEQS}" result2profile "${TARGETDB}" "${TARGETDB}" "${RESULTS}/aln_target_profile" "${RESULTS}/target_profile" 1>&2 || exit 125
lap
"${MMSEQS}" profile2cs "${RESULTS}/target_profile" "${RESULTS}/target_profile_states" 1>&2 || exit 125
lap
"${MMSEQS}" search "${QUERYDB}" "${TARGETDB}" "${RESULTS}/aln_query_target" "${RESULTS}/tmp" -s 2 --e-profile 0.1 -e 0.1 --realign -a 1>&2 || exit 125
lap
"${MMSEQS}" result2profile "${QUERYDB}" "${TARGETDB}" "${RESULTS}/aln_query_target" "${RESULTS}/query_profile" 1>&2 || exit 125
lap
"${MMSEQS}" search "${RESULTS}/query_profile" "${RESULTS}/target_profile_states" "${RESULTS}/results_aln" "${RESULTS}/tmp" --max-seqs 4000 -e 100000.0 -s 4 -k 10  --min-ungapped-score 0 1>&2 || exit 125
lap
"${MMSEQS}" convertalis "${QUERYDB}" "${TARGETDB}" "${RESULTS}/results_aln" "${RESULTS}/results_aln.m8" 1>&2 || exit 125
lap

EVALPREFIX="${RESULTS}/evaluation"
${EVALUATE} "$QUERY" "$DBANNOTATION" "$RESULTS/results_aln.m8" "${EVALPREFIX}_roc5.dat" 4000 1 > "${EVALPREFIX}.log"
AUC=$(grep "^ROC5 AUC:" "${EVALPREFIX}.log" | cut -d" " -f3)
cat "${EVALPREFIX}.log" >&2
echo -e "${NAME}\t${VERSION}\t10\t${AUC}\t$(printf '%s\t' "${TIMERS[@]}")"
