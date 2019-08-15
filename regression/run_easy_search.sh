#!/bin/sh -e
MMSEQS="${1}"
EVALUATE="${2}"
DATADIR="${3}"
RESULTS="${4}"
mkdir -p "${RESULTS}"

QUERY="${DATADIR}/query.fasta"
TARGET="${DATADIR}/targetannotation.fasta"

"${MMSEQS}" easy-search "$QUERY" "$TARGET" "$RESULTS/results_aln.m8" "$RESULTS/tmp" -e 10000 -s 4 --max-seqs 4000 --compressed 1

"${EVALUATE}" "$QUERY" "$TARGET" "$RESULTS/results_aln.m8" "${RESULTS}/evaluation_roc5.dat" 4000 1 | tee "${RESULTS}/evaluation.log"
ACTUAL=$(grep "^ROC5 AUC:" "${RESULTS}/evaluation.log" | cut -d" " -f3)
TARGET="0.238"
awk -v actual="$ACTUAL" -v target="$TARGET" \
    'BEGIN { print (actual >= target) ? "GOOD" : "BAD"; print "Expected: ", target; print "Actual: ", actual; }' \
    > "${RESULTS}/report"
