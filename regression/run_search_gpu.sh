#!/bin/sh -e

head -n 1000 "${DATADIR}/query.fasta" > "${RESULTS}/query_500.fasta"

QUERY="${RESULTS}/query_500.fasta"
QUERYDB="${RESULTS}/query"
"${MMSEQS}" createdb "${QUERY}" "${QUERYDB}"

TARGET="${DATADIR}/targetannotation.fasta"
TARGETDB="${RESULTS}/targetannotation"
"${MMSEQS}" createdb "${TARGET}" "${TARGETDB}_db"
"${MMSEQS}" makepaddedseqdb "${TARGETDB}_db" "${TARGETDB}"

"${MMSEQS}" search "$QUERYDB" "$TARGETDB" "$RESULTS/results_aln" "$RESULTS/tmp" -e 10000 --max-seqs 1000 --gpu 1
"${MMSEQS}" convertalis "$QUERYDB" "$TARGETDB" "$RESULTS/results_aln" "$RESULTS/results_aln.m8"


"${EVALUATE}" "$QUERY" "$TARGET" "$RESULTS/results_aln.m8" "${RESULTS}/evaluation_roc5.dat" 1000 1 | tee "${RESULTS}/evaluation.log"
ACTUAL=$(grep "^ROC5 AUC:" "${RESULTS}/evaluation.log" | cut -d" " -f3)
TARGET="0.454819"
awk -v actual="$ACTUAL" -v target="$TARGET" \
    'BEGIN { print (actual == target) ? "GOOD" : "BAD"; print "Expected: ", target; print "Actual: ", actual; }' \
    > "${RESULTS}.report"
