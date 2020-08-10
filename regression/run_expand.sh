#!/bin/sh -e
QUERY="${DATADIR}/query.fasta"
QUERYDB="${RESULTS}/query"

TARGET="${DATADIR}/targetannotation.fasta"
TARGETDB="${RESULTS}/targetannotation"
"${MMSEQS}" createdb "${QUERY}" "${QUERYDB}"
"${MMSEQS}" createdb "${TARGET}" "${TARGETDB}"
"${MMSEQS}" search "$TARGETDB" "$TARGETDB" "$RESULTS/results_bc" "$RESULTS/tmp" -s 1 -a

"${MMSEQS}" search "$QUERYDB" "$TARGETDB" "$RESULTS/results_ab" "$RESULTS/tmp" -e 10000 -s 2 -a --max-seqs 4000
"${MMSEQS}" expandaln "$QUERYDB" "$TARGETDB" "$RESULTS/results_ab" "$RESULTS/results_bc" "$RESULTS/results_ac" -e 10000
"${MMSEQS}" convertalis "$QUERYDB" "$TARGETDB" "$RESULTS/results_ac" "$RESULTS/results_aln.m8"
"${EVALUATE}" "$QUERY" "$TARGET" "$RESULTS/results_aln.m8" "${RESULTS}/evaluation_roc5.dat" 4000 1 | tee "${RESULTS}/evaluation.log"
ACTUAL1=$(grep "^ROC5 AUC:" "${RESULTS}/evaluation.log" | cut -d" " -f3)
TARGET1="0.16614"

"${MMSEQS}" expand2profile "$QUERYDB" "$TARGETDB" "$RESULTS/results_ab" "$RESULTS/results_bc" "$RESULTS/prof_ac"
"${MMSEQS}" search "$RESULTS/prof_ac" "$TARGETDB" "$RESULTS/results_prof" "$RESULTS/tmp" -e 10000 -s 2 -a --max-seqs 4000
"${MMSEQS}" convertalis "$QUERYDB" "$TARGETDB" "$RESULTS/results_prof" "$RESULTS/results_prof.m8"
"${EVALUATE}" "$QUERY" "$TARGET" "$RESULTS/results_prof.m8" "${RESULTS}/evaluation_roc5.dat" 4000 1 | tee "${RESULTS}/evaluation.log"
ACTUAL2=$(grep "^ROC5 AUC:" "${RESULTS}/evaluation.log" | cut -d" " -f3)
TARGET2="0.171723"

"${MMSEQS}" result2profile "$TARGETDB" "$TARGETDB" "$RESULTS/results_bc" "$RESULTS/prof_bc"
"${MMSEQS}" result2profile "$QUERYDB" "$TARGETDB" "$RESULTS/results_ab" "$RESULTS/prof_ab"
"${MMSEQS}" result2pp "$RESULTS/prof_ab" "$RESULTS/prof_bc" "$RESULTS/results_ab" "$RESULTS/prof_merged"
"${MMSEQS}" search "$RESULTS/prof_merged" "$TARGETDB" "$RESULTS/results_merged" "$RESULTS/tmp" -e 10000 -s 2 -a --max-seqs 4000
"${MMSEQS}" convertalis "$QUERYDB" "$TARGETDB" "$RESULTS/results_merged" "$RESULTS/results_merged.m8"
"${EVALUATE}" "$QUERY" "$TARGET" "$RESULTS/results_merged.m8" "${RESULTS}/evaluation_roc5.dat" 4000 1 | tee "${RESULTS}/evaluation.log"
ACTUAL3=$(grep "^ROC5 AUC:" "${RESULTS}/evaluation.log" | cut -d" " -f3)
TARGET3="0.222982"

awk -v actual1="$ACTUAL1" -v target1="$TARGET1" -v actual2="$ACTUAL2" -v target2="$TARGET2" -v actual3="$ACTUAL3" -v target3="$TARGET3" \
    'BEGIN { print (actual1 >= target1 && actual2 >= target2 && actual3 >= target3) ? "GOOD" : "BAD"; print "Expected: "target1", "target2", "target3; print "Actual: "actual1", "actual2", "actual3; }' \
    > "${RESULTS}.report"
