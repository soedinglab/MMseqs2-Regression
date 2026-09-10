#!/bin/sh -e
fake_pref() {
    QDB="$1"
    TDB="$2"
    RES="$3"
    # create link to data file which contains a list of all targets that should be aligned
    ln -s "${TDB}.index" "${RES}"
    # create new index repeatedly pointing to same entry
    INDEX_SIZE="$(echo $(wc -c < "${TDB}.index"))"
    awk -v size=$INDEX_SIZE '{ print $1"\t0\t"size; }' "${QDB}.index" > "${RES}.index"
    # create dbtype (7)
    awk 'BEGIN { printf("%c%c%c%c",7,0,0,0); exit; }' > "${RES}.dbtype"
}

mkdir -p "${RESULTS}/fwbw" #makes sure index files do not get mixed 

QUERY="${DATADIR}/query.fasta"
QUERYDB="${RESULTS}/fwbw/query_fwbw"
"${MMSEQS}" createdb "${QUERY}" "${QUERYDB}"

TARGET="${DATADIR}/targetannotation.fasta"
TARGETDB="${RESULTS}/fwbw/targetannotation_fwbw"
"${MMSEQS}" createdb "${TARGET}" "${TARGETDB}"

"${MMSEQS}" prefilter "$QUERYDB" "${TARGETDB}" "$RESULTS/fwbw/pref"
awk 'BEGIN { printf("%c%c%c%c",5,0,0,0); exit; }' > "$RESULTS/fwbw/pref.dbtype"
"${MMSEQS}" fwbw "$QUERYDB" "$TARGETDB" "$RESULTS/fwbw/pref" "$RESULTS/fwbw/results_fwbw" -e 10000
"${MMSEQS}" convertalis "$QUERYDB" "$TARGETDB" "$RESULTS/fwbw/results_fwbw" "$RESULTS/fwbw/results_fwbw.m8"
"${EVALUATE}" "$QUERY" "$TARGET" "$RESULTS/fwbw/results_fwbw.m8" "${RESULTS}/fwbw/evaluation_roc5.dat" 4000 1 | tee "${RESULTS}/fwbw/evaluation.log"
ACTUAL=$(grep "^ROC5 AUC:" "${RESULTS}/fwbw/evaluation.log" | cut -d" " -f3)
TARGET="0.207452"
awk -v actual="$ACTUAL" -v target="$TARGET" \
    'BEGIN { print (actual == target) ? "GOOD" : "BAD"; print "Expected: ", target; print "Actual: ", actual; }' \
    > "${RESULTS}.report"
