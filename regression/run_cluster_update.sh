#!/bin/sh -e
MMSEQS="${1}"
DATADIR="${2}"
RESULTS="${3}"
mkdir -p "${RESULTS}"

SEQCLUDB1="${RESULTS}/clu1"
SEQCLUDB2="${RESULTS}/clu2"
awk 'NR%4==0 || NR%4==1{print}' "${DATADIR}/clu.fasta" > "$RESULTS/clu1.fasta"
awk 'NR%4==2 || NR%4==3{print}' "${DATADIR}/clu.fasta" > "$RESULTS/clu2.fasta"
head -n 2 "$RESULTS/clu1.fasta" >> "$RESULTS/clu2.fasta"
cat "$RESULTS/clu1.fasta" "$RESULTS/clu2.fasta" > "$RESULTS/cluCombined.fasta"

"${MMSEQS}" createdb "$RESULTS/clu1.fasta" "${SEQCLUDB1}"
"${MMSEQS}" createdb "$RESULTS/cluCombined.fasta" "${SEQCLUDB2}"

"${MMSEQS}" linclust "${SEQCLUDB1}" "$RESULTS/results_clu" "$RESULTS/tmp" --cov-mode 1 -a -c 0.50 --min-seq-id 0.50
"${MMSEQS}" clusterupdate "${SEQCLUDB1}" "${SEQCLUDB2}" "$RESULTS/results_clu" "$RESULTS/seqdb_update" "$RESULTS/clu_updated" "$RESULTS/tmp"

CLUSTERMEMEBER=$(wc -l "$RESULTS/clu_updated" | awk '{print $1}')
CLUSTER=$(wc -l "$RESULTS/clu_updated.index" | awk '{print $1}')
UPDATEDSEQCNT=$(wc -l "$RESULTS/seqdb_update" | awk '{print $1}')

TARGET="32130 17403 32130"
ACTUAL="$CLUSTERMEMEBER $CLUSTER $UPDATEDSEQCNT"
awk -v actual="$ACTUAL" -v target="$TARGET" 'BEGIN { print (actual == target) ? "GOOD" : "BAD"; \
    print "Expected: ", target; \
    print "Actual:   ", actual; }' > "${RESULTS}/report"
