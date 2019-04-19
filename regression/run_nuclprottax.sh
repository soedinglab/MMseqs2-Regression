#!/bin/sh -e
MMSEQS="${1}"
EVALUATE="${2}"
DATADIR="${3}"
RESULTS="${4}"
mkdir -p "${RESULTS}"

QUERY="${DATADIR}/query.fasta"
QUERYDB="${RESULTS}/query"
"${MMSEQS}" createdb "${QUERY}" "${QUERYDB}"
"${MMSEQS}" translateaa "${QUERYDB}" "${QUERYDB}_nucl"
ln -sf "${QUERYDB}_h" "${QUERYDB}_nucl_h"
ln -sf "${QUERYDB}_h.index" "${QUERYDB}_nucl_h.index"
ln -sf "${QUERYDB}_h.dbtype" "${QUERYDB}_nucl_h.dbtype"

TARGET="${DATADIR}/targetannotation.fasta"
TARGETDB="${RESULTS}/targetannotation"
TARGETDB_MAPPING="${DATADIR}/targetannotation.mapping"
"${MMSEQS}" createdb "${TARGET}" "${TARGETDB}"
"${MMSEQS}" createtaxdb "${TARGETDB}" "${TARGETDB_MAPPING}" "${DATADIR}/ncbitax" "$RESULTS/tmp"  
"${MMSEQS}" taxonomy "${QUERYDB}_nucl" "$TARGETDB" "$RESULTS/results_aln" "$RESULTS/tmp" --blacklist "0" -e 10000 -s 4 --max-seqs 4000
"${MMSEQS}" filtertaxdb "${TARGETDB}" "$RESULTS/results_aln" "$RESULTS/results_aln_bacteria" --taxon-list 2 
"${MMSEQS}" filtertaxdb "${TARGETDB}" "$RESULTS/results_aln" "$RESULTS/results_aln_virus" --taxon-list 10239
"${MMSEQS}" filtertaxdb "${TARGETDB}" "$RESULTS/results_aln" "$RESULTS/results_aln_eukaryota" --taxon-list 2759 

BACTERIA=$(awk '$3 != 1 {print}' "$RESULTS/results_aln_bacteria.index" | wc -l| awk '{print $1}')
VIRUS=$(awk '$3 != 1 {print}' "$RESULTS/results_aln_virus.index" | wc -l| awk '{print $1}')
EUKARYOTA=$(awk '$3 != 1 {print}' "$RESULTS/results_aln_eukaryota.index" | wc -l| awk '{print $1}')

TARGET="2524 259 2713"
ACTUAL="$BACTERIA $VIRUS $EUKARYOTA"
awk -v actual="$ACTUAL" -v target="$TARGET" \
    'BEGIN { print (actual == target) ? "GOOD" : "BAD"; print "Expected: ", target; print "Actual: ", actual; }' \
    > "${RESULTS}/report"
