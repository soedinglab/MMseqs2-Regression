#!/bin/sh -ex
QUERY="${DATADIR}/query.fasta"
TARGET="${DATADIR}/targetannotation.fasta"
TARGETDB="${RESULTS}/targetannotation"
TARGETDB_MAPPING="${DATADIR}/targetannotation.mapping"
"${MMSEQS}" createdb "${TARGET}" "${TARGETDB}"
"${MMSEQS}" createtaxdb "${TARGETDB}" "$RESULTS/tmp" --tax-mapping-file "${TARGETDB_MAPPING}" --ncbi-tax-dump "${DATADIR}/ncbitax"

cat "$QUERY" | "${MMSEQS}" easy-search stdin "${TARGETDB}" "$RESULTS/results_aln_all.m8" "$RESULTS/tmp" --exact-kmer-matching 1 --max-seqs 5 --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage
cat "$QUERY" | "${MMSEQS}" easy-search stdin "${TARGETDB}" "$RESULTS/results_aln_filter.m8" "$RESULTS/tmp" --exact-kmer-matching 1 --max-seqs 5 --taxon-list 40674 --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage

awk '{ print $1"-"$2; }' "$RESULTS/results_aln_filter.m8" | sort > "$RESULTS/results_aln_filter_ids.list"
awk -F'\t' '$15 ~ "c_Mammalia" { print $1"-"$2; }' "$RESULTS/results_aln_all.m8" | sort > "$RESULTS/results_aln_all_ids.list"

ACTUAL=$(
    comm "$RESULTS/results_aln_filter_ids.list"  "$RESULTS/results_aln_all_ids.list"\
        | awk -F'\t' 'BEGIN { l=0; r=0; b=0; } $1 != "" { l+=1; } $2 != "" { r+=1; } $3 != "" { b+=1; } END { print l","r","b}'
)


TARGET="31,0,763"
awk -v actual="$ACTUAL" -v target="$TARGET" \
    'BEGIN { print (actual == target) ? "GOOD" : "BAD"; print "Expected: ", target; print "Actual: ", actual; }' \
    > "${RESULTS}.report"
