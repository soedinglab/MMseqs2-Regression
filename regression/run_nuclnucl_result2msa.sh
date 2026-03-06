#!/bin/sh -e
# Regression test: nucleotide search → result2msa WITHOUT backtrace (-a).
#
# Tests both plain and padded (makepaddedseqdb + createindex) target databases.
#
# This exercises the fix for result2msa where:
#   1. NucleotideMatrix must be used (not SubstitutionMatrix) for NT databases
#   2. Backtraces must be recomputed from alignment records that lack them
#      (e.g. 14-column records from offsetalignment in the blastn workflow)
#
# Without the fix, this produces "DUMMY" all-gap MSA rows or segfaults.

# -- clean previous run --
rm -rf "${RESULTS:?}"/*

# -- create inline nucleotide FASTA files --
cat > "${RESULTS}/target.fasta" <<'EOF'
>seq1 identical_to_query
CTGCAGCTTGCCCTCAGAGACCGATCTCTCAGAGAGGTACATGGAATCGTGTTCCATCCCTGGATAACGGAACTCTCAGTCCTGCAG
>seq2 three_snps
CTGCAGCTTGCCCTCATAGACCGATCTCTCAGAGAGGTACATCGAATCGTGTTCCATCCATGGATAACGGAACTCTCAGTCCTGCAG
>seq3 short_deletion
CTGCAGCTTGCCCTCAGAGACCGATCTCTCAGAGAGCGTGTTCCATCCCTGGATAACGGAACTCTCAGTCCTGCAG
EOF

cat > "${RESULTS}/query.fasta" <<'EOF'
>query test_sequence
CTGCAGCTTGCCCTCAGAGACCGATCTCTCAGAGAGGTACATGGAATCGTGTTCCATCCCTGGATAACGGAACTCTCAGTCCTGCAG
EOF

# -- build databases --
"${MMSEQS}" createdb "${RESULTS}/target.fasta" "${RESULTS}/targetdb" --dbtype 2
"${MMSEQS}" createdb "${RESULTS}/query.fasta"  "${RESULTS}/querydb"  --dbtype 2

# padded + indexed variant
"${MMSEQS}" makepaddedseqdb "${RESULTS}/targetdb" "${RESULTS}/targetdb_padded"
"${MMSEQS}" createindex "${RESULTS}/targetdb_padded" "${RESULTS}/tmp_idx" \
    --remove-tmp-files 1 --split 1 --index-subset 0 --search-type 3

ERR=0

validate_a3m() {
    local label=$1 a3m=$2

    if [ ! -s "$a3m" ]; then
        echo "FAIL [$label]: result.a3m is empty or missing"
        ERR=$((ERR + 1))
    fi

    if grep -q "DUMMY" "$a3m"; then
        echo "FAIL [$label]: result.a3m contains DUMMY placeholder"
        ERR=$((ERR + 1))
    fi

    SEQ_COUNT=$(grep -c "^>" "$a3m" || true)
    if [ "$SEQ_COUNT" -lt 2 ]; then
        echo "FAIL [$label]: expected >=2 sequences, got ${SEQ_COUNT}"
        ERR=$((ERR + 1))
    fi

    if grep -qE "^-+$" "$a3m"; then
        echo "FAIL [$label]: result.a3m contains all-gap rows"
        ERR=$((ERR + 1))
    fi
}

# -- Test 1: plain (non-padded) target DB --
"${MMSEQS}" search "${RESULTS}/querydb" "${RESULTS}/targetdb" \
    "${RESULTS}/result_plain" "${RESULTS}/tmp" \
    --search-type 3

"${MMSEQS}" result2msa "${RESULTS}/querydb" "${RESULTS}/targetdb" \
    "${RESULTS}/result_plain" "${RESULTS}/result_plain.a3m" \
    --msa-format-mode 6

validate_a3m "plain" "${RESULTS}/result_plain.a3m"

# -- Test 2: padded + indexed target DB --
rm -rf "${RESULTS}/tmp"

"${MMSEQS}" search "${RESULTS}/querydb" "${RESULTS}/targetdb_padded" \
    "${RESULTS}/result_padded" "${RESULTS}/tmp" \
    --search-type 3

"${MMSEQS}" result2msa "${RESULTS}/querydb" "${RESULTS}/targetdb_padded" \
    "${RESULTS}/result_padded" "${RESULTS}/result_padded.a3m" \
    --msa-format-mode 6

validate_a3m "padded" "${RESULTS}/result_padded.a3m"

# -- report --
awk -v actual="$ERR" -v target="0" \
    'BEGIN { print (actual == target) ? "GOOD" : "BAD"; print "Expected: ", target, "errors"; print "Actual: ", actual, "errors"; }' \
    > "${RESULTS}.report"
