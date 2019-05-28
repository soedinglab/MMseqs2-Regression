#!/bin/sh -e
MMSEQS="${1}"
DATADIR="${2}"
RESULTS="${3}"
mkdir -p "${RESULTS}"

"${MMSEQS}" createdb "${DATADIR}/query.fasta" "${RESULTS}/query"
"${MMSEQS}" apply "${RESULTS}/query" "$RESULTS/apply" -- wc -c

ACTUAL="$(tr -d '\000' < /Users/mirdita/repositories/mmseqs-benchmark/APPLY/apply | awk '{ l += $1; } END { print l }')"
TARGET="$(grep -v "^>" "${DATADIR}/query.fasta" | wc -c)"
awk -v actual="$ACTUAL" -v target="$TARGET" \
    'BEGIN { print (actual == target) ? "GOOD" : "BAD"; print "Expected: ", target; print "Actual: ", actual; }' \
    > "${RESULTS}/report"

