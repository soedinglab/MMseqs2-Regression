#!/bin/bash -e
BENCHDIR="${1}"
MMSEQS="${2}"
SCRIPTSDIR="${3}"
VERSION="${4}"
RESULTDIR="${5}"
NAME="${6:-regression-test}"
THREADS="${7:-4}"

# rbh test #
# RBHproteinsA.fas has 6 sequences
# RBHproteinsB.fas has 6 sequences
# in each file - one sequence has no match at all
# in each file all other sequences match all other sequences in the other file
# the best matching is:
# seqA1 with seqB1
# seqA2 with seqB2
# seqA2 with seqB2_also_best
# seqA3 with seqB3
# seqA4 with seqB4

APROTEINS="$BENCHDIR/data/RBHproteinsA.fas"
BPROTEINS="$BENCHDIR/data/RBHproteinsB.fas"

RESULTS="${RESULTDIR}/mmseqs-${NAME}-${VERSION}"
mkdir -p "${RESULTS}"

#  1>&2 || exit 125 --> move output to stderr
${MMSEQS} createdb "${APROTEINS}" "${RESULTS}/proteinsA" 1>&2 || exit 125
${MMSEQS} createdb "${BPROTEINS}" "${RESULTS}/proteinsB" 1>&2 || exit 125
${MMSEQS} rbh "${RESULTS}/proteinsA" "${RESULTS}/proteinsB" "${RESULTS}/rbhAB" "${RESULTS}/tmp" 1>&2 || exit 125
${MMSEQS} convertalis "${RESULTS}/proteinsA" "${RESULTS}/proteinsB" "${RESULTS}/rbhAB" "${RESULTS}/rbhAB.m8" 1>&2 || exit 125

# both of these should be 5
TOTAL_NUM_LINES=$(cat "${RESULTS}/rbhAB.m8" | wc -l)
NUM_GOOD_MATCHES=$(grep -P "seqA(\d)(_also_best)?\tseqB\1(_also_best)?\t" "${RESULTS}/rbhAB.m8" | wc -l)
echo -e "${NAME}\t${VERSION}\t23\t${TOTAL_NUM_LINES}"
echo -e "${NAME}\t${VERSION}\t24\t${NUM_GOOD_MATCHES}"

