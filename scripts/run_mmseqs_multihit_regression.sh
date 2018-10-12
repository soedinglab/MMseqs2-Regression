#!/bin/bash -e

BENCHDIR="${1}"
MMSEQS="${2}"
VERSION="${3}"
RESULTDIR="${4}"
NAME="${5:-regression-test}"
THREADS="${6:-4}"

BENCHDB="small-benchmark-db"

mkdir -p "$RESULTDIR" && cd "$RESULTDIR"

##### multihitdb #####
$MMSEQS multihitdb "$BENCHDIR/$BENCHDB/qset_01.fas.gz" "$BENCHDIR/$BENCHDB/qset_02.fas.gz" qsetsdb tmp 1>&2 || exit 125
$MMSEQS multihitdb "$BENCHDIR/$BENCHDB/tset_01.fas.gz" "$BENCHDIR/$BENCHDB/tset_02.fas.gz" tsetsdb tmp 1>&2 || exit 125

##### multihitsearch #####
$MMSEQS multihitsearch qsetsdb tsetsdb result tmp --threads ${THREADS} 1>&2 || exit 125

##### combinepvalperset #####
$MMSEQS combinepvalperset qsetsdb tsetsdb result pval tmp --threads ${THREADS} 1>&2 || exit 125

EVAL1="$(tr -d '\000' < pval | head -n1 | cut -f2)"
EVAL2="$(tr -d '\000' < pval | tail -n1 | cut -f2)"

echo -e "${NAME}\t${VERSION}\t12\t${EVAL1}"
echo -e "${NAME}\t${VERSION}\t13\t${EVAL2}"

