#!/bin/bash -e
BENCHDIR="${1}"
MMSEQS="${2}"
SCRIPTSDIR="${3}"
VERSION="${4}"
RESULTDIR="${5}"
NAME="${6:-regression-test}"
THREADS="${7:-4}"

QUERY="$BENCHDIR/data/dna.fas"

RESULTS="${RESULTDIR}/mmseqs-${NAME}-${VERSION}"
mkdir -p "${RESULTS}"

${MMSEQS} createdb "${QUERY}" "${RESULTS}/dna" 1>&2 || exit 125
${MMSEQS} extractorfs "${RESULTS}/dna" "${RESULTS}/coding_frags_mode_0" --threads ${THREADS} --min-length 1 --orf-start-mode 0 --contig-start-mode 2 --contig-end-mode 2 --forward-frames 1,2,3 --reverse-frames 1,2,3 1>&2 || exit 125
${MMSEQS} extractorfs "${RESULTS}/dna" "${RESULTS}/coding_frags_mode_1" --threads ${THREADS} --min-length 1 --orf-start-mode 1 --contig-start-mode 2 --contig-end-mode 2 --forward-frames 1,2,3 --reverse-frames 1,2,3 1>&2 || exit 125
${MMSEQS} extractorfs "${RESULTS}/dna" "${RESULTS}/coding_frags_mode_2" --threads ${THREADS} --min-length 1 --orf-start-mode 2 --contig-start-mode 2 --contig-end-mode 2 --forward-frames 1,2,3 --reverse-frames 1,2,3 1>&2 || exit 125
perl "${SCRIPTSDIR}/extractorfs.pl" "${QUERY}" "${RESULTS}/perl_coding_frags" 1>&2 || exit 125
RES0="$(perl ${SCRIPTSDIR}/compare_frags.pl ${RESULTS}/perl_coding_frags_mode_0.txt ${RESULTS}/coding_frags_mode_0 1>&2; echo $?)"
RES1="$(perl ${SCRIPTSDIR}/compare_frags.pl ${RESULTS}/perl_coding_frags_mode_1.txt ${RESULTS}/coding_frags_mode_1 1>&2; echo $?)"
RES2="$(perl ${SCRIPTSDIR}/compare_frags.pl ${RESULTS}/perl_coding_frags_mode_2.txt ${RESULTS}/coding_frags_mode_2 1>&2; echo $?)"
echo -e "${NAME}\t${VERSION}\t15\t$((RES0+RES1+RES2))"

