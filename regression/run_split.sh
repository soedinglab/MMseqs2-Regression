#!/bin/sh -e
MMSEQS="${1}"
DATADIR="${2}"
RESULTS="${3}"
SPLIT_MODE="${4}"

mkdir -p "${RESULTS}"
cd "${RESULTS}"
head -n 1000 "${DATADIR}/query.fasta" > "500prots.fasta"
"${MMSEQS}" createdb "500prots.fasta" prots --dbtype 0

NUM_RES_LINES_BASELINE=500

NUM_LINES_ALL_TESTS="${NUM_RES_LINES_BASELINE}"
for NUM_SPLITS in 1 5
do
	# target split:
	OUTPUT_SUBDIR="NS_${NUM_SPLITS}"
	echo "working on: ${OUTPUT_SUBDIR}"
	mkdir -p "${OUTPUT_SUBDIR}"
	cd "${OUTPUT_SUBDIR}"
	"${MMSEQS}" search ../prots ../prots res tmpFolder -s 1 --split-mode "${SPLIT_MODE}" --split "${NUM_SPLITS}"
	"${MMSEQS}" convertalis ../prots ../prots res res.m8
	NUM_LINES_CURR_TEST="$(wc -l < res.m8)"
	if [ "${NUM_LINES_CURR_TEST}" -ne "${NUM_RES_LINES_BASELINE}" ]; then
		NUM_LINES_ALL_TESTS="${NUM_LINES_CURR_TEST}"
	fi
	cd "${RESULTS}"
done

# if ALL tests passed - GOOD, otherwise - BAD
awk -v actual="$NUM_LINES_ALL_TESTS" -v target="$NUM_RES_LINES_BASELINE" \
	'BEGIN { print (actual == target) ? "GOOD" : "BAD"; print "Expected: ", target; print "Actual: ", actual; }' \
	> "${RESULTS}/report"