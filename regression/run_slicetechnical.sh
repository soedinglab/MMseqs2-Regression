#!/bin/sh -e
MMSEQS="${1}"
DATADIR="${2}"
RESULTS="${3}"

mkdir -p "${RESULTS}"
cd "${RESULTS}"

# create profiles db from fasta #
"${MMSEQS}" createdb "${DATADIR}/five_profiles.fasta" prof_5
"${MMSEQS}" mergedbs prof_5 prof_5_fasta prof_5_h prof_5 --prefixes ">"
"${MMSEQS}" msa2profile prof_5_fasta five_profiles --filter-msa 0

# create db from duplicated query #
"${MMSEQS}" createdb "${DATADIR}/query_pow2.fasta" query_pow2

PROFILES_DB="five_profiles"
PROTEINS_DB="query_pow2"

FINAL_COUNTS_AS_SHOULD="512,256,128,64,32"
COUNTS_ALL_TESTS="${FINAL_COUNTS_AS_SHOULD}"

for SPLIT_MODE in 0 1
do
	DISK_SPACE_LIMIT_KB=17
	OUTPUT_SUBDIR="DSL_${DISK_SPACE_LIMIT_KB}K_SPLIT_MODE_${SPLIT_MODE}"
	echo "working on: ${OUTPUT_SUBDIR}"
	mkdir -p "${OUTPUT_SUBDIR}"
	cd "${OUTPUT_SUBDIR}"
	# run slice search #
	"${MMSEQS}" search ../"${PROTEINS_DB}" ../"${PROFILES_DB}" res tmpFolder --slice-search -s 1 --disk-space-limit "${DISK_SPACE_LIMIT_KB}K" --split-mode "${SPLIT_MODE}"
	"${MMSEQS}" convertalis ../"${PROTEINS_DB}" ../"${PROFILES_DB}" res alis
	awk '{print $1}' alis | sort | uniq -c | sort -nr | awk -v ORS=, '{ print $1 }' | sed 's/,$/\n/' > final_counts.txt
	COUNTS_CURR_TEST="$(cat < final_counts.txt)"
	echo "$COUNTS_CURR_TEST"
	if [ "$COUNTS_CURR_TEST" != "$FINAL_COUNTS_AS_SHOULD" ]; then
		COUNTS_ALL_TESTS="$COUNTS_CURR_TEST"
	fi
	cd "${RESULTS}"
done

# if ALL tests passed - GOOD, otherwise - BAD
awk -v actual="$COUNTS_ALL_TESTS" -v target="$FINAL_COUNTS_AS_SHOULD" \
	'BEGIN { print (actual == target) ? "GOOD" : "BAD"; print "Expected: ", target; print "Actual: ", actual; }' \
	> "${RESULTS}/report"
