#!/bin/sh -e
MMSEQS="${1}"
RESULTS="$( (readlink -f -- "$2" 2>/dev/null) || (greadlink -f -- "$2") )"
mkdir -p "${RESULTS}"

BASE="$(dirname "$( (readlink -f -- "$0" 2>/dev/null) || (greadlink -f -- "$0") )")"
cd "${BASE}"

# build the benchmark tools
(mkdir -p build && cd build && cmake -DHAVE_FFINDEX=0 -DCMAKE_BUILD_TYPE=Release .. && make -j4)

DATADIR="${BASE}/data"
SCRIPTS="${BASE}/regression"
EVALUATE="${BASE}/build/evaluate_results"
SAMTOOLS="${BASE}/samtools/samtools.sh"

TESTS=""
# continue on if one test fails
set +e
time "${SCRIPTS}/run_search.sh" "${MMSEQS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/SEARCH"; TESTS="SEARCH ${TESTS}"
time "${SCRIPTS}/run_easy_search.sh" "${MMSEQS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/EASY_SEARCH"; TESTS="EASY_SEARCH ${TESTS}"
time "${SCRIPTS}/run_profile.sh" "${MMSEQS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/PROFILE"; TESTS="PROFILE ${TESTS}"
time "${SCRIPTS}/run_easy_profile.sh" "${MMSEQS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/EASY_PROFILE"; TESTS="EASY_PROFILE ${TESTS}"
time "${SCRIPTS}/run_sliceprofile.sh" "${MMSEQS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/SLICEPROFILE"; TESTS="SLICEPROFILE ${TESTS}"
time "${SCRIPTS}/run_dbprofile.sh" "${MMSEQS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/DBPROFILE"; TESTS="DBPROFILE ${TESTS}"
time "${SCRIPTS}/run_nuclprot.sh" "${MMSEQS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/NUCLPROT_SEARCH"; TESTS="NUCLPROT_SEARCH ${TESTS}"
time "${SCRIPTS}/run_nuclnucl.sh" "${MMSEQS}" "${SAMTOOLS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/NUCLNUCL_SEARCH"; TESTS="NUCLNUCL_SEARCH ${TESTS}"
time "${SCRIPTS}/run_nuclnucl_translated.sh" "${MMSEQS}" "${SAMTOOLS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/NUCLNUCL_TRANS_SEARCH"; TESTS="NUCLNUCL_TRANS_SEARCH ${TESTS}"
time "${SCRIPTS}/run_cluster.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/CLUSTER"; TESTS="CLUSTER ${TESTS}"
time "${SCRIPTS}/run_easy_cluster.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/EASY_CLUSTER"; TESTS="EASY_CLUSTER ${TESTS}"
time "${SCRIPTS}/run_linclust.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/LINCLUST"; TESTS="LINCLUST ${TESTS}"
time "${SCRIPTS}/run_easy_linclust.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/EASY_LINCLUST"; TESTS="EASY_LINCLUST ${TESTS}"
time "${SCRIPTS}/run_protnucl.sh" "${MMSEQS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/PROTNUCL_SEARCH"; TESTS="PROTNUCL_SEARCH ${TESTS}"
time "${SCRIPTS}/run_nuclprottax.sh" "${MMSEQS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/NUCLPROTTAX_SEARCH"; TESTS="NUCLPROTTAX_SEARCH ${TESTS}"
time "${SCRIPTS}/run_dbprofile_index.sh" "${MMSEQS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/DBPROFILE_INDEX"; TESTS="DBPROFILE_INDEX ${TESTS}"
time "${SCRIPTS}/run_nuclnucl_linsearchtranslated.sh" "${MMSEQS}" "${SAMTOOLS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/LINSEARCH_NUCLNUCL_TARNS_SEARCH"; TESTS="LINSEARCH_NUCLNUCL_TARNS_SEARCH ${TESTS}"
time "${SCRIPTS}/run_nuclnucl_linsearch.sh" "${MMSEQS}" "${SAMTOOLS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/LINSEARCH_NUCLNUCL_SEARCH"; TESTS="LINSEARCH_NUCLNUCL_SEARCH ${TESTS}"
time "${SCRIPTS}/run_cluster_update.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/LINCLUST_UPDATE"; TESTS="LINCLUST_UPDATE ${TESTS}"
# time "${SCRIPTS}/run_multihit.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/MULTHIT"; TESTS="MULTHIT ${TESTS}"
time "${SCRIPTS}/run_easynuclnucltax.sh" "${MMSEQS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/EASYNUCLNUCLTAX_SEARCH"; TESTS="EASYNUCLNUCLTAX_SEARCH ${TESTS}"
time "${SCRIPTS}/run_extractorfs.sh" "${MMSEQS}" "${SCRIPTS}/extractorfs.pl" "${SCRIPTS}/compare_frags.pl" "${DATADIR}" "${RESULTS}/EXTRACTORFS"; TESTS="EXTRACTORFS ${TESTS}"
time "${SCRIPTS}/run_rbh.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/RBH"; TESTS="RBH ${TESTS}"
time "${SCRIPTS}/run_apply.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/APPLY"; TESTS="APPLY ${TESTS}"
# time "${SCRIPTS}/run_profilestates.sh" "${MMSEQS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/CSPROFILE"; TESTS="CSPROFILE ${TESTS}"

case "$("${MMSEQS}" version)" in
	*MPI)
		export RUNNER="mpirun -np 1"
		time "${SCRIPTS}/run_split.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/MPI_TARGET_SPLIT_NP1" 0; TESTS="MPI_TARGET_SPLIT_NP1 ${TESTS}"
		time "${SCRIPTS}/run_split.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/MPI_QUERY_SPLIT_NP1" 1; TESTS="MPI_QUERY_SPLIT_NP1 ${TESTS}"
		time "${SCRIPTS}/run_slicetechnical.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/MPI_SLICE_TECH_NP1"; TESTS="MPI_SLICE_TECH_NP1 ${TESTS}"
		
		export RUNNER="mpirun -np 3"
		time "${SCRIPTS}/run_split.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/MPI_TARGET_SPLIT_NP3" 0; TESTS="MPI_TARGET_SPLIT_NP3 ${TESTS}"
		time "${SCRIPTS}/run_split.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/MPI_QUERY_SPLIT_NP3" 1; TESTS="MPI_QUERY_SPLIT_NP3 ${TESTS}"
		time "${SCRIPTS}/run_slicetechnical.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/MPI_SLICE_TECH_NP3"; TESTS="MPI_SLICE_TECH_NP3 ${TESTS}"
		
		unset RUNNER
		;;
	*)
		time "${SCRIPTS}/run_split.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/NOMPI_TARGET_SPLIT" 0; TESTS="NOMPI_TARGET_SPLIT ${TESTS}"
		time "${SCRIPTS}/run_slicetechnical.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/NOMPI_SLICE_TECH"; TESTS="NOMPI_SLICE_TECH ${TESTS}"
esac

set -e

printf "\n"
ERR=0
for i in ${TESTS} ; do
    echo "$i"
    if [ ! -f "${RESULTS}/${i}/report" ]; then
        printf "TEST FAILED (NO REPORT)\n\n"
        ERR=$((ERR+1))
        continue
    fi

    if [ ! -s "${RESULTS}/${i}/report" ]; then
        printf "TEST FAILED (EMPTY REPORT)\n\n"
        ERR=$((ERR+1))
        continue
    fi
    STATUS="$(head -n 1 "${RESULTS}/${i}/report")"
    if [ "$STATUS" != "GOOD" ]; then
        echo "TEST FAILED"
        ERR=$((ERR+1))
    else
        echo "TEST SUCCESS"
    fi
    cat "${RESULTS}/${i}/report"
    printf "\n"
done

exit "$ERR"
