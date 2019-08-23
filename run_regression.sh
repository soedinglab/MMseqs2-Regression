#!/bin/sh -e
abspath() {
    if [ -d "$1" ]; then
        (cd "$1"; pwd)
    elif [ -f "$1" ]; then
        if [ -z "${1##*/*}" ]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    elif [ -d "$(dirname "$1")" ]; then
        echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
    fi
}

MMSEQS="$(abspath "$1")"
RESULTS="$(abspath "$2")"

mkdir -p "${RESULTS}"

BASE="$(dirname "$(abspath "$0")")"
cd "${BASE}"

# build the benchmark tools
(mkdir -p build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make -j4)

DATADIR="${BASE}/data"
SCRIPTS="${BASE}/regression"
EVALUATE="${BASE}/build/evaluate_results"
SAMTOOLS="${BASE}/samtools/samtools.sh"

TESTS=""
run_test() {
	NAME="$1"
	TESTS="${TESTS} ${NAME}"
	shift
	START="$(date +%s)"
	"$@"
	END="$(date +%s)"
    eval "${NAME}_TIME"="$((END-START))"
}

# continue on if one test fail
set +e

run_test SEARCH "${SCRIPTS}/run_search.sh" "${MMSEQS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/SEARCH"
run_test EASY_SEARCH "${SCRIPTS}/run_easy_search.sh" "${MMSEQS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/EASY_SEARCH"
run_test EASY_SEARCH_INDEX_SPLIT "${SCRIPTS}/run_easy_search_index_split.sh" "${MMSEQS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/EASY_SEARCH_INDEX_SPLIT"
run_test PROFILE "${SCRIPTS}/run_profile.sh" "${MMSEQS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/PROFILE"
run_test EASY_PROFILE "${SCRIPTS}/run_easy_profile.sh" "${MMSEQS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/EASY_PROFILE"
run_test SLICEPROFILE "${SCRIPTS}/run_sliceprofile.sh" "${MMSEQS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/SLICEPROFILE"
run_test DBPROFILE "${SCRIPTS}/run_dbprofile.sh" "${MMSEQS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/DBPROFILE"
run_test NUCLPROT_SEARCH "${SCRIPTS}/run_nuclprot.sh" "${MMSEQS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/NUCLPROT_SEARCH"
run_test NUCLNUCL_SEARCH "${SCRIPTS}/run_nuclnucl.sh" "${MMSEQS}" "${SAMTOOLS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/NUCLNUCL_SEARCH"
run_test NUCLNUCL_TRANS_SEARCH "${SCRIPTS}/run_nuclnucl_translated.sh" "${MMSEQS}" "${SAMTOOLS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/NUCLNUCL_TRANS_SEARCH"
run_test CLUSTER "${SCRIPTS}/run_cluster.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/CLUSTER"
run_test EASY_CLUSTER "${SCRIPTS}/run_easy_cluster.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/EASY_CLUSTER"
run_test CLUSTER_REASSIGN "${SCRIPTS}/run_easy_cluster_reassign.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/CLUSTER_REASSIGN"
run_test LINCLUST "${SCRIPTS}/run_linclust.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/LINCLUST"
run_test LINCLUST_SPLIT "${SCRIPTS}/run_linclust_split.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/LINCLUST_SPLIT"
run_test EASY_LINCLUST "${SCRIPTS}/run_easy_linclust.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/EASY_LINCLUST"
run_test PROTNUCL_SEARCH "${SCRIPTS}/run_protnucl.sh" "${MMSEQS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/PROTNUCL_SEARCH"
run_test NUCLPROTTAX_SEARCH "${SCRIPTS}/run_nuclprottax.sh" "${MMSEQS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/NUCLPROTTAX_SEARCH"
run_test DBPROFILE_INDEX "${SCRIPTS}/run_dbprofile_index.sh" "${MMSEQS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/DBPROFILE_INDEX"
run_test LINSEARCH_NUCLNUCL_TARNS_SEARCH "${SCRIPTS}/run_nuclnucl_linsearchtranslated.sh" "${MMSEQS}" "${SAMTOOLS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/LINSEARCH_NUCLNUCL_TARNS_SEARCH"
run_test LINSEARCH_NUCLNUCL_SEARCH "${SCRIPTS}/run_nuclnucl_linsearch.sh" "${MMSEQS}" "${SAMTOOLS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/LINSEARCH_NUCLNUCL_SEARCH"
run_test EASY_LINSEARCH_NUCLNUCL_SEARCH_SPLIT "${SCRIPTS}/run_easy_nuclnucl_linsearch_split.sh" "${MMSEQS}" "${SAMTOOLS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/EASY_LINSEARCH_NUCLNUCL_SEARCH_SPLIT"
run_test LINCLUST_UPDATE "${SCRIPTS}/run_cluster_update.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/LINCLUST_UPDATE"
run_test EASYNUCLNUCLTAX_SEARCH "${SCRIPTS}/run_easynuclnucltax.sh" "${MMSEQS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/EASYNUCLNUCLTAX_SEARCH"
run_test EXTRACTORFS "${SCRIPTS}/run_extractorfs.sh" "${MMSEQS}" "${SCRIPTS}/extractorfs.pl" "${SCRIPTS}/compare_frags.pl" "${DATADIR}" "${RESULTS}/EXTRACTORFS"
run_test RBH "${SCRIPTS}/run_rbh.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/RBH"
run_test APPLY "${SCRIPTS}/run_apply.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/APPLY"
run_test INDEX_COMPATIBLE "${SCRIPTS}/run_index_compatible.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/INDEX_COMPATIBLE"
# run_test MULTHIT "${SCRIPTS}/run_multihit.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/MULTHIT"
# run_test CSPROFILE "${SCRIPTS}/run_profilestates.sh" "${MMSEQS}" "${EVALUATE}" "${DATADIR}" "${RESULTS}/CSPROFILE"

case "$("${MMSEQS}" version)" in
	*MPI)
		export RUNNER="mpirun -np 1"
		run_test MPI_TARGET_SPLIT_NP1 "${SCRIPTS}/run_split.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/MPI_TARGET_SPLIT_NP1" 0
		run_test MPI_QUERY_SPLIT_NP1 "${SCRIPTS}/run_split.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/MPI_QUERY_SPLIT_NP1" 1
		run_test MPI_SLICE_TECH_NP1 "${SCRIPTS}/run_slicetechnical.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/MPI_SLICE_TECH_NP1"
		
		export RUNNER="mpirun -np 3"
		run_test MPI_TARGET_SPLIT_NP3 "${SCRIPTS}/run_split.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/MPI_TARGET_SPLIT_NP3" 0
		run_test MPI_QUERY_SPLIT_NP3 "${SCRIPTS}/run_split.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/MPI_QUERY_SPLIT_NP3" 1
		run_test MPI_SLICE_TECH_NP3 "${SCRIPTS}/run_slicetechnical.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/MPI_SLICE_TECH_NP3"
		
		unset RUNNER
		;;
	*)
		run_test NOMPI_TARGET_SPLIT "${SCRIPTS}/run_split.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/NOMPI_TARGET_SPLIT" 0
		run_test NOMPI_SLICE_TECH "${SCRIPTS}/run_slicetechnical.sh" "${MMSEQS}" "${DATADIR}" "${RESULTS}/NOMPI_SLICE_TECH"
esac

set -e
printf "\n"
ERR=0
for i in ${TESTS} ; do
    VAL="${i}_TIME"
    eval TIME="\$$VAL"
    printf "\033[1m$i (Time: %ss)\033[0m\n" "${TIME}"
    if [ ! -f "${RESULTS}/${i}/report" ]; then
        printf "\033[33mTEST FAILED (NO REPORT)\033[0m\n\n"
        ERR=$((ERR+1))
        continue
    fi

    if [ ! -s "${RESULTS}/${i}/report" ]; then
        printf "\033[33mTEST FAILED (EMPTY REPORT)\033[0m\n\n"
        ERR=$((ERR+1))
        continue
    fi
    STATUS="$(head -n 1 "${RESULTS}/${i}/report")"
    if [ "$STATUS" != "GOOD" ]; then
        printf "\033[31mTEST FAILED\033[0m\n"
        ERR=$((ERR+1))
    else
        printf "\033[32mTEST SUCCESS\033[0m\n"
    fi
    cat "${RESULTS}/${i}/report"
    printf "\n"
done

exit "$ERR"

