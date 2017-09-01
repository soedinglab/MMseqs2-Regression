#!/bin/bash -ex
BASE_DIR="$HOME/clone/regression_test"
MMSEQSSSE="$HOME/clone/build/src/mmseqs"
MMSEQSAVX="$HOME/clone/build_avx2/src/mmseqs"

cd ${BASE_DIR}
# build the benchmark tools
git clone https://bitbucket.org/martin_steinegger/mmseqs-benchmark.git
cd mmseqs-benchmark
git submodule init
git submodule update
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release  ..
make -j 4 VERBOSE=0
cd ..
cd ..

#setup benchmark database
mkdir small-benchmark-db
${MMSEQSSSE} createdb mmseqs-benchmark/data/query.fasta small-benchmark-db/query
${MMSEQSSSE} createdb mmseqs-benchmark/data/targetannotation.fasta small-benchmark-db/db2
cp mmseqs-benchmark/data/query.fasta mmseqs-benchmark/data/targetannotation.fasta small-benchmark-db

# go run it
RUNEVAL="./mmseqs-benchmark/scripts/run_mmseqs_regression.sh"
EVALUATE="./mmseqs-benchmark/build/evaluate_results"
time ${RUNEVAL} . ${MMSEQSSSE} ${EVALUATE} ${CI_COMMIT_ID} results 0 SSE_SEARCH 16 > report-${CI_COMMIT_ID}
time ${RUNEVAL} . ${MMSEQSAVX} ${EVALUATE} ${CI_COMMIT_ID} results 0 AVX2_SEARCH 16 >> report-${CI_COMMIT_ID}
time ${RUNEVAL} . ${MMSEQSSSE} ${EVALUATE} ${CI_COMMIT_ID} results 1 SSE_PROFILE 16 >> report-${CI_COMMIT_ID}
time ${RUNEVAL} . ${MMSEQSAVX} ${EVALUATE} ${CI_COMMIT_ID} results 1 AXX2_PROFILE 16 >> report-${CI_COMMIT_ID}
RUNEVAL="./mmseqs-benchmark/scripts/run_mmseqs_dbprofile_regression.sh"
time ${RUNEVAL} . ${MMSEQSAVX} ${EVALUATE} ${CI_COMMIT_ID} dbprofile-results DBPROFILE 16 >> report-${CI_COMMIT_ID}

cp mmseqs-benchmark/data/clu.fasta mmseqs-benchmark/data/clu-tcov.fasta.gz small-benchmark-db
RUNEVAL="./mmseqs-benchmark/scripts/run_mmseqs_clu_regression.sh"
time ${RUNEVAL} small-benchmark-db/clu.fasta ${MMSEQSAVX} CLU ${CI_COMMIT_ID} clu-results 0 "--cascaded --min-seq-id 0.3 -s 2 --threads 16" \
    >> report-${CI_COMMIT_ID}
time ${RUNEVAL} "small-benchmark-db/query.fasta small-benchmark-db/clu.fasta" ${MMSEQSAVX} LINCLU ${CI_COMMIT_ID} linclu-results 1 "--cov-mode 1 -c 0.90 --min-seq-id 0.50 --threads 16" \
    >> report-${CI_COMMIT_ID}

# fill out the report and fail
cat report-${CI_COMMIT_ID}
curl -F upfile=@report-${CI_COMMIT_ID} https://mmseqs.com/regression.php?secret=${REGRESSIONSECRET}
./mmseqs-benchmark/scripts/regression_report.sh report-${CI_COMMIT_ID} 0.235 0.331 0.22 17670 26896
exit $?
