#!/bin/bash

BASE_DIR="~/clone/regression_test"
MMSEQSSSE="~/clone/build/src/mmseqs"
MMSEQSAVX="~/clone/build_avx2/src/mmseqs"

cd ${BASE_DIR}
# build the benchmark tools
git clone https://bitbucket.org/martin_steinegger/mmseqs-benchmark.git
cd mmseqs-benchmark
git submodule init
git submodule update
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release  ..
make -j 4 VERBOSE=1
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
${RUNEVAL} . ${MMSEQSSSE} ${EVALUATE} ${CI_COMMIT_ID} results 0 SSE_SEARCH 16 > report-${CI_COMMIT_ID}
${RUNEVAL} . ${MMSEQSAVX} ${EVALUATE} ${CI_COMMIT_ID} results 0 AVX2_SEARCH 16 >> report-${CI_COMMIT_ID}
${RUNEVAL} . ${MMSEQSSSE} ${EVALUATE} ${CI_COMMIT_ID} results 1 SSE_PROFILE 16 >> report-${CI_COMMIT_ID}
${RUNEVAL} . ${MMSEQSAVX} ${EVALUATE} ${CI_COMMIT_ID} results 1 AXX2_PROFILE 16 >> report-${CI_COMMIT_ID}

# fill out the report and fail
cat report-${CI_COMMIT_ID}
curl -F upfile=@report-${CI_COMMIT_ID} https://mmseqs.com/regression.php?secret=${REGRESSIONSECRET}
./mmseqs-benchmark/scripts/regression_report.sh report-${CI_COMMIT_ID} 0.236084 0.387362
exit $?