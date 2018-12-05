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

cp mmseqs-benchmark/data/{q,t}set_{01,02}.fas.gz small-benchmark-db

# go run it
RUNEVAL="./mmseqs-benchmark/scripts/run_mmseqs_regression.sh"
EVALUATE="./mmseqs-benchmark/build/evaluate_results"
time ${RUNEVAL} . ${MMSEQSAVX} ${EVALUATE} ${CI_COMMIT_ID} results 0 AVX2_SEARCH 16 >> report-${CI_COMMIT_ID}
time ${RUNEVAL} . ${MMSEQSAVX} ${EVALUATE} ${CI_COMMIT_ID} results 1 AXX2_PROFILE 16 >> report-${CI_COMMIT_ID}
RUNEVAL="./mmseqs-benchmark/scripts/run_mmseqs_easy_regression.sh"
time ${RUNEVAL} . ${MMSEQSSSE} ${EVALUATE} ${CI_COMMIT_ID} results 0 EASY_SSE_SEARCH 16 > report-${CI_COMMIT_ID}
time ${RUNEVAL} . ${MMSEQSSSE} ${EVALUATE} ${CI_COMMIT_ID} results 1 EASY_SSE_PROFILE 16 >> report-${CI_COMMIT_ID}

${MMSEQSAVX} translateaa small-benchmark-db/query small-benchmark-db/nucl
( cd small-benchmark-db && ln -sf query_h nucl_h && ln -sf query_h.index nucl_h.index )
${MMSEQSAVX} translateaa small-benchmark-db/db2 small-benchmark-db/db2nucl
time ${RUNEVAL} . ${MMSEQSAVX} ${EVALUATE} ${CI_COMMIT_ID} results NUCL_SEARCH 16 >> report-${CI_COMMIT_ID}
RUNEVAL="./mmseqs-benchmark/scripts/run_mmseqs_sliceprofile_regression.sh"
time ${RUNEVAL} . ${MMSEQSAVX} ${EVALUATE} ${CI_COMMIT_ID} sliceprofile-results SLICEPROFILE 16 >> report-${CI_COMMIT_ID}
RUNEVAL="./mmseqs-benchmark/scripts/run_mmseqs_dbprofile_regression.sh"
time ${RUNEVAL} . ${MMSEQSAVX} ${EVALUATE} ${CI_COMMIT_ID} dbprofile-results DBPROFILE 16 >> report-${CI_COMMIT_ID}
#RUNEVAL="./mmseqs-benchmark/scripts/run_mmseqs_profilestates_regression.sh"
#time ${RUNEVAL} . ${MMSEQSAVX} ${EVALUATE} ${CI_COMMIT_ID} csprofile-results CSPROFILE 16 >> report-${CI_COMMIT_ID}

cp mmseqs-benchmark/data/clu.fasta mmseqs-benchmark/data/clu-tcov.fasta.gz small-benchmark-db
RUNEVAL="./mmseqs-benchmark/scripts/run_mmseqs_clu_regression.sh"
time ${RUNEVAL} small-benchmark-db/clu.fasta ${MMSEQSAVX} CLU ${CI_COMMIT_ID} clu-results 0 "--min-seq-id 0.3 -s 2 --cluster-steps 3 --threads 16" \
    >> report-${CI_COMMIT_ID}
time ${RUNEVAL} "small-benchmark-db/clu.fasta" ${MMSEQSAVX} LINCLU ${CI_COMMIT_ID} linclu-results 1 "--cov-mode 1 --cluster-mode 0 -c 0.90 --min-seq-id 0.50 --threads 16" \
    >> report-${CI_COMMIT_ID}

RUNEVAL="./mmseqs-benchmark/scripts/run_mmseqs_easy_clu_regression.sh"
time ${RUNEVAL} small-benchmark-db/clu.fasta ${MMSEQSAVX} CLU ${CI_COMMIT_ID} clu-results 0 "--min-seq-id 0.3 -s 2 --cluster-steps 3 --threads 16" \
    >> report-${CI_COMMIT_ID}
time ${RUNEVAL} "small-benchmark-db/clu.fasta" ${MMSEQSAVX} LINCLU ${CI_COMMIT_ID} linclu-results 1 "--cov-mode 1 --cluster-mode 0 -c 0.90 --min-seq-id 0.50 --threads 16" \
    >> report-${CI_COMMIT_ID}

RUNEVAL="./mmseqs-benchmark/scripts/run_mmseqs_multihit_regression.sh"
time ${RUNEVAL} . ${MMSEQSAVX} ${CI_COMMIT_ID} multihit-results MULTHIT 16 >> report-${CI_COMMIT_ID}

RUNEVAL="./mmseqs-benchmark/scripts/run_mmseqs_extractorfs.sh"
time ${RUNEVAL} "./mmseqs-benchmark" ${MMSEQSAVX} "./mmseqs-benchmark/scripts" ${CI_COMMIT_ID} results EXTRACTORFS 16 >> report-${CI_COMMIT_ID}

( cd small-benchmark-db && ln -sf db2_h db2nucl_h && ln -sf db2_h.index db2nucl_h.index && ln -sf db2_h.dbtype db2nucl_h.dbtype )
RUNEVAL="./mmseqs-benchmark/scripts/run_mmseqs_nucl_nucl_regression.sh"
time ${RUNEVAL} . ${MMSEQSAVX} ${EVALUATE} ${CI_COMMIT_ID} results NUCLNUCL_SEARCH 16 >> report-${CI_COMMIT_ID}

# fill out the report and fail
cat report-${CI_COMMIT_ID}
#curl -F upfile=@report-${CI_COMMIT_ID} https://mmseqs.com/regression.php?secret=${REGRESSIONSECRET}
./mmseqs-benchmark/scripts/regression_report.sh report-${CI_COMMIT_ID} 0.235 0.334 0.235 0.142 0.140 0.245 17299 26505 1.112E-202 4.032E-142 0 0.177 
exit $?
