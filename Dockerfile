FROM alpine:latest as benchmark-builder

RUN apk add --no-cache gcc g++ cmake musl-dev vim git ninja

WORKDIR /opt/mmseqs-benchmark
ADD . .

RUN git submodule init && git submodule update

WORKDIR build
RUN cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..
RUN ninja && ninja install

FROM mmseqs2 AS mmseqs
MAINTAINER Milot Mirdita <milot@mirdita.de>

RUN apk add --no-cache coreutils gawk

COPY --from=benchmark-builder /opt/mmseqs-benchmark/build/bin/evaluate_results /usr/local/bin/evaluate_results
COPY --from=benchmark-builder /opt/mmseqs-benchmark/build/bin/evaluate_pos_results /usr/local/bin/evaluate_pos_results
COPY --from=benchmark-builder /opt/mmseqs-benchmark/build/bin/test_main /usr/local/bin/test_main

WORKDIR small-benchmark-db
ADD data/query.fasta query.fasta
ADD data/targetannotation.fasta targetannotation.fasta
ADD data/clu.fasta clu.fasta

RUN mmseqs createdb query.fasta query
RUN mmseqs createdb targetannotation.fasta db2

WORKDIR ..

ADD scripts/run_mmseqs_regression.sh /usr/local/bin/run_mmseqs_regression.sh
ADD scripts/run_mmseqs_dbprofile_regression.sh /usr/local/bin/run_mmseqs_dbprofile_regression.sh
ADD scripts/run_mmseqs_clu_regression.sh /usr/local/bin/run_mmseqs_clu_regression.sh
ADD scripts/regression_report.sh /usr/local/bin/regression_report.sh

RUN time run_mmseqs_regression.sh . mmseqs evaluate_results $(mmseqs | awk '/^MMseqs2? Version:/ {print $3}') results 0 SEARCH $(nproc --all) > report
RUN time run_mmseqs_regression.sh . mmseqs evaluate_results $(mmseqs | awk '/^MMseqs2? Version:/ {print $3}') results 1 PROFILE $(nproc --all) >> report

RUN time run_mmseqs_dbprofile_regression.sh . mmseqs evaluate_results $(mmseqs | awk '/^MMseqs2? Version:/ {print $3}') dbprof-results DBPROFILE $(nproc --all) >> report

RUN time run_mmseqs_clu_regression.sh small-benchmark-db/clu.fasta mmseqs CLU $(mmseqs | awk '/^MMseqs2? Version:/ {print $3}') clu-results 0 "--min-seq-id 0.3 -s 2 --cluster-steps 3 --threads $(nproc --all)" >> report
RUN time run_mmseqs_clu_regression.sh "small-benchmark-db/query.fasta small-benchmark-db/clu.fasta" mmseqs LINCLU $(mmseqs | awk '/^MMseqs2? Version:/ {print $3}') linclu-results 1 "--cov-mode 1 --cluster-mode 0 -c 0.90 --min-seq-id 0.50 --threads $(nproc --all)" >> report

RUN cat report
RUN regression_report.sh report 0.235 0.331 0.22 17285 26819
