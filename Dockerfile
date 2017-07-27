FROM alpine:latest as benchmark-builder

RUN apk add --no-cache gcc g++ cmake musl-dev vim git ninja

WORKDIR /opt/mmseqs-benchmark
ADD . .

RUN git submodule init && git submodule update

WORKDIR build
RUN cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..
RUN ninja && ninja install

FROM soedinglab/mmseqs2 AS mmseqs
MAINTAINER Milot Mirdita <milot@mirdita.de>

RUN apk add --no-cache coreutils

COPY --from=benchmark-builder /opt/mmseqs-benchmark/build/bin/evaluate_results /usr/local/bin/evaluate_results
COPY --from=benchmark-builder /opt/mmseqs-benchmark/build/bin/evaluate_pos_results /usr/local/bin/evaluate_pos_results
COPY --from=benchmark-builder /opt/mmseqs-benchmark/build/bin/test_main /usr/local/bin/test_main

WORKDIR small-benchmark-db
ADD data/query.fasta query.fasta
ADD data/targetannotation.fasta targetannotation.fasta
RUN mmseqs createdb query.fasta query
RUN mmseqs createdb targetannotation.fasta db2

WORKDIR ..

ADD scripts/run_mmseqs_regression.sh /usr/local/bin/run_mmseqs_regression.sh
RUN time run_mmseqs_regression.sh . mmseqs evaluate_results $(mmseqs | awk '/^MMseqs Version:/ {print $3}') results 0 SEARCH $(nproc --all)
RUN time run_mmseqs_regression.sh . mmseqs evaluate_results $(mmseqs | awk '/^MMseqs Version:/ {print $3}') results 1 PROFILE $(nproc --all)

