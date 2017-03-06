# setup commands
cd ~/clone
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DHAVE_MPI=0 -DHAVE_SSE4_1=1 -DBUILD_SHARED_LIBS=OFF -DCMAKE_EXE_LINKER_FLAGS_RELEASE="-static -static-libgcc -static-libstdc++" -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" ..
make -j 4 VERBOSE=1
cd ..
mkdir build_avx2 && cd build_avx2
cmake -DCMAKE_BUILD_TYPE=Release -DHAVE_MPI=0 -DHAVE_AVX2=1 -DBUILD_SHARED_LIBS=OFF -DCMAKE_EXE_LINKER_FLAGS_RELEASE="-static -static-libgcc -static-libstdc++" -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" ..
make -j 4 VERBOSE=1
cd ..

# test pipeline
cd ~/clone
mkdir regression_test && cd regression_test
curl https://bitbucket.org/martin_steinegger/mmseqs-benchmark/raw/master/scripts/run_codeship_pipeline.sh > run_codeship_pipeline.sh
chmod +x run_codeship_pipeline.sh
./run_codeship_pipeline.sh

# deployment
cd ~/clone
mkdir -p ~/clone/mmseqs2.wiki
cd ~/clone/mmseqs2.wiki
curl https://raw.githubusercontent.com/wiki/soedinglab/mmseqs2/Home.md > Home.md
sed '1,/<!--- TOC END -->/d' Home.md > Home-body.md
echo -e "---\ntitle: MMseqs2 User Guide\nauthor:\n- Martin Steinegger\n- Milot Mirdita\n- Lars von den Driesch\n- Clovis Galiez\n- Johannes Söding\n..." > userguide-header.yaml
cat userguide-header.yaml Home-body.md > userguide.md
pandoc -D latex > pdflatex.template
sed -i 's/\[utf8\]/\[utf8x\]/g' pdflatex.template
pandoc userguide.md -o userguide.pdf --template=pdflatex.template --toc
chmod a+r userguide.pdf
scp userguide.pdf codeship@uniclust.mmseqs.com:/home/mirdita/repositories/mmseqs-webserver/latest
cd ~/clone/build
CURR_BUILD="mmseqs2"
mkdir -p ${CURR_BUILD}/bin
mkdir -p ${CURR_BUILD}/util
mkdir -p ${CURR_BUILD}
cp src/mmseqs ${CURR_BUILD}/bin
chmod +x ${CURR_BUILD}/bin/mmseqs
cp ../util/bash-completion.sh ${CURR_BUILD}/util
chmod +x ${CURR_BUILD}/util/bash-completion.sh
cp -r ../LICENCE.md ../README.md ~/clone/mmseqs2.wiki/userguide.pdf ../data ../examples ${CURR_BUILD}
chmod -R g-w,o-w ${CURR_BUILD}
tar czvf mmseqs-static_sse41.tar.gz  ${CURR_BUILD}
scp mmseqs-static_sse41.tar.gz codeship@uniclust.mmseqs.com:/home/mirdita/repositories/mmseqs-webserver/latest
cp ~/clone/build_avx2/src/mmseqs ${CURR_BUILD}/bin
chmod +x ${CURR_BUILD}/bin/mmseqs
tar czvf mmseqs-static_avx2.tar.gz ${CURR_BUILD}
scp mmseqs-static_avx2.tar.gz codeship@uniclust.mmseqs.com:/home/mirdita/repositories/mmseqs-webserver/latest
