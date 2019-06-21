#!/bin/sh -e
git clone https://github.com/samtools/htslib.git
git clone https://github.com/samtools/samtools.git
cd samtools/
autoheader
autoconf -Wno-error
./configure --enable-configure-htslib=yes --without-curses --disable-lzma --disable-bz2 CFLAGS="-Os"
make -j

if [ -e samtools.exe ]; then
    strip --strip-unneeded samtools.exe
    upx samtools.exe
    cp /usr/bin/cygz.dll .
    strip --strip-unneeded cygz.dll
    upx cygz.dll
    cp /usr/bin/cygwin1.dll .
    strip --strip-unneeded cygwin1.dll
    upx cygwin1.dll --force
else
    strip samtools
    upx samtools
fi

