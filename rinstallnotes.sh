#downlaod source
mkdir ~/Downloads
mkdir ~/packages

# install open ssl
cd ~/Downloads
 git clone https://github.com/openssl/openssl.git
 cd ~/Downloads/openssl
 ./config --prefix=$HOME/packages
 make
 make install

# install zlib
cd ~/Downloads
 wget http://zlib.net/zlib-1.2.11.tar.gz
 tar -xzvf zlib-1.2.11.tar.gz
 cd ~/Downloads/zlib-1.2.11
 ./config --prefix=$HOME/packages
 make
 make install

# set flags for future packages
export PATH=$HOME/packages/bin:$PATH
export LD_LIBRARY_PATH=$HOME/packages/lib:/usr/lib64:$LD_LIBRARY_PATH 
export CFLAGS="-I$HOME/packages/include" 
export LDFLAGS="-L$HOME/packages/lib"

# install bzip
cd ~/Downloads
 wget http://www.bzip.org/1.0.6/bzip2-1.0.6.tar.gz
 tar -xzvf bzip2-1.0.6.tar.gz
 cd ~/Downloads/bzip2-1.0.6
 make -f Makefile-libbz2_so
 make clean
 make CFLAGS='-fPIC'
 make -n install PREFIX=$HOME/packages
 make install PREFIX=$HOME/packages

# install xzvf
cd ~/Downloads
 wget http://tukaani.org/xz/xz-5.2.3.tar.gz
 tar -xzvf xz-5.2.3.tar.gz
 cd ~/Downloads/xz-5.2.3
 ./configure --prefix=$HOME/packages
 make -j3
 make install

# install pcre
cd ~/Downloads
 wget ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.39.tar.gz
 tar -xzvf pcre-8.39.tar.gz
 cd ~/Downloads/pcre-8.39
 ./configure --enable-utf8 --prefix=$HOME/packages
 make -j3
 make install

# install curl
cd ~/Downloads
 wget --no-check-certificate https://curl.haxx.se/download/curl-7.54.0.tar.gz
 tar -xzvf curl-7.54.0.tar.gz
 cd ~/Downloads/curl-7.54.0
 CPPFLAGS="-I$HOME/packages/include" LDFLAGS="-L$HOME/packages/lib" ./configure --prefix=$HOME/packages
 make -j3
 make install

# install R
cd ~/Downloads
 wget http://cran.cnr.berkeley.edu/src/base/R-3/R-3.4.0.tar.gz
 tar -xvzf R-3.4.0.tar.gz
 cd ~/Downloads/R-3.4.0
 mkdir builddir
 cd builddir/

 ../configure --prefix=$HOME/packages/R-3.4.0 '--with-cairo' \
  '--with-jpeglib' '--with-readline' '--with-tcltk' \
  '--with-blas' '--with-lapack' '--enable-R-profiling' \
  '--enable-R-shlib' \
  '--enable-memory-profiling'

 make
 make install
