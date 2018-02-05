#!/bin/bash
apt-get update
apt-get install -y build-essential checkinstall g++ gfortran python curl git sharutils cpio g++-multilib less

PREFPATH=/usr/local
LIBPATH=$PREFPATH/lib
INCPATH=$PREFPATH/include
proc=3

# install gsl
cd
FDOWN=gsl-latest.tar.gz
curl -O  http://ftp.igh.cnrs.fr/pub/gnu/gsl/$FDOWN
tar -zxf $FDOWN
cd gsl*/
./configure --prefix=$PREFPATH
make all
make install
cd
rm -rf gsl*/
rm -f $FDOWN

# Install blas
FDOWN=blas-3.7.0.tgz
cd
curl -O http://www.netlib.org/blas/$FDOWN
tar -zxf $FDOWN
cd BLAS*/
make -j$proc all
cp blas_LINUX.a $LIBPATH/libblas.a
ranlib $LIBPATH/libblas.a
cd
rm -rf BLAS*
rm -f $FDOWN

# Install cblas
cd
FDOWN=cblas.tgz
curl -O http://www.netlib.org/blas/blast-forum/$FDOWN
tar -zxf $FDOWN
cd CBLAS*/
cp Makefile.LINUX Makefile.in
ln -s $LIBPATH/libblas.a src/.
ln -s $LIBPATH/libblas.a testing/.
make all
cp lib/cblas_LINUX.a $LIBPATH/libcblas.a
ranlib $LIBPATH/libcblas.a
cp include/*.h $INCPATH
cd
rm -rf CBLAS*
rm -f $FDOWN

# Install Sparse BLAS form NIST
FDOWN=spblastk.shar
cd
curl -O ftp://gams.nist.gov/pub/karin/spblas/${FDOWN}.Z
gzip -d ${FDOWN}.Z
unshar $FDOWN
cd spblas*/
make -j$proc library
cp lib/libsbtk.a $LIBPATH
ranlib $LIBPATH/libsbtk.a
cp include/*.h $INCPATH
cd
rm -f $FDOWN
rm -rf spblas*/

# Insall lapack/lapacke
FDOWN=lapack
cd
git clone https://github.com/Reference-LAPACK/lapack.git $FDOWN
cd $FDOWN

## make LAPACK
cp make.inc.example make.inc
ln -s $LIBPATH/libblas.a librefblas.a
ln -s $LIBPATH/libcblas.a libcblas.a
ulimit -s 65000
# see explanations at http://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=2&t=1596
make -j${proc} all
cp liblapack.a libtmglib.a $LIBPATH
ranlib $LIBPATH/liblapack.a
ranlib $LIBPATH/libtmglib.a

## make LAPACKE
make -j$proc lapackelib
cp LAPACKE/include/*.h $INCPATH/.
cp liblapacke.a $LIBPATH
ranlib $LIBPATH/liblapacke.a
cd
rm -rf $FDOWN

# install mkl
FDOWN=l_mkl_2017.2.174
curl -O http://registrationcenter-download.intel.com/akdlm/irc_nas/tec/11306/${FDOWN}.tgz
tar -zxf ${FDOWN}.tgz
cp silent_online.cfg ${FDOWN}/silent.cfg
cd ${FDOWN}
./install.sh --silent ./silent.cfg
cd
rm -rf ${FDOWN} ${FDOWN}.tgz

## remove unnecessary packages
apt-get purge -y git curl sharutils
