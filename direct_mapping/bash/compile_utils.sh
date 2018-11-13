#!/bin/bash

# get source directory directories
RDIR=$1
if [ ! -d $RDIR ]
then
  echo Directory doesn\'t exist $RDIR
  exit
fi

# get main
SRC1=$2
SRC1PATH=$(python -c "import os.path; print os.path.join(\"$RDIR\",\"$SRC1\")")

if [ ! -f $SRC1PATH ]
then
  echo File doesn\'t exist $SRC1PATH
  exit
fi

# other source files
SRC2='utils_gsl.cpp sparse_matrix.cpp linalg.cpp gem_contacts_kij.cpp'

# method used for linear algebra
#KEY1=gsl
KEY1=lapack
#KEY1=mkl

FFAC=ffactor_gauss
EXEC=prog

# redirect output to log file
LOG=$(python -c "import os.path; b=\"$0\".replace(\".sh\",\".log\");print os.path.basename(b)")
LOG=${RDIR}/$LOG
exec 1<&-
exec 2<&-
exec 1<>$LOG
exec 2>&1

# start compilation
INCLUDE=$RDIR/include
if [ ! -d $INCLUDE ]
then
  echo "directory $INCLUDE does not exist!"
  exit
fi

# set compilation variables
CDEF=$(python -c "print \"-d${KEY1} -d${FFAC}\".upper()")
CCINC="-std=c++11 -O3 ${CDEF}"
OBJ=$(python -c "a=\"${SRC2}\".replace(\".cpp\",\".o\"); print a")

if [ ${KEY1} == gsl ]
then
  CFLAGS='-lgsl -lgslcblas -lm'
elif [ ${KEY1} == lapack ]
then
  CFLAGS='-lgsl -lsbtk -llapacke -lcblas -llapack -lblas -lgfortran -lm'
elif [ ${KEY1} == mkl ]
then
  source /opt/intel/mkl/bin/mklvars.sh intel64
  CFLAGS='-lgsl -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm'
else
  echo "wrong key: ${KEY1}"
  exit 0
fi

# compile source
for cc in ${SRC2}
do
  g++ ${CCINC} -c ${INCLUDE}/${cc} ${CFLAGS}
done

g++ ${CCINC} ${RDIR}/${SRC1} ${OBJ} ${CFLAGS} -o ${EXEC}
ls -lh
pwd

# clean
rm -f ${OBJ}
