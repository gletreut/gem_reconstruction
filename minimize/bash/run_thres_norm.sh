#!/bin/bash

# get directories
if [ ! -z $1 ]
then
	RDIR=$(python -c "import os.path; print os.path.realpath(\"$1\")")
else
	RDIR=$(pwd)
fi
cd $RDIR

# redirect output to log file
LOG=$(python -c "import os.path; b=\"$0\".replace(\".sh\",\".log\");print os.path.basename(b)")
LOG=${RDIR}/$LOG
exec 1<&-
exec 2<&-
exec 1<>$LOG
exec 2>&1

# parameters
EXEC=prog                   # executable
argsstr="argsstr"
datastr="datastr"

# check existence of args and data file
if [ ! -e $EXEC ]
then
  echo "missing executable file $EXEC"
  exit 0
fi

if [ ! -f $argsstr ]
then
  echo "missing arguments argument file $argsstr"
  exit 0
fi

if [ ! -f $datastr ]
then
  echo "missing data argument file $datastr"
  exit 0
fi

# get arguments
#ARGS=$(cat $argsstr)
CMAP=$(sed -n 1p $argsstr)
Z=$(sed -n 2p $argsstr)
N=$(sed -n 3p $argsstr)
THRESMIN=$(sed -n 4p $argsstr)
THRESMAX=$(sed -n 5p $argsstr)
DTHRES=$(sed -n 6p $argsstr)
DATA=$(cat $datastr)

echo "CMAP: $CMAP"
echo "Z: $Z"
echo "N: $N"
echo "THRESMIN: $THRESMIN"
echo "THRESMAX: $THRESMAX"
echo "DTHRES: $DTHRES"
echo "DATA: $DATA"

if [ -z $CMAP ] || [ -z $Z ] || [ -z $THRESMIN ] || [ -z $THRESMAX ] [ -z $DTHRES ] || [ -z $N ]
then
  echo invalid arguments in arguments file: $ARGS
  exit 1
fi

if [ ! -f $CMAP ]
then
	echo file does not exist: ${CMAP}!
fi

# source mkl variables if necessary
if [ -z ${MKLROOT} ]
then
  source /opt/intel/mkl/bin/mklvars.sh intel64
fi

# execute program
time ./$EXEC $CMAP $Z $N $THRESMIN $THRESMAX $DTHRES < $DATA
