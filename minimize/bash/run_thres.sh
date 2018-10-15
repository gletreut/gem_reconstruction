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
#LOG=$(python -c "import os.path; b=\"$0\".replace(\".sh\",\".log\");print os.path.basename(b)")
#LOG=${RDIR}/$LOG
#exec 1<&-
#exec 2<&-
#exec 1<>$LOG
#exec 2>&1

# parameters
EXEC=prog                   # executable
argsstr="argsstr"
configstr="configstr"

# check existence of args and config file
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

if [ ! -f $configstr ]
then
  echo "missing config argument file $configstr"
  exit 0
fi

# get arguments
#ARGS=$(cat $argsstr)
CMAP=$(sed -n 1p $argsstr)
N=$(sed -n 2p $argsstr)
THRESMIN=$(sed -n 3p $argsstr)
THRESMAX=$(sed -n 4p $argsstr)
DTHRES=$(sed -n 5p $argsstr)
CONFIG=$(cat $configstr)

echo "CMAP: $CMAP"
echo "N: $N"
echo "THRESMIN: $THRESMIN"
echo "THRESMAX: $THRESMAX"
echo "DTHRES: $DTHRES"
echo "CONFIG: $CONFIG"

if [ -z $CMAP ] || [ -z $THRESMIN ] || [ -z $THRESMAX ] || [ -z $DTHRES ] || [ -z $N ]
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
time ./$EXEC $CMAP $N $THRESMIN $THRESMAX $DTHRES < $CONFIG
