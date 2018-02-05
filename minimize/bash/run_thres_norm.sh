#!/bin/bash

if [ -z $1 ]
then
  echo you must give argments
  exit 0
fi

# get arguments
EXEC=$(python -c "import os.path; print os.path.realpath(\"$1\")")
N=$2
THRESMIN=$3
THRESMAX=$4
DTHRES=$5
CMAP=$6
Z=$7
CONFIG=$8

## redirect output to log file
#RDIR=$(python -c "import os.path; b=\"$0\".replace(\".sh\",\".log\");print os.path.basename(b)")
#LOG=$(python -c "import os.path; b=\"$0\".replace(\".sh\",\".log\");print os.path.basename(b)")
#LOG=${RDIR}/$LOG
#exec 1<&-
#exec 2<&-
#exec 1<>$LOG
#exec 2>&1

if [ -z $CMAP ] || [ -z $Z ] || [ -z $THRESMIN ] || [ -z $THRESMAX ] || [ -z $DTHRES ] || [ -z $N ]
then
  echo invalid arguments in arguments!
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
$EXEC $CMAP $Z $N $THRESMIN $THRESMAX $DTHRES < $CONFIG
