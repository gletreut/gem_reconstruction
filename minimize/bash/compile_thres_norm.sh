#!/bin/bash

# get target directory
TDIR=$(python -c "import os.path; print os.path.realpath(\"$1\")")
if [ ! -d $TDIR ]
then
  echo Directory doesn\'t exist $TDIR
  exit
fi
cd $TDIR

# source files
SRC1='minimize_thres_norm.cpp'

# get source directory directories
BDIR=$(python -c "import os.path; print os.path.dirname(os.path.realpath(\"$0\"))")
RDIR=$(python -c "import os.path; print os.path.realpath(os.path.join(\"$BDIR\",\"..\"))")

# redirect output to log file
LOG=$(python -c "import os.path; b=\"$0\".replace(\".sh\",\".log\");print os.path.basename(b)")
LOG=${RDIR}/$LOG
exec 1<&-
exec 2<&-
exec 1<>$LOG
exec 2>&1

if [ ! -d $BDIR ]
then
  echo Directory doesn\'t exist $BDIR
  exit
fi

if [ ! -d $RDIR ]
then
  echo Directory doesn\'t exist $RDIR
  exit
fi

COMPILEUTILS=$(python -c "import os.path; print os.path.join(\"$BDIR\",\"compile_utils.sh\")")

# compile
bash $COMPILEUTILS $RDIR ${SRC1}

echo "finished compilation"
