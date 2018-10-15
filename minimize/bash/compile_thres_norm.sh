#!/bin/bash

# source files
SRC1='minimize_thres_norm.cpp'

# get source directory directories
BDIR=$(python -c "import os.path; print os.path.dirname(os.path.realpath(\"$0\"))")
RDIR=$(python -c "import os.path; print os.path.realpath(os.path.join(\"$BDIR\",\"..\"))")

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

# compile
bash $COMPILEUTILS $RDIR ${SRC1}
