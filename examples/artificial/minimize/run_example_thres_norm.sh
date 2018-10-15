#!/bin/bash

# parameters
CMAT='cmat_bd_nconf100_thres1.5.txt' # file containing contact probabilities computed using 100 samples configurations of a GEM simulated by Brownian Dynamics.
N=99  # size of the square matrix imported from the CMAP file. Here restrict to the 100x100 upper left square.
THRESMIN=1.0 # minimum threshold used when using the GEM mapping between contact proba and couplings.
THRESMAX=2.0 # maximum threshold used when using the GEM mapping between contact proba and couplings.
DTHRES=0.1 # threshold increment
ZNORM=2.0
CONFIG='config.txt' # configuration file for the minimization algorithm.
YAML1='make_matrix_plots_cmaps.yml' # configuration file for plot
YAML2='make_matrix_plots_kmaps.yml' # configuration file for plot
argsstr="argsstr"
configstr="configstr"

# compile code
bash ../../../minimize/bash/compile_thres_norm.sh

# create argument files
## ARGUMENT FILES
echo $CMAT > $argsstr
echo $N >> $argsstr
echo $THRESMIN >> $argsstr
echo $THRESMAX >> $argsstr
echo $DTHRES >> $argsstr
echo $ZNORM >> $argsstr
echo ${CONFIG} > $configstr

# execute
bash ../../../minimize/bash/run_thres_norm.sh
rm -f prog

# plot result
python ../../../plots/make_matrix_plots.py ${YAML1}  .

# plot reconstructed couplings versus original couplings
python ../../../plots/make_matrix_plots.py ${YAML2}  .

# LSE versus threshold value
python ../../../plots/make_distance_plots.py make_distance_plots.yml distances.dat

