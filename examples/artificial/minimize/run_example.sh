#!/bin/bash

# parameters
CMAT='cmat_bd_nconf100_thres1.5.txt' # file containing contact probabilities computed using 100 samples configurations of a GEM simulated by Brownian Dynamics.
N=99  # size of the square matrix imported from the CMAP file. Here restrict to the 100x100 upper left square.
THRES=1.5 # threshold used when using the GEM mapping between contact proba and couplings.
CONFIG='config.txt' # configuration file for the minimization algorithm.
YAML1='make_matrix_plots_cmaps.yml' # configuration file for plot
YAML2='make_matrix_plots_kmaps.yml' # configuration file for plot
argsstr="argsstr"
configstr="configstr"

# compile code
bash ../../../minimize/bash/compile.sh ../../../minimize

# create argument files
## ARGUMENT FILES
echo $N > $argsstr
echo $THRES >> $argsstr
echo $CMAT >> $argsstr
echo $THMAX >> $argsstr
echo ${CONFIG} > $configstr

# execute
#bash ../../../minimize/bash/run.sh prog ${CMAT} ${N} ${THRES} ${CONFIG}
bash ../../../minimize/bash/run.sh
rm -f prog

# plot result
python ../../../plots/make_matrix_plots.py ${YAML1}  .

# plot reconstructed couplings versus original couplings
python ../../../plots/make_matrix_plots.py ${YAML2}  .
