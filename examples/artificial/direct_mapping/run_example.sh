#!/bin/bash

# parameters
CMAT='cmat_bd_nconf10000_thres2.0.txt' # file containing contact probabilities computed using 10000 samples configurations of a GEM simulated by Brownian Dynamics.
N=99  # size of the square matrix imported from the CMAP file. Here restrict to the 100x100 upper left square.
THRES=2.0 # threshold used when using the GEM mapping between contact proba and couplings.
CMAPOUT=cmat_out.dat
KMAPOUT=kmat_out.dat
YAML1='make_matrix_plots_cmaps.yml' # configuration file for plot
YAML2='make_matrix_plots_kmaps.yml' # configuration file for plot
argsstr="argsstr"

# compile code
bash ../../../direct_mapping/bash/compile.sh

# create argument files
## ARGUMENT FILES
echo $CMAT > $argsstr
echo $N >> $argsstr
echo $THRES >> $argsstr
echo $CMAPOUT >> $argsstr
echo $KMAPOUT >> $argsstr

# execute
bash ../../../direct_mapping/bash/run.sh
rm -f prog

# plot result
python ../../../plots/make_matrix_plots.py ${YAML1}  .

# plot reconstructed couplings versus original couplings
python ../../../plots/make_matrix_plots.py ${YAML2}  .

