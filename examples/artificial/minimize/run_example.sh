#!/bin/bash

# parameters
CMAT='cmat_bd_nconf100.txt' # file containing contact probabilities computed using 100 samples configurations of a GEM simulated by Brownian Dynamics.
N=99  # size of the square matrix imported from the CMAP file. Here restrict to the 100x100 upper left square.
THRES=1.5 # threshold used when using the GEM mapping between contact proba and couplings.
CONFIG='config.txt' # configuration file for the minimization algorithm.
YAML1='make_matrix_plots_cmaps.yml' # configuration file for plot
YAML2='make_matrix_plots_kmaps.yml' # configuration file for plot

# compile code
bash ../../../minimize/bash/compile.sh ../../../minimize

# execute
#bash ../../../minimize/bash/run.sh prog ${CMAT} ${N} ${THRES} ${CONFIG}
#or
./prog ${CMAT} ${N} ${THRES} < ${CONFIG}
rm prog

# plot result
python ../../../plots/make_matrix_plots.py ${YAML1}  .

# plot reconstructed couplings versus original couplings
python ../../../plots/make_matrix_plots.py ${YAML2}  .
