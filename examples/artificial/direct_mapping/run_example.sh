#!/bin/bash

# parameters
CMAT='cmat_bd_nconf10000_thres2.0.txt' # file containing contact probabilities computed using 10000 samples configurations of a GEM simulated by Brownian Dynamics.
N=99  # size of the square matrix imported from the CMAP file. Here restrict to the 100x100 upper left square.
THRES=2.0 # threshold used when using the GEM mapping between contact proba and couplings.
EXEC=prog
YAML1='make_matrix_plots_cmaps.yml' # configuration file for plot
YAML2='make_matrix_plots_kmaps.yml' # configuration file for plot

# compile code
bash ../../../direct_mapping/bash/compile.sh ../../../direct_mapping

# execute

# execute
#or
bash ../../../direct_mapping/bash/run.sh prog ${N} ${THRES} ${CMAT} cmat_out.dat kmat_out.dat
rm prog

# plot result
python ../../../plots/make_matrix_plots.py ${YAML1}  .

# plot reconstructed couplings versus original couplings
python ../../../plots/make_matrix_plots.py ${YAML2}  .
