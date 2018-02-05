#!/bin/bash

# parameters
CMAT='cmat_exp.txt' # file containing the experimental contact probabilities.
N=100  # size of the square matrix imported from the CMAP file. Here restrict to the 100x100 upper left square.
THRES=1.0 # threshold used when using the GEM mapping between contact proba and couplings.
CONFIG='config.txt' # configuration file for the minimization algorithm.
YAML='make_matrix_plots_cmaps.yml' # configuration file for plot

# compile code
bash ../../minimize/bash/compile.sh ../../minimize

# execute
#bash ../../minimize/bash/run.sh prog ${CMAT} ${N} ${THRES} ${CONFIG}
#or
./prog ${CMAT} ${N} ${THRES} < ${CONFIG}
rm prog

# plot result
python ../../plots/make_matrix_plots.py ${YAML}  .
