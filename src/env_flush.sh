#!/bin/bash

module purge
module load jet
module load omfit
module load intel-compilers/2020.3.110
export PYTHONPATH=/usr/local/depot/flush-2.4.2/lib/python3:$PYTHONPATH
export LD_LIBRARY_PATH=/usr/local/depot/flush-2.4.2/lib/Linux_64_gfortran:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/local/depot/NAGMark24/fll6i24dc/lib:$LD_LIBRARY_PATH
