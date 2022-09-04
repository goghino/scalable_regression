#!/bin/bash

# USER CONFIGURATION
RUN_MATLAB=1 #0 or 1
LOAD_FILES=1 #0 or 1, effective only if RUN_MATLAB=0
SOLVER=$FMINCON #IPOPT or FMINCON, effective only if RUN_MATLAB=1

## load list of benchmarks to run
. config.sh

######################################
## specify tool for memory tracking
MEMORY_WRAPPER="/usr/bin/time -v"
TIMEOUT_WRAPPER="timeout --signal=KILL --verbose --foreground 5h"
MATLAB="matlab"

export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8
########################################
if [[ $RUN_MATLAB -eq 1 ]]
then
    LOAD_FILES=1
fi


if [[ $LOAD_FILES -eq 1 ]]
then

   rm -f optimize_espa
   ln -s optimize_espa_file optimize_espa
    
    Nx=${#gridsX[@]}
    Nb=${#gridsB[@]}
    echo "Using $Nx Xbeta and $Nb beta_eLR files..."
############################################
else
   rm -f optimize_espa
   ln -s optimize_espa_generator optimize_espa
    
    declare -a gridsX=(1024 2048 4096 8192 16384 32768 65536 131072 262144 524288 1048576)
    
    Nx=${#gridsX[@]}
    echo "Using $Nx sizes for the benchmarks"
fi
############################################
echo "Running $Nx benchmarks..."

## now loop through the cases and solvers
#for OPFcase in "${grids[@]}"; do
for i in $(seq 0 $((Nx-1))); do

#========================================================
#           C++ 
#========================================================
if [[ $RUN_MATLAB -ne 1 ]]
then
    echo "RUNNING C++ Code"

    if [[ $LOAD_FILES -eq 1 ]]
    then
    echo "Running benchmark with loading the data files:"
    echo "data/"${gridsX[$i]}
    echo "data/"${gridsB[$i]}
    # Execute C++ code: ./optimize_espa data/Xbeta data/beta_eLR file.out
    ${MEMORY_WRAPPER} ${TIMEOUT_WRAPPER} ./optimize_espa "data/"${gridsX[$i]} "data/"${gridsB[$i]} ${gridsX[$i]}".sol" 2>&1 | tee ${gridsX[$i]}".out"
    else
    echo "Running benchmark with generating data files on the fly:"
    echo "dim="${gridsX[$i]}
    # Execute C++ code: ./optimize_espa dim T reg_param file.out
    ${MEMORY_WRAPPER} ${TIMEOUT_WRAPPER} ./optimize_espa ${gridsX[$i]} 1000 1e-3 ${gridsX[$i]}".sol" 2>&1 | tee ${gridsX[$i]}".out"
    fi
fi

#========================================================
#            MATLAB
#========================================================
if [[ $RUN_MATLAB -eq 1 ]]
then
    # Execute Matlab code: matlab -r test_bench
    echo "Running MATLAB benchmark with loading the data files:"
    echo "data/"${gridsX[$i]}
    echo "data/"${gridsB[$i]}
    ${MEMORY_WRAPPER} ${MATLAB} -singleCompThread -nodisplay -nosplash -nodesktop -r "try test_bench('"data/"${gridsX[$i]}', '"data/"${gridsB[$i]}', $SOLVER); catch; end; quit" 2>&1 |  tee ${gridsX[$i]}".matlab.out"
fi


done
