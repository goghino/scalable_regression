## Synthetic scalable regression model

This repository contains source code of the synthetic scalable regression model
used in the "On cheap entropy-sparsified regression learning" manuscript.

The code can be run using three solvers
* MathWorks Fmincon
* Ipopt
* Panua-Ipopt

MathWorks Fmincon is distributed together with [MATLAB](https://ch.mathworks.com/products/matlab.html "MATLAB's Homepage"),
Ipopt can be accessed at [this repository](https://github.com/coin-or/Ipopt "Ipopt's Homepage") and
Panua-Ipopt can be accessed via [Panua Technologies](http://panua.ch/ "Panua Technologies"). This repository contains
binary files of Ipopt and Panua-Ipopt solvers, located at `Ipopt\` and `PanuaIpopt\` folders respectively.

#### Running the benchmarks

The main benchmark file is `runBench.sh` which configures the regression model to be run and solvers and executes all the runs.
In order to configure the benchmarks, follow the steps below:

1. Select the solver:
```sh
# USER CONFIGURATION
RUN_MATLAB=1 #0 or 1
SOLVER=$FMINCON #IPOPT or FMINCON, effective only if RUN_MATLAB=1
```
The MathWorks' Fmincon and Ipopt solvers are run using the MATLAB implementation of the regression model.
Thus, the variable `RUN_MATLAB` needs to be set to `1`. Then the Fmincon or Ipopt solver is selected using
the variable `SOLVER=$FMINCON` or `SOLVER=$IPOPT`.

On the other hand, Panua-Ipopt codebase is completely C++ based, thus in order to execute the code set `RUN_MATLAB` variable to `0`.

2. Select the regression models:
```sh
LOAD_FILES=1 #0 or 1, effective only if RUN_MATLAB=0
## load list of benchmarks to run
. config.sh 
```
The regression models can be initialized by loading the data files from folder `data/` or they can be generated on the fly.
In order to load the models from the file set `LOAD_FILES=1`. The individual models are then listed in the `config.sh` file, e.g.
```
declare -a gridsX=(
    "Xbeta_1024_1000.csv"
    "Xbeta_2048_1000.csv"
    "Xbeta_4096_1000.csv"
    )

declare -a gridsB=(
    "beta_eLR_1024_1000.csv"
    "beta_eLR_2048_1000.csv"
    "beta_eLR_4096_1000.csv"
    )
```
In order to generate larger models use script `generator.m` in the `data/` folder (they are not included in this repository
due to their excessive size). In order to generate the models on the fly set `LOAD_FILES=0`. Note that this option is currently
supported only in the C++ code and Panua-Ipopt solver.

3. Processing the results:

Each benchmark produces two output files, one with the solution and one with the solver output. In order to extract the memory and timing information,
you can search the solver output logs, using e.g. for Ipopt and Panua-Ipopt logs
```
grep "Maximum resident set size (kbytes)" <output_log_file>
grep "OverallAlgorithm" <output_log_file>
```

or for Fmincon logs use


```
grep "Maximum resident set size (kbytes)" <output_log_file>
grep "Elapsed time" <output_log_file>
```
