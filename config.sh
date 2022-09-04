#!/bin/bash

#########################################################
## declare constants (DO NOT CHANGE THE NUMBERS BELOW!!!)

#solvers
MIPS=1;
MIPSsc=2;
FMINCON=3;
IPOPT=4;
KNITRO=5;
BELTISTOSopf=6;
MIPSscPardiso=7;
IPOPTHSL=8;
BELTISTOSmpopf=9;
BELTISTOSmem=10;

#OPFstart
FLAT=1;
MPC=2;
PF=3;

#OPFvoltage
POLAR=0;
CARTESIAN=1;

#OPFbalance
POWER=0;
CURRENT=1;

## end of constant declaration (DO NOT CHANGE THE NUMBERS ABOVE!!!)
###################################################################


### declare an array of the OPF benchmarks
declare -a gridsX=(
    "Xbeta_1024_1000.csv"
    "Xbeta_2048_1000.csv"
    "Xbeta_4096_1000.csv"
    "Xbeta_8192_1000.csv"
    "Xbeta_16384_1000.csv"
    "Xbeta_32768_1000.csv"
    "Xbeta_65536_1000.csv"
    "Xbeta_131072_1000.csv"
    "Xbeta_262144_1000.csv"
    ) 


## declare an array of the benchmarks
declare -a gridsB=(
    "beta_eLR_1024_1000.csv"
    "beta_eLR_2048_1000.csv"
    "beta_eLR_4096_1000.csv"
    "beta_eLR_8192_1000.csv"
    "beta_eLR_16384_1000.csv"
    "beta_eLR_32768_1000.csv"
    "beta_eLR_65536_1000.csv"
    "beta_eLR_131072_1000.csv"
    "beta_eLR_262144_1000.csv"
    ) 

