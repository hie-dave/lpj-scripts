#!/bin/bash

## JK: this script is an updated version of initial_setup.sh which also runs
#      'submit_to_gadi.sh'


## specify paths
homedir="/home/599/jk8585/LPJ_run"
outdir="/scratch/hw83/jk8585/LPJ_runs"
codedir="/home/599/jk8585/LPJ_code/trunk_r8538"
submit_script="${homedir}/submit_to_gadi.sh"
experiment="test_run"

## run settings
gridfile="${homedir}/gridlist_oz_cf.txt"


## create links and folders
mkdir -p ${outdir}/${experiment}
mkdir -p ${outdir}/${experiment}/output
mkdir -p ${outdir}/${experiment}/logs


## change to outdir and link to executable
cd ${outdir}/${experiment}
ln -sf ${codedir}/guess guess


## copy submit script and other files to outdir (right now a requirement of 'submit_to_gadi.sh')
cp ${gridfile} .
cp ${submit_script} .
cp ${homedir}/global.ins .
cp ${homedir}/global_cf.ins .
cp ${homedir}/global_soiln.ins .



## make changes to ${submit_script}




## run ${submit_script}
./$(basename ${submit_script})
