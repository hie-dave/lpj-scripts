#!/bin/bash

outdir="/scratch/pt17/jk8585/LPJ_runs/output"
logdir="/scratch/pt17/jk8585/LPJ_runs/logs"
codedir="/home/599/jk8585/LPJ_code/trunk_r8538"

if [ -L "exe" ]; then
    rm $exe
fi

ln -s ${codedir}/guess guess

if [ -d ${outdir} ]; then
    rm -rf $outdir
    mkdir $outdir
fi

if [ -d ${logdir} ]; then
    rm -rf $logdir
    mkdir $logdir
fi
