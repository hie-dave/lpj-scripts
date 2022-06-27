#!/bin/bash

## Script to efficiently transfer files from one location on gadi to another location on gadi
# In this case applied to transfer NARCliM files to scratch

## PBS settings
#PBS -N download_clim
#PBS -P hw83
#PBS -q copyq
#PBS -l walltime=02:59:00
#PBS -l mem=192GB
#PBS -l ncpus=1
# #PBS -l jobfs=1GB
#PBS -l storage=gdata/at43
#PBS -r y
##PBS -l wd
#PBS -j oe
#PBS -S /bin/bash

basepath_in="/g/data/at43/output/NARCliMi/UNSW"
basepath_out="/scratch/hw83/jk8585/NARCliM1.5_orig"

globmods="CCCma-CanESM2 CSIRO-BOM-ACCESS1-0 CSIRO-BOM-ACCESS1-3" # global models
regmods="UNSW-WRF360J UNSW-WRF360K"  # regional (downscaled) climate models 
scenarios="historical rcp45 rcp85"   # climate scenarios
vars="hurs pr rsds sfcWind tas tasmin-bc tasmax-bc"  # variables
version="v1"
freq="day"

for globmod in $globmods ; do
    for regmod in $regmods ; do
	for scenario in $scenarios ; do
	    for var in $vars ; do
		# create target directory
		outpath=${basepath_out}/${scenario}/${var}
		mkdir -p ${outpath}

		# list files to copy
		cp_files=$(ls ${basepath_in}/${globmod}/${scenario}/r1i1p1/${regmod}/${version}/${freq}/${var}/*.nc)

		# copy files
		for cp_file in $cp_files ; do
                    echo $(basename $cp_file)
		    cp -u $cp_file ${outpath}/$(basename $cp_file)
		done
	    done
	done
    done
done
    
