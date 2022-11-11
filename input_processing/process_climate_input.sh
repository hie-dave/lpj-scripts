#!/bin/bash

## Script to process climate input data for LPJ-GUESS.
## Note that this applies to NARCliM1.5 data only. Also note that
## is highly dataset-specific and requirements might completely different for
## another data set.

## Juergen Knauer
## June 2022

## PBS settings
#PBS -N process_clim_input
#PBS -P hw83
#PBS -q normal
#PBS -l walltime=04:59:00
#PBS -l mem=192GB
#PBS -l ncpus=1
# #PBS -l jobfs=1GB
#PBS -l storage=scratch/hw83
#PBS -r y
#PBS -l wd
#PBS -j oe
#PBS -S /bin/bash


## Settings
cf_path_in="/home/drew/Downloads/narclim"
cf_path_out="/home/drew/Downloads/narclim-out"

globmods="CSIRO-BOM-ACCESS1-0" # global models
regmods="UNSW-WRF360J"  # regional (downscaled) climate models 
scenarios="historical"   # climate scenarios
vars="tasmax-bc tasmin-bc"      # variables
version="v1"
freq="day"
lonlatbox="130.0,155.0,-40.0,-20.0"  # format: "lon1,lon2,lat1,lat2"
# output name format: var_scenario_globmod_regmod_freq_year.nc

### Start processing
for var in $vars ; do
    for scenario in $scenarios ; do
        if [[ "${scenario}" == "historical" ]] ; then
			startyear=1951
			endyear=2005
		else
			startyear=2006
			endyear=2099
		fi
		mkdir -p ${cf_path_out}/${scenario}/${var}
		for globmod in $globmods ; do
			for regmod in $regmods ; do

				# 1) bring to annual resolution (if necessary)
				if [[ "${globmod}" == "CSIRO-BOM-ACCESS1-3" && "${scenario}" == "historical" ]] ; then
					for ((year=${startyear};year<=${endyear};year++)) ; do
						if [[ $(echo "${year}%5" | bc) -eq 1 ]] ; then
							syear=${year}
							eyear=$(echo "${year} + 4" | bc)
							cf_file=${cf_path_in}/${scenario}/${var}/${var}_NARCliMi_${globmod}_${scenario}_r1i1p1_${regmod}_${version}_${freq}_${syear}0101-${eyear}1231.nc
							cdo splityear $cf_file ${cf_path_out}/${scenario}/${var}/tmp_
						fi
					done
				fi

				#                                            hurs_historical_CSIRO-BOM-ACCESS1-0_UNSW-WRF360J_daily_1951_2005.nc
				out_pfx="${cf_path_out}/${scenario}/${var}/${var}_${scenario}_${globmod}_${regmod}_${freq}"

				for ((year=${startyear};year<=${endyear};year++)) ; do
					# get filename
					if [[ "${globmod}" == "CSIRO-BOM-ACCESS1-3" && "${scenario}" == "historical" ]] ; then
						cf_file=${cf_path_out}/${scenario}/${var}/tmp_${year}.nc
					else
						cf_file=${cf_path_in}/${scenario}/${var}/${var}_NARCliMi_${globmod}_${scenario}_r1i1p1_${regmod}_${version}_${freq}_${year}0101-${year}1231.nc
					fi

					# 1) crop file
					cdo -s sellonlatbox,${lonlatbox} ${cf_file} tmp.nc

					# 2) correct units and change names
					case "$var" in
					hurs)
						cdo -s divc,100 tmp.nc tmp1.nc
						ncatted -O -a units,${var},o,c,"1" tmp1.nc
						mv tmp1.nc tmp.nc
					;;
					#pr)
					# nothing to be done for pr
					#;;
					#tas)
					# nothing to be done for tas
					#;;
					#rsds)
					# nothing to be done for rsds
					#;;
					#sfcWind)
					# nothing to be done for sfcWind
					#;;
					esac

					# 3) rename file and save at final location
					mv tmp.nc "${out_pfx}_${year}.nc"

					# cleanup
					if [[ -f ${cf_path_out}/${scenario}/${var}/tmp_${year}.nc ]] ; then
						rm ${cf_path_out}/${scenario}/${var}/tmp_${year}.nc
					fi
				done

				# Merge annual files into a single file containing entire timeseries.
				cdo mergetime ${out_pfx}*.nc ${out_pfx}_${startyear}_${endyear}.nc

				# Remove intermediate files.
				for ((year=${startyear};year<=${endyear};year++))
				do
					rm "${out_pfx}_${year}.nc"
				done

			done
		done
	done
done

# cleanup
if [[ -f tmp.nc ]] ; then
    rm tmp.nc
fi
