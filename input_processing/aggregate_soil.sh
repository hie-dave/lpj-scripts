#!/bin/bash

## Script to aggregate soil parameters from 90m to a target resolution
## as given by gridfile.
## Juergen Knauer, May 2022
## Input:
## Output:

## PBS settings
#PBS -N pbs_aggregate_soil
#PBS -P hw83
#PBS -q hugemem
#PBS -l walltime=08:30:00
#PBS -l mem=384GB
#PBS -l ncpus=1
# #PBS -l jobfs=1GB
#PBS -l storage=gdata/hw83
#PBS -r y
#PBS -l wd
#PBS -j oe
#PBS -S /bin/bash

## settings
basepath="/g/data/hw83/SLG"
gridfile="/scratch/hw83/jk8585/Narclim_grid.nc"
vars="CLY SLT SND pHc CN SOC"
slg_version="20140801"
outfile_tag="0.1deg"  # descriptor added to output file
layer_depths="5 15 30 60 100"
#layer_depths="5 15 30 60 100 200"  # cumulative depth in cm

## load modules
module load gdal/3.0.2
module load cdo/1.9.10


## Start calculations (no changes needed beyond that point)
total_depth=$(echo "$layer_depths" | rev | cut -d" " -f1 | rev)
for var in $vars ; do	

  for layer in $layer_depths ; do
	
    ## get layer information (could be moved outside of loop, but ok)
    if [[ -z $dmax ]] ; then dmin=0 ; else dmin=$dmax ; fi
    dmax=$layer
    ld=$(echo "$dmax - $dmin" | bc)  # layer depth
    fd=$(echo "scale=3; $ld / $total_depth" | bc) # fraction of layer depth to total depth

    ## convert layer depths into XXX format
    sdmin=$dmin
    sdmax=$dmax
    if [[ ${#sdmin} < 2 ]] ; then sdmin=0${sdmin} ; fi
    if [[ ${#sdmin} < 3 ]] ; then sdmin=0${sdmin} ; fi				  
    if [[ ${#sdmax} < 2 ]] ; then sdmax=0${sdmax} ; fi
    if [[ ${#sdmax} < 3 ]] ; then sdmax=0${sdmax} ; fi

    ## convert geotiff to netcdf
    if [[ "$var" == "CN" ]] ; then
	gdal_translate -of NetCDF ${basepath}/SOC_${sdmin}_${sdmax}_EV_N_P_AU_NAT_C_${slg_version}.tif tmp_C.nc
	gdal_translate -of NetCDF ${basepath}/NTO_${sdmin}_${sdmax}_EV_N_P_AU_NAT_C_${slg_version}.tif tmp_N.nc
	cdo div tmp_C.nc tmp_N.nc tmp.nc
    else
        gdal_translate -of NetCDF ${basepath}/${var}_${sdmin}_${sdmax}_EV_N_P_AU_NAT_C_${slg_version}.tif tmp.nc
    fi
    
    ## weighting by layer depth 
    cdo -O -s mulc,${fd} tmp.nc tmp_${layer}.nc
  done
  unset dmax
  
  ## sum up all (weighted) soil layers
  cdo -O -s enssum tmp_*.nc tmp2.nc

  ## make copies of clay, silt and sand to be used later
  if [[ "$var" == "SND" ]] ; then cp tmp2.nc snd_tmp.nc ; fi
  if [[ "$var" == "SLT" ]] ; then cp tmp2.nc slt_tmp.nc ; fi
  if [[ "$var" == "CLY" ]] ; then cp tmp2.nc cly_tmp.nc ; fi

  
  if [[ "$var" != "SND" && "$var" != "SLT" && "$var" != "CLY" ]] ; then

    ## aggregate to target resolution (as in gridfile)
    cdo -O -s remapcon,${gridfile} tmp2.nc ${var}_${outfile_tag}.nc

    ## rename variable and change attributes
    case $var in
      pHc)
	  varname="pH"
	  ;;
      SOC)
	  varname="orgC"
	  ;;
      *)
	  varname="$var"
          ;;
    esac
  
     ncrename -v Band1,${varname} ${var}_${outfile_tag}.nc
  fi
    
done


## ensure that sand, silt and clay add up to 100% at the original resolution
cdo -s enssum snd_tmp.nc slt_tmp.nc cly_tmp.nc sum_tmp.nc
cdo -s subc,100 sum_tmp.nc resid_tmp.nc
cdo -s div snd_tmp.nc sum_tmp.nc sndfrac_tmp.nc
cdo -s div slt_tmp.nc sum_tmp.nc sltfrac_tmp.nc
cdo -s div cly_tmp.nc sum_tmp.nc clyfrac_tmp.nc

cdo -s mulc,-1 -mul sndfrac_tmp.nc resid_tmp.nc addsnd_tmp.nc
cdo -s mulc,-1 -mul sltfrac_tmp.nc resid_tmp.nc addslt_tmp.nc
cdo -s mulc,-1 -mul clyfrac_tmp.nc resid_tmp.nc addcly_tmp.nc
     
cdo -s -divc,100 -add addsnd_tmp.nc snd_tmp.nc sndnew_tmp.nc
cdo -s -divc,100 -add addslt_tmp.nc slt_tmp.nc sltnew_tmp.nc
cdo -s -divc,100 -add addcly_tmp.nc cly_tmp.nc clynew_tmp.nc

cdo -O -s remapcon,${gridfile} sndnew_tmp.nc SND_${outfile_tag}.nc
cdo -O -s remapcon,${gridfile} sltnew_tmp.nc SLT_${outfile_tag}.nc
cdo -O -s remapcon,${gridfile} clynew_tmp.nc CLY_${outfile_tag}.nc

# rename variables
ncrename -v Band1,sand SND_${outfile_tag}.nc
ncrename -v Band1,silt SLT_${outfile_tag}.nc
ncrename -v Band1,clay CLY_${outfile_tag}.nc


##  merge all variables to one netcdf file
if [[ -f soilvars_${outfile_tag}.nc ]] ; then rm soilvars_${outfile_tag}.nc ; fi
cdo -O merge *_${outfile_tag}.nc soilvars_${outfile_tag}.nc


## cleanup
rm *tmp*.nc
