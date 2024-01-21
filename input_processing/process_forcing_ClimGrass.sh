#!/bin/bash

#PBS -N preprocess_forcing
#PBS -P pt17
#PBS -q normal
#PBS -l walltime=10:30:00
#PBS -l mem=32GB
#PBS -l ncpus=1
#PBS -l storage=gdata/pr09
#PBS -l software=netCDF
#PBS -r y
#PBS -l wd
#PBS -j oe
#PBS -S /bin/bash

# process forcing data as provided by ClimGrass MIP.

# paths and files
forcing_path="/g/data/pr09/ClimGrass/forcing"
files=$(ls ${forcing_path}/*_met.nc)

cd $forcing_path

# constants (note: constants not used below because not sure how to do it)
#Kelvin=273.15
#ts=1800     # seconds per timestep
#Pa2kPa=0.001

for file in ${files} ; do

  echo ${file}

  # Note: variable names do not need to be changed because they can be set in the settings file.  
  # 1) convert units
  cdo -O aexpr,'Tair=Tair-273.15;' ${file} tmp.nc    # K    -> degC
  cdo -O aexpr,'Precip=Precip*1800;' tmp.nc tmp1.nc  # mm/s -> mm/timestep
  ncatted -O -a units,Tair,o,c,"degC" tmp1.nc	
  ncatted -O -a units,Precip,o,c,"mm" tmp1.nc	

  # 2) calculate VPD from RH and Tair (from bigleaf R package function rH.to.VPD)
  cdo -O aexpr,'Esat=0.001*611.2*exp((17.62*Tair)/(243.04+Tair));' tmp1.nc tmp.nc
  cdo -O -z zip_6 aexpr,'VPD=Esat - RH/100 * Esat;' tmp.nc tmp1.nc
  ncatted -O -a units,Esat,c,c,"kPa" tmp1.nc
  ncatted -O -a long_name,Esat,c,c,"Saturation vapour pressure" tmp1.nc
  ncatted -O -a units,VPD,c,c,"kPa" tmp1.nc
  ncatted -O -a long_name,VPD,c,c,"Vapour pressure deficit" tmp1.nc

  # 3) delete redundant dimensions in latitude and longitude and rename
  ncwa -O -v longitude -a y tmp1.nc lon.nc
  ncwa -O -v latitude -a x tmp1.nc lat.nc

  ncks -O -x -v longitude,latitude tmp1.nc tmp.nc 
  ncks -A lon.nc tmp.nc
  ncks -A lat.nc tmp.nc

  ncrename -O -d x,longitude tmp.nc
  ncrename -O -d y,latitude tmp.nc

  # 4) reshuffle dimensions (better for LPJ-GUESS, but does not seem to work in this case)
  ncpdq -O -a latitude,longitude,time tmp.nc tmp1.nc

  # 5) reset name
  mv tmp1.nc ${file}

done


# 6) cleanup
rm tmp*.nc lat.nc lon.nc
chmod 775 *.nc