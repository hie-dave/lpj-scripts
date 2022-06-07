#!/usr/bin/env Rscript

## Script that creates a gridlist file for LPJ-GUESS
# Juergen Knauer, May 2022
# Inputs/Settings:
# 1) Ncdf file containing all soil parameters
# 2) Spatial extent c(lonmin,lonmax,latmin,latmax) and resolution in deg


library(ncdf4)    # make sure package is installed and path linked in .bashrc file
options(digits=9)

## Settings
# Soil parameter file (netcdf file containing all required soil input variables)
soil_file <- "/g/data/hw83/SLG/soilvars_0.1deg.nc"
outname   <- "LPJ_soilparams.dat"  # name of output file (with extension)

## bounding box
lonmin <- 140.8
lonmax <- 150
latmin <- -39.2
latmax <- -33.9

soilvars      <- c("clay","silt","sand","orgC","CN","pH")
prec_lonlat   <- 2  # precision of lon and lat in output file
prec_soilvars <- 3  # precision of soil variables in output file


## get lonlat and create gridlist with all possible lonlat combinations
sf      <- nc_open(soil_file)
lat     <- ncvar_get(sf,"lat")
lon     <- ncvar_get(sf,"lon")
gridlist_all <- round(expand.grid(lon,lat),prec_lonlat) 

## get soilvars and add them to gridllist
for (varname in soilvars){
   var          <- ncvar_get(sf,varname)
   gridlist_all <- cbind(gridlist_all,round(c(var),prec_soilvars))
}

# name columns
colnames(gridlist_all) <- c("lon","lat",soilvars)

# set all lats and lons outside of bounding box to NAs
gridlist_all[gridlist_all[,"lon"] > lonmax | gridlist_all[,"lon"] < lonmin & !is.na(gridlist_all[,"lon"]),-c(1,2)] <- NA
gridlist_all[gridlist_all[,"lat"] > latmax | gridlist_all[,"lat"] < latmin & !is.na(gridlist_all[,"lat"]),-c(1,2)] <- NA

# exclude all lats and lons with NAs (assumed to be not land)
# assumes that all soil variables have the same land extent (maybe add check)
gridlist <- gridlist_all[-which(is.na(gridlist_all[,soilvars[1]])),]

# account for rounding errors and ensure sum of silt, sand, and clay add up to 100%
gridlist[,'silt'] <- gridlist[,'silt'] - (apply(gridlist[,c('sand','silt','clay')],1,sum,na.rm=T) - 100)

# write to file
write.table(gridlist,file=outname,quote=F,row.names=F)







  