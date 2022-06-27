The folder 'input_processing' contains scripts for preparing LPJ-GUESS inputs.


aggregate_soil.sh:

aggregates soil variables as given by the Australian Soil and Landscape Grid at 90m to a
target resolution as specified in gridfile. The script requires input soil maps in tiff format
which are first converted to netcdf and then aggregated over soil layers (weighted by depth)
and then spatially using conservative remapping. The script ensures that sand, silt, and clay
fractions add up to 1 at the target resolution. The output is a netcdf file containing all
selected soil variables.


create_gridlist_LPJ.R

the script creates two files:
1) a soil parameter file as required by LPJ-GUESS.
2) a 'gridlist' file defining the grid cells run by LPJ-GUESS
The input is a netcdf file as created by running 'aggregate_soil.sh' (see above). The input file contains
all required soil parameters and also acts as a 'land mask', i.e. only grid cells containing
land according to the input map are added to the list. The soil params are converted from 2D to 1D and
additionally cropped to a bounding box specified in the script. 


copy_files.sh

script to efficiently copy files (such as climate forcing) from one gadi location to another gadi location.
Can be used as a template for similar jobs in the future.


process_climate_input.sh

script to process climate forcing files into a format directly useable by LPJ-GUESS.
At the moment specific to NARCliM data (which don't require much modifications). But
can be extended for different data sets, whcih may require regrdidding, unit conversions,
time axis adjustments, renaming etc.
