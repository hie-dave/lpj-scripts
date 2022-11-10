#!/usr/bin/env bash
# Script to efficiently transfer files from one location on gadi to another
# location on gadi. In this case applied to transfer NARCliM files to scratch.

# Fail immediately with error message if any individual command fails.
set -euo pipefail

basepath_in="/g/data/at43/output/NARCliMi/UNSW"

dave="/data/hiestorage/WorkingData/MEDLYN_GROUP/PROJECTS/dynamics_simulations"
dest_dir="${dave}/narclim1.5/forcing-raw"

globmods="CCCma-CanESM2 CSIRO-BOM-ACCESS1-0 CSIRO-BOM-ACCESS1-3" # global models
regmods="UNSW-WRF360J UNSW-WRF360K"  # regional (downscaled) climate models 
scenarios="historical rcp45 rcp85"   # climate scenarios
vars="tasmin-bc tasmax-bc"  # variables
version="v1"
freq="day"

gadi_user=jk8585
hies_user=u30044953

src_server="${gadi_user}@gadi-dm.nci.org.au"
dest_server="${hies_user}@hie-storage.sstars.ws"
cp_command="rsync -a --partial --info=progress2 --no-i-r"

if ! ssh "${src_server}" echo >/dev/null 2>&1
then
	echo "Error: Unable to login to src server ${src_server}."
	exit 1
fi

if ! ssh "${src_server}" ssh "${dest_server}" echo >/dev/null 2>&1
then
	echo "Error: Cannot connect from ${src_server} to ${dest_server}"
	exit 1
fi

# # The new method, which recreates the old directory structure, but has the
# # advantage of showing overall progress.
# ssh "${src_server}" ssh "${dest_server}" mkdir -p "${dest_dir}"
# ssh "${src_server}" "${cp_command}" "${basepath_in}/*" "${dest_server}:${dest_dir}/"

# The older method, which creates a nice directory structure in dest. However,
# this doesn't report overall progress (only per-file).
for globmod in $globmods
do
    for regmod in $regmods
	do
		for scenario in $scenarios
		do
			for var in $vars
			do
				# create target directory
				outpath=${dest_dir}/${scenario}/${var}
				#mkdir -p ${outpath}

				# list files to copy
				cp_files=$(echo ${basepath_in}/${globmod}/${scenario}/r1i1p1/${regmod}/${version}/${freq}/${var}/*.nc)

				# copy files
				for src_file in $cp_files
				do
					# Get raw filename without path.
					dest_file="$(basename "${src_file}")"

					# Destination path.
					dest_path="${dest_dir}/${outpath}"
					dest="${dest_server}:${dest_path}/${dest_file}"

					# Create destination directory if it doesn't already exist.
					ssh "${dest_server}" mkdir -p "${dest_path}"

					# Copy the file.
					ssh "${src_server}" "${cp_command}" "${src_file}" "${dest}"
				done
			done
		done
    done
done
