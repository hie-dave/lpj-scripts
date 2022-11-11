#!/usr/bin/env bash
# Script to efficiently transfer files from one location on gadi to another
# location on gadi. In this case applied to transfer NARCliM files to scratch.

# Fail immediately with error message if any individual command fails.
set -euo pipefail

# If environment variable DEBUG is set to 1, run in debug mode.
if [ "${DEBUG:-}" = 1 ]
then
  echo "Running in debug mode. Set DEBUG to 0 to disable."
  set -x
fi

basepath_in="/g/data/at43/output/NARCliMi/UNSW"

dave="/data/hiestorage/WorkingData/MEDLYN_GROUP/PROJECTS/dynamics_simulations"
dest_dir="${dave}/narclim1.5/forcing-raw"

globmods="CCCma-CanESM2 CSIRO-BOM-ACCESS1-0 CSIRO-BOM-ACCESS1-3" # global models
regmods="UNSW-WRF360J UNSW-WRF360K"  # regional (downscaled) climate models 
scenarios="historical rcp45 rcp85"   # climate scenarios
vars="tasmin-bc tasmax-bc"  # variables
version="v1"
freq="day"

# Set to 0 to disable parallel mode.
PARALLEL=1

# Max number of parallel jobs.
MAX_NUM_JOBS=10

# If enabled, count the files to be downlaoded by don't actually download them.
DRY_RUN=0

# Name of the directory used to synchronise progress reporting between workers.
LOCK_FILE="$(mktemp --dry-run --directory .copy_files_XXX.lock)"

# Name of the file used to store number of completed downloads.
PROGRESS_FILE="$(mktemp .copy_files_XXX.progress)"

gadi_user=dh7190
hies_user=u30062639

src_server="${gadi_user}@gadi-dm.nci.org.au"
dest_server="${hies_user}@hie-storage.sstars.ws"
cp_command="rsync -a --partial"

# Initialise progress file.
printf "%d" 0 >"${PROGRESS_FILE}"

if [ ${DRY_RUN} -eq 1 ]
then
	echo "Running in dry-run mode. No files will be downloaded."
fi

# Check if server connection is possible.
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

echo "Connection to server established."

# Function to sleep until number of background jobs is less than MAX_NUM_JOBS.
function wait_for_jobs() {
	while [ $(jobs | wc -l) -ge ${MAX_NUM_JOBS} ]
	do
		sleep 5
	done
}

# Total number of files to de downloaded. This is used for progress reporting.
# It could be determined dynamically, but I leave that as an exercise to the
# reader.
total_files=2756

# Start time of script execution.
START_TIME=$(date +%s)

# # The new method, which recreates the old directory structure, but has the
# # advantage of showing overall progress.
# ssh "${src_server}" ssh "${dest_server}" mkdir -p "${dest_dir}"
# ssh "${src_server}" "${cp_command}" "${basepath_in}/*" "${dest_server}:${dest_dir}/"

# Acquire mutex lock (busy wait until lock is acquired).
function lock() {
	while ! mkdir "${LOCK_FILE}" >/dev/null 2>&1
	do
		sleep 0.5
	done
}

# Release mutex lock.
function unlock() {
	rmdir "${LOCK_FILE}"
}

# Print a time in seconds in hh:mm:ss format.
function get_time_str() {
	TIME=${1}
	HOURS=$(echo "${TIME} / 3600" | bc)
	MINUTES=$(echo "(${TIME} - ${HOURS} * 3600) / 60" | bc)
	SECONDS=$(echo "${TIME} - ${HOURS} * 3600 - ${MINUTES} * 60" | bc)
	printf "%02d:%02d:%02d" ${HOURS} ${MINUTES} ${SECONDS}
}

# Called when a download is completed.
function download_completed() {
	lock

	num_completed=$(echo `cat "${PROGRESS_FILE}"` + 1 | bc)

	# Update progress file.
	printf "%d" ${num_completed} >"${PROGRESS_FILE}"


	progress=$(echo "${num_completed} / ${total_files}" | bc -l)
	progress_percent=$(echo "100 * ${progress}" | bc -l)

	# Get elapsed time in seconds.
	ELAPSED_TIME=$(echo "$(date +%s) - ${START_TIME}" | bc)

	# Calculate remaining time (in seconds).
	TOTAL_TIME=$(echo "${ELAPSED_TIME} / ${progress}" | bc -l)
	REMAINING_TIME=$(echo "${TOTAL_TIME} - ${ELAPSED_TIME}" | bc -l)

	# Report progress to user.
	printf "Working: %.2f%%; elapsed: %s; remaining: %s\r" \
	${progress_percent} $(get_time_str ${ELAPSED_TIME}) \
	$(get_time_str ${REMAINING_TIME})

	# Release mutex lock.
	unlock
}

# The older method, which creates a nice directory structure in dest. However,
# this doesn't report overall progress (only per-file).
total_downloaded=0
for globmod in $globmods
do
    for regmod in $regmods
	do
		for scenario in $scenarios
		do
			for var in $vars
			do
				# Create destination directory if it doesn't exist.
				dest_path=${dest_dir}/${scenario}/${var}
				mkdir -p ${dest_path}

				# List files on the source server to be copied.
				cp_files=$(ssh "${src_server}" echo ${basepath_in}/${globmod}/${scenario}/r1i1p1/${regmod}/${version}/${freq}/${var}/*.nc)

				for src_file in $cp_files
				do
					# Destination path.
					dest="${dest_path}/$(basename "${src_file}")"

					if [ ${DRY_RUN} -ne 1 ]
					then
						# Create destination directory if it doesn't already exist.
						mkdir -p "${dest_path}"

						# Need to fork the command if parallel mode enabled.
						PARA=
						if [ ${PARALLEL} -eq 1 ]
						then
							PARA=&
						fi

						# Sleep until number of ongoing downloads is less than max.
						wait_for_jobs

						# Copy the file.
						eval "${cp_command}" "${src_server}:${src_file}" "${dest}" \
						&& download_completed ${PARA}
					fi
					total_downloaded=$(echo "${total_downloaded} + 1" | bc)
				done
			done
		done
    done
done

printf "\n"
echo "Successfully downloaded ${total_downloaded} files."
