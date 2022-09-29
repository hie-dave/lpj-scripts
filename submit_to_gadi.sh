#!/usr/bin/env bash
#
# Portable bash script to run LPJ-Guess as a parallel job using PBS on Gadi.
#
# Adapted by Juergen Knauer and Drew Holzworth from older scripts.
#
# Usage:
#
# submit_to_gadi.sh -s <config-file> [-h] [-d]
#
# -s  Path to the configuration file, which provides run settings such as
#     number of CPUs, instruction file path, etc.
# -h  Display this help info
# -d  Dry run. Create the run directory but don't actually submit the
#     job to PBS.
#

# Exit immediately if:
# - Any command fails
# - Any part of a pipeline fails
# - We reference any unbound variables
set -euo pipefail

# If environment variable DEBUG is set to 1, run in debug mode.
if [ "${DEBUG:-}" = 1 ]
then
  echo "Running in debug mode. Set DEBUG to 0 to disable."
  set -x
fi

# Get directory containing this script.
SUBMIT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd -P )"

# Parse the command line arguments
CONF=
DRY_RUN=0
USAGE="Usage: ${0} -s <config-file> [-h] [-d]

-s  Path to the configuration file, which provides run settings such as
    number of CPUs, instruction file path, etc.
-h  Display this help info
-d  Dry run. Create the run directory but don't actually submit the job
    to PBS."

while getopts ":s:hd" opt; do
    case $opt in
      s ) source "${OPTARG}"; CONF=1 CONFIG_FILE="${OPTARG}" ;;
      h ) echo "${USAGE}"; exit 0 ;;
      d ) DRY_RUN=1
    esac
done

# Exit if no config file provided.
if [ -z "${CONF:-}" ];
then
  echo "Error: no config file supplied."
  echo "${USAGE}"
  exit 1
fi

# This function will print the absolute path to the given filename.
# We could use readlink, but it's a Linux tool - not available on MacOS (afaik).
get_absolute_path() {
  echo $( cd "$(dirname "$1")"; pwd -P )/$(basename "$1")
}

# Function which prints the name of all instruction files imported,
# recursively, by the specified instruction file.
function get_all_insfiles() {
	echo $1
	# import "file.ins"
	local rx='^[ \t]*import[ \t]*"([^"]+)"\r?'
	local referenced_files="$(sed -rn "s/${rx}/\1/gmp" $1)"
	for file in ${referenced_files}
	do
		if [ -n "${file}" ]
		then
			get_all_insfiles "${file}"
		fi
	done
}

# Function which prints the name of the gridlist file in the specified ins file.
function get_gridlist() {
	local gridlist_variable_name='file_gridlist_cf'
	# param "file_gridlist_cf" (str "CFA_gridlist.txt")
	local gridlist_regex="[ \t]*param[ \t]*\"${gridlist_variable_name}\"[ \t]*\(str[ \t]*\"([^\"]+)\"[ \t]*\)\r?"

	# Note: this will not work if any of the instruction file names contain
	# a newline character.
	while IFS=$'\n' read ins_file
	do
		local gridlist="$(sed -rn "s/${gridlist_regex}/\1/gmp" "${ins_file}")"
		if [ -n "${gridlist}" ]
		then
			echo "${gridlist}"
			return 0
		fi
	done < <(get_all_insfiles $1)

	# Gridlist not found.
  echo "Unable to parse gridlist file from .ins file ${ins_file}"
	return 1
}

# Read gridlist file name from the .ins file.
GRIDLIST="$(get_gridlist "${INSFILE}")"

# Convert to absolute paths.
INSFILE="$(get_absolute_path "${INSFILE}")"
BINARY="$(get_absolute_path "${BINARY}")"
OUT_DIR="$(get_absolute_path "${OUT_DIR}")"
GRIDLIST="$(get_absolute_path "${GRIDLIST}")"

function print_settings() {
  echo OUT_DIR=${OUT_DIR}
  echo BINARY=${BINARY}
  echo EXPERIMENT=${EXPERIMENT}
  echo
  echo NPROCESS=${NPROCESS}
  echo WALLTIME=${WALLTIME}
  echo MEMORY=${MEMORY}
  echo QUEUE=${QUEUE}
  echo PROJECT=${PROJECT}
  echo EMAIL=${EMAIL}
  echo EMAIL_NOTIFICATIONS=${EMAIL_NOTIFICATIONS}
  echo JOB_NAME=${JOB_NAME}
  echo
  echo INSFILE=${INSFILE}
  echo INPUT_MODULE=${INPUT_MODULE}
  echo GRIDLIST=${GRIDLIST}
  echo
}
print_settings

# Input checking/validation.

# This function will perform the specified check on the given file, and if the
# check fails, will print the specified message and exit.
#
# Arguments:
# - An argument to the test command - e.g. -f to check if file exists.
# - The file path
# - Error message
check_permission() {
  if [ ! $1 "$2" ]
  then
    echo $3
    exit 1
  fi
}

# Check if a file exists. If it does not, print an error message and return
# non-zero exit code.
#
# Arguments:
# - The file to be checked
check_exists() {
  check_permission -f $1 "Error: $1 does not exist or is not a file"
}

# Check if a file or directory has read permission for the current user. If it
# doesn't, print an error message and return a non-zero exit code.
#
# Arguments:
# - The file to be checked
check_read() {
  check_permission -r $1 "Error: $1 does not have read permission"
}

# Check if a file or directory has execute permission for the current user. If
# it doesn't, print an error message and return a non-zero exit code.
#
# Arguments:
# - The file to be checked
check_execute() {
  check_permission -x $1 "Error: $1 does not have execute permission"
}

# Ensure the input files exist and have the required permissions.
check_exists "${BINARY}"
check_exists "${INSFILE}"
check_exists "${GRIDLIST}"
check_read "${BINARY}"
check_read "${INSFILE}"
check_read "${GRIDLIST}"
check_execute "${BINARY}"

# Output path for this run.
RUN_OUT_DIR="${OUT_DIR}/${EXPERIMENT}"

# Path to which individual run directories will be created.
RUNS_DIR="${RUN_OUT_DIR}/runs"

# Path to which output will be aggregated.
OUTPUT_DIR="${RUN_OUT_DIR}/output"

# Path to which log files will be aggregated.
LOGS_DIR="${RUN_OUT_DIR}/logs"

# Create links and directories.
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${LOGS_DIR}"
mkdir -p "${RUNS_DIR}"

# Get the file name (without path) of the gridlist file.
GRIDLIST_FILENAME=$(basename ${GRIDLIST})

# This function creates the gridlist files for each run by splitting
# the original gridlist file into approximately equal parts.
function split_gridlist {
    # Create empty gridlists first to make sure each run gets one.
    for ((a=1; a <= NPROCESS ; a++))
    do
        mkdir -p "${RUNS_DIR}/run${a}"
        echo > "${RUNS_DIR}/run${a}/${GRIDLIST_FILENAME}"
    done

    # Use the split command to split the files into temporary files
    # Splitting using r/N mode for round-robin distribution.
    local tmp_prefix=tmpSPLITGRID_

    # By default, split will use a suffix of 2 digits. This will
    # obviously be insufficient if we have >99 CPUs. We could just use
    # an excessively large suffix (ie 30 - we're never going to have
    # 10^30 CPUs), but I think we can be a bit more elegant than that.
    # Here we compute the number of digits required for the suffix as
    # log10(NCPU) + 1. The only complication is that there is no builtin
    # log10 function in bash, so I've implemented a simple integer log
    # function here.
    log(){ local x=$1 n=2 l=-1;if [ "$2" != "" ];then n=$x;x=$2;fi;while((x));do let l+=1 x/=n;done;echo $l; }
    local num_digits_required=$(echo "`log 10 ${NPROCESS}` + 1" | bc)

    split -a${num_digits_required} -dn r/${NPROCESS} ${GRIDLIST} "${tmp_prefix}"

    # Move the temporary files into the runX-directories.
    local files="${tmp_prefix}*"
    local i=1
    for file in ${files}
    do
      mv ${file} "${RUNS_DIR}/run${i}/${GRIDLIST_FILENAME}"
      i=$((i+1))
    done
}

# Create header of progress.sh script
# The progress script will print the progress of each run/process.
progress_sh="${RUN_OUT_DIR}/check_progress"
cat <<EOF > "${progress_sh}"
#!/usr/bin/env bash
set -euo pipefail

# Change to the runs directory containing this script.
cd "${RUNS_DIR}"

if [ ! -d run1 ]
then
  echo "Error: Runs directory is empty. Submit script didn't run or is buggy."
  exit 1
fi

# Get total number of jobs.
num_jobs=\$(ls -ld run* | wc -l)

# Regex for parsing guess logfiles, which look like this:
#   6% complete, 00:07:37 elapsed, 01:49:57 remaining
elapsed_regex="^[^[0-9]]*([0-9]+)% complete.*\$"

# Now loop through the files and scan their logfiles to get percentage complete values.
percent_total=0

# Get the percentage complete for the first node. This should in theory
# be similar to the progress of the other nodes.
expected_percent=\$(tail -n10 run1/guess.log | sed -rn "s/\$elapsed_regex/\1/p" | tail -n1)

echo "Job is probably around \${expected_percent}% complete"

i=0
for r in run*
do
  file="\${r}/guess.log"
  if grep "Finished" "\${file}" >/dev/null
  then
    run_percent=100
  else
    # - get last few lines of file
    # - scan for "x% complete"
    # - get last match only (ie most recent progress report)
    run_percent=\$(tail -n10 "\${file}" | sed -rn "s/\$elapsed_regex/\1/p" | tail -n1)
  fi

  # Add this run's percentage complete to the total.
  percent_total=\$(echo "\${percent_total} + \${run_percent}" | bc)

  i=\$(echo "\${i} + 1" | bc)
  progress=\$(echo "100.0 * \${i} / \${num_jobs}" | bc -l)
  printf "\rParsing guess logfiles: %.2f%%" \${progress}
done
printf "\n"

# Overal progress = sum / count
percent_total=\$(echo "\${percent_total} / \${num_jobs}" | bc)

echo "\${percent_total}% complete"
EOF
chmod a+x "${progress_sh}"

# Copy the config and .ins files into the run directory, and also create
# a readme file there. This is not really necessary, but it's useful to
# have this sort of metadata readily accessible for archival purposes.
INPUTS_DIR="${RUN_OUT_DIR}/inputs"
mkdir -p "${INPUTS_DIR}"
cp "${CONFIG_FILE}" "${INPUTS_DIR}/"

# This won't work if any of the .ins file names contain a newline.
while IFS=$'\n' read ins_file
do
  cp "${ins_file}" "${INPUTS_DIR}/"
done < <(get_all_insfiles "${INSFILE}")

README_FILE="${RUN_OUT_DIR}/README.md"
echo "Job submitted $(date) with the following settings:" >"${README_FILE}"
print_settings >>"${README_FILE}"

# Create a run subdirectory for each process, clean up any existing files, and
# append code to the progress script which reports progress of this run.
for ((a=1; a <= NPROCESS ; a++))
do
    # Ensure the directory exists for this run.
    mkdir -p "${RUNS_DIR}/run${a}"

    # Delete any existing log files.
    cd "${RUNS_DIR}/run${a}"; rm -f guess.log; rm -f ${GRIDLIST_FILENAME}; cd ..
done

# Split the grid list into equally sized chunks, 1 chunk per CPU.
split_gridlist

# Get PBS options for email notifications.
# Check the qsub man page for more details.
if [ ${EMAIL_NOTIFICATIONS} -eq 1 ]
then
  EMAIL_OPT=abe
else
  EMAIL_OPT=n
fi

# Create PBS script to request place in queue
guess_cmd="${RUN_OUT_DIR}/guess.cmd"
cat <<EOF > "${guess_cmd}"
#!/bin/bash
#PBS -l ncpus=${NPROCESS}
#PBS -l walltime=${WALLTIME}
#PBS -l mem=${MEMORY}
#PBS -q ${QUEUE}
#PBS -P ${PROJECT}
#PBS -m ${EMAIL_OPT}
#PBS -M ${EMAIL}
#PBS -j oe
#PBS -l wd
#PBS -W umask=0022
#PBS -S /bin/bash
#PBS -l storage=gdata/${PROJECT}+scratch/${PROJECT}+scratch/hw83
#PBS -l jobfs=40GB
#PBS -o ${RUN_OUT_DIR}/${JOB_NAME}
set -e
module purge
module load openmpi
module load netcdf
umask 022
cd "${RUNS_DIR}"
mpirun -np ${NPROCESS} ${BINARY} -parallel -input ${INPUT_MODULE} ${INSFILE}
EOF
chmod a+x "${guess_cmd}"

# Create PBS script to generate combined output files. This will be run after
# the main job has finished running.
APPEND_JOB_NAME="${JOB_NAME}_append"
append_cmd="${RUN_OUT_DIR}/append.cmd"
cat <<EOF > "${append_cmd}"
#!/usr/bin/env bash
#PBS -l ncpus=1
#PBS -l walltime=${WALLTIME}
#PBS -l mem=8GB
#PBS -q normal
#PBS -l wd
#PBS -j oe
#PBS -m ${EMAIL_OPT}
#PBS -M ${EMAIL}
#PBS -o ${RUN_OUT_DIR}/${APPEND_JOB_NAME}.log
#PBS -W umask=0022
#PBS -l storage=gdata/${PROJECT}+scratch/${PROJECT}
set -euo pipefail
cd "${RUNS_DIR}"
function append_files {
    local number_of_jobs=\$1
    local file=\$2
    local out_dir="\$3"
    local out_file="\${out_dir}/\$file"
    cp run1/\$file "\${out_file}"
    local i=""
    for ((i=2; i <= number_of_jobs; i++))
    do
      if [ -f run\$i/\$file ]; then
        cat run\$i/\$file | awk 'NR!=1 || NF==0 || \$1 == \$1+0 { print \$0 }' >>"\${out_file}"
      fi
    done
}
pushd run1 &> /dev/null
outfiles_unexpanded='${OUTFILES}'
outfiles_expanded=\$(echo \$outfiles_unexpanded)
popd &> /dev/null
for file in \$outfiles_expanded
do
  append_files ${NPROCESS} "\$file" "${OUTPUT_DIR}"
done

# Now combine the log files.
cat run*/guess.log > "${LOGS_DIR}/guess.log"
EOF
chmod a+x "${append_cmd}"

if [ ${DRY_RUN} -eq 1 ]
then
  echo "Dry run completed (job not submitted). Job run directory created at:"
  echo "${RUN_OUT_DIR}"
  exit 0
fi

# Submit guess job
JOB_ID=$(qsub -N "${JOB_NAME}" "${guess_cmd}")
echo JOB_ID=${JOB_ID}

# Submit append job
APPEND_JOB_ID=$(qsub -W depend=afterok:${JOB_ID} -N "${APPEND_JOB_NAME}" "${append_cmd}")
echo APPEND_JOB_ID=${APPEND_JOB_ID}
