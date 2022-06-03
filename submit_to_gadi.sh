#!/usr/bin/env bash
#
# Portable bash script to run LPJ-Guess as a parallel job using PBS on Gadi.
#
# Adapted by Juergen Knauer and Drew Holzworth from older scripts.
#
# Usage:
#
# submit_to_gadi.sh [-n <name>] [-s <file>] [-i <ins-file>]
#
# All arguments are optional and interpreted as:
#
#   name     = the name of the job (shown in PBS queue)
#   file     = filename of a file which can override the input variables. If
#              set, this file should contain lines of the form VARIABLE=Value
#   ins-file = instruction file to use, overrides the hardcoded INSFILE
#

# Exit immediately if:
# - Any command fails
# - Any part of a pipeline fails
# - We reference any unbound variables
set -euo pipefail

# Get directory containing this script.
SUBMIT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd -P )"

# ------------------------------------------------------------------------------
# Input variables - change these as required.
# ------------------------------------------------------------------------------

# File paths
OUT_DIR="/scratch/hw83/jk8585/LPJ_runs" # Outut path.
BINARY="/home/599/jk8585/LPJ_code/trunk_r8538/guess" # Path to LPJ-Guess binary.
EXPERIMENT="test_run"               # Experiment name.

# PBS/Cluster Settings
NPROCESS=48                         # Number of processes in parallel job.
WALLTIME=22:30:00                   # maximum wall (real) time for job hh:mm:ss
MEMORY=96GB                         # Total amount of memory available to job.
QUEUE=normal                        # PBS queue priority(?).
PROJECT=hw83                        # Project name.
EMAIL=j.knauer@westernsydney.edu.au # Email address of job owner.
JOB_NAME=guess                      # Job name (overriden by -n CLI option).

# LPJ-Guess Settings
INSFILE="${SUBMIT_DIR}/global_cf.ins" # path to ins file from run directory.
INPUT_MODULE=cf                     # Input module to use.
OUTFILES='*.out'                    # List of LPJ-Guess output files.

# ------------------------------------------------------------------------------
# End input variables section.
# ------------------------------------------------------------------------------

# Parse the command line arguments
while getopts ":n:s:i:" opt; do
    case $opt in
	n ) JOB_NAME=$OPTARG ;;
	s ) source $OPTARG ;;
	i ) INSFILE=$OPTARG ;;
    esac
done

# This function will print the absolute path to the given filename.
# We could use readlink, but it's a Linux tool - not available on MacOS (afaik).
get_absolute_path() {
  echo $( cd "$(dirname "$1")"; pwd -P )/$(basename "$1")
}

# Convert INSFILE to an absolute path since we will be starting the
# guess instances from different directories.
INSFILE="$(get_absolute_path "${INSFILE}")"

# Read gridlist file name from the .ins file.
GRIDLIST="$(sed -E -n -e 's/^.*"file_gridlist_cf" \(str "([^"]+)"\)\r?$/\1/p' "${INSFILE}")"
GRIDLIST="$(get_absolute_path "${GRIDLIST}")"

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
echo JOB_NAME=${JOB_NAME}
echo
echo INSFILE=${INSFILE}
echo INPUT_MODULE=${INPUT_MODULE}
echo GRIDLIST=${GRIDLIST}
echo

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

check_exists() { check_permission -f $1 "Error: $1 does not exist"; }
check_read() { check_permission -r $1 "Error: $1 does not have read permission"; }
check_execute() { check_permission -x $1 "Error: $1 does not have execute permission"; }

check_exists "${BINARY}"
check_exists "${INSFILE}"
check_exists "${GRIDLIST}"
check_read "${BINARY}"
check_read "${INSFILE}"
check_read "${GRIDLIST}"
check_execute "${BINARY}"

# Output path for this run.
RUN_OUT_DIR="${OUT_DIR}/${EXPERIMENT}"

# Create links and directories.
mkdir -p "${RUN_OUT_DIR}/output"
mkdir -p "${RUN_OUT_DIR}/logs"

# Get the file name (without path) of the gridlist file.
GRIDLIST_FILENAME=$(basename ${GRIDLIST})

# This function creates the gridlist files for each run by splitting
# the original gridlist file into approximately equal parts.
function split_gridlist {
    # Create empty gridlists first to make sure each run gets one.
    for ((a=1; a <= NPROCESS ; a++))
    do
        mkdir -p "${RUN_OUT_DIR}/run${a}"
        echo > "${RUN_OUT_DIR}/run${a}/${GRIDLIST_FILENAME}"
    done

    # Figure out suitable number of lines per gridlist, get the number of
    # lines in original gridlist file, divide by NPROCESS and round up.
    local lines_per_run=$(wc -l ${GRIDLIST} | \
	awk '{ x = $1/'${NPROCESS}'; d = (x == int(x)) ? x : int(x)+1; print d}')

    # Use the split command to split the files into temporary files
    local tmp_prefix=tmpSPLITGRID_
    split --suffix-length=4 --lines ${lines_per_run} ${GRIDLIST} "${tmp_prefix}"

    # Move the temporary files into the runX-directories.
    local files="${tmp_prefix}*"
    local i=1
    for file in ${files}
    do
      mv ${file} "${RUN_OUT_DIR}/run${i}/${GRIDLIST_FILENAME}"
      i=$((i+1))
    done
}

# Create header of progress.sh script
# The progress script will print the progress of each run/process.
progress_sh="${RUN_OUT_DIR}/progress.sh"
echo "##############################################################" > "${progress_sh}"
echo "# PROGRESS.SH" >> "${progress_sh}"
echo "# Upload current guess.log files from local nodes and check" >> "${progress_sh}"
echo "# Usage: sh progress.sh" >> "${progress_sh}"
echo >> "${progress_sh}"
echo 'DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd -P )"' >> "${progress_sh}"
echo >> "${progress_sh}"

# Create a run subdirectory for each process, clean up any existing files, and
# append code to the progress script which reports progress of this run.
for ((a=1; a <= NPROCESS ; a++))
do
    # Ensure the directory exists for this run.
    mkdir -p "${RUN_OUT_DIR}/run${a}"

    # Delete any existing log files.
    cd "${RUN_OUT_DIR}/run${a}"; rm -f guess.log; rm -f ${GRIDLIST_FILENAME}; cd ..
    echo "echo '********** Last few lines of ./run${a}/guess.log: **********'" >> "${progress_sh}"
    echo "tail \${DIR}/run${a}/guess.log" >> "${progress_sh}"
done

split_gridlist

# Create PBS script to request place in queue
guess_cmd="${RUN_OUT_DIR}/guess.cmd"
cat <<EOF > "${guess_cmd}"
#!/bin/bash
#PBS -l ncpus=${NPROCESS}
#PBS -l walltime=${WALLTIME}
#PBS -l mem=${MEMORY}
#PBS -q ${QUEUE}
#PBS -P ${PROJECT}
#PBS -m n
#PBS -M ${EMAIL}
#PBS -j oe
#PBS -l wd
#PBS -W umask=0022
#PBS -S /bin/bash
#PBS -l storage=scratch/pt17
set -e
module purge
module load intel-mpi/2019.6.166
umask 022
mpirun ${BINARY} -parallel -input ${INPUT_MODULE} ${INSFILE}
EOF

append_cmd="${RUN_OUT_DIR}/append.cmd"
cat <<EOF > "${append_cmd}"
#!/bin/bash
#PBS -l ncpus=1
#PBS -l walltime=${WALLTIME}
#PBS -l mem=${MEMORY}
#PBS -q normal
#PBS -l wd
#PBS -j oe
#PBS -m n
#PBS -W umask=0022
set -e
function append_files {
    local number_of_jobs=\$1
    local file=\$2
    cp run1/\$file \$file
    local i=""
    for ((i=2; i <= number_of_jobs; i++))
    do
      if [ -f run\$i/\$file ]; then
        cat run\$i/\$file | awk 'NR!=1 || NF==0 || \$1 == \$1+0 { print \$0 }' >> \$file
      fi
    done
}
pushd run1 &> /dev/null
outfiles_unexpanded='${OUTFILES}'
outfiles_expanded=\$(echo \$outfiles_unexpanded)
popd &> /dev/null
for file in \$outfiles_expanded
do
  append_files ${NPROCESS} \$file
done
cat run*/guess.log > guess.log
EOF

# Submit guess job
append_dependency=$(qsub -N "${JOB_NAME}" "${guess_cmd}")
echo JOB_ID=${append_dependency}

# Submit append job
qsub -W depend=afterok:${append_dependency} -N "${JOB_NAME}_append" "${append_cmd}"
