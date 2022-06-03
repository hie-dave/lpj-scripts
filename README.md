# lpj-scripts

Scripts for running LPJ-Guess on Gadi.

## Usage

`submit_to_gadi.sh` is the main job submission script. Running this script will submit the LPJ-Guess job using PBS. The submit script supports a few optional arguments:

- `-n <job-name>`: Job name
- `-s <config-file>`: Path to a file which contains job configuration
- `-i <ins-file>`: Path to the .ins file to be used by LPJ guess(optional - may be set in the config file)

## Configuration

The submit script contains several parameters which may be overriden for a particular job. These may be modified in the script, but the recommended way of setting them is to create a config file which sets the parameter values, and pass this config file to the submit script via the `-s` command-line argument.

```bash
./submit_to_gadi.sh -s settings.ini
```

The config file will be sourced by the submit script, so the variables in the config file should be set using standard bash syntax (e.g. `NAME=Value`). An example configuration file is provided in this repository (`example.conf`).
