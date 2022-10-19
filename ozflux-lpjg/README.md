## ozflux processing scripts

This directory contains scripts for processing ozflux met data:

### Downloading the data

```bash
python -m pip install -r requirements.txt
python ./download-ozflux <out-dir> # or just ./download-ozflux <out-dir>
```

Run the script with 1 argument: path to a directory. All files will be
downloaded into this directory. It will be created if it doesn't exist.

### Converting to format suitable for lpj-guess (lsminput)

```bash
python -m pip install -r requirements.txt
python ./met-processing --help # or just ./met-processing --help
```

This script has a few more arguments. Run with `--help` to view them.

- Using the `--parallel` option requires 1.6.1 of the netcdf4 python module
- Note: I think `--parallel` is not working correctly. Will revisit this
  when possible, but for now, I'd advise against using it.

### Generating the benchmarks

```bash
python -m pip install -r requirements.txt
python ./create-benchmarks --help # or just ./create-benchmarks --help
```

This script requires a number of arguments. Run with `--help` to view them.
