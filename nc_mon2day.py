#!/usr/bin/env python3

from argparse import ArgumentParser
from attr import dataclass
import xarray as xr

@dataclass(frozen = True)
class Options:
    in_file: str
    out_file: str
    ref_file: str

def parse_args() -> Options:
    parser = ArgumentParser(description="Convert monthly ps to daily using daily time axis.")
    parser.add_argument("--in-file", required=True, help="Input monthly file.")
    parser.add_argument("--out-file", required=True, help="Output daily file.")
    parser.add_argument("--ref-file", required=True, help="Reference daily file.")
    args = parser.parse_args()
    return Options(
        in_file=args.in_file,
        out_file=args.out_file,
        ref_file=args.ref_file,
    )

def main(opts: Options):
    ds_in = xr.open_mfdataset(opts.in_file, combine="by_coords")
    ds_ref = xr.open_mfdataset(opts.ref_file, combine="by_coords")

    ds_daily = ds_in.interp(time=ds_ref.time).ffill("time").bfill("time")

    ds_daily.attrs["history"] = (
        ds_daily.attrs.get("history", "")
        + "\nMonthly ps linearly interpolated to daily using tasmax daily time axis."
    )

    ds_daily.to_netcdf(opts.out_file)

if __name__ == "__main__":
    opts = parse_args()
    main(opts)
