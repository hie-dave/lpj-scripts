#!/usr/bin/env python3

from argparse import ArgumentParser
from attr import dataclass
import xarray as xr
import pandas as pd

@dataclass(frozen = True)
class Options:
    in_file: str
    out_file: str

def parse_args() -> Options:
    parser = ArgumentParser(description="Convert monthly ps to daily using daily time axis.")
    parser.add_argument("--in-file", required=True, help="Input monthly file.")
    parser.add_argument("--out-file", required=True, help="Output daily file.")
    args = parser.parse_args()
    return Options(in_file=args.in_file, out_file=args.out_file)

def main(opts: Options):
    ds_in = xr.open_mfdataset(opts.in_file, combine="by_coords")

    # Decode first/last monthly timestamps
    t0 = pd.Timestamp(ds_in.time.values[0])
    t1 = pd.Timestamp(ds_in.time.values[-1])

    # Daily axis from start of first month to end of last month
    start = t0.to_period("M").start_time
    end = t1.to_period("M").end_time.floor("D")

    daily_time = pd.date_range(start=start, end=end, freq="D")

    ds_daily = ds_in.interp(time=daily_time).ffill("time").bfill("time")

    ds_daily.attrs["history"] = (
        ds_daily.attrs.get("history", "")
        + "\nMonthly ps linearly interpolated to daily using tasmax daily time axis."
    )

    ds_daily.to_netcdf(opts.out_file)

if __name__ == "__main__":
    opts = parse_args()
    main(opts)
