#!/usr/bin/env python3

import pandas as pd
import xarray as xr
import numpy as np
from pathlib import Path
import datetime as dt
import math
from argparse import ArgumentParser
from typing import Dict, List, Tuple, Optional, Union

class Options:
    def __init__(self, input_dir: Path, grid_file: Path, output_file: Path):
        self.input_dir = input_dir
        self.grid_file = grid_file
        self.output_file = output_file

def parse_args() -> Options:
    """Parse command line arguments.

    Returns:
        Options: Parsed and validated command line arguments
    """
    parser = ArgumentParser(
        description="Convert biomass CSV files to a single NetCDF file.",
        formatter_class=ArgumentParser.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "-i", "--input-dir",
        type=Path,
        required=True,
        help="Directory containing biomass CSV files"
    )

    parser.add_argument(
        "-g", "--grid-file",
        type=Path,
        required=True,
        help="Path to grid file containing site coordinates"
    )

    parser.add_argument(
        "-o", "--output-file",
        type=Path,
        required=True,
        help="Output NetCDF file path"
    )

    args = parser.parse_args()

    # Validate paths
    if not args.input_dir.is_dir():
        parser.error(f"Input directory does not exist: {args.input_dir}")

    if not args.grid_file.is_file():
        parser.error(f"Grid file does not exist: {args.grid_file}")

    # Create output directory if it doesn't exist
    args.output_file.parent.mkdir(parents=True, exist_ok=True)

    return Options(args.input_dir, args.grid_file, args.output_file)

def load_coordinates(grid_file: Path) -> Dict[str, Tuple[float, float]]:
    """Load site coordinates from grid file.
    
    Args:
        grid_file: Path to grid file containing site coordinates
        
    Returns:
        Dictionary mapping site names to (longitude, latitude) tuples
    """
    site_coords = {}
    with open(grid_file) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                site, lon, lat = line.strip().split()
                site_coords[site] = (float(lon), float(lat))
    return site_coords

def haversine_distance(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
    """Calculate the great circle distance between two points on earth in kilometers.
    
    Args:
        lon1: Longitude of first point in degrees
        lat1: Latitude of first point in degrees
        lon2: Longitude of second point in degrees
        lat2: Latitude of second point in degrees
        
    Returns:
        Distance between points in kilometers
    """
    R = 6371  # Earth's radius in kilometers

    # Convert to radians
    lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])

    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
    c = 2 * math.asin(math.sqrt(a))
    return R * c

def find_closest_site(target_lon: float, target_lat: float, 
                     site_coords: Dict[str, Tuple[float, float]], 
                     max_distance: float = 200) -> Tuple[Optional[str], Optional[float]]:
    """Find the closest site geographically within max_distance (km).
    
    Args:
        target_lon: Target longitude in degrees
        target_lat: Target latitude in degrees
        site_coords: Dictionary mapping site names to (longitude, latitude) tuples
        max_distance: Maximum allowed distance in kilometers
        
    Returns:
        Tuple of (site_name, distance) if a site is found within max_distance,
        otherwise (None, None)
    """
    if not site_coords:
        return None, None

    closest_site = None
    min_distance = float("inf")

    for site, (lon, lat) in site_coords.items():
        distance = haversine_distance(target_lon, target_lat, lon, lat)
        if distance < min_distance:
            min_distance = distance
            closest_site = site

    # Only return if within maximum distance
    if min_distance <= max_distance:
        return closest_site, min_distance
    return None, None

def process_biomass_files(biomass_dir: Path, grid_file: Path) -> xr.Dataset:
    """Process all biomass CSV files and create NetCDF.
    
    Args:
        biomass_dir: Directory containing biomass CSV files
        grid_file: Path to grid file containing site coordinates
        
    Returns:
        xarray Dataset containing merged biomass data
    """
    # Load coordinates
    site_coords = load_coordinates(grid_file)

    # Initialize lists to store data
    lons = []
    lats = []
    times = []
    biomass = []

    # Process each CSV file
    for csv_file in Path(biomass_dir).glob("*.csv"):
        if csv_file.name == "print_csv_columns.py":
            continue

        df = pd.read_csv(csv_file)

        # Extract coordinates from the CSV
        try:
            csv_lon = float(df["longitude"].iloc[0])
            csv_lat = float(df["latitude"].iloc[0])
        except (KeyError, ValueError) as e:
            print(f"Warning: Could not extract coordinates from {csv_file.name}: {str(e)}")
            continue

        # Find closest site based on coordinates
        closest_site, distance = find_closest_site(csv_lon, csv_lat, site_coords)
        if closest_site:
            print(f"Found site "{closest_site}" ({distance:.1f}km away) for coordinates ({csv_lon}, {csv_lat}) in {csv_file.name}")
            site_name = closest_site
        else:
            print(f"Warning: No nearby sites found within 200km for coordinates ({csv_lon}, {csv_lat}) in {csv_file.name}")
            continue

        # Get coordinates - constant for each site
        lon, lat = site_coords[site_name]

        # Process all time points in the file
        for _, row in df.iterrows():
            # Convert time string to datetime and ensure it's numpy datetime64
            time = pd.to_datetime(row["phenomenonTime"]).to_datetime64()

            # Convert biomass from T/ha to kg/m2 (multiply by 0.1)
            biomass_value = row["standBiomass_tonnesPerHectare"] * 0.1

            # Append to lists
            lons.append(lon)
            lats.append(lat)
            times.append(time)
            biomass.append(biomass_value)

    # Create xarray dataset
    ds = xr.Dataset(
        data_vars={
            "biomass": (("time", "latitude", "longitude"), np.zeros((len(set(times)), len(set(lats)), len(set(lons))))),
        },
        coords={
            "longitude": sorted(set(lons)),
            "latitude": sorted(set(lats), reverse=True),  # Reverse to go from north to south
            "time": sorted(set(times)),  # Use set to remove duplicates
        }
    )

    # Fill in biomass data
    for t, lat, lon, bio in zip(times, lats, lons, biomass):
        try:
            t_idx = ds.time.values.tolist().index(t)
            lat_idx = ds.latitude.values.tolist().index(lat)
            lon_idx = ds.longitude.values.tolist().index(lon)
            ds.biomass[t_idx, lat_idx, lon_idx] = bio
        except ValueError as e:
            print(f"Warning: Could not find index for time {t}, coordinates ({lon}, {lat}): {str(e)}")
            continue

    # Add metadata
    ds.biomass.attrs = {
        "units": "kg/m2",
        "long_name": "Stand Biomass",
        "description": "Total stand biomass"
    }

    ds.longitude.attrs = {
        "units": "degrees_east",
        "long_name": "Longitude"
    }

    ds.latitude.attrs = {
        "units": "degrees_north",
        "long_name": "Latitude"
    }

    ds.time.attrs = {
        "long_name": "Time",
        "axis": "T"
    }

    # Add global attributes
    ds.attrs = {
        "title": "Australian Site Biomass Data",
        "creation_date": dt.datetime.now().strftime("%Y-%m-%d"),
        "source": "Merged from individual site CSV files",
        "conversion": "Original biomass values in T/ha converted to kg/m2"
    }

    return ds

def main(args: Options) -> None:
    """Main application logic.
    
    Args:
        args: Parsed command line arguments
    """
    # Load coordinates and process biomass files
    ds = process_biomass_files(args.input_dir, args.grid_file)

    # Save to NetCDF file
    print(f"Saving dataset to {args.output_file}...")
    ds.to_netcdf(args.output_file)
    print("Done!")

if __name__ == "__main__":
    main(parse_args())
