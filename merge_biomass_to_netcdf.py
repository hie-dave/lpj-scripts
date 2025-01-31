#!/usr/bin/env python3

import pandas as pd
import xarray as xr
import numpy as np
from pathlib import Path
import datetime as dt
import math

def load_coordinates(grid_file):
    """Load site coordinates from grid file."""
    coords = {}
    with open(grid_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 3:
                lon, lat, site = parts
                coords[site] = (float(lon), float(lat))
    return coords

def haversine_distance(lon1, lat1, lon2, lat2):
    """Calculate the great circle distance between two points on earth in kilometers."""
    R = 6371  # Earth's radius in kilometers

    # Convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])

    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
    c = 2 * math.asin(math.sqrt(a))
    return R * c

def find_closest_site(target_lon, target_lat, site_coords, max_distance=200):
    """Find the closest site geographically within max_distance (km)."""
    if not site_coords:
        return None
        
    closest_site = None
    min_distance = float('inf')
    
    for site, (lon, lat) in site_coords.items():
        distance = haversine_distance(target_lon, target_lat, lon, lat)
        if distance < min_distance:
            min_distance = distance
            closest_site = site
    
    # Only return if within maximum distance
    if min_distance <= max_distance:
        return closest_site, min_distance
    return None, None

def clean_site_name(name):
    """Convert CSV site names to match grid file format."""
    # Remove common words and special characters
    replacements = {
        'Mulga': '',
        'Savanna': '',
        'Woodland': '',
        'Woodlands': '',
        'Plain': '',
        'Wet': '',
        'Dry': '',
        'Tall': '',
        'Stringbark': '',
        '_': '',
        '-': '',
        ' ': '',
    }
    clean = name
    for old, new in replacements.items():
        clean = clean.replace(old, new)
    return clean

def process_biomass_files(biomass_dir, grid_file):
    """Process all biomass CSV files and create NetCDF."""
    # Load coordinates
    site_coords = load_coordinates(grid_file)
    
    # Initialize lists to store data
    lons = []
    lats = []
    times = []
    biomass = []
    
    # Process each CSV file
    for csv_file in Path(biomass_dir).glob('*.csv'):
        if csv_file.name == 'print_csv_columns.py':
            continue
            
        df = pd.read_csv(csv_file)
        
        # Extract site name and coordinates from the CSV
        site_name = clean_site_name(df['siteId'].iloc[0])
        try:
            csv_lon = float(df['lon'].iloc[0])
            csv_lat = float(df['lat'].iloc[0])
        except (KeyError, ValueError) as e:
            print(f"Warning: Could not extract coordinates from {csv_file.name}: {str(e)}")
            continue
        
        # Look up coordinates from grid file
        if site_name not in site_coords:
            closest_site, distance = find_closest_site(csv_lon, csv_lat, site_coords)
            if closest_site:
                print(f"Warning: Using closest site '{closest_site}' ({distance:.1f}km away) for coordinates in {csv_file.name}")
                site_name = closest_site
            else:
                print(f"Warning: No nearby sites found within 200km for coordinates ({csv_lon}, {csv_lat}) in {csv_file.name}")
                continue
            
        # Get coordinates - constant for each site
        lon, lat = site_coords[site_name]
        
        # Process all time points in the file
        for _, row in df.iterrows():
            # Convert time string to datetime
            time = pd.to_datetime(row['phenomenonTime'])
            
            # Convert biomass from T/ha to kg/m2 (multiply by 0.1)
            biomass_value = row['standBiomass_tonnesPerHectare'] * 0.1
            
            # Append to lists
            lons.append(lon)
            lats.append(lat)
            times.append(time)
            biomass.append(biomass_value)
    
    # Create xarray dataset
    ds = xr.Dataset(
        data_vars={
            'biomass': (('time', 'latitude', 'longitude'), np.zeros((len(times), len(set(lats)), len(set(lons))))),
        },
        coords={
            'longitude': sorted(set(lons)),
            'latitude': sorted(set(lats), reverse=True),  # Reverse to go from north to south
            'time': sorted(times),
        }
    )
    
    # Fill in biomass data
    for t, lat, lon, bio in zip(times, lats, lons, biomass):
        t_idx = ds.time.values.tolist().index(np.datetime64(t))
        lat_idx = ds.latitude.values.tolist().index(lat)
        lon_idx = ds.longitude.values.tolist().index(lon)
        ds.biomass[t_idx, lat_idx, lon_idx] = bio
    
    # Add metadata
    ds.biomass.attrs = {
        'units': 'kg/m2',
        'long_name': 'Stand Biomass',
        'description': 'Total stand biomass converted from tonnes per hectare'
    }
    
    ds.longitude.attrs = {
        'units': 'degrees_east',
        'long_name': 'Longitude'
    }
    
    ds.latitude.attrs = {
        'units': 'degrees_north',
        'long_name': 'Latitude'
    }
    
    ds.time.attrs = {
        'long_name': 'Time',
        'axis': 'T'
    }
    
    # Add global attributes
    ds.attrs = {
        'title': 'Australian Site Biomass Data',
        'creation_date': dt.datetime.now().strftime('%Y-%m-%d'),
        'source': 'Merged from individual site CSV files',
        'conversion': 'Original biomass values in T/ha converted to kg/m2'
    }
    
    return ds

if __name__ == "__main__":
    biomass_dir = Path('/home/drew/code/lpj-guess/scripts/ozflux-lpjg/data/biomass/new')
    grid_file = Path('/home/drew/code/lpj-guess/output-analysis/inst/extdata/ozflux.grid')
    
    # Process files and create dataset
    ds = process_biomass_files(biomass_dir, grid_file)
    
    # Save to NetCDF
    output_file = biomass_dir / 'merged_biomass.nc'
    ds.to_netcdf(output_file)
    print(f"Created NetCDF file: {output_file}")
