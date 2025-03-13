import netCDF4 as nc
import numpy as np

# Create a 1-degree global grid
nx, ny = 360, 180
lon = np.linspace(-179.5, 179.5, nx)
lat = np.linspace(-89.5, 89.5, ny)

# Create the NetCDF file
with nc.Dataset("1deg.grid", "w", format="NETCDF4") as ds:
    # Create dimensions
    ds.createDimension("lon", nx)
    ds.createDimension("lat", ny)
    
    # Create variables
    lon_var = ds.createVariable("lon", "f8", ("lon",))
    lat_var = ds.createVariable("lat", "f8", ("lat",))
    
    # Create a dummy variable to make CDO happy
    data_var = ds.createVariable("topo", "f8", ("lat", "lon"))
    
    # Assign values
    lon_var[:] = lon
    lat_var[:] = lat
    data_var[:] = np.zeros((ny, nx))
    
    # Add attributes
    lon_var.units = "degrees_east"
    lat_var.units = "degrees_north"
    lon_var.standard_name = "longitude"
    lat_var.standard_name = "latitude"
    lon_var.axis = "X"
    lat_var.axis = "Y"
    
    # Add global attributes
    ds.Conventions = "CF-1.6"
    ds.gridtype = "lonlat"
