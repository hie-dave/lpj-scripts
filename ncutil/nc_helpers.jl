using NCDatasets

# Name of the "standard_name" attribute, as per the CF conventions.
const ATTR_STD_NAME = "standard_name"

# Standard name for latitude, as per the CF conventions.
const STD_LAT = "latitude"

# Standard name for longitude, as per the CF conventions.
const STD_LON = "longitude"

# Standard name for time, as per the CF conventions.
const STD_TIME = "time"

# Encapsulates the dimension order of a 3D variable. Note that this is Julia
# order, not "real" (C) order, and is also 1-indexed.
struct DimensionOrder
    index_lon::Int
    index_lat::Int
    index_time::Int
end

# Get a variable from the dataset by its "standard_name" attribute.
# @param nc The NetCDF dataset to search.
# @param std_name The value of the "standard_name" attribute to look for.
function var_from_std_name(nc::NCDataset,
                           std_name::String)::NCDatasets.CFVariable
    for (_, var) in values(nc)
        if haskey(var.attrib, ATTR_STD_NAME) &&
                  var.attrib[ATTR_STD_NAME] == std_name
            return var
        end
    end
    error("No variable with 'standard_name' attribute '$std_name' found in the dataset.")
end

# Get the latitude variable from a CF-compliant dataset.
function get_lat_var(nc::NCDataset)::NCDatasets.CFVariable
    return var_from_std_name(nc, STD_LAT)
end

# Get the longitude variable from a CF-compliant dataset.
function get_lon_var(nc::NCDataset)::NCDatasets.CFVariable
    return var_from_std_name(nc, STD_LON)
end

# Get the time variable from a CF-compliant dataset.
function get_time_var(nc::NCDataset)::NCDatasets.CFVariable
    return var_from_std_name(nc, STD_TIME)
end

# Determine the dimension order of a 3D variable.
function get_dim_order(var::NCDatasets.CFVariable)::DimensionOrder
    if length(dimnames(var)) != 3
        error("Variable $(name(var)) is not 3D: $(dimnames(var))")
    end

    var_lon = get_lon_var(var.var.ds)
    var_lat = get_lat_var(var.var.ds)
    var_time = get_time_var(var.var.ds)

    if ndims(var_lat) != 1
        error("Latitude variable $(name(var_lat)) is not 1D.")
    end
    if ndims(var_lon) != 1
        error("Longitude variable $(name(var_lon)) is not 1D.")
    end
    if ndims(var_time) != 1
        error("Time variable $(name(var_time)) is not 1D.")
    end

    dim_lon = dimnames(var_lon)[1]
    dim_lat = dimnames(var_lat)[1]
    dim_time = dimnames(var_time)[1]

    index_lon = findfirst(==(dim_lon), dimnames(var))
    index_lat = findfirst(==(dim_lat), dimnames(var))
    index_time = findfirst(==(dim_time), dimnames(var))

    return DimensionOrder(index_lon, index_lat, index_time)
end

# Read a time series for a specific lat/lon index from a 3D variable.
function read_timeseries(var::NCDatasets.CFVariable, ilat::Int,
                         ilon::Int)::AbstractVector{AbstractFloat}
    dims = get_dim_order(var)
    hyperslab = selectdim(var, dims.index_lat, ilat)
    hyperslab = selectdim(hyperslab, dims.index_lon, ilon)
    return hyperslab[:]
end

# Read a time series for the closest grid point to a specified lat/lon from a 3D
# variable. Returns a tuple of the data, as well as the actual latitude and
# longitude used.
function read_timeseries(var::NCDatasets.CFVariable, lat::AbstractFloat,
                         lon::AbstractFloat)::Tuple{AbstractVector{AbstractFloat}, AbstractFloat, AbstractFloat}
    lat_var = get_lat_var(var.var.ds)
    lon_var = get_lon_var(var.var.ds)

    # Find the indices of the closest grid point to the specified lat/lon.
    lat_idx = findmin(abs.(lat_var[:] .- lat))[2]
    lon_idx = findmin(abs.(lon_var[:] .- lon))[2]

    data = read_timeseries(var, lat_idx, lon_idx)
    return data, lat_var[lat_idx], lon_var[lon_idx]
end
