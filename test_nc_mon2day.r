library(ncdf4)
library("data.table")
library("dplyr")
library("CFtime")
library("ggplot2")

dir <- file.path(getwd(), "ozflux-lpjg/data/mon2day")

file_day <- file.path(dir, "ps_day.nc")
file_mon <- file.path(dir, "ps_Amon_EC-Earth3_ssp245_r1i1p1f1_gr_210001-210012.nc")

index_of_closest <- function(needle, haystack) {
    return(which.min(abs(haystack - needle)))
}

read_nc <- function(file, var = "ps", lon = NULL, lat = NULL) {
    nc <- nc_open(file)

    calendar <- ncatt_get(nc, "time", "calendar")$value
    tunits <- ncatt_get(nc, "time", "units")$value

    values <- ncvar_get(nc, var)
    times <- ncvar_get(nc, "time")

    dim_names <- sapply(nc$var[[var]]$dim, function(dim) dim$name)
    dim_values <- lapply(dim_names, function(dim_name) ncvar_get(nc, dim_name))
    names(dim_values) <- dim_names

    nc_close(nc)

    dates <- as.Date(as_timestamp(CFtime(tunits, calendar, times)))
    if ("time" %in% names(dim_values)) {
        dim_values$time <- dates
    }

    # Flatten
    dt <- data.table(expand.grid(dim_values),
                     value = as.vector(values))

    if (!is.null(lon) && !is.null(lat)) {
        target_lon <- dt$lon[index_of_closest(lon, dt$lon)]
        target_lat <- dt$lat[index_of_closest(lat, dt$lat)]

        dt <- dt[lon == target_lon & lat == target_lat]
    }

    return(dt)
}

lon <- 150
lat <- -33

dt_day <- read_nc(file_day, lon = lon, lat = lat)
dt_mon <- read_nc(file_mon, lon = lon, lat = lat)

ggplot() +
    geom_line(data = dt_day, aes(x = time, y = value), linewidth = 0.4) +
    geom_point(data = dt_mon, aes(x = time, y = value), size = 2) +
    labs(x = "Date", y = "Surface pressure", title = "Daily and monthly surface pressure") +
    theme_classic()
