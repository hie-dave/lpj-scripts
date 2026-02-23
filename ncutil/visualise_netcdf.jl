#!/usr/bin/env -S julia --project=@.

using ArgParse
using Dates
using NCDatasets
using PlotlyJS
using Logging

include("nc_helpers.jl")

struct Options
    filename::String
    variable::String
    lat::Float64
    lon::Float64
    log_level::Logging.LogLevel
    output::Union{Nothing, String}
    width::Int
    height::Int
end

function parse_log_level(level::Int)::Logging.LogLevel
    if level == 0
        return Logging.Error
    elseif level == 1
        return Logging.Warn
    elseif level == 2
        return Logging.Info
    elseif level == 3
        return Logging.Debug
    else
        error("Invalid log level: $level. Must be between 0 and 3.")
    end
end

function parse_cli()::Options
    parser = ArgParseSettings()
    @add_arg_table parser begin
        "--filename", "-f"
            help = "Path to the NetCDF file"
            required = true
            arg_type = String
        "--variable", "-v"
            help = "Name of the variable to visualize"
            required = true
            arg_type = String
        "--lat", "--latitude"
            help = "Latitude for the point to visualize"
            required = true
            arg_type = Float64
        "--lon", "--longitude"
            help = "Longitude for the point to visualize"
            required = true
            arg_type = Float64
        "-l", "--log-level"
            help = "Logging level (e.g., 0=error, 1=warn, 2=info, 3=debug)"
            required = false
            arg_type = Int
            default = 2
        "-o", "--output"
            help = "Write plot to file (e.g. out.png, out.pdf). If omitted, show plot interactively."
            required = false
            arg_type = String
        "--width"
            help = "Width of the output plot (if saving to a raster) in pixels (default: 800)"
            required = false
            arg_type = Int
            default = 800
        "--height"
            help = "Height of the output plot (if saving to a raster) in pixels (default: 600)"
            required = false
            arg_type = Int
            default = 600
    end
    p = parse_args(parser)
    level = parse_log_level(p["log-level"])
    return Options(p["filename"], p["variable"], p["lat"], p["lon"], level,
                   get(p, "output", nothing), p["width"], p["height"])
end

function plot_timeseries(ds::NCDataset, variable::String, lat::Float64,
                         lon::Float64)
    var = ds[variable]
    var_time = get_time_var(ds)
    times = reinterpret.(DateTime, var_time[:])
    data, alat, alon = read_timeseries(var, lat, lon)
    xlab = "Date"
    units = haskey(var.attrib, "units") ? var.attrib["units"] : "unknown units"
    ylab = "$variable ($units)"
    title = "Timeseries of $variable at (lat: $alat, lon: $alon)"
    return plot(times, data, Layout(xaxis_title = xlab,
                                             yaxis_title = ylab,
                                             title = title))
end

function main(opts::Options)
    logger = ConsoleLogger(stderr, opts.log_level)
    global_logger(logger)
    NCDataset(opts.filename) do ds
        plt = plot_timeseries(ds, opts.variable, opts.lat, opts.lon)
        if isnothing(opts.output)
            @info "Displaying plot interactively"
            display(plt)
            println("Press enter to exit...")
            readline()
        else
            @info "Saving plot to $(opts.output)"
            savefig(plt, opts.output, width = opts.width, height = opts.height)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    opts = parse_cli()
    main(opts)
end
