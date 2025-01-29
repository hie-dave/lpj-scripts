using CommandLine;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using ClimateProcessing.Units;

namespace ClimateProcessing.Models;

public class ProcessingConfig
{
    [Option('i', "input-dir", Required = true, HelpText = "Input directory containing climate data files")]
    public string InputDirectory { get; set; } = string.Empty;

    [Option('o', "output-dir", Required = false, HelpText = "Output directory for processed files")]
    public string OutputDirectory { get; set; } = string.Empty;

    [Option('d', "dataset-type", Required = true, HelpText = "Type of climate dataset (e.g., narclim2)")]
    public string DatasetType { get; set; } = string.Empty;

    [Option('g', "grid-file", Required = false, HelpText = "Path to target grid file for spatial remapping")]
    public string? GridFile { get; set; }

    [Option('c', "chunk-size", Required = false, HelpText = "Chunk size for time variable")]
    public int ChunkSizeTime { get; set; } = 8760;

    [Option("chunk-size-spatial", Required = false, HelpText = "Chunk size for spatial variables. Default is 1, and this should always be used unless you really know what you're doing or enjoy inefficiency.")]
    public int ChunkSizeSpatial { get; set; } = 1;

    [Option('q', "queue", Required = false, Default = "normal", HelpText = "PBS queue to which job will be submitted")]
    public string Queue { get; set; } = "normal";

    [Option('w', "walltime", Required = false, Default = "02:00:00", HelpText = "Walltime limit (HH:MM:SS)")]
    public string Walltime { get; set; } = "12:00:00";

    [Option('m', "memory", Required = false, Default = 16, HelpText = "Memory limit in GiB")]
    public int Memory { get; set; } = 16;

    [Option('p', "project", Required = true, HelpText = "Project for job submission")]
    public string Project { get; set; } = string.Empty;

    [Option('n', "job-name", Required = false, Default = "process_climate", HelpText = "Job name")]
    public string JobName { get; set; } = "process_climate";

    [Option("ncpus", Required = false, Default = 4, HelpText = "Number of CPUs to request")]
    public int Ncpus { get; set; } = 4;

    [Option("compress", Default = true, HelpText = "Enable netCDF compression")]
    public bool CompressOutput { get; set; } = true;

    [Option("compression-level", Default = 5, HelpText = "Compression level (1-9)")]
    public int CompressionLevel { get; set; } = 5;

    [Option("vpd-method", Default = VPDMethod.Magnus, HelpText = "Method to calculate VPD: Magnus (default), Buck1981, AlduchovEskridge1996, AllenFAO1998, or Sonntag1990")]
    public VPDMethod VPDMethod { get; set; } = VPDMethod.Magnus;

    private TimeStep inputTimeStep = TimeStep.Hourly;
    private TimeStep outputTimeStep = TimeStep.ThreeHourly;

    [Option("input-timestep", HelpText = "Input time step in hours")]
    public int InputTimeStepHours
    {
        get => inputTimeStep.Hours;
        set => inputTimeStep = new TimeStep(value);
    }

    [Option("output-timestep", HelpText = "Output time step in hours")]
    public int OutputTimeStepHours
    {
        get => outputTimeStep.Hours;
        set => outputTimeStep = new TimeStep(value);
    }

    public TimeStep InputTimeStep => inputTimeStep;
    public TimeStep OutputTimeStep => outputTimeStep;

    private static readonly Dictionary<ClimateVariable, (string units, AggregationMethod aggregation)> DefaultVariableConfig = new()
    {
        { ClimateVariable.Temperature, ("degC", AggregationMethod.Mean) },
        { ClimateVariable.Precipitation, ("mm", AggregationMethod.Sum) },
        { ClimateVariable.SpecificHumidity, ("1", AggregationMethod.Mean) },
        { ClimateVariable.SurfacePressure, ("Pa", AggregationMethod.Mean) },
        { ClimateVariable.ShortwaveRadiation, ("W m-2", AggregationMethod.Mean) },
        { ClimateVariable.WindSpeed, ("m s-1", AggregationMethod.Mean) }
    };

    public string GetTargetUnits(ClimateVariable variable)
    {
        if (!DefaultVariableConfig.TryGetValue(variable, out var config))
        {
            throw new ArgumentException($"No configuration found for variable {variable}");
        }
        return config.units;
    }

    public AggregationMethod GetAggregationMethod(ClimateVariable variable)
    {
        if (!DefaultVariableConfig.TryGetValue(variable, out var config))
        {
            throw new ArgumentException($"No configuration found for variable {variable}");
        }
        return config.aggregation;
    }

    public IEnumerable<PBSStorageDirective> GetRequiredStorageDirectives()
    {
        var paths = new List<string> { InputDirectory };
        if (!string.IsNullOrEmpty(OutputDirectory))
        {
            paths.Add(OutputDirectory);
        }
        return PBSStorageHelper.GetStorageDirectives(paths);
    }

    public void Validate()
    {
        if (!Directory.Exists(InputDirectory))
            throw new ArgumentException($"Input directory does not exist: '{InputDirectory}'");

        // Create output directory if it doesn't already exist.
        if (!string.IsNullOrEmpty(OutputDirectory))
            Directory.CreateDirectory(OutputDirectory);

        if (string.IsNullOrEmpty(Project))
            throw new ArgumentException("Project must be specified");

        if (CompressOutput && (CompressionLevel < 1 || CompressionLevel > 9))
            throw new ArgumentException("Compression level must be between 1 and 9");

        if (Memory < 0)
            throw new ArgumentException("Memory must be greater than 0");

        if (ChunkSizeTime < 1)
            throw new ArgumentException("Chunk size must be greater than 0");

        if (ChunkSizeSpatial < 1)
            throw new ArgumentException("Chunk size for spatial variables must be greater than 0");

        if (!string.IsNullOrEmpty(GridFile) && !File.Exists(GridFile))
            throw new ArgumentException($"Grid file does not exist: {GridFile}");
    }
}
