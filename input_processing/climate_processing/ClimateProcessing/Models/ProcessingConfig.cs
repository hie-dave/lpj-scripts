using CommandLine;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using ClimateProcessing.Units;
using ClimateProcessing.Services;

namespace ClimateProcessing.Models;

public abstract class ProcessingConfig
{
    [Option('i', "input-dir", Required = true, HelpText = "Input directory containing climate data files")]
    public string InputDirectory { get; set; } = string.Empty;

    [Option('o', "output-dir", Required = false, HelpText = "Output directory for processed files")]
    public string OutputDirectory { get; set; } = string.Empty;

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

    [Option("jobfs", Required = false, Default = 128, HelpText = "PBS JobFS size in GiB (recommended: 100+)")]
    public int JobFS { get; set; } = 128;

    [Option('p', "project", Required = true, HelpText = "Project for job submission")]
    public string Project { get; set; } = string.Empty;

    [Option("ncpus", Required = false, Default = 1, HelpText = "Number of CPUs to request")]
    public int Ncpus { get; set; } = 1;

    [Option("compress", Default = true, HelpText = "Enable netCDF compression")]
    public bool CompressOutput { get; set; } = true;

    [Option("compression-level", Default = 5, HelpText = "Compression level (1-9)")]
    public int CompressionLevel { get; set; } = 5;

    [Option("email", Required = false, HelpText = "Email address for job notifications (optional)")]
    public string Email { get; set; } = string.Empty;

    [Option("vpd-method", Default = VPDMethod.Magnus, HelpText = "Method to calculate VPD: Magnus (default), Buck1981, AlduchovEskridge1996, AllenFAO1998, or Sonntag1990")]
    public VPDMethod VPDMethod { get; set; } = VPDMethod.Magnus;

    [Option("version", Default = ModelVersion.Dave, HelpText = "The version of LPJ-Guess by which the data is going to be used.")]
    public ModelVersion Version { get; set; } = ModelVersion.Dave;

    private TimeStep inputTimeStep = new TimeStep(0);
    private TimeStep outputTimeStep = new TimeStep(0);

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

    [Option("dry-run", Default = false, HelpText = "If set, scripts will be generated but not submitted to the PBS queue.")]
    public bool DryRun { get; set; } = false;

    public TimeStep InputTimeStep => inputTimeStep;
    public TimeStep OutputTimeStep => outputTimeStep;

    public virtual void Validate()
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

        if (Version == ModelVersion.Dave)
        {
            if (InputTimeStepHours == 0)
                throw new ArgumentException("Input timestep must be specified when using dave version");
            if (OutputTimeStepHours == 0)
                throw new ArgumentException("Output timestep must be specified when using dave version");
        }
        else if (Version == ModelVersion.Trunk)
        {
            if (InputTimeStepHours != 0)
                throw new ArgumentException("Input timestep must not be specified when using trunk version (it is ignored, as daily will always be used).");
            if (OutputTimeStepHours != 0)
                throw new ArgumentException("Output timestep must not be specified when using trunk version (it is ignored, as daily will always be used).");

            // Always use a daily timestep.
            inputTimeStep = TimeStep.Daily;
            outputTimeStep = TimeStep.Daily;
        }

        if (InputTimeStepHours > OutputTimeStepHours)
            throw new ArgumentException("Input timestep cannot be coarser than the output timestep.");
    }

    /// <summary>
    /// Creates the datasets that will be used to generate the script.
    /// </summary>
    public abstract IEnumerable<NarClim2Dataset> CreateDatasets();

    /// <summary>
    /// Creates the script generator used to generate the script.
    /// </summary>
    public abstract ScriptGenerator CreateScriptGenerator();
}
