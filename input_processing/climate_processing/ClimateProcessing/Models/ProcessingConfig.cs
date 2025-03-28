using CommandLine;
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

    [Option("email-notification", Required = false, Default = EmailNotificationType.None, HelpText = "Email notification types (can be combined): none, before, after, error")]
    public EmailNotificationType EmailNotifications { get; set; } = EmailNotificationType.None;

    [Option("vpd-method", Default = VPDMethod.Magnus, HelpText = "Method to calculate VPD: Magnus (default), Buck1981, AlduchovEskridge1996, AllenFAO1998, or Sonntag1990")]
    public VPDMethod VPDMethod { get; set; } = VPDMethod.Magnus;

    [Option("version", Default = ModelVersion.Dave, HelpText = "The version of LPJ-Guess by which the data is going to be used.")]
    public ModelVersion Version { get; set; } = ModelVersion.Dave;

    [Option("input-timestep", HelpText = "Input time step in hours")]
    public int InputTimeStepHours { get; set; }

    [Option("output-timestep", HelpText = "Output time step in hours")]
    public int OutputTimeStepHours { get; set; }

    [Option("dry-run", Default = false, HelpText = "If set, scripts will be generated but not submitted to the PBS queue.")]
    public bool DryRun { get; set; } = false;

    public TimeStep InputTimeStep => new TimeStep(InputTimeStepHours);
    public TimeStep OutputTimeStep => new TimeStep(OutputTimeStepHours);

    public virtual void Validate()
    {
        ValidateDirectories();
        ValidateBasicParameters();
        ValidateTimeStepSettings();
        ValidateEmailSettings();
    }

    internal virtual void ValidateDirectories()
    {
        if (!Directory.Exists(InputDirectory))
            throw new ArgumentException($"Input directory does not exist: '{InputDirectory}'");

        if (!string.IsNullOrEmpty(GridFile) && !File.Exists(GridFile))
            throw new ArgumentException($"Grid file does not exist: '{GridFile}'");
    }

    internal virtual void ValidateBasicParameters()
    {
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
    }

    internal virtual void ValidateTimeStepSettings()
    {
        if (Version == ModelVersion.Dave)
            ValidateDaveTimeStepSettings();
        else if (Version == ModelVersion.Trunk)
            ValidateTrunkTimeStepSettings();
        else
            throw new ArgumentException($"Invalid version: {Version}");

        if (InputTimeStepHours > OutputTimeStepHours)
            throw new ArgumentException("Input timestep cannot be coarser than the output timestep.");
    }

    internal virtual void ValidateTrunkTimeStepSettings()
    {
        if (InputTimeStepHours != 24 && InputTimeStepHours != 0)
            throw new ArgumentException("Input timestep must be daily (24 hours) or not specified when processing for trunk.");
        if (OutputTimeStepHours != 24 && OutputTimeStepHours != 0)
            throw new ArgumentException("Output timestep must be daily (24 hours) or not specified when processing for trunk.");

        // Always use a daily timestep.
        InputTimeStepHours = 24;
        OutputTimeStepHours = 24;
    }

    internal virtual void ValidateDaveTimeStepSettings()
    {
        if (OutputTimeStepHours == 24)
            throw new ArgumentException("Output timestep must be subdaily when processing for Dave.");
        if (InputTimeStepHours == 0)
            throw new ArgumentException("Input timestep must be specified when processing for Dave.");
        if (OutputTimeStepHours == 0)
            throw new ArgumentException("Output timestep must be specified when processing for Dave.");
    }

    internal void ValidateEmailSettings()
    {
        if (EmailNotifications.HasFlag(EmailNotificationType.None))
        {
            // Ensure that "none" is not combined with other values.
            if ((EmailNotifications & EmailNotificationType.None) != EmailNotificationType.None)
                throw new ArgumentException("Cannot combine \"none\" and other values for email notifications.");

            // No email notifications are required, so no further validation is
            // necessary.
            return;
        }

        // Some email notifications are required, so ensure that the email
        // address is set.
        if (string.IsNullOrEmpty(Email))
            throw new ArgumentException("Must specify an email address when email notifications are enabled.");
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
