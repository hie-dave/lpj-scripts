namespace ClimateProcessing.Services;

/// <summary>
/// Handles management of paths used in PBS scripts.
/// </summary>
public class PathManager : IPathManager
{
    /// <summary>
    /// Subdirectory of the output directory into which the scripts are written.
    /// </summary>
    private const string scriptDirectoryName = "scripts";

    /// <summary>
    /// Subdirectory of the output directory into which the logs are written.
    /// </summary>
    private const string logDirectoryName = "logs";

    /// <summary>
    /// Subdirectory of the output directory into which stdout will be streamed.
    /// </summary>
    private const string streamDirectoryName = "streams";

    /// <summary>
    /// Name of the directory into which all output files will be written.
    /// </summary>
    private const string outputDirectoryName = "output";

    /// <summary>
    /// Name of the working directory, into which intermediate output files are
    /// written.
    /// </summary>
    private const string workDirectoryName = "tmp";

    /// <summary>
    /// Name of the checksum file.
    /// </summary>
    private const string checksumFilename = "sha512sums.txt";

    /// <summary>
    /// The output base directory.
    /// </summary>
    private readonly string outputDirectory;

    /// <summary>
    /// Creates a new instance of the <see cref="PathManager"/> class.
    /// </summary>
    /// <param name="outputPath">The output base directory.</param>
    public PathManager(string outputPath)
    {
        this.outputDirectory = outputPath;
    }

    /// <summary>
    /// Get the name of the log file for a given job name.
    /// </summary>
    /// <param name="jobName">The name of the job.</param>
    /// <returns>The name of the log file.</returns>
    private string GetLogFileName(string jobName) => $"{jobName}.log";

    /// <inheritdoc />
    public string GetLogFilePath(string jobName)
    {
        string logFileName = GetLogFileName(jobName);
        string logFile = Path.Combine(GetLogPath(), logFileName);
        return logFile;
    }

    /// <inheritdoc />
    public string GetStreamFilePath(string jobName)
    {
        string logFileName = GetLogFileName(jobName);
        string streamFile = Path.Combine(GetStreamPath(), logFileName);
        return streamFile;
    }

    /// <inheritdoc /> 
    public string GetScriptFilePath(string jobName)
    {
        string scriptName = jobName;
        string scriptFile = Path.Combine(GetScriptPath(), scriptName);
        return scriptFile;
    }

    /// <inheritdoc />
    public string GetOutputFilePath(string outputFileName)
    {
        string outputPath = GetOutputPath();
        string outputFile = Path.Combine(outputPath, outputFileName);
        return outputFile;
    }

    /// <inheritdoc />
    public string GetTempFilePath(string outputFileName)
    {
        string outputPath = GetWorkingPath();
        string outputFile = Path.Combine(outputPath, outputFileName);
        return outputFile;
    }

    /// <inheritdoc />
    public string GetChecksumFilePath()
    {
        return Path.Combine(GetOutputPath(), checksumFilename);
    }

    /// <inheritdoc />
    public void CreateDirectoryTree()
    {
        Directory.CreateDirectory(GetLogPath());
        Directory.CreateDirectory(GetScriptPath());
        Directory.CreateDirectory(GetStreamPath());
        Directory.CreateDirectory(GetWorkingPath());
        Directory.CreateDirectory(GetOutputPath());
    }

    /// <inheritdoc />
    public string GetWorkingPath()
    {
        return Path.Combine(outputDirectory, workDirectoryName);
    }

    /// <inheritdoc />
    public string GetOutputPath()
    {
        return Path.Combine(outputDirectory, outputDirectoryName);
    }

    /// <summary>
    /// Get the path to the directory in which the scripts will be stored, and
    /// create it if it doesn't exist.
    /// </summary>
    /// <returns>The path to the script directory.</returns>
    private string GetScriptPath()
    {
        return Path.Combine(outputDirectory, scriptDirectoryName);
    }

    /// <summary>
    /// Get the path to the directory in which stdout/stderr will be streamed,
    /// and create it if it doesn't exist.
    /// </summary>
    /// <returns>The path to the stream directory.</returns>
    private string GetStreamPath()
    {
        return Path.Combine(outputDirectory, streamDirectoryName);
    }

    /// <summary>
    /// Get the path to the directory in which log files will be stored, and
    /// create it if it doesn't exist.
    /// </summary>
    /// <returns>The path to the log directory.</returns>
    private string GetLogPath()
    {
        return Path.Combine(outputDirectory, logDirectoryName);
    }
}
