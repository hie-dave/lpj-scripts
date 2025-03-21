namespace ClimateProcessing.Services;

/// <summary>
/// Handles management of paths used in PBS scripts.
/// </summary>
public interface IPathManager
{
    /// <summary>
    /// Get the path to the log file for the specified job.
    /// </summary>
    /// <param name="jobName">The name of the job.</param>
    /// <returns>The path to the log file.</returns>
    string GetLogFilePath(string jobName);

    /// <summary>
    /// Get the path to the stream file for the specified job.
    /// </summary>
    /// <param name="jobName">The name of the job.</param>
    /// <returns>The path to the stream file.</returns>
    string GetStreamFilePath(string jobName);

    /// <summary>
    /// Get the path to the script file for the specified job.
    /// </summary>
    /// <param name="jobName">The name of the job.</param>
    /// <returns>The path to the script file.</returns>
    string GetScriptFilePath(string jobName);

    /// <summary>
    /// Get the path to the output file for the specified output file name.
    /// </summary>
    /// <param name="outputFileName">The name of the output file.</param>
    /// <returns>The path to the output file.</returns>
    string GetOutputFilePath(string outputFileName);

    /// <summary>
    /// Get the path to the temporary file for the specified output file name.
    /// </summary>
    /// <param name="outputFileName">The name of the output file.</param>
    /// <returns>The path to the temporary file.</returns>
    string GetTempFilePath(string outputFileName);

    /// <summary>
    /// Get the path to the checksum file.
    /// </summary>
    /// <returns>The path to the checksum file.</returns>
    string GetChecksumFilePath();

    /// <summary>
    /// Get the path to the working directory for a dataset.
    /// </summary>
    /// <param name="dataset">The dataset.</param>
    /// <returns>The working directory path.</returns>
    string GetWorkingPath();


    /// <summary>
    /// Get the directory path into which output file tree will be generated.
    /// </summary>
    /// <returns>The output directory path.</returns>
    string GetOutputPath();

    /// <summary>
    /// Create the directory tree for the job.
    /// </summary>
    void CreateDirectoryTree();
}

