using ClimateProcessing.Models;

namespace ClimateProcessing.Services;

/// <summary>
/// Handles management of paths used in PBS scripts.
/// </summary>
public interface IPathManager
{
    // /// <summary>
    // /// Get the path to the log file for the specified job.
    // /// </summary>
    // /// <param name="jobName">The name of the job.</param>
    // /// <returns>The path to the log file.</returns>
    // string GetLogFilePath(string jobName);

    // /// <summary>
    // /// Get the path to the stream file for the specified job.
    // /// </summary>
    // /// <param name="jobName">The name of the job.</param>
    // /// <returns>The path to the stream file.</returns>
    // string GetStreamFilePath(string jobName);

    // /// <summary>
    // /// Get the path to the script file for the specified job.
    // /// </summary>
    // /// <param name="jobName">The name of the job.</param>
    // /// <returns>The path to the script file.</returns>
    // string GetScriptFilePath(string jobName);

    // /// <summary>
    // /// Get the path to the output file for the specified output file name.
    // /// </summary>
    // /// <param name="outputFileName">The name of the output file.</param>
    // /// <returns>The path to the output file.</returns>
    // string GetOutputFilePath(string outputFileName);

    // /// <summary>
    // /// Get the path to the temporary file for the specified output file name.
    // /// </summary>
    // /// <param name="outputFileName">The name of the output file.</param>
    // /// <returns>The path to the temporary file.</returns>
    // string GetTempFilePath(string outputFileName);

    /// <summary>
    /// Get the path to the checksum file.
    /// </summary>
    /// <returns>The path to the checksum file.</returns>
    string GetChecksumFilePath();

    /// <summary>
    /// Create the required output directory tree for the specified dataset.
    /// </summary>
    /// <param name="dataset">The dataset.</param>
    void CreateDirectoryTree(IClimateDataset dataset);

    /// <summary>
    /// Get a dataset-level path for the specified path type, and create the
    /// directory if it doesn't exist.
    /// </summary>
    /// <param name="dataset">The dataset.</param>
    /// <param name="type">The type of path.</param>
    /// <returns>The path.</returns>
    string GetDatasetPath(IClimateDataset dataset, PathType type);

    /// <summary>
    /// Get a dataset-level file name for the specified variable and path type.
    /// </summary>
    /// <param name="dataset">The dataset.</param>
    /// <param name="variable">The variable.</param>
    /// <param name="type">The type of path.</param>
    /// <returns>The file name.</returns>
    string GetDatasetFileName(IClimateDataset dataset, ClimateVariable variable, PathType type);

    /// <summary>
    /// Gets the base directory for the specified path type.
    /// </summary>
    /// <param name="pathType">The path type.</param>
    /// <returns>The base directory.</returns>
    string GetBasePath(PathType pathType);
}
