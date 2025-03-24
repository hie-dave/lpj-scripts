using ClimateProcessing.Models;

namespace ClimateProcessing.Services;

/// <summary>
/// Handles management of paths used in PBS scripts.
/// </summary>
/// <remarks>
/// There may be some redundant directory creation going on here. It's useful
/// for some of the unit tests though, so it can stay for now.
/// </remarks
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
        outputDirectory = outputPath;
    }

    /// <inheritdoc />
    public string GetChecksumFilePath()
    {
        return Path.Combine(GetBasePath(PathType.Output), checksumFilename);
    }

    /// <inheritdoc />
    public void CreateDirectoryTree(IClimateDataset dataset)
    {
        // Create top-level directories.
        Directory.CreateDirectory(GetBasePath(PathType.Log));
        Directory.CreateDirectory(GetBasePath(PathType.Script));
        Directory.CreateDirectory(GetBasePath(PathType.Stream));
        Directory.CreateDirectory(GetBasePath(PathType.Working));
        Directory.CreateDirectory(GetBasePath(PathType.Output));

        // Create dataset-level directories.
        Directory.CreateDirectory(GetDatasetPath(dataset, PathType.Output));
        Directory.CreateDirectory(GetDatasetPath(dataset, PathType.Working));
    }

    /// <inheritdoc />
    public string GetDatasetPath(IClimateDataset dataset, PathType pathType)
    {
        if (pathType == PathType.Log || pathType == PathType.Script || pathType == PathType.Stream)
            throw new ArgumentException($"Path type {pathType} is not valid at the dataset-level. Only output/working paths make sense here.", nameof(pathType));

        string basePath = GetBasePath(pathType);
        string relativePath = dataset.GetOutputDirectory();
        string fullPath = Path.Combine(basePath, relativePath);

        // Create directory if needed
        Directory.CreateDirectory(fullPath);

        return fullPath;
    }

    /// <inheritdoc />
    public string GetDatasetFileName(IClimateDataset dataset, ClimateVariable variable, PathType pathType)
    {
        string directory = GetDatasetPath(dataset, pathType);
        Directory.CreateDirectory(directory);
        string fileName = dataset.GenerateOutputFileName(variable);
        return Path.Combine(directory, fileName);
    }

    /// <inheritdoc /> 
    public string GetBasePath(PathType pathType)
    {
        string directoryName = pathType switch
        {
            PathType.Output => outputDirectoryName,
            PathType.Working => workDirectoryName,
            PathType.Script => scriptDirectoryName,
            PathType.Log => logDirectoryName,
            PathType.Stream => streamDirectoryName,
            _ => throw new ArgumentException($"Unsupported path type: {pathType}", nameof(pathType))
        };
        string directory = Path.Combine(outputDirectory, directoryName);
        Directory.CreateDirectory(directory);
        return directory;
    }
}
