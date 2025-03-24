namespace ClimateProcessing.Models;

public record VariableInfo(string Name, string Units);

/// <summary>
/// An interface to a climate dataset.
/// </summary>
public interface IClimateDataset
{
    /// <summary>
    /// Name of the dataset.
    /// </summary>
    string DatasetName { get; }

    /// <summary>
    /// Get the directory containing the input files for a given variable.
    /// </summary>
    /// <param name="variable">The variable.</param>
    string GetInputFilesDirectory(ClimateVariable variable);

    /// <summary>
    /// Get the input files for a given variable, sorted by date ascending.
    /// </summary>
    /// <param name="variable">The variable.</param>
    IEnumerable<string> GetInputFiles(ClimateVariable variable);

    /// <summary>
    /// Get metadata for a given variable as it exists in this dataset.
    /// </summary>
    /// <param name="variable">The variable.</param>
    /// <returns>The variable info.</returns>
    VariableInfo GetVariableInfo(ClimateVariable variable);

    /// <summary>
    /// Get the path to the output file for a given variable.
    /// </summary>
    /// <param name="variable">The variable.</param>
    /// <returns>The path to the output file.</returns>
    string GenerateOutputFileName(ClimateVariable variable);

    /// <summary>
    /// Get a non-rooted (ie relative) path to the output directory. This will
    /// be appended to the base output path for all output files generated from
    /// this dataset.
    /// </summary>
    /// <returns>The output directory.</returns>
    public string GetOutputDirectory();
}
