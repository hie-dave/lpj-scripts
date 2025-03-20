namespace ClimateProcessing.Models;

public record VariableInfo(string Name, string Units);

public interface IClimateDataset
{
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
    VariableInfo GetVariableInfo(ClimateVariable variable);
    string GenerateOutputFileName(ClimateVariable variable);

    /// <summary>
    /// Get a non-rooted (ie relative) path to the output directory. This will
    /// be appended to the base output path for all output files generated from
    /// this dataset.
    /// </summary>
    /// <remarks>
    /// This is a poor design choice. Need to rethink this - the onus for this
    /// probably shouldn't be on the dataset, but I need to work out how this
    /// should be encapsulated. For now it can live here.
    /// </remarks>
    /// <returns>The output directory.</returns>
    string GetOutputDirectory();
}
