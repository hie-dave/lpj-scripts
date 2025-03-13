using ClimateProcessing.Models;

namespace ClimateProcessing.Services;

/// <summary>
/// Factory for creating dataset-specific script generators.
/// </summary>
public interface IScriptGeneratorFactory
{
    /// <summary>
    /// Generate processing scripts for the specified dataset.
    /// </summary>
    /// <param name="dataset">The dataset to process.</param>
    /// <returns>The path to the top-level script.</returns>
    Task<string> GenerateScriptsAsync(IClimateDataset dataset);
}
