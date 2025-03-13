using ClimateProcessing.Models;

namespace ClimateProcessing.Services;

/// <summary>
/// Interface for dataset-specific script generators.
/// </summary>
/// <typeparam name="T">The type of climate dataset this generator handles.</typeparam>
public interface IScriptGenerator<T> where T : IClimateDataset
{
    /// <summary>
    /// Generate processing scripts for the dataset.
    /// </summary>
    /// <param name="dataset">The dataset to process.</param>
    /// <returns>The path to the top-level script.</returns>
    Task<string> GenerateScriptsAsync(T dataset);
}
