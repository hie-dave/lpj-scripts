using ClimateProcessing.Models;
using ClimateProcessing.Services;

namespace ClimateProcessing.Tests.Mocks;

/// <summary>
/// This class extends <see cref="NarClim2ScriptGenerator"/> and exposes some
/// of its private methods for testing purposes.
/// </summary>
/// <remarks>
/// None of the base class' behaviour is changed.
/// </remarks>
internal class TestNarClim2ScriptGenerator : NarClim2ScriptGenerator
{
    /// <summary>
    /// Create a new <see cref="NarClim2ScriptGenerator" /> instance.
    /// </summary>
    /// <param name="config">The configuration to use.</param>
    public TestNarClim2ScriptGenerator(NarClim2Config config) : base(config)
    {
    }

    /// <summary>
    /// Call WritePreMerge().
    /// </summary>
    /// <param name="writer">Any text writer.</param>
    /// <param name="dataset">A climate dataset.</param>
    /// <param name="variable">A climate variable.</param>
    public async Task PreMergeAsync(TextWriter writer, IClimateDataset dataset, ClimateVariable variable)
    {
        await base.WritePreMerge(writer, dataset, variable);
    }

    /// <summary>
    /// Call GetRlonValuesFile().
    /// </summary>
    /// <param name="dataset"></param>
    /// <returns></returns>
    public string ReadRlonValuesFile(NarClim2Dataset dataset)
    {
        return GetRlonValuesFile(dataset);
    }
}
