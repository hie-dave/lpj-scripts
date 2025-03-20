using ClimateProcessing.Models;

namespace ClimateProcessing.Tests.Mocks;

/// <summary>
/// A simple mock dataset that returns the same variable name and units for all variables.
/// Useful for basic tests that don't need specific variable information.
/// </summary>
internal class StaticMockDataset : IClimateDataset
{
    private readonly string basePath;
    private readonly string varName;
    private readonly string varUnits;

    public StaticMockDataset(
        string basePath,
        string varName = "tas",
        string varUnits = "K")
    {
        this.basePath = basePath;
        this.varName = varName;
        this.varUnits = varUnits;
    }

    public string DatasetName => "mock";

    public string GetInputFilesDirectory(ClimateVariable variable)
        => basePath;

    public IEnumerable<string> GetInputFiles(ClimateVariable variable)
        => new[] { Path.Combine(basePath, "input.nc") };

    public string GetOutputDirectory() => "static_mock";

    public VariableInfo GetVariableInfo(ClimateVariable variable)
        => new(varName, varUnits);

    public string GenerateOutputFileName(ClimateVariable variable)
        => "output.nc";
}
