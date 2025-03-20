using ClimateProcessing.Models;

namespace ClimateProcessing.Tests.Mocks;

/// <summary>
/// A more sophisticated mock dataset that maintains a mapping of climate variables to their names and units.
/// Useful for integration-style tests that need to work with multiple variables.
/// </summary>
internal class DynamicMockDataset : IClimateDataset
{
    private readonly string inputPath;
    private readonly string outputPath;
    private static readonly Dictionary<ClimateVariable, (string Name, string Units)> _variableInfo = new()
    {
        [ClimateVariable.Temperature] = ("tas", "K"),
        [ClimateVariable.SpecificHumidity] = ("huss", "1"),
        [ClimateVariable.SurfacePressure] = ("ps", "Pa"),
        [ClimateVariable.Precipitation] = ("pr", "kg m-2 s-1"),
        [ClimateVariable.ShortwaveRadiation] = ("rsds", "W m-2"),
        [ClimateVariable.WindSpeed] = ("sfcWind", "m s-1"),
        [ClimateVariable.MaxTemperature] = ("tasmax", "K"),
        [ClimateVariable.MinTemperature] = ("tasmin", "K")
    };

    public string DatasetName => GetType().Name;

    public DynamicMockDataset(string inputPath, string outputPath)
    {
        this.inputPath = inputPath;
        this.outputPath = outputPath;
    }

    public void SetVariableInfo(ClimateVariable variable, string name, string units)
    {
        _variableInfo[variable] = (name, units);
    }

    // IClimateDataset implementation.

    public string GetInputFilesDirectory(ClimateVariable variable)
    {
        return Path.Combine(inputPath, variable.ToString());
    }

    public IEnumerable<string> GetInputFiles(ClimateVariable variable)
    {
        string inputDirectory = GetInputFilesDirectory(variable);
        return new[] { "input1.nc", "input2.nc" }
            .Select(f => Path.Combine(inputDirectory, f));
    }

    VariableInfo IClimateDataset.GetVariableInfo(ClimateVariable variable)
    {
        (string name, string units) = _variableInfo[variable];
        return new VariableInfo(name, units);
    }

    public string GenerateOutputFileName(ClimateVariable variable) =>
        $"{_variableInfo[variable].Name}_output.nc";

    public string GetOutputDirectory() => "dynamic_mock";
}
