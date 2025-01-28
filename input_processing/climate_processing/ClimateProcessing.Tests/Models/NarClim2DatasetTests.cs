using Xunit;
using ClimateProcessing.Models;
using System.IO;

namespace ClimateProcessing.Tests.Models;

public class NarClim2DatasetTests : IDisposable
{
    private readonly string _testDir;
    private readonly NarClim2Dataset _dataset;

    public NarClim2DatasetTests()
    {
        // Create a temporary directory for test files
        _testDir = Path.Combine(Path.GetTempPath(), Path.GetRandomFileName());
        Directory.CreateDirectory(_testDir);

        // Create some test files
        CreateTestFile("tas_AUS18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_mon_195101-195112.nc");
        CreateTestFile("tas_AUS18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_mon_195201-195212.nc");
        CreateTestFile("pr_AUS18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_mon_195101-195112.nc");

        _dataset = new NarClim2Dataset(_testDir);
    }

    public void Dispose()
    {
        // Clean up test directory
        if (Directory.Exists(_testDir))
        {
            Directory.Delete(_testDir, true);
        }
    }

    private void CreateTestFile(string filename)
    {
        File.WriteAllText(Path.Combine(_testDir, filename), "");
    }

    [Fact]
    public void GetOutputFilePattern_ForTemperature_ReturnsCorrectPattern()
    {
        var pattern = _dataset.GetOutputFilePattern(ClimateVariable.Temperature);
        Assert.Equal("tas_AUS18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_mon_XXXXXXXX-XXXXXXXX.nc", pattern);
    }

    [Fact]
    public void GetOutputFilePattern_ForPrecipitation_ReturnsCorrectPattern()
    {
        var pattern = _dataset.GetOutputFilePattern(ClimateVariable.Precipitation);
        Assert.Equal("pr_AUS18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_mon_XXXXXXXX-XXXXXXXX.nc", pattern);
    }

    [Fact]
    public void GetOutputFilePattern_ForMissingVariable_ThrowsException()
    {
        var ex = Assert.Throws<InvalidOperationException>(() => 
            _dataset.GetOutputFilePattern(ClimateVariable.WindSpeed));
        Assert.Contains("No input files found for variable", ex.Message);
    }

    [Fact]
    public void GetOutputFilePattern_WithInvalidFiles_ThrowsException()
    {
        // Create a dataset with an empty directory
        var emptyDir = Path.Combine(_testDir, "empty");
        Directory.CreateDirectory(emptyDir);
        var emptyDataset = new NarClim2Dataset(emptyDir);

        var ex = Assert.Throws<InvalidOperationException>(() => 
            emptyDataset.GetOutputFilePattern(ClimateVariable.Temperature));
        Assert.Contains("No input files found for variable", ex.Message);
    }
}
