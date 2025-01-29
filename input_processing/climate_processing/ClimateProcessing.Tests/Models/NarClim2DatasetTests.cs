using Xunit;
using ClimateProcessing.Models;
using System.IO;

namespace ClimateProcessing.Tests.Models;

public class NarClim2DatasetTests : IDisposable
{
    private readonly string _testDir;
    private readonly NarClim2Dataset _dataset;

    private readonly NarClim2Domain _domain;
    private readonly NarClim2GCM _gcm;
    private readonly NarClim2Experiment _experiment;
    private readonly NarClim2RCM _rcm;
    private readonly NarClim2Frequency _frequency;

    public NarClim2DatasetTests()
    {
        // TOOD: test other parameters?
        _testDir = Path.Combine(Path.GetTempPath(), Path.GetRandomFileName());
        _domain = NarClim2Domain.AUS18;
        _gcm = NarClim2GCM.AccessEsm15;
        _rcm = NarClim2RCM.WRF412R3;
        _experiment = NarClim2Experiment.Historical;
        _frequency = NarClim2Frequency.Day;

        string baseDir = Path.Combine(
            _testDir,
            NarClim2Constants.Paths.MipEra,
            NarClim2Constants.Paths.ActivityId,
            NarClim2Constants.DomainNames.ToString(_domain),
            NarClim2Constants.Paths.Institution,
            NarClim2Constants.GCMNames.ToString(_gcm),
            NarClim2Constants.ExperimentNames.ToString(_experiment),
            NarClim2Constants.VariantLabels.GetVariantLabel(_gcm),
            NarClim2Constants.RCMNames.ToString(_rcm),
            NarClim2Constants.Paths.Version,
            NarClim2Constants.FrequencyNames.ToString(_frequency));

        // Create directories for each variable
        foreach (string var in new[] { "tas", "pr" })
        {
            string varDir = Path.Combine(baseDir, var, NarClim2Constants.Paths.LatestVersion);
            Directory.CreateDirectory(varDir);

            // Create test files for each variable
            CreateTestFile(Path.Combine(varDir, $"{var}_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_mon_195101-195112.nc"));
            CreateTestFile(Path.Combine(varDir, $"{var}_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_mon_195201-195212.nc"));
        }

        _dataset = new NarClim2Dataset(
            _testDir,
            domain: _domain,
            gcm: _gcm,
            rcm: _rcm,
            frequency: _frequency,
            experiment: _experiment);
    }

    public void Dispose()
    {
        if (Directory.Exists(_testDir))
            Directory.Delete(_testDir, true);
    }

    private void CreateTestFile(string filepath)
    {
        Directory.CreateDirectory(Path.GetDirectoryName(filepath)!);
        File.WriteAllText(filepath, "");
    }

    [Theory]
    [InlineData(ClimateVariable.Temperature, 2)]
    [InlineData(ClimateVariable.Precipitation, 2)]
    public void GetInputFiles_ReturnsCorrectFiles(ClimateVariable variable, int nfile)
    {
        List<string> files = _dataset.GetInputFiles(variable).ToList();

        Assert.Equal(nfile, files.Count);

        string name = _dataset.GetVariableInfo(variable).Name;
        Assert.Contains(files, f => Path.GetFileName(f).StartsWith($"{name}_") && f.Contains("195101-195112"));
        Assert.Contains(files, f => Path.GetFileName(f).StartsWith($"{name}_") && f.Contains("195201-195212"));
    }

    [Theory]
    [InlineData(ClimateVariable.Temperature, "tas", "K")]
    [InlineData(ClimateVariable.Precipitation, "pr", "kg m-2 s-1")]
    [InlineData(ClimateVariable.SpecificHumidity, "huss", "1")]
    [InlineData(ClimateVariable.SurfacePressure, "ps", "Pa")]
    [InlineData(ClimateVariable.ShortwaveRadiation, "rsds", "W m-2")]
    [InlineData(ClimateVariable.WindSpeed, "sfcWind", "m s-1")]
    public void GetVariableInfo_ReturnsCorrectInfo(ClimateVariable variable, string expectedName, string expectedUnits)
    {
        var info = _dataset.GetVariableInfo(variable);
        Assert.Equal(expectedName, info.Name);
        Assert.Equal(expectedUnits, info.Units);
    }

    [Fact]
    public void GetVariableInfo_ForInvalidVariable_ThrowsException()
    {
        var invalidVariable = (ClimateVariable)999;
        var ex = Assert.Throws<ArgumentException>(() => _dataset.GetVariableInfo(invalidVariable));
        Assert.Contains("not supported in NARCliM2 dataset", ex.Message);
    }

    [Theory]
    [InlineData(ClimateVariable.Temperature)]
    [InlineData(ClimateVariable.Precipitation)]
    public void GenerateOutputFileName_ReturnsCorrectPattern(ClimateVariable variable)
    {
        var filename = _dataset.GenerateOutputFileName(variable);
        string name = _dataset.GetVariableInfo(variable).Name;

        string expected = $"{name}_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_mon_195101-195201.nc";
        Assert.Equal(expected, filename);
    }

    [Fact]
    public void GenerateOutputFileName_ForMissingVariable_ThrowsException()
    {
        // Create dataset pointing to empty directory
        var emptyDir = Path.Combine(_testDir, "empty");
        Directory.CreateDirectory(emptyDir);
        var emptyDataset = new NarClim2Dataset(emptyDir);

        var ex = Assert.Throws<InvalidOperationException>(() => 
            emptyDataset.GenerateOutputFileName(ClimateVariable.Temperature));
        Assert.Contains("No input files found for variable", ex.Message);
    }

    [Fact]
    public void GetMetadata_ReturnsCorrectValues()
    {
        var metadata = _dataset.GetMetadata();
        Assert.Equal("NARCliM2.0", metadata["source"]);
        Assert.Equal(NarClim2Constants.DomainNames.ToString(_domain), metadata["domain"]);
        Assert.Equal(NarClim2Constants.GCMNames.ToString(_gcm), metadata["gcm"]);
        Assert.Equal(NarClim2Constants.ExperimentNames.ToString(_experiment), metadata["experiment"]);
        Assert.Equal(NarClim2Constants.RCMNames.ToString(_rcm), metadata["rcm"]);
        Assert.Equal(NarClim2Constants.FrequencyNames.ToString(_frequency), metadata["frequency"]);
        Assert.Equal(NarClim2Constants.Paths.Version, metadata["version"]);
    }
}
