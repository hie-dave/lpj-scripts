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

        string expected = $"{name}_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_mon_195101-195212.nc";
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

    [Theory]
    [InlineData("tas_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_1hr_1951010100-1951123123.nc", NarClim2Frequency.Hour1)]
    [InlineData("tas_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_3hr_1951010100-1951123121.nc", NarClim2Frequency.Hour3)]
    [InlineData("tas_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_day_19510101-19511231.nc", NarClim2Frequency.Day)]
    [InlineData("tas_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_mon_195101-195112.nc", NarClim2Frequency.Month)]
    [InlineData("/some/path/to/tas_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_mon_195101-195112.nc", NarClim2Frequency.Month)]
    public void GetFrequencyFromFilename_ReturnsCorrectFrequency(string filename, NarClim2Frequency expected)
    {
        var result = NarClim2Dataset.GetFrequencyFromFilename(filename);
        Assert.Equal(expected, result);
    }

    [Theory]
    [InlineData("invalid_filename.nc")]
    [InlineData("tas_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1.nc")]
    [InlineData("tas_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_invalid_19510101-19511231.nc")]
    public void GetFrequencyFromFilename_ThrowsForInvalidFilename(string filename)
    {
        Assert.ThrowsAny<ArgumentException>(() => NarClim2Dataset.GetFrequencyFromFilename(filename));
    }

    [Theory]
    [InlineData("tas_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_1hr_1951010100-1951123123.nc", true, 1951, 1, 1, 0)]
    [InlineData("tas_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_1hr_1951010100-1951123123.nc", false, 1951, 12, 31, 23)]
    [InlineData("tas_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_3hr_1951010100-1951123121.nc", true, 1951, 1, 1, 0)]
    [InlineData("tas_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_3hr_1951010100-1951123121.nc", false, 1951, 12, 31, 21)]
    [InlineData("tas_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_day_19510101-19511231.nc", true, 1951, 1, 1, 0)]
    [InlineData("tas_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_day_19510101-19511231.nc", false, 1951, 12, 31, 0)]
    [InlineData("tas_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_mon_195101-195112.nc", true, 1951, 1, 1, 0)]
    [InlineData("tas_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_mon_195101-195112.nc", false, 1951, 12, 1, 0)]
    public void GetDateFromFilename_ReturnsCorrectDate(string filename, bool start, int year, int month, int day, int hour)
    {
        var expected = new DateTime(year, month, day, hour, 0, 0);
        var result = NarClim2Dataset.GetDateFromFilename(filename, start);
        Assert.Equal(expected, result);
    }

    [Theory]
    [InlineData("invalid_filename.nc")]
    [InlineData("tas_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1.nc")]
    [InlineData("tas_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_mon.nc")]
    [InlineData("tas_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_mon_invalid.nc")]
    public void GetDateFromFilename_ThrowsForInvalidFilename(string filename)
    {
        var ex = Assert.ThrowsAny<ArgumentException>(() => NarClim2Dataset.GetDateFromFilename(filename, true));
        Assert.Contains("Unable to determine frequency from filename. Invalid filename format", ex.Message);
    }

    [Theory]
    [InlineData("tas_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_1hr_1951010199-1951123199.nc")]
    [InlineData("tas_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_mon_195113-195112.nc")]
    [InlineData("tas_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_day_19510132-19511231.nc")]
    public void GetDateFromFilename_ThrowsForInvalidDates(string filename)
    {
        Assert.ThrowsAny<FormatException>(() => NarClim2Dataset.GetDateFromFilename(filename, true));
    }

    [Theory]
    [InlineData("invalid.filename.nc")]
    [InlineData("tas_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_day_f-19511231.nc")]
    [InlineData("tas_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_day_19510132-g.nc")]
    public void GetDateFromFileName_ThrowsForInvalidFileNames(string filename)
    {
        Assert.ThrowsAny<ArgumentException>(() => NarClim2Dataset.GetDateFromFilename(filename, true));
    }

    [Theory]
    [InlineData(NarClim2Domain.AUS18, NarClim2GCM.AccessEsm15, NarClim2Experiment.Historical, NarClim2RCM.WRF412R3, NarClim2Frequency.Hour1, "NARCliM2.0_AUS-18_ACCESS-ESM1-5_historical_NARCliM2-0-WRF412R3")]
    [InlineData(NarClim2Domain.SEAus04, NarClim2GCM.Ukesm10Ll, NarClim2Experiment.SSP370, NarClim2RCM.WRF412R5, NarClim2Frequency.Day, "NARCliM2.0_NARCliM2-0-SEAus-04_UKESM1-0-LL_ssp370_NARCliM2-0-WRF412R5")]
    public void TestNarClim2DatasetName(
        NarClim2Domain domain,
        NarClim2GCM gcm,
        NarClim2Experiment experiment,
        NarClim2RCM rcm,
        NarClim2Frequency frequency,
        string expected)
    {
        NarClim2Dataset dataset = new NarClim2Dataset(
            "/input",
            domain,
            gcm,
            experiment,
            rcm,
            frequency);
        Assert.Equal(expected, dataset.DatasetName);
    }

    [Theory]
    [InlineData(ClimateVariable.MaxTemperature, NarClim2Frequency.Hour1)]
    [InlineData(ClimateVariable.MinTemperature, NarClim2Frequency.Hour1)]
    [InlineData(ClimateVariable.MaxTemperature, NarClim2Frequency.Hour3)]
    [InlineData(ClimateVariable.MinTemperature, NarClim2Frequency.Hour3)]
    public void GetInputFiles_InvalidFrequency(ClimateVariable variable, NarClim2Frequency frequency)
    {
        NarClim2Dataset dataset = new NarClim2Dataset(
            "/input",
            NarClim2Domain.AUS18,
            NarClim2GCM.AccessEsm15,
            NarClim2Experiment.Historical,
            NarClim2RCM.WRF412R3,
            frequency);
        Assert.ThrowsAny<ArgumentException>(() => dataset.GetInputFiles(variable));
    }

    [Theory]
    [InlineData(NarClim2Domain.AUS18, NarClim2GCM.EcEarth3Veg, NarClim2Experiment.SSP126, NarClim2RCM.WRF412R5, NarClim2Frequency.Hour3, "EC-Earth3-Veg/ssp126/NARCliM2-0-WRF412R5")]
    [InlineData(NarClim2Domain.SEAus04, NarClim2GCM.NorEsm2Mm, NarClim2Experiment.Historical, NarClim2RCM.WRF412R3, NarClim2Frequency.Hour1, "NorESM2-MM/historical/NARCliM2-0-WRF412R3")]
    public void TestGetOuptutDirectory(
        NarClim2Domain domain,
        NarClim2GCM gcm,
        NarClim2Experiment experiment,
        NarClim2RCM rcm,
        NarClim2Frequency frequency,
        string expected)
    {
        NarClim2Dataset dataset = new NarClim2Dataset(
            "/input",
            domain,
            gcm,
            experiment,
            rcm,
            frequency);
        Assert.Equal(expected, dataset.GetOutputDirectory());
    }

    [Theory]
    [InlineData(NarClim2Domain.AUS18)]
    [InlineData(NarClim2Domain.SEAus04)]
    public void TestGetDomain(NarClim2Domain domain)
    {
        NarClim2Dataset dataset = new NarClim2Dataset(
            "/input",
            domain,
            NarClim2GCM.AccessEsm15,
            NarClim2Experiment.Historical,
            NarClim2RCM.WRF412R3,
            NarClim2Frequency.Hour1);
        Assert.Equal(domain, dataset.Domain);
    }

    [Theory]
    [InlineData("/input")]
    [InlineData("a/b/c/d.txt")]
    public void TestGetInputPath(string inputPath)
    {
        NarClim2Dataset dataset = new NarClim2Dataset(
            inputPath,
            NarClim2Domain.AUS18,
            NarClim2GCM.AccessEsm15,
            NarClim2Experiment.Historical,
            NarClim2RCM.WRF412R3,
            NarClim2Frequency.Hour1);
        Assert.Equal(inputPath, dataset.BasePath);
    }
}
