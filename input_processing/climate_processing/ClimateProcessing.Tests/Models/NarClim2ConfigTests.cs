using ClimateProcessing.Models;
using ClimateProcessing.Services;
using Xunit;

namespace ClimateProcessing.Tests.Models;

public class NarClim2ConfigTests
{
    private readonly string _validInputDir;
    private readonly string _validProject;

    public NarClim2ConfigTests()
    {
        _validInputDir = Path.Combine(Path.GetTempPath(), "test_input");
        Directory.CreateDirectory(_validInputDir);
        _validProject = "test_project";
    }

    private NarClim2Config CreateValidConfig()
    {
        return new NarClim2Config
        {
            InputDirectory = _validInputDir,
            Project = _validProject,
            InputTimeStepHours = 24,
            OutputTimeStepHours = 24,
            Version = ModelVersion.Trunk
        };
    }

    [Fact]
    public void CreateDatasets_WithNoFilters_ReturnsAllCombinations()
    {
        var config = CreateValidConfig();

        var datasets = config.CreateDatasets().ToList();

        var expectedCount = Enum.GetValues<NarClim2Domain>().Length *
                            Enum.GetValues<NarClim2GCM>().Length *
                            Enum.GetValues<NarClim2Experiment>().Length *
                            Enum.GetValues<NarClim2RCM>().Length;
        Assert.Equal(expectedCount, datasets.Count);

        foreach (NarClim2Domain domain in Enum.GetValues<NarClim2Domain>())
        foreach (NarClim2GCM gcm in Enum.GetValues<NarClim2GCM>())
        foreach (NarClim2Experiment experiment in Enum.GetValues<NarClim2Experiment>())
        foreach (NarClim2RCM rcm in Enum.GetValues<NarClim2RCM>())
            Assert.Contains(datasets, ds =>
                ds.Domain == domain &&
                ds.GCM == gcm &&
                ds.Experiment == experiment &&
                ds.RCM == rcm);
    }

    [Theory]
    [InlineData(NarClim2Domain.AUS18)]
    [InlineData(NarClim2Domain.SEAus04)]
    public void CreateDatasets_WithSingleDomain_FiltersCorrectly(NarClim2Domain domain)
    {
        NarClim2Config config = CreateValidConfig();
        config.Domains = [NarClim2Constants.DomainNames.ToString(domain)];

        List<NarClim2Dataset> datasets = config.CreateDatasets().ToList();

        Assert.NotEmpty(datasets);
        Assert.All(datasets, ds => Assert.Equal(domain, ds.Domain));
    }

    [Theory]
    [InlineData(NarClim2GCM.AccessEsm15)]
    [InlineData(NarClim2GCM.EcEarth3Veg)]
    [InlineData(NarClim2GCM.MpiEsm12Hr)]
    [InlineData(NarClim2GCM.NorEsm2Mm)]
    [InlineData(NarClim2GCM.Ukesm10Ll)]
    public void CreateDatasets_WithSingleGCM_FiltersCorrectly(NarClim2GCM gcm)
    {
        NarClim2Config config = CreateValidConfig();
        config.GCMs = [NarClim2Constants.GCMNames.ToString(gcm)];

        List<NarClim2Dataset> datasets = config.CreateDatasets().ToList();

        Assert.NotEmpty(datasets);
        Assert.All(datasets, ds => Assert.Equal(gcm, ds.GCM));
    }

    [Theory]
    [InlineData(NarClim2Experiment.Historical)]
    [InlineData(NarClim2Experiment.SSP126)]
    [InlineData(NarClim2Experiment.SSP370)]
    public void CreateDatasets_WithSingleExperiment_FiltersCorrectly(NarClim2Experiment experiment)
    {
        NarClim2Config config = CreateValidConfig();
        config.Experiments = [NarClim2Constants.ExperimentNames.ToString(experiment)];

        List<NarClim2Dataset> datasets = config.CreateDatasets().ToList();

        Assert.NotEmpty(datasets);
        Assert.All(datasets, ds => Assert.Equal(experiment, ds.Experiment));
    }

    [Theory]
    [InlineData(NarClim2RCM.WRF412R3)]
    [InlineData(NarClim2RCM.WRF412R5)]
    public void CreateDatasets_WithSingleRCM_FiltersCorrectly(NarClim2RCM rcm)
    {
        NarClim2Config config = CreateValidConfig();
        config.RCMs = [NarClim2Constants.RCMNames.ToString(rcm)];

        List<NarClim2Dataset> datasets = config.CreateDatasets().ToList();

        Assert.NotEmpty(datasets);
        Assert.All(datasets, ds => Assert.Equal(rcm, ds.RCM));
    }

    [Theory]
    [InlineData(NarClim2Domain.AUS18, NarClim2GCM.AccessEsm15, NarClim2Experiment.Historical, NarClim2RCM.WRF412R3)]
    [InlineData(NarClim2Domain.SEAus04, NarClim2GCM.Ukesm10Ll, NarClim2Experiment.SSP370, NarClim2RCM.WRF412R5)]
    public void CreateDatasets_WithMultipleFilters_AppliesAllFilters(
        NarClim2Domain domain,
        NarClim2GCM gcm,
        NarClim2Experiment experiment,
        NarClim2RCM rcm)
    {
        NarClim2Config config = CreateValidConfig();
        config.Domains = [NarClim2Constants.DomainNames.ToString(domain)];
        config.GCMs = [NarClim2Constants.GCMNames.ToString(gcm)];
        config.Experiments = [NarClim2Constants.ExperimentNames.ToString(experiment)];
        config.RCMs = [NarClim2Constants.RCMNames.ToString(rcm)];

        List<NarClim2Dataset> datasets = config.CreateDatasets().ToList();

        Assert.Single(datasets);
        NarClim2Dataset dataset = datasets.First();
        Assert.Equal(domain, dataset.Domain);
        Assert.Equal(gcm, dataset.GCM);
        Assert.Equal(experiment, dataset.Experiment);
        Assert.Equal(rcm, dataset.RCM);
    }

    [Theory]
    [InlineData("invalid-domain")]
    public void CreateDatasets_WithInvalidDomain_ThrowsException(string invalidDomain)
    {
        var config = CreateValidConfig();
        config.Domains = new[] { invalidDomain };

        Assert.Throws<ArgumentException>(() => config.CreateDatasets().ToList());
    }

    [Theory]
    [InlineData("invalid-gcm")]
    public void CreateDatasets_WithInvalidGCM_ThrowsException(string invalidGCM)
    {
        var config = CreateValidConfig();
        config.GCMs = new[] { invalidGCM };

        Assert.Throws<ArgumentException>(() => config.CreateDatasets().ToList());
    }

    [Theory]
    [InlineData("invalid-experiment")]
    public void CreateDatasets_WithInvalidExperiment_ThrowsException(string invalidExperiment)
    {
        var config = CreateValidConfig();
        config.Experiments = new[] { invalidExperiment };

        Assert.Throws<ArgumentException>(() => config.CreateDatasets().ToList());
    }

    [Theory]
    [InlineData("invalid-rcm")]
    public void CreateDatasets_WithInvalidRCM_ThrowsException(string invalidRCM)
    {
        var config = CreateValidConfig();
        config.RCMs = new[] { invalidRCM };

        Assert.Throws<ArgumentException>(() => config.CreateDatasets().ToList());
    }

    [Fact]
    public void CreateScriptGenerator_ReturnsNarClim2ScriptGenerator()
    {
        var config = CreateValidConfig();

        var generator = config.CreateScriptGenerator();

        Assert.IsType<NarClim2ScriptGenerator>(generator);
    }

    [Fact]
    public void Validate_WithValidConfig_DoesNotThrow()
    {
        var config = CreateValidConfig();

        var exception = Record.Exception(() => config.Validate());
        Assert.Null(exception);
    }

    [Fact]
    public void Validate_WithInvalidInputDirectory_ThrowsException()
    {
        var config = CreateValidConfig();
        config.InputDirectory = "non-existent-directory";

        Assert.False(Directory.Exists(config.InputDirectory));
        Assert.Throws<ArgumentException>(() => config.Validate());
    }
}
