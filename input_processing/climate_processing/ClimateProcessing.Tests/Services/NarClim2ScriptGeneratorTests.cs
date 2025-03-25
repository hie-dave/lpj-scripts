using Xunit;
using ClimateProcessing.Models;
using Moq;
using ClimateProcessing.Services;
using ClimateProcessing.Tests.Mocks;

namespace ClimateProcessing.Tests.Services;

public class NarClim2ScriptGeneratorTests : IDisposable
{
    private const string outputDirectoryPrefix = "narclim2_script_generator_tests_output";
    private readonly string outputDirectory;
    private readonly NarClim2Config config;
    private readonly NarClim2ScriptGenerator generator;

    public NarClim2ScriptGeneratorTests()
    {
        outputDirectory = Path.Combine(Path.GetTempPath(), $"{outputDirectoryPrefix}_{Guid.NewGuid()}");
        Directory.CreateDirectory(outputDirectory);

        config = new NarClim2Config
        {
            OutputDirectory = outputDirectory,
            ChunkSizeSpatial = 10,
            ChunkSizeTime = 5,
            CompressOutput = true,
            CompressionLevel = 4,
            InputTimeStepHours = 3,
            OutputTimeStepHours = 3,
            Version = ModelVersion.Dave
        };

        generator = new NarClim2ScriptGenerator(config);
    }

    public void Dispose()
    {
        if (Directory.Exists(outputDirectory))
            Directory.Delete(outputDirectory, true);
    }

    /// <summary>
    /// Create a NarClim2Dataset which doesn't require the existence of a full
    /// dataset on the filesystem.
    /// </summary>
    /// <param name="basePath">The base path of the dataset.</param>
    /// <param name="domain">The domain of the dataset.</param>
    /// <param name="gcm">The GCM of the dataset.</param>
    /// <param name="experiment">The experiment of the dataset.</param>
    /// <param name="rcm">The RCM of the dataset.</param>
    /// <param name="frequency">The frequency of the dataset.</param>
    /// <param name="outputFileName">The name of the output file returned by GenerateOutputFileName() for any climate variable.</param>
    /// <returns>A mocked narclim2 dataset.</returns>
    private Mock<NarClim2Dataset> CreateMockDataset(
        string basePath = "/mock/path",
        NarClim2Domain domain = NarClim2Domain.AUS18,
        NarClim2GCM gcm = NarClim2GCM.AccessEsm15,
        NarClim2Experiment experiment = NarClim2Experiment.Historical,
        NarClim2RCM rcm = NarClim2RCM.WRF412R3,
        NarClim2Frequency frequency = NarClim2Frequency.Month,
        string outputFileName = "/mock/output/file.nc"
    )
    {
        var mockDataset = new Mock<NarClim2Dataset>(
            basePath,
            domain,
            gcm,
            experiment,
            rcm,
            frequency);
        mockDataset.CallBase = true;
        mockDataset.Setup(x => x.GenerateOutputFileName(It.IsAny<ClimateVariable>()))
            .Returns(outputFileName);
        return mockDataset;
    }

    [Theory]
    [InlineData(NarClim2Domain.AUS18, NarClim2Constants.Files.RlonValuesFileAUS18)]
    [InlineData(NarClim2Domain.SEAus04, NarClim2Constants.Files.RlonValuesFileSEAus04)]
    public async Task GenerateVariableMergeScript_UsesCorrectPath(NarClim2Domain domain, string expectedFileName)
    {
        Mock<NarClim2Dataset> mockDataset = CreateMockDataset(domain: domain);
        string file = await generator.GenerateVariableMergeScript(mockDataset.Object, ClimateVariable.Temperature);
        string script = await File.ReadAllTextAsync(file);

        // The script should use the correct rlon values file for this domain.
        Assert.Contains("setvar.py", script);
        Assert.Contains("--var rlon", script);
        Assert.Contains(expectedFileName, script);
    }

    [Fact]
    public async Task GenerateVariableMergeScript_WithNonNarClim2Dataset_ThrowsArgumentException()
    {
        IClimateDataset mockDataset = new StaticMockDataset("/input");

        await Assert.ThrowsAsync<ArgumentException>(() => 
            generator.GenerateVariableMergeScript(mockDataset, ClimateVariable.ShortwaveRadiation));
    }

    [Fact]
    public async Task GenerateScriptsAsync_WithValidNarClim2Dataset_GeneratesCompleteScript()
    {
        NarClim2Dataset dataset = new NarClim2Dataset("/path/to/narclim2");

        // Attempting to generate scripts without setting up an appropriate
        // file tree is going to result in an exception.
        // TODO: is it worth setting up a suitable file tree?
        InvalidOperationException ex = await Assert.ThrowsAsync<InvalidOperationException>(() => generator.GenerateScriptsAsync(dataset));

        Assert.NotNull(ex);
        Assert.Contains("No input files found for variable", ex.Message);
    }
}
