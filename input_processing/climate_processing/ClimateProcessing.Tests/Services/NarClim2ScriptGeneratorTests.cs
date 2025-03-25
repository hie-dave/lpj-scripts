using Xunit;
using ClimateProcessing.Models;
using ClimateProcessing.Tests.Mocks;
using Moq;

namespace ClimateProcessing.Tests.Services;

public class NarClim2ScriptGeneratorTests : IDisposable
{
    private const string outputDirectoryPrefix = "narclim2_script_generator_tests_output";
    private readonly string outputDirectory;
    private readonly NarClim2Config config;
    private readonly TestNarClim2ScriptGenerator generator;

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
            CompressionLevel = 4
        };

        generator = new TestNarClim2ScriptGenerator(config);
    }

    public void Dispose()
    {
        if (Directory.Exists(outputDirectory))
            Directory.Delete(outputDirectory, true);
    }

    [Fact]
    public async Task WritePreMerge_WithValidNarClim2Dataset_GeneratesCorrectRlonCorrectionScript()
    {
        NarClim2Dataset dataset = new NarClim2Dataset("/path/to/narclim2");
        using StringWriter writer = new();

        await generator.PreMergeAsync(writer, dataset, ClimateVariable.MaxTemperature);
        string script = writer.ToString();

        Assert.Contains("setvar.py", script);
        Assert.Contains("--var rlon", script);
    }

    [Fact]
    public async Task WritePreMerge_WithNonNarClim2Dataset_ThrowsArgumentException()
    {
        // Arrange
        Mock<IClimateDataset> mockDataset = new();
        using StringWriter writer = new();

        // Act & Assert
        await Assert.ThrowsAsync<ArgumentException>(() => 
            generator.PreMergeAsync(writer, mockDataset.Object, ClimateVariable.MaxTemperature));
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

    [Theory]
    [InlineData(NarClim2Domain.AUS18, NarClim2Constants.Files.RlonValuesFileAUS18)]
    [InlineData(NarClim2Domain.SEAus04, NarClim2Constants.Files.RlonValuesFileSEAus04)]
    public void ReadRlonValuesFile_ReturnsCorrectPath(NarClim2Domain domain, string expectedFileName)
    {
        string basePath = "/path/to/narclim2";
        NarClim2Dataset dataset = new NarClim2Dataset(basePath, domain);

        string result = generator.ReadRlonValuesFile(dataset);

        Assert.NotNull(result);
        Assert.EndsWith(expectedFileName, result);
        Assert.StartsWith(basePath, result);
    }
}
