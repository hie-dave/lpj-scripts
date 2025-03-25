using ClimateProcessing.Models;
using ClimateProcessing.Tests.Mocks;
using Xunit;

namespace ClimateProcessing.Tests.Models;

public class ProcessingConfigTests : IDisposable
{
    private readonly TestProcessingConfig _config;
    private readonly string _tempDir;

    public ProcessingConfigTests()
    {
        _tempDir = Path.Combine(Path.GetTempPath(), Path.GetRandomFileName());
        Directory.CreateDirectory(_tempDir);
        _config = new TestProcessingConfig
        {
            InputDirectory = _tempDir,
            Project = "TestProject",
            Memory = 1000,
            ChunkSizeTime = 10,
            ChunkSizeSpatial = 5,
            CompressionLevel = 5,
            Version = ModelVersion.Trunk
        };
    }

    [Fact]
    public void ValidateDirectories_WithNonExistentInputDirectory_ThrowsArgumentException()
    {
        _config.InputDirectory = Path.Combine(_tempDir, "nonexistent");
        var ex = Assert.Throws<ArgumentException>(() => _config.ValidateDirectories());
        Assert.Contains("Input directory does not exist", ex.Message);
    }

    [Fact]
    public void ValidateDirectories_WithNonExistentGridFile_ThrowsArgumentException()
    {
        _config.GridFile = Path.Combine(_tempDir, "nonexistent.grid");
        var ex = Assert.Throws<ArgumentException>(() => _config.ValidateDirectories());
        Assert.Contains("Grid file does not exist", ex.Message);
    }

    [Fact]
    public void ValidateBasicParameters_WithEmptyProject_ThrowsArgumentException()
    {
        _config.Project = "";
        var ex = Assert.Throws<ArgumentException>(() => _config.ValidateBasicParameters());
        Assert.Equal("Project must be specified", ex.Message);
    }

    [Theory]
    [InlineData(0)]
    [InlineData(10)]
    public void ValidateBasicParameters_WithInvalidCompressionLevel_ThrowsArgumentException(int level)
    {
        _config.CompressOutput = true;
        _config.CompressionLevel = level;
        var ex = Assert.Throws<ArgumentException>(() => _config.ValidateBasicParameters());
        Assert.Equal("Compression level must be between 1 and 9", ex.Message);
    }

    [Fact]
    public void ValidateBasicParameters_WithNegativeMemory_ThrowsArgumentException()
    {
        _config.Memory = -1;
        var ex = Assert.Throws<ArgumentException>(() => _config.ValidateBasicParameters());
        Assert.Equal("Memory must be greater than 0", ex.Message);
    }

    [Fact]
    public void ValidateBasicParameters_WithInvalidChunkSizeTime_ThrowsArgumentException()
    {
        _config.ChunkSizeTime = 0;
        var ex = Assert.Throws<ArgumentException>(() => _config.ValidateBasicParameters());
        Assert.Equal("Chunk size must be greater than 0", ex.Message);
    }

    [Fact]
    public void ValidateBasicParameters_WithInvalidChunkSizeSpatial_ThrowsArgumentException()
    {
        _config.ChunkSizeSpatial = 0;
        var ex = Assert.Throws<ArgumentException>(() => _config.ValidateBasicParameters());
        Assert.Equal("Chunk size for spatial variables must be greater than 0", ex.Message);
    }

    [Theory]
    [InlineData(12)]
    [InlineData(1)]
    public void ValidateTrunkTimeStepSettings_WithNonDailyInputTimeStep_ThrowsArgumentException(int hours)
    {
        _config.Version = ModelVersion.Trunk;
        _config.InputTimeStepHours = hours;
        var ex = Assert.Throws<ArgumentException>(() => _config.ValidateTrunkTimeStepSettings());
        Assert.Contains("Input timestep must be daily (24 hours)", ex.Message);
    }

    [Theory]
    [InlineData(12)]
    [InlineData(8)]
    public void ValidateTrunkTimeStepSettings_WithNonDailyOutputTimeStep_ThrowsArgumentException(int hours)
    {
        _config.Version = ModelVersion.Trunk;
        _config.OutputTimeStepHours = hours;
        var ex = Assert.Throws<ArgumentException>(() => _config.ValidateTrunkTimeStepSettings());
        Assert.Contains("Output timestep must be daily (24 hours)", ex.Message);
    }

    [Fact]
    public void ValidateDaveTimeStepSettings_WithDailyOutputTimeStep_ThrowsArgumentException()
    {
        _config.Version = ModelVersion.Dave;
        _config.OutputTimeStepHours = 24;
        var ex = Assert.Throws<ArgumentException>(() => _config.ValidateDaveTimeStepSettings());
        Assert.Equal("Output timestep must be subdaily when processing for Dave.", ex.Message);
    }

    [Fact]
    public void ValidateDaveTimeStepSettings_WithUnspecifiedInputTimeStep_ThrowsArgumentException()
    {
        _config.Version = ModelVersion.Dave;
        _config.InputTimeStepHours = 0;
        _config.OutputTimeStepHours = 3;
        var ex = Assert.Throws<ArgumentException>(() => _config.ValidateDaveTimeStepSettings());
        Assert.Equal("Input timestep must be specified when processing for Dave.", ex.Message);
    }

    [Fact]
    public void ValidateDaveTimeStepSettings_WithUnspecifiedOutputTimeStep_ThrowsArgumentException()
    {
        _config.Version = ModelVersion.Dave;
        _config.OutputTimeStepHours = 0;
        _config.InputTimeStepHours = 3;
        var ex = Assert.Throws<ArgumentException>(() => _config.ValidateDaveTimeStepSettings());
        Assert.Equal("Output timestep must be specified when processing for Dave.", ex.Message);
    }

    [Fact]
    public void ValidateTimeStepSettings_WithCoarserInputThanOutput_ThrowsArgumentException()
    {
        _config.Version = ModelVersion.Dave;
        _config.InputTimeStepHours = 6;
        _config.OutputTimeStepHours = 12;
        _config.ValidateTimeStepSettings();  // Should not throw

        _config.InputTimeStepHours = 12;
        _config.OutputTimeStepHours = 6;
        var ex = Assert.Throws<ArgumentException>(() => _config.ValidateTimeStepSettings());
        Assert.Equal("Input timestep cannot be coarser than the output timestep.", ex.Message);
    }

    [Fact]
    public void Validate_WithValidConfig_DoesNotThrow()
    {
        var tempGridFile = Path.Combine(_tempDir, "test.grid");
        File.WriteAllText(tempGridFile, "test");
        
        _config.GridFile = tempGridFile;
        _config.Version = ModelVersion.Trunk;
        _config.InputTimeStepHours = 24;
        _config.OutputTimeStepHours = 24;

        var exception = Record.Exception(() => _config.Validate());
        Assert.Null(exception);
    }

    public void Dispose()
    {
        try
        {
            if (Directory.Exists(_tempDir))
                Directory.Delete(_tempDir, true);
        }
        catch
        {
            // Ignore cleanup errors
        }
    }
}
