using ClimateProcessing.Models;
using ClimateProcessing.Services;
using Moq;
using Xunit;

namespace ClimateProcessing.Tests.Services;

public class PathManagerTests : IDisposable
{
    private readonly string _testOutputPath;
    private readonly PathManager _pathManager;
    private readonly Mock<IClimateDataset> _mockDataset;

    public PathManagerTests()
    {
        _testOutputPath = Path.Combine(Path.GetTempPath(), Path.GetRandomFileName());
        _pathManager = new PathManager(_testOutputPath);
        _mockDataset = new Mock<IClimateDataset>();
        Directory.CreateDirectory(_testOutputPath);
    }

    public void Dispose()
    {
        if (Directory.Exists(_testOutputPath))
        {
            Directory.Delete(_testOutputPath, true);
        }
    }

    [Theory]
    [InlineData(PathType.Log)]
    [InlineData(PathType.Script)]
    [InlineData(PathType.Stream)]
    public void GetDatasetPath_ThrowsForInvalidPathTypes(PathType pathType)
    {
        _mockDataset.Setup(d => d.GetOutputDirectory()).Returns("test_dataset");
        var ex = Assert.Throws<ArgumentException>(() => _pathManager.GetDatasetPath(_mockDataset.Object, pathType));
        Assert.Contains("not valid at the dataset-level", ex.Message);
    }

    [Theory]
    [InlineData(PathType.Output)]
    [InlineData(PathType.Working)]
    public void GetDatasetPath_ReturnsCorrectPathForValidTypes(PathType pathType)
    {
        const string datasetDir = "test_dataset";
        _mockDataset.Setup(d => d.GetOutputDirectory()).Returns(datasetDir);

        string result = _pathManager.GetDatasetPath(_mockDataset.Object, pathType);

        string expectedBaseDir = pathType == PathType.Output ? "output" : "tmp";
        string expectedPath = Path.Combine(_testOutputPath, expectedBaseDir, datasetDir);
        Assert.Equal(expectedPath, result);
        Assert.True(Directory.Exists(result));
    }

    [Fact]
    public void GetDatasetFileName_ReturnsCorrectPath()
    {
        const string datasetDir = "test_dataset";
        const string fileName = "test_file.nc";
        _mockDataset.Setup(d => d.GetOutputDirectory()).Returns(datasetDir);
        _mockDataset.Setup(d => d.GenerateOutputFileName(It.IsAny<ClimateVariable>()))
            .Returns(fileName);

        string result = _pathManager.GetDatasetFileName(_mockDataset.Object, ClimateVariable.Precipitation, PathType.Output);

        string expectedPath = Path.Combine(_testOutputPath, "output", datasetDir, fileName);
        Assert.Equal(expectedPath, result);
        Assert.True(Directory.Exists(Path.GetDirectoryName(result)));
    }

    [Fact]
    public void GetChecksumFilePath_ReturnsCorrectPath()
    {
        string result = _pathManager.GetChecksumFilePath();

        string expectedPath = Path.Combine(_testOutputPath, "output", "sha512sums.txt");
        Assert.Equal(expectedPath, result);
    }

    [Fact]
    public void CreateDirectoryTree_CreatesAllRequiredDirectories()
    {
        const string datasetDir = "test_dataset";
        _mockDataset.Setup(d => d.GetOutputDirectory()).Returns(datasetDir);

        _pathManager.CreateDirectoryTree(_mockDataset.Object);

        // Check top-level directories
        Assert.True(Directory.Exists(Path.Combine(_testOutputPath, "logs")));
        Assert.True(Directory.Exists(Path.Combine(_testOutputPath, "scripts")));
        Assert.True(Directory.Exists(Path.Combine(_testOutputPath, "streams")));
        Assert.True(Directory.Exists(Path.Combine(_testOutputPath, "output")));
        Assert.True(Directory.Exists(Path.Combine(_testOutputPath, "tmp")));

        // Check dataset-specific directories
        Assert.True(Directory.Exists(Path.Combine(_testOutputPath, "output", datasetDir)));
        Assert.True(Directory.Exists(Path.Combine(_testOutputPath, "tmp", datasetDir)));
    }
}
