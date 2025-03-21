using Xunit;
using ClimateProcessing.Services;
using ClimateProcessing.Models;
using ClimateProcessing.Units;
using System.Text.RegularExpressions;
using Microsoft.VisualStudio.TestPlatform.ObjectModel.Engine.ClientProtocol;
using ClimateProcessing.Tests.Mocks;
using System.Reflection;
using ClimateProcessing.Configuration;

namespace ClimateProcessing.Tests.Services;

public class PBSWriterTests : IDisposable
{
    private const string outputDirectoryPrefix = "climate_processing_pbs_writer_tests";
    private readonly string outputDirectory;
    private readonly PathManager pathManager;

    public PBSWriterTests()
    {
        outputDirectory = CreateOutputDirectory();
        pathManager = new PathManager(outputDirectory);
    }

    /// <summary>
    /// Teardown method - delete the temporary output directory.
    /// </summary>
    public void Dispose()
    {
        try
        {
            Directory.Delete(outputDirectory, true);
        }
        catch (IOException error)
        {
            // Log an error but don't throw an exception.
            Console.Error.WriteLine($"Warning: could not delete temporary output directory used by {GetType().Name}: {error}");
        }
    }

    private static string CreateOutputDirectory()
    {
        return Directory.CreateTempSubdirectory(outputDirectoryPrefix).FullName;
    }

    [Theory]
    [InlineData(1, 8, 100, "normal", "01:00:00", "test_job", "p123", "test@example.com")]
    [InlineData(2, 192, 1024, "hugemem", "48:00:00", "asdf_name", "x987", null)]
    public async Task WritesPBSHeader(int ncpu, int memory, int jobfs, string queue, string walltime, string jobName, string project, string? email)
    {
        // Arrange
        StringWriter writer = new();
        PBSConfig config = new(
            queue,
            ncpu,
            memory,
            jobfs,
            project,
            TimeSpan.ParseExact(walltime, @"hh\:mm\:ss", null),
            email);
        PathManager pathManager = new PathManager(outputDirectory);
        PBSWriter generator = new(config, pathManager);

        // Act
        await generator.WritePBSHeader(writer, jobName, Array.Empty<PBSStorageDirective>());
        string result = writer.ToString();

        // Assert
        Assert.Contains($"#PBS -N {jobName}", result);
        Assert.Contains($"#PBS -P {project}", result);
        Assert.Contains($"#PBS -q {queue}", result);
        Assert.Contains($"#PBS -l walltime={walltime}", result);
        Assert.Contains($"#PBS -l ncpus={ncpu}", result);
        Assert.Contains($"#PBS -l mem={memory}GB", result);
        Assert.Contains($"#PBS -l jobfs={jobfs}GB", result);
        Assert.Contains("#PBS -j oe", result);

        if (email != null)
        {
            Assert.Contains($"#PBS -M {email}", result);
            Assert.Contains("#PBS -m abe", result);
        }
        else
        {
            Assert.DoesNotContain("#PBS -M", result);
            Assert.DoesNotContain("#PBS -m", result);
        }

        // Verify proper line endings.
        Assert.DoesNotContain(result, "\r");

        // Verify proper line endings and no empty lines with whitespace.
        string[] lines = result.Split("\n");

        // Get index of last line starting with #PBS.
        string lastPBSLine = lines.Last(line => line.StartsWith("#PBS "));
        int lastPBSLineIndex = Array.LastIndexOf(lines, lastPBSLine);

        // Assert that all lines before the last #PBS line are not empty.
        Assert.All(lines[0..lastPBSLineIndex], Assert.NotEmpty);
    }

    [Theory]
    [InlineData(0, "/home/user/input", "/home/user/grid.txt", "/home/user/output")]  // No storage directives needed
    [InlineData(1, "/g/data/proj1/input", "/home/user/grid.txt", "/home/user/output")]  // Single gdata directive
    [InlineData(1, "/home/user/input", "/scratch/test/grid.txt", "/home/user/output")]  // Single scratch directive
    [InlineData(2, "/g/data/proj1/input", "/home/user/grid.txt", "/scratch/test/output")]  // Both gdata and scratch
    public async Task WritePBSHeader_GeneratesCorrectHeader(
        int expectedDirectives,
        params string[] filePaths)
    {
        // Arrange
        StringWriter writer = new();
        PBSConfig config = new(
            "normal",
            2,
            8,
            100,
            "p123",
            TimeSpan.FromHours(1)
        );
        PathManager pathManager = new PathManager(outputDirectory);
        PBSWriter generator = new(config, pathManager);
        IEnumerable<PBSStorageDirective> directives =
            PBSStorageHelper.GetStorageDirectives(filePaths);

        // Act
        await generator.WritePBSHeader(writer, "test_job", directives);
        string result = writer.ToString();

        // Assert
        // Verify basic header elements are present
        Assert.Contains("#PBS -P p123", result);
        Assert.Contains("#PBS -q normal", result);
        Assert.Contains("#PBS -l walltime=01:00:00", result);
        Assert.Contains("#PBS -l ncpus=2", result);
        Assert.Contains("#PBS -l mem=8GB", result);
        Assert.Contains("#PBS -N test_job", result);

        // Verify storage directive presence and format
        string prefix = "#PBS -l storage=";
        var storageLines = result.Split('\n').Where(l => l.StartsWith(prefix)).ToList();
        var needsStorage = expectedDirectives > 0;

        if (needsStorage)
        {
            Assert.Single(storageLines);
            string storageLine = storageLines[0];
            string[] parts = storageLine.Replace(prefix, string.Empty).Split('+');
            Assert.Equal(expectedDirectives, parts.Length);
        }
        else
        {
            Assert.Empty(storageLines);
        }
    }
}
