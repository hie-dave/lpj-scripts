using ClimateProcessing.Models;
using System.Linq;
using Xunit;

namespace ClimateProcessing.Tests.Models;

public class PBSStorageHelperTests
{
    [Theory]
    [InlineData("/g/data/asdf/x.y", "gdata", "asdf")]
    [InlineData("/g/data/p123/some/path", "gdata", "p123")]
    [InlineData("/g/data/wx789/data.nc", "gdata", "wx789")]
    public void GetStorageDirective_ValidGdataPath_ReturnsCorrectDirective(string path, string expectedKind, string expectedId)
    {
        var directive = PBSStorageHelper.GetStorageDirective(path);
        
        Assert.NotNull(directive);
        Assert.Equal(path, directive.Path);
        Assert.Equal(expectedKind, directive.Kind);
        Assert.Equal(expectedId, directive.StorageIdentifier);
    }

    [Theory]
    [InlineData("/scratch/proj123/work", "scratch", "proj123")]
    [InlineData("/scratch/test/output.txt", "scratch", "test")]
    [InlineData("/scratch/xyz/nested/path", "scratch", "xyz")]
    public void GetStorageDirective_ValidScratchPath_ReturnsCorrectDirective(string path, string expectedKind, string expectedId)
    {
        var directive = PBSStorageHelper.GetStorageDirective(path);
        
        Assert.NotNull(directive);
        Assert.Equal(path, directive.Path);
        Assert.Equal(expectedKind, directive.Kind);
        Assert.Equal(expectedId, directive.StorageIdentifier);
    }

    [Theory]
    [InlineData("/data")]
    [InlineData("/home/user/file.txt")]
    [InlineData("/g/data")] // Missing project ID
    [InlineData("/scratch")] // Missing project ID
    [InlineData("/g/data/")] // Missing project ID and path
    [InlineData("/scratch/")] // Missing project ID and path
    public void GetStorageDirective_InvalidPath_ReturnsNull(string path)
    {
        var directive = PBSStorageHelper.GetStorageDirective(path);
        Assert.Null(directive);
    }

    [Fact]
    public void GetStorageDirectives_MultiplePathsSameProject_ReturnsSingleDirective()
    {
        var paths = new[]
        {
            "/g/data/proj1/file1.txt",
            "/g/data/proj1/file2.txt",
            "/g/data/proj1/subdir/file3.txt"
        };

        var directives = PBSStorageHelper.GetStorageDirectives(paths).ToList();
        
        Assert.Single(directives);
        Assert.Equal("gdata", directives[0].Kind);
        Assert.Equal("proj1", directives[0].StorageIdentifier);
    }

    [Fact]
    public void GetStorageDirectives_MultipleProjects_ReturnsUniqueDirectives()
    {
        var paths = new[]
        {
            "/g/data/proj1/file1.txt",
            "/g/data/proj2/file2.txt",
            "/scratch/proj3/file3.txt"
        };

        var directives = PBSStorageHelper.GetStorageDirectives(paths).ToList();
        
        Assert.Equal(3, directives.Count);
        Assert.Contains(directives, d => d.Kind == "gdata" && d.StorageIdentifier == "proj1");
        Assert.Contains(directives, d => d.Kind == "gdata" && d.StorageIdentifier == "proj2");
        Assert.Contains(directives, d => d.Kind == "scratch" && d.StorageIdentifier == "proj3");
    }

    [Fact]
    public void FormatStorageDirectives_MultipleDirectives_FormatsCorrectly()
    {
        var directives = new[]
        {
            new PBSStorageDirective("/g/data/proj1/file.txt", "gdata", "proj1"),
            new PBSStorageDirective("/scratch/test/file.txt", "scratch", "test")
        };

        var formatted = PBSStorageHelper.FormatStorageDirectives(directives);
        
        Assert.Equal("#PBS -l storage=gdata/proj1+scratch/test", formatted);
    }
}
