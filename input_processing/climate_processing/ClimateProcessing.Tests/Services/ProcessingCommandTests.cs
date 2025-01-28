using Xunit;
using ClimateProcessing.Services;

namespace ClimateProcessing.Tests.Services;

public class ProcessingCommandTests
{
    [Fact]
    public void Skip_CreatesSymlink()
    {
        var command = ProcessingCommand.Skip("input.nc", "output.nc");
        
        Assert.Equal("ln -sf input.nc output.nc", command.Command);
        Assert.False(command.RequiresProcessing);
    }

    [Fact]
    public void Process_CreatesProcessingCommand()
    {
        var command = ProcessingCommand.Process("cdo -O expr,'x=y*2' input.nc output.nc");
        
        Assert.Equal("cdo -O expr,'x=y*2' input.nc output.nc", command.Command);
        Assert.True(command.RequiresProcessing);
    }

    [Theory]
    [InlineData("input.nc", "output.nc")]
    [InlineData("/path/to/input.nc", "../output.nc")]
    [InlineData("./input.nc", "/tmp/output.nc")]
    public void Skip_HandlesVariousPathFormats(string input, string output)
    {
        var command = ProcessingCommand.Skip(input, output);
        Assert.Equal($"ln -sf {input} {output}", command.Command);
    }
}
