using Xunit;
using ClimateProcessing.Services;
using ClimateProcessing.Models;
using ClimateProcessing.Units;

namespace ClimateProcessing.Tests.Services;

public class ScriptGeneratorTests
{
    private readonly ProcessingConfig _config = new()
    {
        JobName = "test",
        Project = "test",
        Queue = "normal",
        Walltime = "01:00:00",
        Ncpus = 1,
        Memory = 4,
        OutputDirectory = "/output",
        InputTimeStep = TimeStep.Hourly,
        OutputTimeStep = TimeStep.Daily
    };

    private readonly ScriptGenerator _generator;

    public ScriptGeneratorTests()
    {
        _generator = new ScriptGenerator(_config);
    }

    [Theory]
    [InlineData("K", "K", "temp", "temp", false)]  // No conversion needed
    [InlineData("K", "degC", "temp", "temp", true)]  // Conversion needed
    [InlineData("W/m2", "W m-2", "rad", "rad", true)]  // Only renaming needed
    public void GenerateUnitConversionCommand_DetectsProcessingNeeds(
        string inputUnits,
        string targetUnits,
        string inputVar,
        string outputVar,
        bool requiresProcessing)
    {
        var command = _generator.GenerateUnitConversionCommand(
            inputVar,
            outputVar,
            inputUnits,
            targetUnits,
            TimeStep.Hourly,
            "input.nc",
            "output.nc");

        Assert.Equal(requiresProcessing, command.RequiresProcessing);
    }

    [Theory]
    [InlineData(1, 24, true)]   // Hourly to daily
    [InlineData(24, 24, false)] // Daily to daily
    [InlineData(3, 24, true)]   // 3-hourly to daily
    public void GenerateTimeAggregationCommand_DetectsProcessingNeeds(
        int inputHours,
        int outputHours,
        bool requiresProcessing)
    {
        var config = _config with
        {
            InputTimeStep = new TimeStep(inputHours),
            OutputTimeStep = new TimeStep(outputHours)
        };
        var generator = new ScriptGenerator(config);

        var command = generator.GenerateTimeAggregationCommand(
            ClimateVariable.Temperature,
            "input.nc",
            "output.nc");

        Assert.Equal(requiresProcessing, command.RequiresProcessing);
    }

    [Theory]
    [InlineData("tas", "tas", false)]      // Same name
    [InlineData("temp", "tas", true)]      // Different name
    [InlineData("pr", "prec", true)]       // Different name
    public void GenerateVariableRenameCommand_DetectsProcessingNeeds(
        string inputVar,
        string outputVar,
        bool requiresProcessing)
    {
        var command = _generator.GenerateVariableRenameCommand(
            inputVar,
            outputVar,
            "input.nc",
            "output.nc");

        Assert.Equal(requiresProcessing, command.RequiresProcessing);
    }
}
