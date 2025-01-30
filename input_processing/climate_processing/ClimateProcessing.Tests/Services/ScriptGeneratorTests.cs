using Xunit;
using ClimateProcessing.Services;
using ClimateProcessing.Models;
using ClimateProcessing.Units;

namespace ClimateProcessing.Tests.Services;

public class ScriptGeneratorTests
{
    private readonly ProcessingConfig _config = new()
    {
        Project = "test",
        Queue = "normal",
        Walltime = "01:00:00",
        Ncpus = 1,
        Memory = 4,
        OutputDirectory = "/output",
        InputTimeStepHours = 1,
        OutputTimeStepHours = 24
    };

    private readonly ScriptGenerator _generator;

    public ScriptGeneratorTests()
    {
        _generator = new ScriptGenerator(_config);
    }

    [Theory]
    [InlineData("K", "K", "temp", false, false)]  // No conversion needed
    [InlineData("K", "degC", "temp", true, true)]  // Conversion needed
    [InlineData("W/m2", "W m-2", "rad", false, true)]  // Only renaming needed
    public void GenerateUnitConversionOperators_GeneratesCorrectOperators(
        string inputUnits,
        string targetUnits,
        string outputVar,
        bool expectsConversion,
        bool expectsRenaming)
    {
        var operators = _generator.GenerateUnitConversionOperators(
            outputVar,
            inputUnits,
            targetUnits,
            TimeStep.Hourly).ToList();

        if (expectsConversion)
        {
            Assert.Contains(operators, op => op.StartsWith("-expr"));
        }
        if (expectsRenaming)
        {
            Assert.Contains(operators, op => op.StartsWith("-setattribute"));
        }
        if (!expectsConversion && !expectsRenaming)
        {
            Assert.Empty(operators);
        }
    }

    [Theory]
    [InlineData(1, 24, true)]   // Hourly to daily
    [InlineData(24, 24, false)] // Daily to daily
    [InlineData(3, 24, true)]   // 3-hourly to daily
    public void GenerateTimeAggregationOperators_GeneratesCorrectOperators(
        int inputHours,
        int outputHours,
        bool requiresAggregation)
    {
        var config = new ProcessingConfig
        {
            Project = _config.Project,
            Queue = _config.Queue,
            Walltime = _config.Walltime,
            Ncpus = _config.Ncpus,
            Memory = _config.Memory,
            OutputDirectory = _config.OutputDirectory,
            InputTimeStepHours = inputHours,
            OutputTimeStepHours = outputHours
        };
        var generator = new ScriptGenerator(config);

        string @operator = generator.GenerateTimeAggregationOperator(
            ClimateVariable.Temperature);

        if (requiresAggregation)
        {
            if (outputHours == 24)
                Assert.StartsWith("-daymean", @operator);
            else
                Assert.StartsWith("-timemean", @operator);
        }
        else
        {
            Assert.Empty(@operator);
        }
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
        string op = _generator.GenerateRenameOperator(inputVar, outputVar);
        Assert.Equal(requiresProcessing, !string.IsNullOrEmpty(op));
    }
}
