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
                Assert.StartsWith("-timselmean", @operator);
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

    [Theory]
    [InlineData("kg m-2 s-1", true)]  // Standard notation
    [InlineData("kg/m2/s", true)]     // Division notation
    [InlineData("kg.m-2.s-1", true)]  // Dot notation
    [InlineData("kg m^-2 s^-1", true)] // Caret notation
    [InlineData("KG M-2 S-1", true)]   // Case insensitive
    [InlineData("kg  m-2  s-1", true)] // Extra spaces
    [InlineData("kgm-2s-1", true)]     // No separators
    [InlineData("W", false)]           // No per-area units
    [InlineData("kg s-1", false)]      // Time only, no area
    [InlineData("", false)]            // Empty string
    [InlineData("m2", false)]          // Area but not per-area
    [InlineData("kg/s/m2", true)]      // Different order
    public void HasPerAreaUnits_DetectsUnitsCorrectly(string units, bool expectedResult)
    {
        // Now we can test HasPerAreaUnits directly since it's internal
        var result = ScriptGenerator.HasPerAreaUnits(units);
        Assert.Equal(expectedResult, result);
    }

    [Theory]
    [InlineData(ClimateVariable.Temperature)]
    [InlineData(ClimateVariable.SpecificHumidity)]
    public void GetRemapAlgorithm_NonAreaSensitiveVariables_AlwaysReturnsBilinear(ClimateVariable variable)
    {
        var info = new VariableInfo(Enum.GetName(variable)!, "any_units");
        var result = _generator.GetInterpolationAlgorithm(info, variable);
        
        Assert.Equal(InterpolationAlgorithm.Bilinear, result);
    }

    [Theory]
    [InlineData(ClimateVariable.Precipitation)]
    [InlineData(ClimateVariable.ShortwaveRadiation)]
    public void GetRemapAlgorithm_AreaSensitiveVariables_DependsOnUnits(ClimateVariable variable)
    {
        // Should use conservative remapping when not per-area
        var nonAreaInfo = new VariableInfo(Enum.GetName(variable)!, "W");
        var nonAreaResult = _generator.GetInterpolationAlgorithm(nonAreaInfo, variable);
        Assert.Equal(InterpolationAlgorithm.Conservative, nonAreaResult);

        // Should use bilinear remapping when per-area
        var areaInfo = new VariableInfo(Enum.GetName(variable)!, "W m-2");
        var areaResult = _generator.GetInterpolationAlgorithm(areaInfo, variable);
        Assert.Equal(InterpolationAlgorithm.Bilinear, areaResult);
    }
}
