using Xunit;
using ClimateProcessing.Units;

namespace ClimateProcessing.Tests.Units;

public class UnitConverterTests
{
    [Theory]
    [InlineData("W/m2", "W/m2", false, false)]  // Exact match
    [InlineData("W/m2", "W m-2", false, true)]  // Different notation, same meaning
    [InlineData("W m-2", "W/m2", false, true)]  // Different notation, same meaning (reverse)
    [InlineData("kg/kg", "1", false, true)]     // Different notation, same meaning
    [InlineData("K", "degC", true, true)]       // Requires conversion
    public void AnalyzeConversion_HandlesVariousUnitCombinations(
        string inputUnits,
        string targetUnits,
        bool expectedRequiresConversion,
        bool expectedRequiresRenaming)
    {
        var result = UnitConverter.AnalyzeConversion(inputUnits, targetUnits);
        
        Assert.Equal(expectedRequiresConversion, result.RequiresConversion);
        Assert.Equal(expectedRequiresRenaming, result.RequiresRenaming);
    }

    [Theory]
    [InlineData("kg kg-1", "1")]        // Space-separated with hyphen
    [InlineData("kg/kg", "1")]          // Forward slash notation
    [InlineData("W/m^2", "W m-2")]      // Power notation to hyphen
    [InlineData("W/m2", "W m-2")]       // Forward slash to hyphen
    public void AnalyzeConversion_RecognizesEquivalentUnits(string units1, string units2)
    {
        var result = UnitConverter.AnalyzeConversion(units1, units2);
        Assert.False(result.RequiresConversion);
        Assert.True(result.RequiresRenaming);
    }

    [Theory]
    [InlineData("K", "degC", "test_in", "test_out", null, "test_out=test_in-273.15")]
    [InlineData("kg m-2 s-1", "mm", "pr_in", "pr_out", 24, "pr_out=pr_in*86400")] // Daily accumulation
    [InlineData("kg m-2 s-1", "mm", "pr_in", "pr_out", 3, "pr_out=pr_in*10800")]  // 3-hourly accumulation
    [InlineData("kg m-2 s-1", "mm", "pr_in", "pr_out", 1, "pr_out=pr_in*3600")]   // Hourly accumulation
    public void GenerateConversionExpression_GeneratesCorrectExpressions(
        string inputUnits,
        string targetUnits,
        string inputVar,
        string outputVar,
        int? hours,
        string expectedExpression)
    {
        var timeStep = hours.HasValue ? new TimeStep(hours.Value) : null;
        var expression = UnitConverter.GenerateConversionExpression(inputVar, outputVar, inputUnits, targetUnits, timeStep);
        Assert.Equal(expectedExpression, expression);
    }

    [Theory]
    [InlineData("J m-2", "W m-2")]  // Invalid conversion
    public void AnalyzeConversion_ThrowsForUnsupportedConversion(
        string inputUnits,
        string targetUnits
    )
    {
        Assert.Throws<ArgumentException>(() => 
            UnitConverter.AnalyzeConversion(inputUnits, targetUnits));
    }

    [Theory]
    [InlineData("W m-2", "W m-2")]  // Exact match
    [InlineData("kg kg-1", "kg kg-1")]  // Exact match with spaces
    public void AnalyzeConversion_HandlesExactMatches(string input, string target)
    {
        var result = UnitConverter.AnalyzeConversion(input, target);
        Assert.False(result.RequiresConversion);
        Assert.False(result.RequiresRenaming);
    }

    [Theory]
    [InlineData(0)]    // Too small
    [InlineData(25)]   // Too large
    [InlineData(7)]    // Doesn't divide 24 evenly
    public void TimeStep_ThrowsForInvalidHours(int hours)
    {
        Assert.Throws<ArgumentException>(() => new TimeStep(hours));
    }

    [Theory]
    [InlineData(1)]    // Hourly
    [InlineData(2)]    // 2-hourly
    [InlineData(3)]    // 3-hourly
    [InlineData(4)]    // 4-hourly
    [InlineData(6)]    // 6-hourly
    [InlineData(8)]    // 8-hourly
    [InlineData(12)]   // 12-hourly
    [InlineData(24)]   // Daily
    public void TimeStep_AcceptsValidHours(int hours)
    {
        var timeStep = new TimeStep(hours);
        Assert.Equal(hours * 3600, timeStep.GetSecondsInPeriod());
    }

    [Fact]
    public void TimeStep_StaticProperties_HaveCorrectValues()
    {
        Assert.Equal(1, TimeStep.Hourly.Hours);
        Assert.Equal(3, TimeStep.ThreeHourly.Hours);
        Assert.Equal(24, TimeStep.Daily.Hours);
    }
}
