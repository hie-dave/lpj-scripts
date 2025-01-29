using Xunit;
using ClimateProcessing.Units;

namespace ClimateProcessing.Tests.Units;

public class UnitConverterTests
{
    [Theory]
    [InlineData("W/m2", "W/m2", false, false, false)]  // Exact match
    [InlineData("W/m2", "W m-2", false, true, false)]  // Different notation, same meaning
    [InlineData("K", "degC", true, true, false)]       // Requires conversion, no timestep
    [InlineData("kg m-2 s-1", "mm", true, true, true)] // Requires conversion with timestep
    public void AnalyzeConversion_HandlesVariousUnitCombinations(
        string inputUnits,
        string targetUnits,
        bool expectedRequiresConversion,
        bool expectedRequiresRenaming,
        bool expectedRequiresTimeStep)
    {
        var result = UnitConverter.AnalyzeConversion(inputUnits, targetUnits);
        
        Assert.Equal(expectedRequiresConversion, result.RequiresConversion);
        Assert.Equal(expectedRequiresRenaming, result.RequiresRenaming);
        Assert.Equal(expectedRequiresTimeStep, result.RequiresTimeStep);
        
        if (expectedRequiresConversion)
        {
            Assert.NotNull(result.ConversionExpression);
        }
    }

    [Theory]
    [InlineData("K", "degC", 1, "subc,273.15")]
    [InlineData("kg m-2 s-1", "mm", 24, "mulc,86400")] // Daily accumulation
    [InlineData("kg m-2 s-1", "mm", 3, "mulc,10800")]  // 3-hourly accumulation
    [InlineData("kg m-2 s-1", "mm", 1, "mulc,3600")]   // Hourly accumulation
    public void GenerateConversionExpression_GeneratesCorrectExpressions(
        string inputUnits,
        string targetUnits,
        int hours,
        string expectedExpression)
    {
        TimeStep timeStep = new TimeStep(hours);
        string expression = UnitConverter.GenerateConversionExpression(inputUnits, targetUnits, timeStep);
        Assert.Equal(expectedExpression, expression);
    }

    [Theory]
    [InlineData("J m-2", "W m-2")]      // Invalid conversion
    [InlineData("unknown", "W m-2")]     // Unknown unit
    public void AnalyzeConversion_ThrowsForUnsupportedConversion(
        string inputUnits,
        string targetUnits
    )
    {
        Assert.Throws<ArgumentException>(() => 
            UnitConverter.AnalyzeConversion(inputUnits, targetUnits));
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
    [InlineData(1, 3600)]       // Hourly
    [InlineData(3, 10800)]      // 3-hourly
    [InlineData(24, 86400)]     // Daily
    public void TimeStep_GetSecondsInPeriod_ReturnsCorrectValue(int hours, int expectedSeconds)
    {
        var timeStep = new TimeStep(hours);
        Assert.Equal(expectedSeconds, timeStep.GetSecondsInPeriod());
    }

    [Fact]
    public void TimeStep_StaticProperties_HaveCorrectValues()
    {
        Assert.Equal(1, TimeStep.Hourly.Hours);
        Assert.Equal(24, TimeStep.Daily.Hours);
    }
}
