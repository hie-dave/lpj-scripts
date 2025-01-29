using Xunit;
using ClimateProcessing.Units;

namespace ClimateProcessing.Tests.Units;

public class AggregationMethodTests
{
    [Theory]
    [InlineData(AggregationMethod.Mean, 24, "daymean")]
    [InlineData(AggregationMethod.Mean, 1, "timemean")]
    [InlineData(AggregationMethod.Mean, 3, "timemean")]
    [InlineData(AggregationMethod.Sum, 24, "daysum")]
    [InlineData(AggregationMethod.Sum, 1, "timesum")]
    [InlineData(AggregationMethod.Sum, 3, "timesum")]
    [InlineData(AggregationMethod.Maximum, 24, "daymax")]
    [InlineData(AggregationMethod.Maximum, 1, "timemax")]
    [InlineData(AggregationMethod.Maximum, 3, "timemax")]
    [InlineData(AggregationMethod.Minimum, 24, "daymin")]
    [InlineData(AggregationMethod.Minimum, 1, "timemin")]
    [InlineData(AggregationMethod.Minimum, 3, "timemin")]
    public void ToCdoOperator_ReturnsCorrectOperator(
        AggregationMethod method,
        int hours,
        string expectedOperator)
    {
        TimeStep timeStep = new TimeStep(hours);
        string result = method.ToCdoOperator(timeStep);
        Assert.Equal(expectedOperator, result);
    }

    [Fact]
    public void ToCdoOperator_ThrowsForInvalidMethod()
    {
        // Cast an invalid value to AggregationMethod
        var invalidMethod = (AggregationMethod)999;
        var timeStep = TimeStep.Daily;

        var ex = Assert.Throws<ArgumentException>(() =>
            invalidMethod.ToCdoOperator(timeStep));
        Assert.Contains("Unknown aggregation method", ex.Message);
    }
}
