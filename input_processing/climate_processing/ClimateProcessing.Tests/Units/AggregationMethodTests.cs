using Xunit;
using ClimateProcessing.Units;

namespace ClimateProcessing.Tests.Units;

public class AggregationMethodTests
{
    [Theory]
    [InlineData(AggregationMethod.Mean, 24, "daymean")]
    [InlineData(AggregationMethod.Mean, 1, "timmean")]
    [InlineData(AggregationMethod.Mean, 3, "timmean")]
    [InlineData(AggregationMethod.Sum, 24, "daysum")]
    [InlineData(AggregationMethod.Sum, 1, "timsum")]
    [InlineData(AggregationMethod.Sum, 3, "timsum")]
    [InlineData(AggregationMethod.Maximum, 24, "daymax")]
    [InlineData(AggregationMethod.Maximum, 1, "timmax")]
    [InlineData(AggregationMethod.Maximum, 3, "timmax")]
    [InlineData(AggregationMethod.Minimum, 24, "daymin")]
    [InlineData(AggregationMethod.Minimum, 1, "timmin")]
    [InlineData(AggregationMethod.Minimum, 3, "timmin")]
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
