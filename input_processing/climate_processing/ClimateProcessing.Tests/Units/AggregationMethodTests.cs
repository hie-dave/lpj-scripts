using Xunit;
using ClimateProcessing.Units;

namespace ClimateProcessing.Tests.Units;

public class AggregationMethodTests
{
    [Theory]
    [InlineData(AggregationMethod.Mean, 24, "daymean")]
    [InlineData(AggregationMethod.Mean, 1, "timselmean")]
    [InlineData(AggregationMethod.Mean, 3, "timselmean")]
    [InlineData(AggregationMethod.Sum, 24, "daysum")]
    [InlineData(AggregationMethod.Sum, 1, "timselsum")]
    [InlineData(AggregationMethod.Sum, 3, "timselsum")]
    [InlineData(AggregationMethod.Maximum, 24, "daymax")]
    [InlineData(AggregationMethod.Maximum, 1, "timselmax")]
    [InlineData(AggregationMethod.Maximum, 3, "timselmax")]
    [InlineData(AggregationMethod.Minimum, 24, "daymin")]
    [InlineData(AggregationMethod.Minimum, 1, "timselmin")]
    [InlineData(AggregationMethod.Minimum, 3, "timselmin")]
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
