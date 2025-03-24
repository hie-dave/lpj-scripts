using Xunit;
using ClimateProcessing.Units;

namespace ClimateProcessing.Tests.Units;

public class TimeStepTests
{
    [Fact]
    public void TestStaticConstructors()
    {
        TimeStep hourly = TimeStep.Hourly;
        Assert.Equal(1, hourly.Hours);

        TimeStep threeHour = TimeStep.ThreeHourly;
        Assert.Equal(3, threeHour.Hours);

        TimeStep oneDay = TimeStep.Daily;
        Assert.Equal(24, oneDay.Hours);
    }

    [Theory]
    [InlineData(1, "1hour")]
    [InlineData(3, "3hour")]
    [InlineData(12, "12hour")]
    [InlineData(24, "day")]
    public void TestToString(
        int hours,
        string expected)
    {
        TimeStep timeStep = new TimeStep(hours);
        Assert.Equal(expected, timeStep.ToString());
    }
}
