using Xunit;
using ClimateProcessing.Models;

namespace ClimateProcessing.Tests.Models;

public class NarClim2ConstantsTests
{
    [Theory]
    [InlineData("1hr", NarClim2Frequency.Hour1)]
    [InlineData("3hr", NarClim2Frequency.Hour3)]
    [InlineData("day", NarClim2Frequency.Day)]
    [InlineData("mon", NarClim2Frequency.Month)]
    public void FrequencyNames_FromString_ReturnsCorrectEnum(string input, NarClim2Frequency expected)
    {
        var result = NarClim2Constants.FrequencyNames.FromString(input);
        Assert.Equal(expected, result);
    }

    [Theory]
    [InlineData("")]
    [InlineData("invalid")]
    [InlineData("hr")]
    [InlineData("month")]
    public void FrequencyNames_FromString_ThrowsForInvalidInput(string input)
    {
        var ex = Assert.Throws<ArgumentException>(() =>
            NarClim2Constants.FrequencyNames.FromString(input));
        Assert.Contains("Unknown frequency string", ex.Message);
    }

    [Theory]
    [InlineData(NarClim2Frequency.Hour1, "1hr")]
    [InlineData(NarClim2Frequency.Hour3, "3hr")]
    [InlineData(NarClim2Frequency.Day, "day")]
    [InlineData(NarClim2Frequency.Month, "mon")]
    public void FrequencyNames_ToString_ReturnsCorrectString(NarClim2Frequency input, string expected)
    {
        var result = NarClim2Constants.FrequencyNames.ToString(input);
        Assert.Equal(expected, result);
    }
}
