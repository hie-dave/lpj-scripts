using ClimateProcessing.Models;
using Xunit;

namespace ClimateProcessing.Tests.Models;

public class PBSWalltimeTests
{
    [Theory]
    [InlineData(0, 0, 0)]
    [InlineData(1, 0, 0)]
    [InlineData(0, 1, 0)]
    [InlineData(0, 0, 1)]
    [InlineData(48, 0, 0)]
    [InlineData(1, 59, 59)]
    public void Constructor_ValidValues_CreatesInstance(int hours, int minutes, int seconds)
    {
        var walltime = new PBSWalltime(hours, minutes, seconds);
        Assert.Equal(hours, walltime.Hours);
        Assert.Equal(minutes, walltime.Minutes);
        Assert.Equal(seconds, walltime.Seconds);
    }

    [Theory]
    [InlineData(-1, 0, 0, "Invalid time components")]
    [InlineData(0, -1, 0, "Invalid time components")]
    [InlineData(0, 0, -1, "Invalid time components")]
    [InlineData(0, 60, 0, "Invalid time components")]
    [InlineData(0, 0, 60, "Invalid time components")]
    [InlineData(49, 0, 0, "PBS walltime cannot exceed 48 hours")]
    [InlineData(48, 1, 0, "PBS walltime cannot exceed 48 hours")]
    [InlineData(48, 0, 1, "PBS walltime cannot exceed 48 hours")]
    public void Constructor_InvalidValues_ThrowsArgumentException(int hours, int minutes, int seconds, string expectedMessage)
    {
        var ex = Assert.Throws<ArgumentException>(() => new PBSWalltime(hours, minutes, seconds));
        Assert.Equal(expectedMessage, ex.Message);
    }

    [Theory]
    [InlineData("00:00:00", 0, 0, 0)]
    [InlineData("01:00:00", 1, 0, 0)]
    [InlineData("00:01:00", 0, 1, 0)]
    [InlineData("00:00:01", 0, 0, 1)]
    [InlineData("03:02:01", 3, 2, 1)]
    [InlineData("48:00:00", 48, 0, 0)]
    [InlineData("01:59:59", 1, 59, 59)]
    public void Parse_ValidStrings_ReturnsCorrectInstance(string input, int expectedHours, int expectedMinutes, int expectedSeconds)
    {
        var walltime = PBSWalltime.Parse(input);
        Assert.Equal(expectedHours, walltime.Hours);
        Assert.Equal(expectedMinutes, walltime.Minutes);
        Assert.Equal(expectedSeconds, walltime.Seconds);
    }

    [Theory]
    [InlineData("")]
    [InlineData("1:00:00")]
    [InlineData("01:0:00")]
    [InlineData("01:00:0")]
    [InlineData("1:0:0")]
    [InlineData("abc")]
    [InlineData("01:00")]
    [InlineData("01:00:00:00")]
    public void Parse_InvalidFormat_ThrowsFormatException(string input)
    {
        Assert.Throws<FormatException>(() => PBSWalltime.Parse(input));
    }

    [Theory]
    [InlineData("49:00:00")]
    [InlineData("48:01:00")]
    [InlineData("48:00:01")]
    [InlineData("-1:00:00")]
    [InlineData("00:-1:00")]
    [InlineData("00:60:00")]
    [InlineData("00:00:60")]
    public void Parse_InvalidValues_ThrowsArgumentException(string input)
    {
        Assert.Throws<ArgumentException>(() => PBSWalltime.Parse(input));
    }

    [Theory]
    [InlineData(0, 0, 0, "00:00:00")]
    [InlineData(1, 0, 0, "01:00:00")]
    [InlineData(0, 1, 0, "00:01:00")]
    [InlineData(0, 0, 1, "00:00:01")]
    [InlineData(48, 0, 0, "48:00:00")]
    [InlineData(1, 59, 59, "01:59:59")]
    public void ToString_ReturnsCorrectFormat(int hours, int minutes, int seconds, string expected)
    {
        var walltime = new PBSWalltime(hours, minutes, seconds);
        Assert.Equal(expected, walltime.ToString());
    }

    [Theory]
    [InlineData(0, 0, 0)]
    [InlineData(1, 0, 0)]
    [InlineData(0, 1, 0)]
    [InlineData(0, 0, 1)]
    [InlineData(48, 0, 0)]
    [InlineData(1, 59, 59)]
    public void TimeSpanConversion_Roundtrip_PreservesValues(int hours, int minutes, int seconds)
    {
        var original = new PBSWalltime(hours, minutes, seconds);
        TimeSpan timespan = original;
        PBSWalltime converted = (PBSWalltime)timespan;

        Assert.Equal(original, converted);
    }

    [Theory]
    [InlineData(49, 0, 0)]
    [InlineData(48, 1, 0)]
    [InlineData(48, 0, 1)]
    public void TimeSpanConversion_InvalidValues_ThrowsArgumentException(int hours, int minutes, int seconds)
    {
        var timespan = new TimeSpan(hours, minutes, seconds);
        Assert.Throws<ArgumentException>(() => (PBSWalltime)timespan);
    }

    [Fact]
    public void MaxValue_IsFortyEightHours()
    {
        Assert.Equal(48, PBSWalltime.MaxValue.Hours);
        Assert.Equal(0, PBSWalltime.MaxValue.Minutes);
        Assert.Equal(0, PBSWalltime.MaxValue.Seconds);
    }

    [Fact]
    public void Equals_SameValues_ReturnsTrue()
    {
        var walltime1 = new PBSWalltime(1, 2, 3);
        var walltime2 = new PBSWalltime(1, 2, 3);

        Assert.True(walltime1.Equals(walltime2));
        Assert.True(walltime1.Equals((object)walltime2));
        Assert.True(walltime1 == walltime2);
        Assert.False(walltime1 != walltime2);
    }

    [Fact]
    public void Equals_DifferentValues_ReturnsFalse()
    {
        var walltime1 = new PBSWalltime(1, 2, 3);
        var walltime2 = new PBSWalltime(1, 2, 4);

        Assert.False(walltime1.Equals(walltime2));
        Assert.False(walltime1.Equals((object)walltime2));
        Assert.False(walltime1 == walltime2);
        Assert.True(walltime1 != walltime2);
    }

    [Fact]
    public void Equals_Null_ReturnsFalse()
    {
        var walltime = new PBSWalltime(1, 2, 3);
        Assert.False(walltime.Equals(null));
    }

    [Fact]
    public void GetHashCode_SameValues_ReturnsSameHashCode()
    {
        var walltime1 = new PBSWalltime(1, 2, 3);
        var walltime2 = new PBSWalltime(1, 2, 3);

        Assert.Equal(walltime1.GetHashCode(), walltime2.GetHashCode());
    }
}
