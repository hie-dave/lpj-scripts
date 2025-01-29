using ClimateProcessing.Extensions;
using Xunit;

namespace ClimateProcessing.Tests.Extensions;

public class StringExtensionsTests
{
    [Fact]
    public void ReplaceFirst_WhenTextIsNull_ThrowsArgumentNullException()
    {
        string? text = null;
        Assert.Throws<ArgumentNullException>(() => text!.ReplaceFirst("search", "replace"));
    }

    [Fact]
    public void ReplaceFirst_WhenSearchIsNull_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() => "text".ReplaceFirst(null!, "replace"));
    }

    [Fact]
    public void ReplaceFirst_WhenReplaceIsNull_ReplacesWithEmptyString()
    {
        string result = "hello world".ReplaceFirst("world", string.Empty);
        Assert.Equal("hello ", result);
    }

    [Fact]
    public void ReplaceFirst_WhenSearchNotFound_ReturnsOriginalString()
    {
        string original = "hello world";
        string result = original.ReplaceFirst("xyz", "replace");
        Assert.Equal(original, result);
    }

    [Fact]
    public void ReplaceFirst_WhenSearchIsEmpty_ReturnsOriginalString()
    {
        string original = "hello world";
        string result = original.ReplaceFirst("", "replace");
        Assert.Equal(original, result);
    }

    [Theory]
    [InlineData("hello world", "hello", "hi", "hi world")]
    [InlineData("hello hello world", "hello", "hi", "hi hello world")]
    [InlineData("temperature_data.nc", "temperature_", "vpd_", "vpd_data.nc")]
    [InlineData("prefix_prefix_suffix", "prefix_", "new_", "new_prefix_suffix")]
    public void ReplaceFirst_ReplacesOnlyFirstOccurrence(string input, string search, string replace, string expected)
    {
        string result = input.ReplaceFirst(search, replace);
        Assert.Equal(expected, result);
    }
}
