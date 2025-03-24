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

    [Theory]
    [InlineData(NarClim2Domain.AUS18, NarClim2Constants.Files.RlonValuesFileAUS18)]
    [InlineData(NarClim2Domain.SEAus04, NarClim2Constants.Files.RlonValuesFileSEAus04)]
    public void GetRlonValueFile_ReturnsCorrectFile(NarClim2Domain domain, string expected)
    {
        string actual = NarClim2Constants.Files.GetRlonValuesFile(domain);
        Assert.Equal(expected, actual);
    }

    [Theory]
    [InlineData(NarClim2Domain.AUS18)]
    [InlineData(NarClim2Domain.SEAus04)]
    public void DomainNames_TestRoundTrip(NarClim2Domain domain)
    {
        string str = NarClim2Constants.DomainNames.ToString(domain);
        Assert.Equal(domain, NarClim2Constants.DomainNames.FromString(str));
    }

    [Theory]
    [InlineData(NarClim2GCM.AccessEsm15)]
    [InlineData(NarClim2GCM.EcEarth3Veg)]
    [InlineData(NarClim2GCM.MpiEsm12Hr)]
    [InlineData(NarClim2GCM.NorEsm2Mm)]
    [InlineData(NarClim2GCM.Ukesm10Ll)]
    public void GCMNames_TestRoundTrip(NarClim2GCM gcm)
    {
        string str = NarClim2Constants.GCMNames.ToString(gcm);
        Assert.Equal(gcm, NarClim2Constants.GCMNames.FromString(str));
    }

    [Theory]
    [InlineData(NarClim2Experiment.Historical)]
    [InlineData(NarClim2Experiment.SSP126)]
    [InlineData(NarClim2Experiment.SSP370)]
    public void ExperimentNames_TestRoundTrip(NarClim2Experiment experiment)
    {
        string str = NarClim2Constants.ExperimentNames.ToString(experiment);
        Assert.Equal(experiment, NarClim2Constants.ExperimentNames.FromString(str));
    }

    [Theory]
    [InlineData(NarClim2RCM.WRF412R3)]
    [InlineData(NarClim2RCM.WRF412R5)]
    public void RCMNames_TestRoundTrip(NarClim2RCM rcm)
    {
        string str = NarClim2Constants.RCMNames.ToString(rcm);
        Assert.Equal(rcm, NarClim2Constants.RCMNames.FromString(str));
    }

    [Theory]
    [InlineData(NarClim2Frequency.Hour1)]
    [InlineData(NarClim2Frequency.Hour3)]
    [InlineData(NarClim2Frequency.Day)]
    [InlineData(NarClim2Frequency.Month)]
    public void FrequencyNames_TestRoundTrip(NarClim2Frequency frequency)
    {
        string str = NarClim2Constants.FrequencyNames.ToString(frequency);
        Assert.Equal(frequency, NarClim2Constants.FrequencyNames.FromString(str));
    }

    [Theory]
    [InlineData(1, NarClim2Frequency.Hour1)]
    [InlineData(3, NarClim2Frequency.Hour3)]
    [InlineData(24, NarClim2Frequency.Day)]
    public void ParseFrequency_ValidInput_ReturnsExpected(int input, NarClim2Frequency expected)
    {
        NarClim2Frequency result = NarClim2Constants.ParseFrequency(input);
        Assert.Equal(expected, result);
    }

    [Theory]
    [InlineData(6)]
    [InlineData(12)]
    [InlineData(720)]
    public void ParseFrequency_InvalidInput_ThrowsException(int frequency)
    {
        Assert.ThrowsAny<ArgumentException>(() => NarClim2Constants.ParseFrequency(frequency));
    }
}
