namespace ClimateProcessing.Models;

/// <summary>
/// Domain options for NARCliM2.0 dataset
/// </summary>
public enum NarClim2Domain
{
    /// <summary>
    /// 20km resolution domain covering Australasia
    /// </summary>
    AUS18,

    /// <summary>
    /// 4km resolution domain covering south-eastern Australia
    /// </summary>
    SEAus04
}

/// <summary>
/// Global Climate Models (GCMs) available in NARCliM2.0
/// </summary>
public enum NarClim2GCM
{
    AccessEsm15,
    EcEarth3Veg,
    MpiEsm12Hr,
    NorEsm2Mm,
    Ukesm10Ll
}

/// <summary>
/// Experiments available in NARCliM2.0
/// </summary>
public enum NarClim2Experiment
{
    Historical,
    SSP126,
    SSP245,
    SSP370
}

/// <summary>
/// Regional Climate Models (RCMs) available in NARCliM2.0
/// </summary>
public enum NarClim2RCM
{
    /// <summary>
    /// WRF v4.1.2 with MYNN2 PBL physics
    /// </summary>
    WRF412R3,

    /// <summary>
    /// WRF v4.1.2 with ACM2 PBL physics
    /// </summary>
    WRF412R5
}

/// <summary>
/// Time frequency options for NARCliM2.0 data
/// </summary>
public enum NarClim2Frequency
{
    Hour1,
    Hour3,
    Day,
    Month
}

/// <summary>
/// Constants and string mappings for NARCliM2.0 dataset
/// </summary>
public static class NarClim2Constants
{
    public static class Paths
    {
        public const string MipEra = "CMIP6";
        public const string ActivityId = "DD";
        public const string Institution = "NSW-Government";
        public const string Version = "v1-r1";
        public const string LatestVersion = "latest";
    }

    public static class DomainNames
    {
        public const string AUS18 = "AUS18";
        public const string SEAus04 = "NARCliM2-0-SEAus-04";

        public static string ToString(NarClim2Domain domain) => domain switch
        {
            NarClim2Domain.AUS18 => AUS18,
            NarClim2Domain.SEAus04 => SEAus04,
            _ => throw new ArgumentException($"Unknown domain: {domain}")
        };
    }

    public static class GCMNames
    {
        public const string AccessEsm15 = "ACCESS-ESM1-5";
        public const string EcEarth3Veg = "EC-Earth3-Veg";
        public const string MpiEsm12Hr = "MPI-ESM1-2-HR";
        public const string NorEsm2Mm = "NorESM2-MM";
        public const string Ukesm10Ll = "UKESM1-0-LL";

        public static string ToString(NarClim2GCM gcm) => gcm switch
        {
            NarClim2GCM.AccessEsm15 => AccessEsm15,
            NarClim2GCM.EcEarth3Veg => EcEarth3Veg,
            NarClim2GCM.MpiEsm12Hr => MpiEsm12Hr,
            NarClim2GCM.NorEsm2Mm => NorEsm2Mm,
            NarClim2GCM.Ukesm10Ll => Ukesm10Ll,
            _ => throw new ArgumentException($"Unknown GCM: {gcm}")
        };
    }

    public static class ExperimentNames
    {
        public const string Historical = "historical";
        public const string SSP126 = "ssp126";
        public const string SSP245 = "ssp245";
        public const string SSP370 = "ssp370";

        public static string ToString(NarClim2Experiment experiment) => experiment switch
        {
            NarClim2Experiment.Historical => Historical,
            NarClim2Experiment.SSP126 => SSP126,
            NarClim2Experiment.SSP245 => SSP245,
            NarClim2Experiment.SSP370 => SSP370,
            _ => throw new ArgumentException($"Unknown experiment: {experiment}")
        };
    }

    public static class RCMNames
    {
        public const string WRF412R3 = "NARCliM2-0-WRF412R3";
        public const string WRF412R5 = "NARCliM2-0-WRF412R5";

        public static string ToString(NarClim2RCM rcm) => rcm switch
        {
            NarClim2RCM.WRF412R3 => WRF412R3,
            NarClim2RCM.WRF412R5 => WRF412R5,
            _ => throw new ArgumentException($"Unknown RCM: {rcm}")
        };
    }

    public static class FrequencyNames
    {
        public const string Hour1 = "1hr";
        public const string Hour3 = "3hr";
        public const string Day = "day";
        public const string Month = "mon";

        public static string ToString(NarClim2Frequency frequency) => frequency switch
        {
            NarClim2Frequency.Hour1 => Hour1,
            NarClim2Frequency.Hour3 => Hour3,
            NarClim2Frequency.Day => Day,
            NarClim2Frequency.Month => Month,
            _ => throw new ArgumentException($"Unknown frequency: {frequency}")
        };
    }

    public static class VariantLabels
    {
        public const string AccessEsm15 = "r6i1p1f1";
        public const string EcEarth3Veg = "r1i1p1f1";
        public const string MpiEsm12Hr = "r1i1p1f1";
        public const string NorEsm2Mm = "r1i1p1f1";
        public const string Ukesm10Ll = "r1i1p1f2";

        public static string GetVariantLabel(NarClim2GCM gcm) => gcm switch
        {
            NarClim2GCM.AccessEsm15 => AccessEsm15,
            NarClim2GCM.EcEarth3Veg => EcEarth3Veg,
            NarClim2GCM.MpiEsm12Hr => MpiEsm12Hr,
            NarClim2GCM.NorEsm2Mm => NorEsm2Mm,
            NarClim2GCM.Ukesm10Ll => Ukesm10Ll,
            _ => throw new ArgumentException($"Unknown GCM: {gcm}")
        };
    }
}
