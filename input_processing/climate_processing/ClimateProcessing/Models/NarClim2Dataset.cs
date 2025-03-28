using System.Text.RegularExpressions;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Microsoft.VisualBasic;
using System.Globalization;

namespace ClimateProcessing.Models;

public class NarClim2Dataset : IClimateDataset
{
    private readonly string _inputPath;
    private readonly NarClim2Domain _domain;
    private readonly NarClim2GCM _gcm;
    private readonly NarClim2Experiment _experiment;
    private readonly NarClim2RCM _rcm;

    /// <summary>
    /// Timestep of input data.
    /// </summary>
    private readonly NarClim2Frequency _frequency;

    /// <summary>
    /// The domain of the data.
    /// </summary>
    public NarClim2Domain Domain => _domain;

    /// <summary>
    /// The GCM of the data.
    /// </summary>
    public NarClim2GCM GCM => _gcm;

    /// <summary>
    /// The experiment of the data.
    /// </summary>
    public NarClim2Experiment Experiment => _experiment;

    /// <summary>
    /// The RCM of the data.
    /// </summary>
    public NarClim2RCM RCM => _rcm;

    /// <summary>
    /// The base path for all files in the dataset.
    /// </summary>
    public string BasePath { get => _inputPath; }

    // Variable names and units as they exist in the NARCliM2 dataset.
    private static readonly Dictionary<ClimateVariable, (string name, string units)> _variableMap = new()
    {
        // Specific humidity
        { ClimateVariable.SpecificHumidity, ("huss", "1") },

        // Surface pressure
        { ClimateVariable.SurfacePressure, ("ps", "Pa") },

        // Surface downwelling shortwave radiation
        { ClimateVariable.ShortwaveRadiation, ("rsds", "W m-2") },

        // Near-surface wind speed
        { ClimateVariable.WindSpeed, ("sfcWind", "m s-1") },

        // Near-surface air temperature
        { ClimateVariable.Temperature, ("tas", "K") },

        // Precipitation
        { ClimateVariable.Precipitation, ("pr", "kg m-2 s-1") },

        // Daily maximum near-surface air temperature
        { ClimateVariable.MaxTemperature, ("tasmax", "K") },

        // Daily minimum near-surface air temperature
        { ClimateVariable.MinTemperature, ("tasmin", "K") },
    };

    /// <summary>
    /// Create a new instance of the NARCliM2 dataset.
    /// </summary>
    /// Create a new instance of the NARCliM2 dataset.
    /// </summary>
    /// <param name="inputPath">Path to the NARCliM2 dataset.</param>
    /// <param name="domain">Domain of the dataset.</param>
    /// <param name="gcm">Global Climate Model.</param>
    /// <param name="experiment">The scenario.</param>
    /// <param name="rcm">The Regional Climate Model.</param>
    /// <param name="frequency">Timestep of input data.</param>
    public NarClim2Dataset(
        string inputPath,
        NarClim2Domain domain = NarClim2Domain.AUS18,
        NarClim2GCM gcm = NarClim2GCM.AccessEsm15,
        NarClim2Experiment experiment = NarClim2Experiment.Historical,
        NarClim2RCM rcm = NarClim2RCM.WRF412R3,
        NarClim2Frequency frequency = NarClim2Frequency.Month)
    {
        _inputPath = inputPath;
        _domain = domain;
        _gcm = gcm;
        _experiment = experiment;
        _rcm = rcm;
        _frequency = frequency;
    }

    public string DatasetName =>
        $"NARCliM2.0_{NarClim2Constants.DomainNames.ToString(_domain)}_{NarClim2Constants.GCMNames.ToString(_gcm)}_{NarClim2Constants.ExperimentNames.ToString(_experiment)}_{NarClim2Constants.RCMNames.ToString(_rcm)}";

    public string GetInputFilesDirectory(ClimateVariable variable)
    {
        string baseDir = Path.Combine(_inputPath,
            NarClim2Constants.Paths.MipEra,
            NarClim2Constants.Paths.ActivityId,
            NarClim2Constants.DomainNames.ToString(_domain),
            NarClim2Constants.Paths.Institution,
            NarClim2Constants.GCMNames.ToString(_gcm),
            NarClim2Constants.ExperimentNames.ToString(_experiment),
            NarClim2Constants.VariantLabels.GetVariantLabel(_gcm),
            NarClim2Constants.RCMNames.ToString(_rcm),
            NarClim2Constants.Paths.Version,
            NarClim2Constants.FrequencyNames.ToString(_frequency));

        string varName = GetVariableInfo(variable).Name;
        return Path.Combine(baseDir, varName, NarClim2Constants.Paths.LatestVersion);
    }

    public IEnumerable<string> GetInputFiles(ClimateVariable variable)
    {
        if ( (variable == ClimateVariable.MaxTemperature || variable == ClimateVariable.MinTemperature) && _frequency != NarClim2Frequency.Day)
            throw new ArgumentException($"Variable {variable} not supported with frequency {_frequency}");

        string dir = GetInputFilesDirectory(variable);
        if (!Directory.Exists(dir))
            return Enumerable.Empty<string>();

        return Directory.GetFiles(dir, "*.nc").OrderBy(f => GetDateFromFilename(f, true));
    }

    public VariableInfo GetVariableInfo(ClimateVariable variable)
    {
        if (!_variableMap.TryGetValue(variable, out var info))
            throw new ArgumentException($"Variable {variable} not supported in NARCliM2 dataset");
        return new VariableInfo(info.name, info.units);
    }

    /// <summary>
    /// Get the regular expression and managed DateTime format string which
    /// can be used to parse NARCliM2 filenames for files on the specified
    /// frequency.
    /// </summary>
    /// <param name="frequency">The frequency of the filename.</param>
    /// <returns>Tuple of regular expression and DateTime format string.</returns>
    /// <exception cref="ArgumentException">Thrown when an invalid frequency is specified. Should never happen in practice unless narclim add more timestep frequencies in the future (e.g. 6hr).</exception>
    internal static (string regex, string format) GetDateRegexFormat(NarClim2Frequency frequency)
    {
        const string regexMonth = @"_(\d{6})-(\d{6})\.nc$";
        const string regexDay = @"_(\d{8})-(\d{8})\.nc$";
        const string regexHour = @"_(\d{10})-(\d{10})\.nc$";

        const string formatHour = "yyyyMMddHH";
        const string formatDay = "yyyyMMdd";
        const string formatMonth = "yyyyMM";

        return frequency switch
        {
            NarClim2Frequency.Hour1 => (regexHour, formatHour),
            NarClim2Frequency.Hour3 => (regexHour, formatHour),
            NarClim2Frequency.Day => (regexDay, formatDay),
            NarClim2Frequency.Month => (regexMonth, formatMonth),
            _ => throw new ArgumentException($"Unknown frequency: {frequency}")
        };
    }

    /// <summary>
    /// Gets the start or end date from a NARCliM2 filename.
    /// </summary>
    /// <param name="filename">The filename to parse.</param>
    /// <param name="start">If true, get the start date. If false, get the end date.</param>
    /// <returns>The date.</returns>
    /// <exception cref="ArgumentException">If the date cannot be determined from the filename.</exception>
    internal static DateTime GetDateFromFilename(string filename, bool start)
    {
        // Example filename: tas_AUS18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_mon_195101-195112.nc
        // Remove the directory component, if one is present.
        filename = Path.GetFileName(filename);

        // Determine if the file is subdaily and use the correct regex.
        // Note: 1hr and 3hr files contain yyyyMMddHH, but day, mon, ..., etc
        //       use yyyyMMdd format (ie not time component).
        NarClim2Frequency frequency = GetFrequencyFromFilename(filename);
        (string regex, string fmt) = GetDateRegexFormat(frequency);

        // Parse the filename.
        Match match = Regex.Match(filename, regex);
        if (!match.Success)
            throw new ArgumentException($"Unable to get date from filename. Invalid filename format: {filename}");

        // Choose the correct group and format, depending on whether we are
        // looking for the start or end date.
        string dateStr = match.Groups[start ? 1 : 2].Value;

        // Parse the date.
        return DateTime.ParseExact(dateStr, fmt, CultureInfo.InvariantCulture);
    }

    /// <summary>
    /// Parse the frequency from a NARCliM2 filename.
    /// </summary>
    /// <param name="filename">The filename to parse.</param>
    /// <returns>The frequency enum value.</returns>
    /// <exception cref="ArgumentException">If the frequency cannot be determined from the filename.</exception>
    internal static NarClim2Frequency GetFrequencyFromFilename(string filename)
    {
        // Example: tas_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_mon_195101-195112.nc
        //                                                                                              ^^^
        // Remove the directory component, if one is present
        filename = Path.GetFileName(filename);

        // Build regex pattern to match any of the frequency strings followed by underscore and digits
        var freqStrings = new[]
        {
            NarClim2Constants.FrequencyNames.Hour1,
            NarClim2Constants.FrequencyNames.Hour3,
            NarClim2Constants.FrequencyNames.Day,
            NarClim2Constants.FrequencyNames.Month
        };
        string pattern = $@"_({string.Join("|", freqStrings)})_\d+";

        Match match = Regex.Match(filename, pattern);
        if (!match.Success)
            throw new ArgumentException($"Unable to determine frequency from filename. Invalid filename format: {filename}");

        string freqStr = match.Groups[1].Value;
        return NarClim2Constants.FrequencyNames.FromString(freqStr);
    }

    /// <inheritdoc />
    /// <remarks>
    /// Marked as virtual for testing.
    /// </remarks>
    public virtual string GenerateOutputFileName(ClimateVariable variable)
    {
        // Get the first input file for this variable to use as a pattern
        IEnumerable<string> inputFiles = GetInputFiles(variable);

        if (!inputFiles.Any())
            throw new InvalidOperationException($"No input files found for variable {GetVariableInfo(variable).Name} in domain {_domain}, GCM {_gcm}, experiment {_experiment}, RCM {_rcm}.");

        DateTime startDate = GetDateFromFilename(inputFiles.First(), true);
        DateTime endDate = GetDateFromFilename(inputFiles.Last(), false);

        string firstFile = Path.GetFileName(inputFiles.First());

        // Extract the pattern before the date range
        // pr_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_mon_198301-198312.nc
        string prefix = string.Join("_", firstFile.Split('_').TakeWhile(p => !p.Contains(".nc")));

        // Add the date range and extension
        return $"{prefix}_{startDate:yyyyMM}-{endDate:yyyyMM}.nc";
    }

    /// <inheritdoc /> 
    public string GetOutputDirectory()
    {
        string gcm = NarClim2Constants.GCMNames.ToString(_gcm);
        string experiment = NarClim2Constants.ExperimentNames.ToString(_experiment);
        string rcm = NarClim2Constants.RCMNames.ToString(_rcm);
        return Path.Combine(gcm, experiment, rcm);
    }
}
