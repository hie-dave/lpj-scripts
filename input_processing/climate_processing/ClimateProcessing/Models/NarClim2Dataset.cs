using System.Text.RegularExpressions;
using System.Collections.Generic;
using System.IO;
using System.Linq;

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
        { ClimateVariable.Precipitation, ("pr", "kg m-2 s-1") }
    };

    /// <summary>
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

    public static NarClim2Dataset Create(ProcessingConfig config)
    {
        return new NarClim2Dataset(
            inputPath: config.InputDirectory,
            frequency: NarClim2Constants.ParseFrequency(config.InputTimeStep.Hours));
    }

    public string DatasetName =>
        $"NARCliM2.0_{NarClim2Constants.GCMNames.ToString(_gcm)}_{NarClim2Constants.ExperimentNames.ToString(_experiment)}_{NarClim2Constants.RCMNames.ToString(_rcm)}";

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
        string dir = GetInputFilesDirectory(variable);
        if (!Directory.Exists(dir))
            return Enumerable.Empty<string>();

        return Directory.GetFiles(dir, "*.nc").OrderBy(GetDateFromFilename);
    }

    public VariableInfo GetVariableInfo(ClimateVariable variable)
    {
        if (!_variableMap.TryGetValue(variable, out var info))
            throw new ArgumentException($"Variable {variable} not supported in NARCliM2 dataset");
        return new VariableInfo(info.name, info.units);
    }

    private static DateTime GetDateFromFilename(string filename)
    {
        // Example filename: tas_AUS18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_mon_195101-195112.nc
        Match match = Regex.Match(filename, @"_(\d{6})-\d{6}\.nc$");
        if (!match.Success)
            throw new ArgumentException($"Invalid filename format: {filename}");

        string dateStr = match.Groups[1].Value;
        return DateTime.ParseExact(dateStr, "yyyyMM", null);
    }

    public string GenerateOutputFileName(ClimateVariable variable)
    {
        // Get the first input file for this variable to use as a pattern
        IEnumerable<string> inputFiles = GetInputFiles(variable);

        if (!inputFiles.Any())
            throw new InvalidOperationException($"No input files found for variable {GetVariableInfo(variable).Name}");

        DateTime startDate = GetDateFromFilename(inputFiles.First());
        DateTime endDate = GetDateFromFilename(inputFiles.Last());

        string firstFile = Path.GetFileName(inputFiles.First());

        // Extract the pattern before the date range
        // pr_AUS-18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_mon_198301-198312.nc
        string prefix = string.Join("_", firstFile.Split('_').TakeWhile(p => !p.Contains(".nc")));

        // Add the date range and extension
        return $"{prefix}_{startDate:yyyyMM}-{endDate:yyyyMM}.nc";
    }

    public string GetOutputDirectory()
    {
        string gcm = NarClim2Constants.GCMNames.ToString(_gcm);
        string experiment = NarClim2Constants.ExperimentNames.ToString(_experiment);
        string rcm = NarClim2Constants.RCMNames.ToString(_rcm);
        return Path.Combine(gcm, experiment, rcm);
    }
}
