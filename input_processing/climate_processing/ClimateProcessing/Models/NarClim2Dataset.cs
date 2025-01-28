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
    private readonly NarClim2Frequency _frequency;

    private static readonly Dictionary<ClimateVariable, (string name, string units)> _variableMap = new()
    {
        { ClimateVariable.SpecificHumidity, ("huss", "1") },        // Specific humidity
        { ClimateVariable.SurfacePressure, ("ps", "Pa") },         // Surface pressure
        { ClimateVariable.ShortwaveRadiation, ("rsds", "W m-2") }, // Surface downwelling shortwave radiation
        { ClimateVariable.WindSpeed, ("sfcwind", "m s-1") },       // Near-surface wind speed
        { ClimateVariable.Temperature, ("tas", "K") }              // Near-surface air temperature
    };

    private readonly Dictionary<string, string> _metadata;

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
        _metadata = GetMetadata();
    }

    public string DatasetName => 
        $"NARCliM2.0_{NarClim2Constants.GCMNames.ToString(_gcm)}_{NarClim2Constants.ExperimentNames.ToString(_experiment)}_{NarClim2Constants.RCMNames.ToString(_rcm)}";

    public IEnumerable<string> GetInputFiles()
    {
        var baseDir = Path.Combine(_inputPath,
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

        var files = new List<string>();
        foreach (var variable in _variableMap.Values)
        {
            var varDir = Path.Combine(baseDir, variable.name, NarClim2Constants.Paths.LatestVersion);
            if (Directory.Exists(varDir))
            {
                files.AddRange(Directory.GetFiles(varDir, "*.nc")
                    .OrderBy(f => GetDateFromFilename(f)));
            }
        }

        return files;
    }

    public VariableInfo GetVariableInfo(ClimateVariable variable)
    {
        if (!_variableMap.TryGetValue(variable, out var info))
        {
            throw new ArgumentException($"Variable {variable} not supported in NARCliM2 dataset");
        }
        return new VariableInfo(info.name, info.units);
    }

    public Dictionary<string, string> GetMetadata() => new()
    {
        { "source", "NARCliM2.0" },
        { "domain", NarClim2Constants.DomainNames.ToString(_domain) },
        { "gcm", NarClim2Constants.GCMNames.ToString(_gcm) },
        { "experiment", NarClim2Constants.ExperimentNames.ToString(_experiment) },
        { "rcm", NarClim2Constants.RCMNames.ToString(_rcm) },
        { "frequency", NarClim2Constants.FrequencyNames.ToString(_frequency) },
        { "version", NarClim2Constants.Paths.Version }
    };

    public string GetOutputFilePattern(ClimateVariable variable)
    {
        // Get the first input file for this variable to use as a pattern
        var inputFiles = Directory.GetFiles(_inputPath,
            "*.nc")
            .Where(f => Path.GetFileName(f).Contains(_variableMap[variable].name))
            .OrderBy(f => f);

        if (!inputFiles.Any())
        {
            throw new InvalidOperationException($"No input files found for variable {variable}");
        }

        var firstFile = Path.GetFileName(inputFiles.First());
        
        // Extract the pattern before the date range
        var prefix = string.Join("_", firstFile.Split('_').TakeWhile(p => !p.Contains("-")));
        
        // Add placeholder for full date range and extension
        return $"{prefix}_XXXXXXXX-XXXXXXXX.nc";
    }

    private static DateTime GetDateFromFilename(string filename)
    {
        // Example filename: tas_AUS18_ACCESS-ESM1-5_historical_r6i1p1f1_NSW-Government_NARCliM2-0-WRF412R3_v1-r1_mon_195101-195112.nc
        var match = Regex.Match(filename, @"_(\d{6})-\d{6}\.nc$");
        if (!match.Success)
            throw new ArgumentException($"Invalid filename format: {filename}");

        var dateStr = match.Groups[1].Value;
        return DateTime.ParseExact(dateStr, "yyyyMM", null);
    }
}
