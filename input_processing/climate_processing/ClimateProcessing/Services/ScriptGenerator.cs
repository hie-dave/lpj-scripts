using System.Text;
using ClimateProcessing.Models;
using System.IO;
using ClimateProcessing.Units;
using System.Runtime.CompilerServices;

[assembly: InternalsVisibleTo("ClimateProcessing.Tests")]

namespace ClimateProcessing.Services;

public class ScriptGenerator
{
    private readonly ProcessingConfig _config;

    // List of standard variables and their output names and units.
    private static readonly Dictionary<ClimateVariable, (string outName, string outUnits)> _standardVariables = new()
    {
        { ClimateVariable.SpecificHumidity, ("huss", "1") },
        { ClimateVariable.SurfacePressure, ("ps", "Pa") },
        { ClimateVariable.ShortwaveRadiation, ("rsds", "W m-2") },
        { ClimateVariable.WindSpeed, ("sfcWind", "m s-1") },
        { ClimateVariable.Temperature, ("tas", "degC") },
        { ClimateVariable.Precipitation, ("pr", "mm") }
    };

    private (string outName, string outUnits) GetStandardConfig(ClimateVariable variable)
    {
        if (!_standardVariables.TryGetValue(variable, out var config))
            throw new ArgumentException($"No configuration found for variable {variable}");
        return config;
    }

    public ScriptGenerator(ProcessingConfig config)
    {
        _config = config;
    }

    internal string GenerateRenameOperator(string inName, string outName)
    {
        if (inName == outName)
            return string.Empty;
        return $"-chname,'{inName}','{outName}'";
    }

    internal IEnumerable<string> GenerateUnitConversionOperators(
        string outputVar,
        string inputUnits,
        string targetUnits,
        TimeStep timeStep)
    {
        var result = UnitConverter.AnalyzeConversion(inputUnits, targetUnits);

        List<string> operators = [];

        if (result.RequiresConversion)
        {
            string expression = UnitConverter.GenerateConversionExpression(
                inputUnits,
                targetUnits,
                timeStep);
            operators.Add($"-expr,'{expression}'");
        }

        if (result.RequiresRenaming)
            operators.Add($"-setattribute,'{outputVar}@units={targetUnits}'");

        return operators;
    }

    internal string GenerateTimeAggregationOperator(
        ClimateVariable variable)
    {
        // Only aggregate if input and output timesteps differ
        if (_config.InputTimeStep == _config.OutputTimeStep)
            return string.Empty;

        // Calculate the number of timesteps to aggregate
        int stepsToAggregate = _config.OutputTimeStep.Hours / _config.InputTimeStep.Hours;

        var aggregationMethod = _config.GetAggregationMethod(variable);
        var @operator = aggregationMethod.ToCdoOperator(_config.OutputTimeStep);

        return $"-{@operator},{stepsToAggregate}";
    }

    internal ProcessingCommand GenerateVariableRenameCommand(
        string inputVar,
        string outputVar,
        string inputFile,
        string outputFile)
    {
        if (inputVar == outputVar)
        {
            return ProcessingCommand.Skip(inputFile, outputFile);
        }

        return ProcessingCommand.Process(
            $"cdo -O chname,{inputVar},{outputVar} {inputFile} {outputFile}");
    }

    private string GetVPDEquations(VPDMethod method)
    {
        // All methods follow the same general pattern:
        // 1. Calculate saturation vapor pressure (_esat)
        // 2. Calculate actual vapor pressure (_e)
        // 3. Calculate VPD as (_esat - _e) / 1000 to convert to kPa

        var esatEquation = method switch
        {
            // Magnus equation (default)
            VPDMethod.Magnus => "_esat=0.611*exp((17.27*tas)/(tas+237.3))*1000",

            // Buck (1981)
            // Buck's equation for temperatures above 0°C
            VPDMethod.Buck1981 => "_esat=0.61121*exp((18.678-tas/234.5)*(tas/(257.14+tas)))*1000",

            // Alduchov and Eskridge (1996)
            // More accurate coefficients for the Magnus equation
            VPDMethod.AlduchovEskridge1996 => "_esat=0.61094*exp((17.625*tas)/(tas+243.04))*1000",

            // Allen et al. (1998) FAO
            // Tetens equation with FAO coefficients
            VPDMethod.AllenFAO1998 => "_esat=0.6108*exp((17.27*tas)/(tas+237.3))*1000",

            // Sonntag (1990)
            // Based on ITS-90 temperature scale
            VPDMethod.Sonntag1990 => "_esat=0.61078*exp((17.08085*tas)/(234.175+tas))*1000",

            _ => throw new ArgumentException($"Unsupported VPD calculation method: {method}")
        };

        return $@"# Saturation vapor pressure (Pa)
{esatEquation}
# Actual vapor pressure (Pa)
_e=(huss*ps)/(0.622+0.378*huss)
# VPD (kPa)
vpd=(_esat-_e)/1000";
    }

    private string GenerateVPDOperator()
    {
        var sb = new StringBuilder();
        sb.AppendLine("# Calculate VPD using equation file");

        // Create equation file with selected method
        sb.AppendLine("cat > \"${TMPDIR}/vpd_equations.txt\" <<EOF");
        sb.AppendLine(GetVPDEquations(_config.VPDMethod));
        sb.AppendLine("EOF");
        sb.AppendLine();

        // Calculate VPD using the equation file
        return $"exprf,\"${TMPDIR}/vpd_equations.txt\";
    }

    private void WritePBSHeader(TextWriter writer)
    {
        writer.WriteLine("#!/bin/bash");
        writer.WriteLine($"#PBS -N {_config.JobName}");
        writer.WriteLine($"#PBS -P {_config.Project}");
        writer.WriteLine($"#PBS -q {_config.Queue}");
        writer.WriteLine($"#PBS -l walltime={_config.Walltime}");
        writer.WriteLine($"#PBS -l ncpus={_config.Ncpus}");
        writer.WriteLine($"#PBS -l mem={_config.Memory}GB");
        writer.WriteLine($"#PBS -j oe");

        // Add storage directives if required
        var storageDirectives = _config.GetRequiredStorageDirectives();
        if (storageDirectives.Any())
            writer.WriteLine(PBSStorageHelper.FormatStorageDirectives(storageDirectives));

        writer.WriteLine("");  // Add blank line after header
    }

    public string GenerateProcessingScript(IClimateDataset dataset)
    {
        var sb = new StringBuilder();

        // Add PBS header
        using (TextWriter writer = new StringWriter(sb))
            WritePBSHeader(writer);

        sb.AppendLine("set -euo pipefail");
        sb.AppendLine();

        // Load required modules
        sb.AppendLine("module load cdo");
        sb.AppendLine("module load nco");
        sb.AppendLine();

        // Create working directory
        sb.AppendLine("WORKDIR=\"$PBS_O_WORKDIR\"");
        sb.AppendLine("cd \"$WORKDIR\"");
        sb.AppendLine();

        // Create temporary directory and cd into it
        sb.AppendLine("TMPDIR=$(mktemp -d)");
        sb.AppendLine("cd \"$TMPDIR\"");
        sb.AppendLine("trap 'cd \"$WORKDIR\"; rm -rf \"$TMPDIR\"' EXIT");
        sb.AppendLine();

        // Ensure output directory exists
        sb.AppendLine($"mkdir -p \"{_config.OutputDirectory}\"");
        sb.AppendLine();

        // Process each variable
        foreach (ClimateVariable variable in Enum.GetValues<ClimateVariable>())
        {
            VariableInfo varInfo = dataset.GetVariableInfo(variable);
            (string outVar, string targetUnits) = GetStandardConfig(variable);

            sb.AppendLine($"echo \"Processing {variable}...\"");

            string rename = GenerateRenameOperator(varInfo.Name, outVar);
            string conversion = string.Join(" ", GenerateUnitConversionOperators(outVar, varInfo.Units, targetUnits, _config.InputTimeStep));
            string aggregation = GenerateTimeAggregationOperator(variable);
            string unpack = "-unpack";
            string remap = $"-remapcon,{_config.GridFile}";
            string operators = $"{aggregation} {conversion} {rename} {unpack} {remap}";

            string outFileName = dataset.GenerateOutputFileName(variable);
            string tmpFile = Path.Combine("\"${WORKDIR}\"", outFileName);

            // TODO: quote file names.
            // TODO: use variable name for output file?

            // Merge files and perform all operations in a single step.
            string inFiles = string.Join(" ", dataset.GetInputFiles(variable, false));
            sb.AppendLine($"cdo -O mergetime {operators} {inFiles} {tmpFile}");
            sb.AppendLine();

            // Reorder dimensions, improve chunking, and enable compression
            string ordering = "-a lat,lon,time";
            string chunking = $"--cnk_dmn lat,{_config.ChunkSizeSpatial} --cnk_dmn lon,{_config.ChunkSizeSpatial} --cnk_dmn time,{_config.ChunkSizeTime}";
            string compression = _config.CompressOutput ? $"-L{_config.CompressionLevel}" : "";
            string outFile = Path.Combine(_config.OutputDirectory, outFileName);

            sb.AppendLine($"ncpdq -O {ordering} {chunking} {compression} \"{tmpFile}\" \"{outFile}\"");
            sb.AppendLine($"rm -f \"{tmpFile}\"");
            sb.AppendLine();
        }

        return sb.ToString();
    }

    public string GenerateSubmissionScript(string processingScriptPath)
    {
        var sb = new StringBuilder();

        sb.AppendLine("#!/bin/bash");
        sb.AppendLine();

        // Submit the processing script
        sb.AppendLine($"qsub {processingScriptPath}");

        return sb.ToString();
    }
}
