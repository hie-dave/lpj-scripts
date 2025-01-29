using System.Text;
using ClimateProcessing.Models;
using System.IO;
using ClimateProcessing.Units;
using System.Runtime.CompilerServices;
using ClimateProcessing.Extensions;

[assembly: InternalsVisibleTo("ClimateProcessing.Tests")]

namespace ClimateProcessing.Services;

public class ScriptGenerator
{
    // Subdirectory of the output directory into which the scripts are written.
    private const string scriptDirectory = "scripts";

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
        // TBI
        var sb = new StringBuilder();
        sb.AppendLine("# Calculate VPD using equation file");

        // Create equation file with selected method
        sb.AppendLine("cat > \"${TMPDIR}/vpd_equations.txt\" <<EOF");
        sb.AppendLine(GetVPDEquations(_config.VPDMethod));
        sb.AppendLine("EOF");
        sb.AppendLine();

        // Calculate VPD using the equation file
        return "exprf,\"${TMPDIR}/vpd_equations.txt\"";
    }

    private string GenerateVPDScript(IClimateDataset dataset)
    {
        string script = CreateScript("calc_vpd");
        using TextWriter writer = new StreamWriter(script);
        WritePBSHeader(writer);
        writer.WriteLine();

        writer.WriteLine("# Calculate VPD using equation file");

        // Create equation file with selected method.
        string eqnFile = "${TMPDIR}/vpd_equations.txt";
        writer.WriteLine($"cat >\"{eqnFile}\" <<EOF");
        writer.WriteLine(GetVPDEquations(_config.VPDMethod));
        writer.WriteLine("EOF");
        writer.WriteLine();

        string humidityFile = GetOutputFilePath(dataset, ClimateVariable.SpecificHumidity);
        string pressureFile = GetOutputFilePath(dataset, ClimateVariable.SurfacePressure);
        string temperatureFile = GetOutputFilePath(dataset, ClimateVariable.Temperature);
        string inFiles = $"\"{humidityFile}\" \"{pressureFile}\" \"{temperatureFile}\"";

        // Generate an output file name.
        string fileName = Path.GetFileName(temperatureFile);
        string tempName = dataset.GetVariableInfo(ClimateVariable.Temperature).Name;
        string baseName = fileName.ReplaceFirst($"{tempName}_", "vpd_");
        string outFile = Path.Combine(Path.GetDirectoryName(temperatureFile)!, baseName);

        // Calculate VPD using the equation file.
        writer.WriteLine($"cdo -O exprf,{eqnFile} {inFiles} {outFile}");
        writer.WriteLine();

        // Return the path to the generated script.
        return script;
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

        // Add blank line after header
        writer.WriteLine("");

        // Error handling.
        writer.WriteLine("set -euo pipefail");
        writer.WriteLine();

        // Load required modules.
        writer.WriteLine("module load cdo");
        writer.WriteLine("module load nco");
        writer.WriteLine();

        // Create working directory.
        writer.WriteLine("WORKDIR=\"${PBS_O_WORKDIR}\"");
        writer.WriteLine("cd \"${WORKDIR}\"");
        writer.WriteLine();

        // Create temporary directory and cd into it.
        writer.WriteLine("TMPDIR=\"$(mktemp -d)\"");
        writer.WriteLine("cd \"${TMPDIR}\"");
        writer.WriteLine("trap 'cd \"${WORKDIR}\"; rm -rf \"${TMPDIR}\"' EXIT");
        writer.WriteLine();
    }

    private string GetScriptPath()
    {
        string scriptPath = Path.Combine(_config.OutputDirectory, scriptDirectory);
        Directory.CreateDirectory(scriptPath);
        return scriptPath;
    }

    // Generate processing scripts, and return the path to the top-level script.
    public string GenerateScripts(IClimateDataset dataset)
    {
        string scriptFile = CreateScript($"process_{dataset.DatasetName}");
        using (TextWriter writer = new StreamWriter(scriptFile))
            GenerateProcessingScript(writer, dataset);

        return scriptFile;
    }

    // Create an empty script file and set execute permissions.
    private string CreateScript(string scriptName)
    {
        string script = Path.Combine(GetScriptPath(), scriptName);
        File.WriteAllText(script, string.Empty);
        if (OperatingSystem.IsLinux() || OperatingSystem.IsMacOS())
            File.SetUnixFileMode(script, UnixFileMode.UserRead | UnixFileMode.UserWrite | UnixFileMode.UserExecute);
        else
            // Windows not supported.
            throw new PlatformNotSupportedException();
        return script;
    }

    private string GetOutputFilePath(IClimateDataset dataset, ClimateVariable variable)
    {
        string outFileName = dataset.GenerateOutputFileName(variable);
        return Path.Combine(_config.OutputDirectory, outFileName);
    }

    private string GenerateVariableMergeScript(IClimateDataset dataset, ClimateVariable variable)
    {
        VariableInfo varInfo = dataset.GetVariableInfo(variable);
        (string outVar, string targetUnits) = GetStandardConfig(variable);

        // Create script directory if it doesn't already exist.
        // This should be unnecessary at this point.
        string scriptFile = CreateScript($"mergetime_{varInfo.Name}");
        using TextWriter writer = new StreamWriter(scriptFile);
        WritePBSHeader(writer);

        writer.WriteLine($"echo \"Processing {variable}...\"");

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
        writer.WriteLine($"cdo -O mergetime {operators} {inFiles} {tmpFile}");
        writer.WriteLine();

        // Reorder dimensions, improve chunking, and enable compression.
        string ordering = "-a lat,lon,time";
        string chunking = $"--cnk_dmn lat,{_config.ChunkSizeSpatial} --cnk_dmn lon,{_config.ChunkSizeSpatial} --cnk_dmn time,{_config.ChunkSizeTime}";
        string compression = _config.CompressOutput ? $"-L{_config.CompressionLevel}" : "";
        string outFile = GetOutputFilePath(dataset, variable);

        writer.WriteLine($"ncpdq -O {ordering} {chunking} {compression} \"{tmpFile}\" \"{outFile}\"");
        writer.WriteLine($"rm -f \"{tmpFile}\"");
        writer.WriteLine();

        return scriptFile;
    }

    private void GenerateProcessingScript(TextWriter writer, IClimateDataset dataset)
    {
        // Add PBS header.
        WritePBSHeader(writer);

        // Ensure output directory exists.
        writer.WriteLine($"mkdir -p \"{_config.OutputDirectory}\"");
        writer.WriteLine();

        // Process each variable.
        List<string> vpdDependencies = new List<string>();
        List<string> variableScripts = new List<string>();
        foreach (ClimateVariable variable in Enum.GetValues<ClimateVariable>())
        {
            string subscript = GenerateVariableMergeScript(dataset, variable);
            if (variable == ClimateVariable.SpecificHumidity
                || variable == ClimateVariable.SurfacePressure
                || variable == ClimateVariable.Temperature)
                vpdDependencies.Add(subscript);
            else
                variableScripts.Add(subscript);
        }

        // These checks should be redundant.
        if (vpdDependencies.Count == 0)
            throw new InvalidOperationException("No VPD dependencies were generated.");
        if (variableScripts.Count == 0)
            throw new InvalidOperationException("No scripts were generated.");

        string vpdScript = GenerateVPDScript(dataset);

        // Add job submission logic.
        writer.WriteLine("echo \"Submitting jobs...\"");
        writer.WriteLine();

        writer.WriteLine($"DEPS=\"$(qsub \"{variableScripts[0]}\")\"");
        for (int i = 1; i < variableScripts.Count; i++)
            writer.WriteLine($"DEPS=\"${{DEPS}}:$(qsub \"{variableScripts[i]}\")\"");

        writer.WriteLine($"VPD_DEPS=\"$(qsub \"{vpdDependencies[0]}\")\"");
        for (int i = 1; i < vpdDependencies.Count; i++)
            writer.WriteLine($"VPD_DEPS=\"${{VPD_DEPS}}:$(qsub \"{vpdDependencies[i]}\")\"");

        writer.WriteLine($"DEPS=\"${{DEPS}}:$(qsub -W depend=afterok:\"${{VPD_DEPS}}\" \"{vpdScript}\")\"");
        writer.WriteLine();

        writer.WriteLine("echo \"Job submission complete.\"");
        writer.WriteLine();
    }
}
