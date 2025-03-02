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
    private const string logDirectory = "logs";
    private const string streamDirectory = "streams";

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

    private void WriteVPDEquations(TextWriter writer, VPDMethod method)
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

        writer.WriteLine($@"# Saturation vapor pressure (Pa)");
        writer.WriteLine(esatEquation);
        writer.WriteLine("# Actual vapor pressure (Pa)");
        writer.WriteLine("_e=(huss*ps)/(0.622+0.378*huss)");
        writer.WriteLine("# VPD (kPa)");
        writer.WriteLine("vpd=(_esat-_e)/1000");
    }

    /// <summary>
    /// Standard arguments used for all CDO invocations.
    /// TODO: make verbosity configurable?
    /// </summary>
    private string GetCDOArgs()
    {
        return $"-L -O -v -z zip1";
    }

    private string GenerateVPDScript(IClimateDataset dataset)
    {
        string jobName = $"calc_vpd_{dataset.DatasetName}";
        string script = CreateScript(jobName);
        using TextWriter writer = new StreamWriter(script);
        WritePBSHeader(writer, jobName);

        writer.WriteLine("# File paths.");
        string humidityFile = GetOutputFilePath(dataset, ClimateVariable.SpecificHumidity);
        writer.WriteLine($"HUSS_FILE=\"{humidityFile}\"");

        string pressureFile = GetOutputFilePath(dataset, ClimateVariable.SurfacePressure);
        writer.WriteLine($"PS_FILE=\"{pressureFile}\"");

        string temperatureFile = GetOutputFilePath(dataset, ClimateVariable.Temperature);
        writer.WriteLine($"TAS_FILE=\"{temperatureFile}\"");

        string inFiles = "\"${HUSS_FILE}\" \"${PS_FILE}\" \"${TAS_FILE}\"";

        // Generate an output file name.
        string fileName = Path.GetFileName(temperatureFile);
        string tempName = dataset.GetVariableInfo(ClimateVariable.Temperature).Name;
        string baseName = fileName.ReplaceFirst($"{tempName}_", "vpd_");
        string outFile = Path.Combine(Path.GetDirectoryName(temperatureFile)!, baseName);

        writer.WriteLine($"OUT_FILE=\"{outFile}\"");

        string eqnFile = "${WORK_DIR}/vpd_equations.txt";
        writer.WriteLine($"EQN_FILE=\"{eqnFile}\"");
        writer.WriteLine();

        writer.WriteLine("# Generate equation file.");
        writer.WriteLine("log \"Generating VPD equation file...\"");

        // Create equation file with selected method.
        writer.WriteLine($"cat >\"${{EQN_FILE}}\" <<EOF");
        WriteVPDEquations(writer, _config.VPDMethod);
        writer.WriteLine("EOF");
        writer.WriteLine();

        // Calculate VPD using the equation file.
        writer.WriteLine("# Calculate VPD.");
        writer.WriteLine($"log \"Calculating VPD...\"");
        writer.WriteLine($"cdo {GetCDOArgs()} exprf,\"${{EQN_FILE}}\" {inFiles} \"${{OUT_FILE}}\"");
        writer.WriteLine($"log \"VPD calculation completed successfully.\"");
        writer.WriteLine();

        // Return the path to the generated script.
        return script;
    }

    private void WritePBSHeader(TextWriter writer, string jobName)
    {
        string logFileName = $"{jobName}.log";
        string logFile = Path.Combine(GetLogPath(), logFileName);
        string streamFile = Path.Combine(GetStreamPath(), logFileName);

        writer.WriteLine("#!/usr/bin/env bash");
        writer.WriteLine($"#PBS -N {jobName}");
        writer.WriteLine($"#PBS -o {logFile}");
        writer.WriteLine($"#PBS -P {_config.Project}");
        writer.WriteLine($"#PBS -q {_config.Queue}");
        writer.WriteLine($"#PBS -l walltime={_config.Walltime}");
        writer.WriteLine($"#PBS -l ncpus={_config.Ncpus}");
        writer.WriteLine($"#PBS -l mem={_config.Memory}GB");
        writer.WriteLine($"#PBS -l jobfs={_config.JobFS}GB");
        writer.WriteLine($"#PBS -j oe");
        if (!string.IsNullOrEmpty(_config.Email))
        {
            writer.WriteLine($"#PBS -M {_config.Email}");
            writer.WriteLine($"#PBS -m abe");
        }

        // Add storage directives if required
        var storageDirectives = _config.GetRequiredStorageDirectives();
        if (storageDirectives.Any())
            writer.WriteLine(PBSStorageHelper.FormatStorageDirectives(storageDirectives));

        // Add blank line after header
        writer.WriteLine("");

        // Error handling.
        writer.WriteLine("# Exit immediately if any command fails.");
        writer.WriteLine("set -euo pipefail");
        writer.WriteLine();

        // Load required modules.
        writer.WriteLine("# Load required modules.");
        writer.WriteLine("module purge");
        writer.WriteLine("module load pbs netcdf cdo nco");
        writer.WriteLine();

        // Create temporary directory and cd into it.
        writer.WriteLine("# Create temporary directory and cd into it.");
        writer.WriteLine("WORK_DIR=\"$(mktemp -d -p \"${PBS_JOBFS}\")\"");
        writer.WriteLine("cd \"${WORK_DIR}\"");
        writer.WriteLine();

        // Technically, deleting the temporary directory is unnecessary, because
        // the tempfs on the compute nodes will be deleted when the job
        // finishes. However, it's a good practice to clean up after ourselves.
        writer.WriteLine("# Delete the temporary directory on exit.");
        writer.WriteLine("trap 'cd \"${PBS_JOBFS}\"; rm -rf \"${WORK_DIR}\"' EXIT");
        writer.WriteLine();

        // Set up logging that streams all output into the stream directory.
        writer.WriteLine("# Stream all output to a log file without buffering.");
        writer.WriteLine($"STREAM_FILE=\"{streamFile}\"");
        writer.WriteLine("rm -f \"${STREAM_FILE}\"");
        writer.WriteLine("exec 1> >(tee -a \"${STREAM_FILE}\") 2>&1");
        writer.WriteLine();

        writer.WriteLine("# Print a log message.");
        writer.WriteLine("log() {");
        writer.WriteLine("    echo \"[$(date)] $*\"");
        writer.WriteLine("}");

        // Add blank line after header.
        writer.WriteLine("");
    }

    private string GetScriptPath()
    {
        string scriptPath = Path.Combine(_config.OutputDirectory, scriptDirectory);
        Directory.CreateDirectory(scriptPath);
        return scriptPath;
    }

    private string GetLogPath()
    {
        string logPath = Path.Combine(_config.OutputDirectory, logDirectory);
        Directory.CreateDirectory(logPath);
        return logPath;
    }

    private string GetStreamPath()
    {
        string streamPath = Path.Combine(_config.OutputDirectory, streamDirectory);
        Directory.CreateDirectory(streamPath);
        return streamPath;
    }

    // Generate processing scripts, and return the path to the top-level script.
    public string GenerateScripts(IClimateDataset dataset)
    {
        string jobName = $"submit_{dataset.DatasetName}";
        string scriptFile = CreateScript(jobName);
        using TextWriter writer = new StreamWriter(scriptFile);

        // Add PBS header.
        writer.WriteLine("#!/usr/bin/env bash");
        writer.WriteLine($"# Job submission script for: {dataset.DatasetName}");
        writer.WriteLine();
        writer.WriteLine("# Exit immediately if any command fails.");
        writer.WriteLine("set -euo pipefail");
        writer.WriteLine();

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
        string directory = Path.Combine(_config.OutputDirectory, dataset.GetOutputDirectory());
        Directory.CreateDirectory(directory);
        string outFileName = dataset.GenerateOutputFileName(variable);
        return Path.Combine(directory, outFileName);
    }

    private string GenerateVariableMergeScript(IClimateDataset dataset, ClimateVariable variable)
    {
        VariableInfo varInfo = dataset.GetVariableInfo(variable);
        (string outVar, string targetUnits) = GetStandardConfig(variable);

        // Create script directory if it doesn't already exist.
        // This should be unnecessary at this point.
        string jobName = $"mergetime_{varInfo.Name}_{dataset.DatasetName}";
        string scriptFile = CreateScript(jobName);
        using TextWriter writer = new StreamWriter(scriptFile);
        WritePBSHeader(writer, jobName);

        // File paths.
        string inDir = dataset.GetInputFilesDirectory(variable);
        string outFileName = dataset.GenerateOutputFileName(variable);
        string tmpFile = Path.Combine("${WORK_DIR}", outFileName);
        string outFile = GetOutputFilePath(dataset, variable);
        writer.WriteLine("# File paths.");
        writer.WriteLine($"IN_DIR=\"{inDir}\"");
        writer.WriteLine($"TMP_FILE=\"{tmpFile}\"");
        writer.WriteLine($"OUT_FILE=\"{outFile}\"");
        writer.WriteLine();

        string rename = GenerateRenameOperator(varInfo.Name, outVar);
        string conversion = string.Join(" ", GenerateUnitConversionOperators(outVar, varInfo.Units, targetUnits, _config.InputTimeStep));
        string aggregation = GenerateTimeAggregationOperator(variable);
        string unpack = "-unpack";
        string remap = string.IsNullOrEmpty(_config.GridFile) ? "" : $"-remapcon,{_config.GridFile}";
        string operators = $"{aggregation} {conversion} {rename} {unpack} {remap}";
        operators = operators.Replace("  ", " ");

        // Merge files and perform all operations in a single step.
        writer.WriteLine("log \"Merging files...\"");
        writer.WriteLine($"cdo {GetCDOArgs()} mergetime {operators} \"${{IN_DIR}}\"/*.nc \"${{TMP_FILE}}\"");
        writer.WriteLine("log \"All files merged successfully.\"");
        writer.WriteLine();

        // Reorder dimensions, improve chunking, and enable compression.
        string ordering = "-a lat,lon,time";
        string chunking = $"--cnk_dmn lat,{_config.ChunkSizeSpatial} --cnk_dmn lon,{_config.ChunkSizeSpatial} --cnk_dmn time,{_config.ChunkSizeTime}";
        string compression = _config.CompressOutput ? $"-L{_config.CompressionLevel}" : "";

        writer.WriteLine("log \"Rechunking files...\"");
        writer.WriteLine($"ncpdq -O {ordering} {chunking} {compression} \"${{TMP_FILE}}\" \"${{OUT_FILE}}\"");
        writer.WriteLine("log \"All files rechunked successfully.\"");
        writer.WriteLine();

        writer.WriteLine("# Delete temporary file.");
        writer.WriteLine($"rm -f \"${{TMP_FILE}}\"");

        return scriptFile;
    }
}
