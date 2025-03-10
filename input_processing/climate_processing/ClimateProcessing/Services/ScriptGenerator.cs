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
    /// <summary>
    /// Subdirectory of the output directory into which the scripts are written.
    /// </summary>
    private const string scriptDirectory = "scripts";

    /// <summary>
    /// Subdirectory of the output directory into which the logs are written.
    /// </summary>
    private const string logDirectory = "logs";

    /// <summary>
    /// Subdirectory of the output directory into which stdout will be streamed.
    /// </summary>
    private const string streamDirectory = "streams";

    /// <summary>
    /// CDO's conservative remapping operator.
    /// </summary>
    private const string remapConservative = "-remapcon";

    /// <summary>
    /// CDO's bilinear remapping operator.
    /// </summary>
    private const string remapBilinear = "remapbil";

    /// <summary>
    /// The processing configuration.
    /// </summary>
    private readonly ProcessingConfig _config;

    /// <summary>
    /// List of standard variables and their output names and units.
    /// </summary>
    private static readonly Dictionary<ClimateVariable, (string outName, string outUnits)> _standardVariables = new()
    {
        { ClimateVariable.SpecificHumidity, ("huss", "1") },
        { ClimateVariable.SurfacePressure, ("ps", "Pa") },
        { ClimateVariable.ShortwaveRadiation, ("rsds", "W m-2") },
        { ClimateVariable.WindSpeed, ("sfcWind", "m s-1") },
        { ClimateVariable.Temperature, ("tas", "degC") },
        { ClimateVariable.Precipitation, ("pr", "mm") }
    };

    /// <summary>
    /// Gets the standard configuration for the specified variable.
    /// </summary>
    /// <param name="variable">The variable.</param>
    /// <returns>The standard configuration.</returns>
    private (string outName, string outUnits) GetStandardConfig(ClimateVariable variable)
    {
        if (!_standardVariables.TryGetValue(variable, out var config))
            throw new ArgumentException($"No configuration found for variable {variable}");
        return config;
    }

    /// <summary>
    /// Creates a new script generator.
    /// </summary>
    /// <param name="config">The processing configuration.</param>
    public ScriptGenerator(ProcessingConfig config)
    {
        _config = config;
    }

    /// <summary>
    /// Generates the operator to rename a variable.
    /// </summary>
    /// <param name="inName">The name of the input variable.</param>
    /// <param name="outName">The name of the output variable.</param>
    /// <returns>The CDO operator to use for renaming.</returns>
    internal string GenerateRenameOperator(string inName, string outName)
    {
        if (inName == outName)
            return string.Empty;
        return $"-chname,'{inName}','{outName}'";
    }

    /// <summary>
    /// Generates the operators needed to convert the units of a variable.
    /// </summary>
    /// <param name="outputVar">The name of the output variable.</param>
    /// <param name="inputUnits">The units of the input variable.</param>
    /// <param name="targetUnits">The units of the output variable.</param>
    /// <param name="timeStep">The time step of the variable.</param>
    /// <returns>The CDO operators needed to convert the units.</returns>
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

    /// <summary>
    /// Generate the operator to temporally aggregate the data.
    /// </summary>
    /// <param name="variable">The variable to aggregate.</param>
    /// <returns>The CDO operator to use for temporal aggregation.</returns>
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

    /// <summary>
    /// Write the equations for estimating VPD to a text writer.
    /// </summary>
    /// <param name="writer">The writer to write to.</param>
    /// <param name="method">The VPD estimation method to use.</param>
    /// <exception cref="ArgumentException">If the specified VPD method is not supported.</exception>
    private async Task WriteVPDEquations(TextWriter writer, VPDMethod method)
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

        await writer.WriteLineAsync($@"# Saturation vapor pressure (Pa)");
        await writer.WriteLineAsync(esatEquation);
        await writer.WriteLineAsync("# Actual vapor pressure (Pa)");
        await writer.WriteLineAsync("_e=(huss*ps)/(0.622+0.378*huss)");
        await writer.WriteLineAsync("# VPD (kPa)");
        await writer.WriteLineAsync("vpd=(_esat-_e)/1000");
    }

    /// <summary>
    /// Standard arguments used for all CDO invocations.
    /// TODO: make verbosity configurable?
    /// </summary>
    private string GetCDOArgs()
    {
        return $"-L -O -v -z zip1";
    }

    /// <summary>
    /// Generate a processing script for calculating the VPD for a dataset.
    /// </summary>
    /// <param name="dataset">The dataset.</param>
    /// <returns>The script.</returns>
    private async Task<string> GenerateVPDScript(IClimateDataset dataset)
    {
        string jobName = $"calc_vpd_{dataset.DatasetName}";
        string script = CreateScript(jobName);
        using TextWriter writer = new StreamWriter(script);
        await WritePBSHeader(writer, jobName);

        await writer.WriteLineAsync("# File paths.");
        string humidityFile = GetOutputFilePath(dataset, ClimateVariable.SpecificHumidity);
        await writer.WriteLineAsync($"HUSS_FILE=\"{humidityFile}\"");

        string pressureFile = GetOutputFilePath(dataset, ClimateVariable.SurfacePressure);
        await writer.WriteLineAsync($"PS_FILE=\"{pressureFile}\"");

        string temperatureFile = GetOutputFilePath(dataset, ClimateVariable.Temperature);
        await writer.WriteLineAsync($"TAS_FILE=\"{temperatureFile}\"");

        string inFiles = "\"${HUSS_FILE}\" \"${PS_FILE}\" \"${TAS_FILE}\"";

        // Generate an output file name.
        string fileName = Path.GetFileName(temperatureFile);
        string tempName = dataset.GetVariableInfo(ClimateVariable.Temperature).Name;
        string baseName = fileName.ReplaceFirst($"{tempName}_", "vpd_");
        string outFile = Path.Combine(Path.GetDirectoryName(temperatureFile)!, baseName);

        await writer.WriteLineAsync($"OUT_FILE=\"{outFile}\"");

        string eqnFile = "${WORK_DIR}/vpd_equations.txt";
        await writer.WriteLineAsync($"EQN_FILE=\"{eqnFile}\"");
        await writer.WriteLineAsync();

        await writer.WriteLineAsync("# Generate equation file.");
        await writer.WriteLineAsync("log \"Generating VPD equation file...\"");

        // Create equation file with selected method.
        await writer.WriteLineAsync($"cat >\"${{EQN_FILE}}\" <<EOF");
        await WriteVPDEquations(writer, _config.VPDMethod);
        await writer.WriteLineAsync("EOF");
        await writer.WriteLineAsync();

        // Calculate VPD using the equation file.
        await writer.WriteLineAsync("# Calculate VPD.");
        await writer.WriteLineAsync($"log \"Calculating VPD...\"");
        await writer.WriteLineAsync($"cdo {GetCDOArgs()} exprf,\"${{EQN_FILE}}\" {inFiles} \"${{OUT_FILE}}\"");
        await writer.WriteLineAsync($"log \"VPD calculation completed successfully.\"");
        await writer.WriteLineAsync();

        // Return the path to the generated script.
        return script;
    }

    /// <summary>
    /// Write the PBS header for a job.
    /// </summary>
    /// <param name="writer">The text writer to which the header will be written.</param>
    /// <param name="jobName">The job name.</param>
    private async Task WritePBSHeader(TextWriter writer, string jobName)
    {
        string logFileName = $"{jobName}.log";
        string logFile = Path.Combine(GetLogPath(), logFileName);
        string streamFile = Path.Combine(GetStreamPath(), logFileName);

        await writer.WriteLineAsync("#!/usr/bin/env bash");
        await writer.WriteLineAsync($"#PBS -N {jobName}");
        await writer.WriteLineAsync($"#PBS -o {logFile}");
        await writer.WriteLineAsync($"#PBS -P {_config.Project}");
        await writer.WriteLineAsync($"#PBS -q {_config.Queue}");
        await writer.WriteLineAsync($"#PBS -l walltime={_config.Walltime}");
        await writer.WriteLineAsync($"#PBS -l ncpus={_config.Ncpus}");
        await writer.WriteLineAsync($"#PBS -l mem={_config.Memory}GB");
        await writer.WriteLineAsync($"#PBS -l jobfs={_config.JobFS}GB");
        await writer.WriteLineAsync($"#PBS -j oe");
        if (!string.IsNullOrEmpty(_config.Email))
        {
            await writer.WriteLineAsync($"#PBS -M {_config.Email}");
            await writer.WriteLineAsync($"#PBS -m abe");
        }

        // Add storage directives if required
        var storageDirectives = _config.GetRequiredStorageDirectives();
        if (storageDirectives.Any())
            await writer.WriteLineAsync(PBSStorageHelper.FormatStorageDirectives(storageDirectives));

        // Add blank line after header
        await writer.WriteLineAsync("");

        // Error handling.
        await writer.WriteLineAsync("# Exit immediately if any command fails.");
        await writer.WriteLineAsync("set -euo pipefail");
        await writer.WriteLineAsync();

        // Load required modules.
        await writer.WriteLineAsync("# Load required modules.");
        await writer.WriteLineAsync("module purge");
        await writer.WriteLineAsync("module load pbs netcdf cdo nco");
        await writer.WriteLineAsync();

        // Create temporary directory and cd into it.
        await writer.WriteLineAsync("# Create temporary directory and cd into it.");
        await writer.WriteLineAsync("WORK_DIR=\"$(mktemp -d -p \"${PBS_JOBFS}\")\"");
        await writer.WriteLineAsync("cd \"${WORK_DIR}\"");
        await writer.WriteLineAsync();

        // Technically, deleting the temporary directory is unnecessary, because
        // the tempfs on the compute nodes will be deleted when the job
        // finishes. However, it's a good practice to clean up after ourselves.
        await writer.WriteLineAsync("# Delete the temporary directory on exit.");
        await writer.WriteLineAsync("trap 'cd \"${PBS_JOBFS}\"; rm -rf \"${WORK_DIR}\"' EXIT");
        await writer.WriteLineAsync();

        // Set up logging that streams all output into the stream directory.
        await writer.WriteLineAsync("# Stream all output to a log file without buffering.");
        await writer.WriteLineAsync($"STREAM_FILE=\"{streamFile}\"");
        await writer.WriteLineAsync("rm -f \"${STREAM_FILE}\"");
        await writer.WriteLineAsync("exec 1> >(tee -a \"${STREAM_FILE}\") 2>&1");
        await writer.WriteLineAsync();

        await writer.WriteLineAsync("# Print a log message.");
        await writer.WriteLineAsync("log() {");
        await writer.WriteLineAsync("    echo \"[$(date)] $*\"");
        await writer.WriteLineAsync("}");

        // Add blank line after header.
        await writer.WriteLineAsync("");
    }

    /// <summary>
    /// Get the path to the directory in which the scripts will be stored, and
    /// create it if it doesn't exist.
    /// </summary>
    /// <returns>The path to the script directory.</returns>
    private string GetScriptPath()
    {
        string scriptPath = Path.Combine(_config.OutputDirectory, scriptDirectory);
        Directory.CreateDirectory(scriptPath);
        return scriptPath;
    }

    /// <summary>
    /// Get the path to the directory in which log files will be stored, and
    /// create it if it doesn't exist.
    /// </summary>
    /// <returns>The path to the log directory.</returns>
    private string GetLogPath()
    {
        string logPath = Path.Combine(_config.OutputDirectory, logDirectory);
        Directory.CreateDirectory(logPath);
        return logPath;
    }

    /// <summary>
    /// Get the path to the directory in which stdout/stderr will be streamed,
    /// and create it if it doesn't exist.
    /// </summary>
    /// <returns>The path to the stream directory.</returns>
    private string GetStreamPath()
    {
        string streamPath = Path.Combine(_config.OutputDirectory, streamDirectory);
        Directory.CreateDirectory(streamPath);
        return streamPath;
    }

    /// <summary>
    /// Generate processing scripts, and return the path to the top-level script.
    /// </summary>
    /// <param name="dataset">The dataset.</param>
    /// <returns>The path to the top-level script.</returns>
    public async Task<string> GenerateScriptsAsync(IClimateDataset dataset)
    {
        string jobName = $"submit_{dataset.DatasetName}";
        string scriptFile = CreateScript(jobName);
        using TextWriter writer = new StreamWriter(scriptFile);

        // Add PBS header.
        await writer.WriteLineAsync("#!/usr/bin/env bash");
        await writer.WriteLineAsync($"# Job submission script for: {dataset.DatasetName}");
        await writer.WriteLineAsync();
        await writer.WriteLineAsync("# Exit immediately if any command fails.");
        await writer.WriteLineAsync("set -euo pipefail");
        await writer.WriteLineAsync();

        // Ensure output directory exists.
        await writer.WriteLineAsync($"mkdir -p \"{_config.OutputDirectory}\"");
        await writer.WriteLineAsync();

        // Process each variable.
        List<string> vpdDependencies = new List<string>();
        List<string> variableScripts = new List<string>();
        foreach (ClimateVariable variable in Enum.GetValues<ClimateVariable>())
        {
            string subscript = await GenerateVariableMergeScript(dataset, variable);
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

        string vpdScript = await GenerateVPDScript(dataset);

        // Add job submission logic.
        await writer.WriteLineAsync("echo \"Submitting jobs...\"");
        await writer.WriteLineAsync();

        await writer.WriteLineAsync($"DEPS=\"$(qsub \"{variableScripts[0]}\")\"");
        for (int i = 1; i < variableScripts.Count; i++)
            await writer.WriteLineAsync($"DEPS=\"${{DEPS}}:$(qsub \"{variableScripts[i]}\")\"");

        await writer.WriteLineAsync($"VPD_DEPS=\"$(qsub \"{vpdDependencies[0]}\")\"");
        for (int i = 1; i < vpdDependencies.Count; i++)
            await writer.WriteLineAsync($"VPD_DEPS=\"${{VPD_DEPS}}:$(qsub \"{vpdDependencies[i]}\")\"");

        await writer.WriteLineAsync($"DEPS=\"${{DEPS}}:$(qsub -W depend=afterok:\"${{VPD_DEPS}}\" \"{vpdScript}\")\"");
        await writer.WriteLineAsync();

        await writer.WriteLineAsync("echo \"Job submission complete.\"");
        await writer.WriteLineAsync();

        return scriptFile;
    }

    /// <summary>
    /// Create an empty script file and set execute permissions.
    /// </summary>
    /// <param name="scriptName">The name of the script file to create.</param>
    /// <returns>The path to the script file.</returns>
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

    /// <summary>
    /// Generate a path that will be used as the output file for a particular
    /// variable in a dataset.
    /// </summary>
    /// <param name="dataset">The dataset.</param>
    /// <param name="variable">The variable.</param>
    /// <returns>The path to the output file.</returns>
    private string GetOutputFilePath(IClimateDataset dataset, ClimateVariable variable)
    {
        string directory = Path.Combine(_config.OutputDirectory, dataset.GetOutputDirectory());
        Directory.CreateDirectory(directory);
        string outFileName = dataset.GenerateOutputFileName(variable);
        return Path.Combine(directory, outFileName);
    }

    /// <summary>
    /// Generate a mergetime script for the specified variable, and return the
    /// path to the generated script file.
    /// </summary>
    /// <param name="dataset">The dataset to process.</param>
    /// <param name="variable">The variable to process.</param>
    /// <returns>The path to the generated script file.</returns>
    private async Task<string> GenerateVariableMergeScript(IClimateDataset dataset, ClimateVariable variable)
    {
        VariableInfo varInfo = dataset.GetVariableInfo(variable);
        (string outVar, string targetUnits) = GetStandardConfig(variable);

        // Create script directory if it doesn't already exist.
        // This should be unnecessary at this point.
        string jobName = $"mergetime_{varInfo.Name}_{dataset.DatasetName}";
        string scriptFile = CreateScript(jobName);
        using TextWriter writer = new StreamWriter(scriptFile);
        await WritePBSHeader(writer, jobName);

        // File paths.
        string inDir = dataset.GetInputFilesDirectory(variable);
        string outFileName = dataset.GenerateOutputFileName(variable);
        string tmpFile = Path.Combine("${WORK_DIR}", outFileName);
        string outFile = GetOutputFilePath(dataset, variable);
        await writer.WriteLineAsync("# File paths.");
        await writer.WriteLineAsync($"IN_DIR=\"{inDir}\"");
        await writer.WriteLineAsync($"TMP_FILE=\"{tmpFile}\"");
        await writer.WriteLineAsync($"OUT_FILE=\"{outFile}\"");
        await writer.WriteLineAsync();

        string rename = GenerateRenameOperator(varInfo.Name, outVar);
        string conversion = string.Join(" ", GenerateUnitConversionOperators(outVar, varInfo.Units, targetUnits, _config.InputTimeStep));
        string aggregation = GenerateTimeAggregationOperator(variable);
        string unpack = "-unpack";
        string remapOperator = GetRemapOperator(varInfo, variable);
        string remap = string.IsNullOrEmpty(_config.GridFile) ? "" : $"-{remapOperator},{_config.GridFile}";
        string operators = $"{aggregation} {conversion} {rename} {unpack} {remap}";
        operators = operators.Replace("  ", " ");

        // Merge files and perform all operations in a single step.
        await writer.WriteLineAsync("log \"Merging files...\"");
        await writer.WriteLineAsync($"cdo {GetCDOArgs()} mergetime {operators} \"${{IN_DIR}}\"/*.nc \"${{TMP_FILE}}\"");
        await writer.WriteLineAsync("log \"All files merged successfully.\"");
        await writer.WriteLineAsync();

        // Reorder dimensions, improve chunking, and enable compression.
        string ordering = "-a lat,lon,time";
        string chunking = $"--cnk_dmn lat,{_config.ChunkSizeSpatial} --cnk_dmn lon,{_config.ChunkSizeSpatial} --cnk_dmn time,{_config.ChunkSizeTime}";
        string compression = _config.CompressOutput ? $"-L{_config.CompressionLevel}" : "";

        await writer.WriteLineAsync("log \"Rechunking files...\"");
        await writer.WriteLineAsync($"ncpdq -O {ordering} {chunking} {compression} \"${{TMP_FILE}}\" \"${{OUT_FILE}}\"");
        await writer.WriteLineAsync("log \"All files rechunked successfully.\"");
        await writer.WriteLineAsync();

        await writer.WriteLineAsync("# Delete temporary file.");
        await writer.WriteLineAsync($"rm -f \"${{TMP_FILE}}\"");
        await writer.WriteLineAsync();

        return scriptFile;
    }

    /// <summary>
    /// Check if a variable is expressed on a per-ground-area basis.
    /// </summary>
    /// <param name="units">The units of the variable.</param>
    /// <returns>True iff the variable is expressed on a per-ground-area basis.</returns>
    /// <remarks>This is used to decide whether to perform conservative remapping.</remarks>
    internal static bool HasPerAreaUnits(string units)
    {
        // Convert to lowercase and remove whitespace and periods for consistent matching
        units = units.ToLower().Replace(" ", "").Replace(".", "");
        
        // Match any of these patterns:
        // - m-2 or m^-2 (negative exponent notation)
        // - /m2 (division notation)
        return System.Text.RegularExpressions.Regex.IsMatch(units, 
            @"(m\^?-2|/m2)");
    }

    /// <summary>
    /// Get an interpolation algorithm to be used when remapping the specified
    /// variable.
    /// </summary>
    /// <param name="info">Metadata for the variable in the dataset being processed.</param>
    /// <param name="variable">The variable to remap.</param>
    /// <returns>The interpolation algorithm to use.</returns>
    internal InterpolationAlgorithm GetInterpolationAlgorithm(VariableInfo info, ClimateVariable variable)
    {
        // Precipitation and shortwave radiation may require conservative
        // remapping, if they are NOT expressed on a per-ground-area basis.
        if (variable != ClimateVariable.Precipitation
            && variable != ClimateVariable.ShortwaveRadiation)
            return InterpolationAlgorithm.Bilinear;

        // Check if units are expressed on a per-ground-area basis.
        if (!HasPerAreaUnits(info.Units))
            // E.g. W
            // E.g. kg s-1
            return InterpolationAlgorithm.Conservative;

        // If, for example, precipitation is expressed in kg m-2 s-1, there's
        // no need for conservative remapping.
        return InterpolationAlgorithm.Bilinear;
    }

    /// <summary>
    /// Get the CDO remap operator to be used when remapping the specified
    /// variable.
    /// </summary>
    /// <param name="info">Metadata for the variable in the dataset being processed.</param>
    /// <param name="variable">The variable to remap.</param>
    /// <returns>The CDO remap operator to use.</returns>
    /// <exception cref="ArgumentException"></exception>
    private string GetRemapOperator(VariableInfo info, ClimateVariable variable)
    {
        return GetInterpolationAlgorithm(info, variable) switch
        {
            InterpolationAlgorithm.Bilinear => remapBilinear,
            InterpolationAlgorithm.Conservative => remapConservative,
            _ => throw new ArgumentException($"Unknown remap algorithm: {GetInterpolationAlgorithm(info, variable)}")
        };
    }
}
