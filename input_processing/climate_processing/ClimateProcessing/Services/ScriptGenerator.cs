using System.Text;
using ClimateProcessing.Models;
using System.IO;
using ClimateProcessing.Units;
using System.Runtime.CompilerServices;
using ClimateProcessing.Extensions;
using System.Web;
using System.Text.RegularExpressions;

[assembly: InternalsVisibleTo("ClimateProcessing.Tests")]

namespace ClimateProcessing.Services;

/// <summary>
/// Default implementation of script generator that works with any climate dataset.
/// </summary>
public class ScriptGenerator : IScriptGenerator<IClimateDataset>
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
    protected readonly ProcessingConfig _config;

    /// <summary>
    /// Name of the variable containing the directory holding all input files
    /// used as operands for the remap command.
    /// </summary>
    protected const string inDirVariable = "IN_DIR";

    /// <summary>
    /// Name of the variable containing the directory holding all remapped input
    /// files used as operands for the mergetime command.
    /// </summary>
    protected const string remapDirVariable = "REMAP_DIR";

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
    /// Get the path to the unoptimised VPD output file for a dataset.
    /// </summary>
    /// <param name="dataset">The dataset.</param>
    /// <returns>The output file path.</returns>
    private string GetUnoptimisedVpdOutputFilePath(IClimateDataset dataset)
    {
        string temperatureFile = GetMergetimeOutputPath(dataset, ClimateVariable.Temperature);
        return GetVpdFilePath(dataset, temperatureFile);
    }

    /// <summary>
    /// Get the path to the optimised VPD output file for a dataset.
    /// </summary>
    /// <param name="dataset">The dataset.</param>
    /// <returns>The output file path.</returns>
    private string GetOptimisedVpdOutputFilePath(IClimateDataset dataset)
    {
        string temperatureFile = GetOutputFilePath(dataset, ClimateVariable.Temperature);
        return GetVpdFilePath(dataset, temperatureFile);
    }

    /// <summary>
    /// Given a temperature file, generate a file path for an equivalent file
    /// containing VPD data.
    /// </summary>
    /// <param name="dataset">The climate dataset.</param>
    /// <param name="temperatureFile">The temperature file.</param>
    /// <returns>The equivalent VPD file path.</returns>
    private string GetVpdFilePath(IClimateDataset dataset, string temperatureFile)
    {
        string fileName = Path.GetFileName(temperatureFile);
        string tempName = dataset.GetVariableInfo(ClimateVariable.Temperature).Name;
        string baseName = fileName.ReplaceFirst($"{tempName}_", "vpd_");
        string outFile = Path.Combine(Path.GetDirectoryName(temperatureFile)!, baseName);
        return outFile;
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
        await WritePBSHeader(writer, jobName, true);

        await writer.WriteLineAsync("# File paths.");
        string humidityFile = GetMergetimeOutputPath(dataset, ClimateVariable.SpecificHumidity);
        await writer.WriteLineAsync($"HUSS_FILE=\"{humidityFile}\"");

        string pressureFile = GetMergetimeOutputPath(dataset, ClimateVariable.SurfacePressure);
        await writer.WriteLineAsync($"PS_FILE=\"{pressureFile}\"");

        string temperatureFile = GetMergetimeOutputPath(dataset, ClimateVariable.Temperature);
        await writer.WriteLineAsync($"TAS_FILE=\"{temperatureFile}\"");

        string inFiles = "\"${HUSS_FILE}\" \"${PS_FILE}\" \"${TAS_FILE}\"";

        // Generate an output file name.
        string outFile = GetUnoptimisedVpdOutputFilePath(dataset);

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
    /// Generate a script for rechunking VPD data.
    /// </summary>
    /// <param name="dataset">The dataset.</param>
    /// <returns>The script path.</returns>
    private async Task<string> GenerateVPDRechunkScript(IClimateDataset dataset)
    {
        string jobName = $"rechunk_vpd_{dataset.DatasetName}";
        string inFile = GetUnoptimisedVpdOutputFilePath(dataset);
        string outFile = GetOptimisedVpdOutputFilePath(dataset);
        return await GenerateRechunkScript(jobName, inFile, outFile, true);
    }

    /// <summary>
    /// Write the PBS header for a job.
    /// </summary>
    /// <param name="writer">The text writer to which the header will be written.</param>
    /// <param name="jobName">The job name.</param>
    private async Task WritePBSHeader(TextWriter writer, string jobName, bool lightweight)
    {
        string logFileName = $"{jobName}.log";
        string logFile = Path.Combine(GetLogPath(), logFileName);
        string streamFile = Path.Combine(GetStreamPath(), logFileName);

        string queue = lightweight ? PBSConstants.QueueNormal : _config.Queue;
        int ncpus = lightweight ? PBSConstants.LightweightNcpus : _config.Ncpus;
        int mem = lightweight ? PBSConstants.LightweightMemory : _config.Memory;

        await writer.WriteLineAsync("#!/usr/bin/env bash");
        await writer.WriteLineAsync($"#PBS -N {jobName}");
        await writer.WriteLineAsync($"#PBS -o {logFile}");
        await writer.WriteLineAsync($"#PBS -P {_config.Project}");
        await writer.WriteLineAsync($"#PBS -q {queue}");
        await writer.WriteLineAsync($"#PBS -l walltime={_config.Walltime}");
        await writer.WriteLineAsync($"#PBS -l ncpus={ncpus}");
        await writer.WriteLineAsync($"#PBS -l mem={mem}GB");
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
        await writer.WriteLineAsync("module load pbs netcdf cdo nco python3/3.12.1");
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
    /// <param name="outputDirectory">The output directory.</param>
    /// <returns>The path to the script directory.</returns>
    private static string GetScriptPath(string outputDirectory)
    {
        string scriptPath = Path.Combine(outputDirectory, scriptDirectory);
        Directory.CreateDirectory(scriptPath);
        return scriptPath;
    }

    /// <summary>
    /// Get the path to the directory in which the scripts will be stored, and
    /// create it if it doesn't exist.
    /// </summary>
    /// <returns>The path to the script directory.</returns>
    private string GetScriptPath()
    {
        return GetScriptPath(_config.OutputDirectory);
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
    /// Check if the VPD calculation depends on the specified variable.
    /// </summary>
    /// <param name="variable">The variable.</param>
    /// <returns>True iff the VPD calculation depends on the specified variable.</returns>
    private bool IsVpdDependency(ClimateVariable variable)
    {
        return variable == ClimateVariable.SpecificHumidity
            || variable == ClimateVariable.SurfacePressure
            || variable == ClimateVariable.Temperature;
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
        Dictionary<ClimateVariable, string> mergetimeScripts = new();
        Dictionary<ClimateVariable, string> rechunkScripts = new();
        foreach (ClimateVariable variable in Enum.GetValues<ClimateVariable>())
        {
            string mergetime = await GenerateVariableMergeScript(dataset, variable);
            string rechunk = await GenerateVariableRechunkScript(dataset, variable);

            mergetimeScripts[variable] = mergetime;
            rechunkScripts[variable] = rechunk;
        }

        string vpdScript = await GenerateVPDScript(dataset);
        string vpdRechunkScript = await GenerateVPDRechunkScript(dataset);

        // Add job submission logic.
        await writer.WriteLineAsync("echo \"Submitting jobs...\"");
        await writer.WriteLineAsync();
        bool vpdEmpty = true;

        foreach (ClimateVariable variable in Enum.GetValues<ClimateVariable>())
        {
            // Submit mergetime script.
            await writer.WriteLineAsync($"JOB_ID=\"$(qsub \"{mergetimeScripts[variable]}\")\"");

            // Append this job to the list of VPD dependencies if necessary.
            if (IsVpdDependency(variable))
            {
                if (vpdEmpty)
                    await writer.WriteLineAsync("VPD_DEPS=\"${JOB_ID}\"");
                else
                    await writer.WriteLineAsync($"VPD_DEPS=\"${{VPD_DEPS}}:${{JOB_ID}}\"");
                vpdEmpty = false;
            }

            // Submit rechunk script.
            await writer.WriteLineAsync($"qsub -W depend=afterok:\"${{JOB_ID}}\" \"{rechunkScripts[variable]}\"");
            await writer.WriteLineAsync();
        }

        // Submit VPD scripts.
        await writer.WriteLineAsync($"JOB_ID=\"$(qsub -W depend=afterok:\"${{VPD_DEPS}}\" \"{vpdScript}\")\"");
        await writer.WriteLineAsync($"qsub -W depend=afterok:\"${{JOB_ID}}\" \"{vpdRechunkScript}\"");
        await writer.WriteLineAsync();

        await writer.WriteLineAsync("echo \"Job submission complete.\"");
        await writer.WriteLineAsync();

        return scriptFile;
    }

    /// <summary>
    /// Create an empty script file and set execute permissions.
    /// </summary>
    /// <param name="scriptPath">The path to the directory in which the script will be created.</param>
    /// <param name="scriptName">The name of the script file to create.</param>
    /// <returns>The path to the script file.</returns>
    private static string CreateScript(string scriptPath, string scriptName)
    {
        string script = Path.Combine(scriptPath, scriptName);
        File.WriteAllText(script, string.Empty);
        if (OperatingSystem.IsLinux() || OperatingSystem.IsMacOS())
            File.SetUnixFileMode(script, UnixFileMode.UserRead | UnixFileMode.UserWrite | UnixFileMode.UserExecute);
        else
            // Windows not supported.
            throw new PlatformNotSupportedException();
        return script;
    }

    /// <summary>
    /// Create an empty script file and set execute permissions.
    /// </summary>
    /// <param name="scriptName">The name of the script file to create.</param>
    /// <returns>The path to the script file.</returns>
    private string CreateScript(string scriptName)
    {
        string scriptPath = GetScriptPath();
        return CreateScript(scriptPath, scriptName);
    }

    /// <summary>
    /// Generate a path that will be used as the output file for the mergetime
    /// operation.
    /// </summary>
    /// <param name="dataset">The climate dataset.</param>
    /// <param name="variable">The variable being processed.</param>
    /// <returns>The path to the output file.</returns>
    private string GetMergetimeOutputPath(IClimateDataset dataset, ClimateVariable variable)
    {
        string directory = Path.Combine(_config.OutputDirectory, "tmp", dataset.GetOutputDirectory());
        Directory.CreateDirectory(directory);
        string outFileName = dataset.GenerateOutputFileName(variable);
        return Path.Combine(directory, outFileName);
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
    /// Write the pre-merge commands to the specified writer.
    /// </summary>
    /// <param name="writer">The text writer.</param>
    /// <param name="dataset">The climate dataset being processed.</param>
    /// <param name="variable">The variable of the dataset being processed.</param>
    protected virtual Task WritePreMerge(TextWriter writer, IClimateDataset dataset, ClimateVariable variable)
    {
        return Task.CompletedTask;
    }

    /// <summary>
    /// Generate a job name.
    /// </summary>
    /// <param name="prefix">Job name prefix which provides context about the job.</param>
    /// <param name="info">Metadata for the variable in the dataset being processed.</param>
    /// <param name="dataset">The dataset being processed.</param>
    /// <returns>A job name.</returns>
    private string GetJobName(string prefix, VariableInfo info, IClimateDataset dataset)
    {
        return $"{prefix}_{info.Name}_{dataset.DatasetName}";
    }

    /// <summary>
    /// Generate a mergetime script for the specified variable, and return the
    /// path to the generated script file.
    /// </summary>
    /// <param name="dataset">The dataset to process.</param>
    /// <param name="variable">The variable to process.</param>
    /// <param name="outFile">The path to the output file that should be generated by the script.</param>
    /// <returns>The path to the generated script file.</returns>
    private async Task<string> GenerateVariableMergeScript(IClimateDataset dataset, ClimateVariable variable)
    {
        VariableInfo varInfo = dataset.GetVariableInfo(variable);
        (string outVar, string targetUnits) = GetStandardConfig(variable);

        // Create script directory if it doesn't already exist.
        // This should be unnecessary at this point.
        string jobName = GetJobName("mergetime", varInfo, dataset);
        string scriptFile = CreateScript(jobName);
        using TextWriter writer = new StreamWriter(scriptFile);
        await WritePBSHeader(writer, jobName, lightweight: true);

        // File paths.
        string inDir = dataset.GetInputFilesDirectory(variable);
        string outFile = GetMergetimeOutputPath(dataset, variable);

        await writer.WriteLineAsync("# File paths.");
        await writer.WriteLineAsync($"{inDirVariable}=\"{inDir}\"");
        if (!string.IsNullOrEmpty(_config.GridFile))
            await writer.WriteLineAsync($"{remapDirVariable}=\"${{WORK_DIR}}/remap\"");
        await writer.WriteLineAsync($"OUT_FILE=\"{outFile}\"");
        if (!string.IsNullOrEmpty(_config.GridFile))
            await writer.WriteLineAsync($"GRID_FILE=\"{_config.GridFile}\"");
        await writer.WriteLineAsync();

        if (!string.IsNullOrEmpty(_config.GridFile))
        {
            await writer.WriteLineAsync($"mkdir -p \"${{{remapDirVariable}}}\"");
            await writer.WriteLineAsync();
        }

        await WritePreMerge(writer, dataset, variable);

        string rename = GenerateRenameOperator(varInfo.Name, outVar);
        string conversion = string.Join(" ", GenerateUnitConversionOperators(outVar, varInfo.Units, targetUnits, _config.InputTimeStep));
        string aggregation = GenerateTimeAggregationOperator(variable);
        string unpack = "-unpack";
        string remapOperator = GetRemapOperator(varInfo, variable);
        string remap = string.IsNullOrEmpty(_config.GridFile) ? string.Empty : $"-{remapOperator},\"${{GRID_FILE}}\"";
        string operators = $"{aggregation} {conversion} {rename} {unpack} {remap}";
        operators = Regex.Replace(operators, " +", " ");

        // The above operators all take a single file as input; therefore we
        // must perform them as a separate step to the mergetime.
        if (!string.IsNullOrWhiteSpace(operators))
        {
            // Write description of processing steps.
            await writer.WriteLineAsync("# Perform corrective operations on input files:");
            if (!string.IsNullOrEmpty(remap))
                await writer.WriteLineAsync("# - Remap input files to target grid.");
            if (!string.IsNullOrEmpty(unpack))
                await writer.WriteLineAsync("# - Unpack data.");
            if (!string.IsNullOrEmpty(rename))
                await writer.WriteLineAsync($"# - Rename variable from {varInfo.Name} to {outVar}.");
            if (!string.IsNullOrEmpty(conversion))
                await writer.WriteLineAsync($"# - Convert units from {varInfo.Units} to {targetUnits}.");
            if (!string.IsNullOrEmpty(aggregation))
                await writer.WriteLineAsync($"# - Aggregate data from {_config.InputTimeStep} to {_config.OutputTimeStep}.");

            await writer.WriteLineAsync($"for FILE in \"${{{inDirVariable}}}\"/*.nc");
            await writer.WriteLineAsync($"do");
            await writer.WriteLineAsync($"    cdo {GetCDOArgs()} {operators} \"${{FILE}}\" \"${{{remapDirVariable}}}/$(basename \"${{FILE}}\")\"");
            await writer.WriteLineAsync("done");
            await writer.WriteLineAsync($"{inDirVariable}=\"${{{remapDirVariable}}}\"");
            await writer.WriteLineAsync();
        }

        // Merge files and perform all operations in a single step.
        await writer.WriteLineAsync("log \"Merging files...\"");
        await writer.WriteLineAsync($"cdo {GetCDOArgs()} mergetime \"${{{inDirVariable}}}\"/*.nc \"${{OUT_FILE}}\"");
        await writer.WriteLineAsync("log \"All files merged successfully.\"");
        await writer.WriteLineAsync();

        // Remapped files are in jobfs, and will be automatically deleted upon
        // job completion (or failure).

        return scriptFile;
    }

    /// <summary>
    /// Generate a rechunking script.
    /// </summary>
    /// <param name="dataset">The dataset.</param>
    /// <param name="variable">The variable.</param>
    /// <param name="inputFile">The path to the file emitted by the mergetime script.</param>
    /// <remarks>
    /// This is separate to the mergetime script, because it's more
    /// resource-intensive. We can reduce costs by having the mergetime script
    /// run on a low-memory node, but the rechunking script needs a high-memory
    /// node.</remarks>
    /// <returns>Path to the generated script file.</returns>
    private async Task<string> GenerateVariableRechunkScript(IClimateDataset dataset, ClimateVariable variable)
    {
        VariableInfo varInfo = dataset.GetVariableInfo(variable);

        string jobName = GetJobName("rechunk", varInfo, dataset);
        string inFile = GetMergetimeOutputPath(dataset, variable);
        string outFile = GetOutputFilePath(dataset, variable);
        bool cleanup = !IsVpdDependency(variable);

        return await GenerateRechunkScript(jobName, inFile, outFile, cleanup);
    }

    private async Task<string> GenerateRechunkScript(string jobName, string inFile, string outFile, bool cleanup)
    {
        string scriptFile = CreateScript(jobName);
        using TextWriter writer = new StreamWriter(scriptFile);

        await WritePBSHeader(writer, jobName, lightweight: false);

        // File paths.
        // The output of the mergetime script is the input file for this script.
        await writer.WriteLineAsync($"IN_FILE=\"{inFile}\"");
        await writer.WriteLineAsync($"OUT_FILE=\"{outFile}\"");
        // Reorder dimensions, improve chunking, and enable compression.

        // Note: we could use lon,lat,time but ncview works better if the
        // x-dimension precedes the y-dimension.
        string ordering = "-a lat,lon,time";
        string chunking = $"--cnk_dmn lat,{_config.ChunkSizeSpatial} --cnk_dmn lon,{_config.ChunkSizeSpatial} --cnk_dmn time,{_config.ChunkSizeTime}";
        string compression = _config.CompressOutput ? $"-L{_config.CompressionLevel}" : "";

        await writer.WriteLineAsync("log \"Rechunking files...\"");
        await writer.WriteLineAsync($"ncpdq -O {ordering} {chunking} {compression} \"${{IN_FILE}}\" \"${{OUT_FILE}}\"");
        await writer.WriteLineAsync("log \"All files rechunked successfully.\"");
        await writer.WriteLineAsync();

        // We can now delete the temporary input file, but only if it's not also
        // required for the VPD estimation, which may not have occurred yet.
        if (cleanup)
        {
            await writer.WriteLineAsync("# Delete temporary file.");
            await writer.WriteLineAsync($"rm -f \"${{TMP_FILE}}\"");
            await writer.WriteLineAsync();
        }
        else
            await writer.WriteLineAsync("# Input file cannot (necessarily) be deleted yet, since it is required for VPD estimation.");

        return scriptFile;
    }

    /// <summary>
    /// Generate a wrapper script that executes the given scripts.
    /// </summary>
    /// <param name="scripts">The script files to execute.</param>
    public static async Task<string> GenerateWrapperScript(string outputDirectory, IEnumerable<string> scripts)
    {
        string jobName = "wrapper";
        string scriptFile = CreateScript(GetScriptPath(outputDirectory), jobName);
        using TextWriter writer = new StreamWriter(scriptFile);

        await writer.WriteLineAsync("#!/usr/bin/env bash");
        await writer.WriteLineAsync("# Master-level script which executes all job submission scripts to submit all PBS jobs.");
        await writer.WriteLineAsync();

        await writer.WriteLineAsync("set -euo pipefail");
        await writer.WriteLineAsync();

        // Execute all scripts (making assumptions about file permissions).
        foreach (string script in scripts)
            await writer.WriteLineAsync($"{script}");
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
        return Regex.IsMatch(units, 
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
