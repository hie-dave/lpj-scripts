using System.Text;
using ClimateProcessing.Models;
using System.IO;
using ClimateProcessing.Units;
using System.Runtime.CompilerServices;
using ClimateProcessing.Extensions;
using System.Web;
using System.Text.RegularExpressions;
using ClimateProcessing.Configuration;

[assembly: InternalsVisibleTo("ClimateProcessing.Tests")]

namespace ClimateProcessing.Services;

/// <summary>
/// Default implementation of script generator that works with any climate dataset.
/// </summary>
public class ScriptGenerator : IScriptGenerator<IClimateDataset>
{
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
    /// The path manager service.
    /// </summary>
    protected readonly IPathManager pathManager;

    /// <summary>
    /// The PBS script generator service.
    /// </summary>
    protected readonly PBSWriter pbsHeavyweight;

    /// <summary>
    /// The PBS script generator service.
    /// </summary>
    protected readonly PBSWriter pbsLightweight;

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
    private static readonly Dictionary<ClimateVariable, string> outputNames = new()
    {
        { ClimateVariable.SpecificHumidity, "huss" }, // "1"
        { ClimateVariable.SurfacePressure, "ps" }, // "Pa"
        { ClimateVariable.ShortwaveRadiation, "rsds" }, // "W m-2"
        { ClimateVariable.WindSpeed, "sfcWind" }, // "m s-1"
        { ClimateVariable.Temperature, "tas" }, // "degC"
        { ClimateVariable.Precipitation, "pr" }, // "mm"
        { ClimateVariable.MaxTemperature, "tasmax" }, // "degC"
        { ClimateVariable.MinTemperature, "tasmin" }, // "degC"
    };

    private static readonly Dictionary<ClimateVariable, (string units, AggregationMethod aggregation)> daveVariables = new()
    {
        { ClimateVariable.Temperature, ("degC", AggregationMethod.Mean) },
        { ClimateVariable.Precipitation, ("mm", AggregationMethod.Sum) },
        { ClimateVariable.SpecificHumidity, ("1", AggregationMethod.Mean) },
        { ClimateVariable.SurfacePressure, ("Pa", AggregationMethod.Mean) },
        { ClimateVariable.ShortwaveRadiation, ("W m-2", AggregationMethod.Mean) },
        { ClimateVariable.WindSpeed, ("m s-1", AggregationMethod.Mean) }
    };

    private static readonly Dictionary<ClimateVariable, (string units, AggregationMethod aggregation)> trunkVariables = new()
    {
        { ClimateVariable.Temperature, ("K", AggregationMethod.Mean) },
        { ClimateVariable.Precipitation, ("mm", AggregationMethod.Sum) },
        { ClimateVariable.SpecificHumidity, ("1", AggregationMethod.Mean) },
        { ClimateVariable.SurfacePressure, ("Pa", AggregationMethod.Mean) },
        { ClimateVariable.ShortwaveRadiation, ("W m-2", AggregationMethod.Mean) },
        { ClimateVariable.WindSpeed, ("m s-1", AggregationMethod.Mean) },
        { ClimateVariable.MaxTemperature, ("K", AggregationMethod.Maximum) },
        { ClimateVariable.MinTemperature, ("K", AggregationMethod.Minimum) },
    };

    /// <summary>
    /// Creates a new script generator.
    /// </summary>
    /// <param name="config">The processing configuration.</param>
    public ScriptGenerator(ProcessingConfig config)
    {
        _config = config;
        pathManager = new PathManager(config.OutputDirectory);
        PBSWalltime walltime = PBSWalltime.Parse(config.Walltime);
        PBSConfig pbsConfig = new(
            config.Queue,
            config.Ncpus,
            config.Memory,
            config.JobFS,
            config.Project,
            walltime,
            config.EmailNotifications,
            config.Email
        );
        pbsHeavyweight = new PBSWriter(pbsConfig, pathManager);

        PBSConfig lightweightConfig = PBSConfig.LightWeight(
            config.JobFS,
            config.Project,
            config.EmailNotifications,
            config.Email,
            walltime
        );
        pbsLightweight = new PBSWriter(lightweightConfig, pathManager);
    }

    /// <summary>
    /// Gets the standard configuration for the specified variable.
    /// </summary>
    /// <param name="variable">The variable.</param>
    /// <returns>The standard configuration.</returns>
    private (string outName, string outUnits) GetStandardConfig(ClimateVariable variable)
    {
        if (!outputNames.TryGetValue(variable, out string? outName))
            throw new ArgumentException($"No configuration found for variable {variable}");
        string outUnits = GetTargetUnits(variable);
        return (outName, outUnits);
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
        var result = UnitConverter.AnalyseConversion(inputUnits, targetUnits);

        List<string> operators = [];

        if (result.RequiresConversion)
        {
            string expression = UnitConverter.GenerateConversionExpression(
                inputUnits,
                targetUnits,
                timeStep);
            operators.Add(expression);
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

        var aggregationMethod = GetAggregationMethod(variable);
        var @operator = aggregationMethod.ToCdoOperator(_config.OutputTimeStep);

        return $"-{@operator},{stepsToAggregate}";
    }

    /// <summary>
    /// Write the equations for estimating VPD to a text writer.
    /// </summary>
    /// <param name="writer">The writer to write to.</param>
    /// <param name="method">The VPD estimation method to use.</param>
    /// <exception cref="ArgumentException">If the specified VPD method is not supported.</exception>
    internal async Task WriteVPDEquationsAsync(TextWriter writer, VPDMethod method)
    {
        // All methods follow the same general pattern:
        // 1. Calculate saturation vapor pressure (_esat)
        // 2. Calculate actual vapor pressure (_e)
        // 3. Calculate VPD as (_esat - _e) / 1000 to convert to kPa

        // Other assumptions:
        // - Temperature is assumed to be in ℃. (TODO: dynamic equations based on actual tas output units)
        // - Temperature variable name is assumed to be tas. fixme!
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

        await writer.WriteLineAsync($@"# Saturation vapor pressure (Pa) (tas in degC)");
        await writer.WriteLineAsync($"{esatEquation};");
        await writer.WriteLineAsync("# Actual vapor pressure (Pa)");
        await writer.WriteLineAsync("_e=(huss*ps)/(0.622+0.378*huss);");
        await writer.WriteLineAsync("# VPD (kPa)");
        await writer.WriteLineAsync("vpd=(_esat-_e)/1000;");
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
        string temperatureFile = pathManager.GetDatasetFileName(dataset, ClimateVariable.Temperature, PathType.Working);
        return GetVpdFilePath(dataset, temperatureFile);
    }

    /// <summary>
    /// Get the path to the optimised VPD output file for a dataset.
    /// </summary>
    /// <param name="dataset">The dataset.</param>
    /// <returns>The output file path.</returns>
    private string GetOptimisedVpdOutputFilePath(IClimateDataset dataset)
    {
        string temperatureFile = pathManager.GetDatasetFileName(dataset, ClimateVariable.Temperature, PathType.Output);
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
        // TODO: this could be simplified if VPD were a ClimateVariable.
        // Currently, these are the input variables though - so adding it would
        // require some refactoring, and would probably break the encapsulation
        // offered by this enum type, as most datasets *don't* have VPD as an
        // input variable.
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
    internal async Task<string> GenerateVPDScript(IClimateDataset dataset)
    {
        string jobName = $"calc_vpd_{dataset.DatasetName}";
        string script = await CreateScript(jobName);

        string humidityFile = pathManager.GetDatasetFileName(dataset, ClimateVariable.SpecificHumidity, PathType.Working);
        string pressureFile = pathManager.GetDatasetFileName(dataset, ClimateVariable.SurfacePressure, PathType.Working);
        string temperatureFile = pathManager.GetDatasetFileName(dataset, ClimateVariable.Temperature, PathType.Working);

        // Generate an output file name.
        string outFile = GetUnoptimisedVpdOutputFilePath(dataset);

        // Equation file is written to JobFS, so it will never require a storage
        // directive.
        string[] requiredFiles = [
            humidityFile,
            pressureFile,
            temperatureFile,
            outFile,
        ];
        IEnumerable<PBSStorageDirective> storageDirectives = PBSStorageHelper.GetStorageDirectives(requiredFiles);

        using TextWriter writer = new StreamWriter(script);
        await pbsLightweight.WritePBSHeader(writer, jobName, storageDirectives);

        await writer.WriteLineAsync("# File paths.");
        await writer.WriteLineAsync($"HUSS_FILE=\"{humidityFile}\"");

        await writer.WriteLineAsync($"PS_FILE=\"{pressureFile}\"");

        await writer.WriteLineAsync($"TAS_FILE=\"{temperatureFile}\"");

        string inFiles = "\"${HUSS_FILE}\" \"${PS_FILE}\" \"${TAS_FILE}\"";

        await writer.WriteLineAsync($"OUT_FILE=\"{outFile}\"");

        string eqnFile = "${WORK_DIR}/vpd_equations.txt";
        await writer.WriteLineAsync($"EQN_FILE=\"{eqnFile}\"");
        await writer.WriteLineAsync();

        await writer.WriteLineAsync("# Generate equation file.");
        await writer.WriteLineAsync("log \"Generating VPD equation file...\"");

        // Create equation file with selected method.
        await writer.WriteLineAsync($"cat >\"${{EQN_FILE}}\" <<EOF");
        await WriteVPDEquationsAsync(writer, _config.VPDMethod);
        await writer.WriteLineAsync("EOF");
        await writer.WriteLineAsync();

        // Calculate VPD using the equation file.
        await writer.WriteLineAsync("# Calculate VPD.");
        await writer.WriteLineAsync($"log \"Calculating VPD...\"");
        await writer.WriteLineAsync($"cdo {GetCDOArgs()} exprf,\"${{EQN_FILE}}\" -merge {inFiles} \"${{OUT_FILE}}\"");
        await writer.WriteLineAsync($"log \"VPD calculation completed successfully.\"");

        // We can't delete the intermediate files yet, because they are required
        // by the rechunk_X jobs, which may not have run yet.

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

    internal (string, AggregationMethod) GetTargetConfig(ClimateVariable variable)
    {
        var variables = GetVariables();
        if (!variables.TryGetValue(variable, out (string units, AggregationMethod _) config))
            throw new ArgumentException($"No configuration found for variable {variable}");
        return config;
    }

    /// <summary>
    /// Get the target units for the specified variable.
    /// </summary>
    /// <param name="variable">The variable.</param>
    /// <returns>The target units.</returns>
    /// <exception cref="ArgumentException">If no configuration is found for the specified variable.</exception>
    public string GetTargetUnits(ClimateVariable variable)
    {
        (string units, _) = GetTargetConfig(variable);
        return units;
    }

    /// <summary>
    /// Get the aggregation method required for the processing of the specified variable.
    /// </summary>
    /// <param name="variable">The variable.</param>
    /// <returns>The aggregation method.</returns>
    /// <exception cref="ArgumentException">If no configuration is found for the specified variable.</exception>
    public AggregationMethod GetAggregationMethod(ClimateVariable variable)
    {
        (string _, AggregationMethod aggregation) = GetTargetConfig(variable);
        return aggregation;
    }

    /// <summary>
    /// Generate processing scripts, and return the path to the top-level script.
    /// </summary>
    /// <param name="dataset">The dataset.</param>
    /// <returns>The path to the top-level script.</returns>
    public async Task<string> GenerateScriptsAsync(IClimateDataset dataset)
    {
        pathManager.CreateDirectoryTree(dataset);

        string submitJobName = $"submit_{dataset.DatasetName}";
        string submitScript = await CreateScript(submitJobName);
        using TextWriter writer = new StreamWriter(submitScript);

        // Add PBS header.
        await writer.WriteLineAsync("#!/usr/bin/env bash");
        await writer.WriteLineAsync($"# Job submission script for: {dataset.DatasetName}");
        await writer.WriteLineAsync();
        await WriteAutoGenerateHeader(writer);
        await writer.WriteLineAsync("# Exit immediately if any command fails.");
        await writer.WriteLineAsync("set -euo pipefail");
        await writer.WriteLineAsync();

        // Ensure output directory exists.
        await writer.WriteLineAsync($"mkdir -p \"{_config.OutputDirectory}\"");
        await writer.WriteLineAsync();

        // Process each variable.
        Dictionary<ClimateVariable, string> mergetimeScripts = new();
        Dictionary<ClimateVariable, string> rechunkScripts = new();
        ClimateVariable[] variables = GetVariables().Keys.ToArray();
        foreach (ClimateVariable variable in variables)
        {
            string mergetime = await GenerateVariableMergeScript(dataset, variable);
            string rechunk = await GenerateVariableRechunkScript(dataset, variable);

            mergetimeScripts[variable] = mergetime;
            rechunkScripts[variable] = rechunk;
        }

        string vpdScript = await GenerateVPDScript(dataset);
        string vpdRechunkScript = await GenerateVPDRechunkScript(dataset);
        string cleanupScript = await GenerateCleanupScript(dataset);

        // Add job submission logic.
        bool requiresVpd = _config.Version == ModelVersion.Dave;
        bool vpdEmpty = true;
        bool allDepsEmpty = true;

        foreach (ClimateVariable variable in variables)
        {
            // Submit mergetime script.
            await writer.WriteLineAsync($"JOB_ID=\"$(qsub \"{mergetimeScripts[variable]}\")\"");

            // Append this job to the list of VPD dependencies if necessary.
            if (requiresVpd && IsVpdDependency(variable))
            {
                if (vpdEmpty)
                    await writer.WriteLineAsync("VPD_DEPS=\"${JOB_ID}\"");
                else
                    await writer.WriteLineAsync($"VPD_DEPS=\"${{VPD_DEPS}}:${{JOB_ID}}\"");
                vpdEmpty = false;
            }

            bool variableRequired = true;
            if (variable == ClimateVariable.SpecificHumidity && _config.Version == ModelVersion.Dave)
                variableRequired = false;

            if (variableRequired)
            {
                // DAVE version doesn't require specific humidity. We need to do
                // the mergetime step, because that's used as an input for the
                // VPD computation, but the quantity itself is not used by the
                // model and we therefore don't need to rechunk it.

                // Submit rechunk script.
                await writer.WriteLineAsync($"JOB_ID=\"$(qsub -W depend=afterok:\"${{JOB_ID}}\" \"{rechunkScripts[variable]}\")\"");

                // Append the ID of the rechunk job to the "all jobs" list.
                if (allDepsEmpty)
                {
                    await writer.WriteLineAsync("ALL_JOBS=\"${JOB_ID}\"");
                    allDepsEmpty = false;
                }
                else
                    await writer.WriteLineAsync($"ALL_JOBS=\"${{ALL_JOBS}}:${{JOB_ID}}\"");
            }
            await writer.WriteLineAsync();
        }

        // Submit VPD scripts.
        if (requiresVpd)
        {
            await writer.WriteLineAsync($"JOB_ID=\"$(qsub -W depend=afterok:\"${{VPD_DEPS}}\" \"{vpdScript}\")\"");
            await writer.WriteLineAsync($"JOB_ID=\"$(qsub -W depend=afterok:\"${{JOB_ID}}\" \"{vpdRechunkScript}\")\"");
            await writer.WriteLineAsync($"ALL_JOBS=\"${{ALL_JOBS}}:${{JOB_ID}}\"");
            await writer.WriteLineAsync();
        }

        // Submit cleanup script.
        await writer.WriteLineAsync($"JOB_ID=\"$(qsub -W depend=afterok:\"${{ALL_JOBS}}\" \"{cleanupScript}\")\"");
        await writer.WriteLineAsync();

        await writer.WriteLineAsync($"echo \"Job submission complete for dataset {dataset.DatasetName}. Job ID: ${{JOB_ID}}\"");
        await writer.WriteLineAsync();

        return submitScript;
    }

    /// <summary>
    /// Get the climate variables required by the version of the model specified
    /// by the current configuration.
    /// </summary>
    /// <returns>The climate variables required by the model.</returns>
    /// <exception cref="ArgumentException">Thrown when an invalid version is specified.</exception>
    private IDictionary<ClimateVariable, (string units, AggregationMethod aggregation)> GetVariables()
    {
        if (_config.Version == ModelVersion.Dave)
            return daveVariables;
        if (_config.Version == ModelVersion.Trunk)
            return trunkVariables;
        throw new ArgumentException($"Invalid version: {_config.Version}");
    }

    /// <summary>
    /// Create the specified script file, and give it the required permissions.
    /// </summary>
    /// <param name="script">Path to a script file.</param>
    /// <exception cref="PlatformNotSupportedException">Thrown when the current operating system is not supported.</exception>
    private static async Task InitialiseScript(string script)
    {
        await File.WriteAllTextAsync(script, string.Empty);
        if (OperatingSystem.IsLinux() || OperatingSystem.IsMacOS())
            File.SetUnixFileMode(script, UnixFileMode.UserRead | UnixFileMode.UserWrite | UnixFileMode.UserExecute);
        else
            // Windows not supported.
            throw new PlatformNotSupportedException();
    }

    /// <summary>
    /// Create an empty script file and set execute permissions.
    /// </summary>
    /// <param name="jobName">The name of the job.</param>
    /// <returns>The path to the script file.</returns>
    private async Task<string> CreateScript(string jobName)
    {
        string directory = pathManager.GetBasePath(PathType.Script);

        // Use job name as script file name. This avoids the problems associated
        // with naming a file with a ".sh" (or similar) extension.
        string script = Path.Combine(directory, jobName);

        // Create empty script file, set execute permissions.
        await InitialiseScript(script);

        return script;
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
    /// Sanitise a string to be stored in a bash variable.
    /// </summary>
    /// <param name="input">The string.</param>
    /// <returns>The sanitised string.</returns>
    internal static string SanitiseString(string input)
    {
        return input.Replace("$", "\\$");
    }

    /// <summary>
    /// Generate a mergetime script for the specified variable, and return the
    /// path to the generated script file.
    /// </summary>
    /// <param name="dataset">The dataset to process.</param>
    /// <param name="variable">The variable to process.</param>
    /// <param name="outFile">The path to the output file that should be generated by the script.</param>
    /// <returns>The path to the generated script file.</returns>
    internal async Task<string> GenerateVariableMergeScript(IClimateDataset dataset, ClimateVariable variable)
    {
        VariableInfo varInfo = dataset.GetVariableInfo(variable);
        (string outVar, string targetUnits) = GetStandardConfig(variable);

        // File paths.
        string inDir = dataset.GetInputFilesDirectory(variable);

        string outFile = pathManager.GetDatasetFileName(dataset, variable, PathType.Working);

        // Sanitise - e.g. /tmp/./x -> /tmp/x
        outFile = Path.GetFullPath(outFile);

        List<string> requiredFiles = [
            inDir,
            outFile
        ];
        if (!string.IsNullOrEmpty(_config.GridFile))
            requiredFiles.Add(_config.GridFile);
        IEnumerable<PBSStorageDirective> storageDirectives =
            PBSStorageHelper.GetStorageDirectives(requiredFiles);

        // Create script directory if it doesn't already exist.
        // This should be unnecessary at this point.
        string jobName = GetJobName("mergetime", varInfo, dataset);
        string scriptFile = await CreateScript(jobName);
        using TextWriter writer = new StreamWriter(scriptFile);
        await pbsLightweight.WritePBSHeader(writer, jobName, storageDirectives);

        await writer.WriteLineAsync("# File paths.");
        await writer.WriteLineAsync($"{inDirVariable}=\"{SanitiseString(inDir)}\"");
        if (!string.IsNullOrEmpty(_config.GridFile))
            await writer.WriteLineAsync($"{remapDirVariable}=\"${{WORK_DIR}}/remap\"");
        await writer.WriteLineAsync($"OUT_FILE=\"{SanitiseString(outFile)}\"");
        if (!string.IsNullOrEmpty(_config.GridFile))
            await writer.WriteLineAsync($"GRID_FILE=\"{SanitiseString(_config.GridFile)}\"");
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
    /// Generate a cleanup script for this job.
    /// </summary>
    /// <param name="dataset">The dataset to process.</param>
    /// <returns>The path to the script file.</returns>
    private async Task<string> GenerateCleanupScript(IClimateDataset dataset)
    {
        string jobName = $"cleanup_{dataset.DatasetName}";
        string scriptFile = await CreateScript(jobName);
        string workDir = pathManager.GetDatasetPath(dataset, PathType.Working);
        IEnumerable<PBSStorageDirective> storageDirectives =
            PBSStorageHelper.GetStorageDirectives([workDir]);

        using TextWriter writer = new StreamWriter(scriptFile);

        await pbsLightweight.WritePBSHeader(writer, jobName, storageDirectives);
        await writer.WriteLineAsync("# File paths.");
        await writer.WriteLineAsync($"IN_DIR=\"{workDir}\"");
        await writer.WriteLineAsync("rm -rf \"${IN_DIR}\"");

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
    internal async Task<string> GenerateVariableRechunkScript(IClimateDataset dataset, ClimateVariable variable)
    {
        VariableInfo varInfo = dataset.GetVariableInfo(variable);

        string jobName = GetJobName("rechunk", varInfo, dataset);
        string inFile = pathManager.GetDatasetFileName(dataset, variable, PathType.Working);
        string outFile = pathManager.GetDatasetFileName(dataset, variable, PathType.Output);
        bool cleanup = !IsVpdDependency(variable);

        return await GenerateRechunkScript(jobName, inFile, outFile, cleanup);
    }

    private async Task<string> GenerateRechunkScript(string jobName, string inFile, string outFile, bool cleanup)
    {
        string scriptFile = await CreateScript(jobName);
        IEnumerable<PBSStorageDirective> storageDirectives =
            PBSStorageHelper.GetStorageDirectives([inFile, outFile]);

        using TextWriter writer = new StreamWriter(scriptFile);

        await pbsHeavyweight.WritePBSHeader(writer, jobName, storageDirectives);

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

        // Calculate checksum.
        // Note: we change directory and use a relative file path, to ensure
        // that the checksum file remains portable.

        string outputPath = pathManager.GetBasePath(PathType.Output);
        string checksumFile = pathManager.GetChecksumFilePath();
        string relativePath = Path.GetRelativePath(outputPath, outFile);

        await writer.WriteLineAsync("# Calculate checksum.");
        await writer.WriteLineAsync($"log \"Calculating checksum...\"");
        await writer.WriteLineAsync($"cd \"{outputPath}\"");
        await writer.WriteLineAsync($"REL_PATH=\"{relativePath}\"");
        await writer.WriteLineAsync($"sha512sum \"${{REL_PATH}}\" >>\"{checksumFile}\"");
        await writer.WriteLineAsync("log \"Checksum calculation completed successfully.\"");
        await writer.WriteLineAsync();

        // We can now delete the temporary input file, but only if it's not also
        // required for the VPD estimation, which may not have occurred yet.
        if (cleanup)
        {
            await writer.WriteLineAsync("# Delete temporary file.");
            await writer.WriteLineAsync($"rm -f \"${{IN_FILE}}\"");
            await writer.WriteLineAsync();
        }
        else
            await writer.WriteLineAsync("# Input file cannot (necessarily) be deleted yet, since it is required for VPD estimation.");

        return scriptFile;
    }

    /// <summary>
    /// Write a comment to a script which indicates that it was automatically
    /// generated.
    /// </summary>
    /// <param name="writer">The text writer to which the comment will be written.</param>
    private static async Task WriteAutoGenerateHeader(TextWriter writer)
    {
        await writer.WriteLineAsync("# This script was automatically generated. Do not modify.");
        await writer.WriteLineAsync();
    }

    /// <summary>
    /// Generate a wrapper script that executes the given scripts.
    /// </summary>
    /// <param name="scripts">The script files to execute.</param>
    public static async Task<string> GenerateWrapperScript(string outputDirectory, IEnumerable<string> scripts)
    {
        PathManager pathManager = new(outputDirectory);

        string jobName = "wrapper";
        string scriptPath = pathManager.GetBasePath(PathType.Script);
        string scriptFile = Path.Combine(scriptPath, jobName);

        await InitialiseScript(scriptFile);
        using TextWriter writer = new StreamWriter(scriptFile);

        await writer.WriteLineAsync("#!/usr/bin/env bash");
        await writer.WriteLineAsync("# Master-level script which executes all job submission scripts to submit all PBS jobs.");
        await writer.WriteLineAsync();

        await WriteAutoGenerateHeader(writer);

        await writer.WriteLineAsync("set -euo pipefail");
        await writer.WriteLineAsync();

        // Execute all scripts (making assumptions about file permissions).
        foreach (string script in scripts)
            await writer.WriteLineAsync(script);
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
