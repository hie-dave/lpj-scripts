using System.Text;
using ClimateProcessing.Models;
using System.IO;
using ClimateProcessing.Units;

namespace ClimateProcessing.Services;

public class ScriptGenerator
{
    private readonly ProcessingConfig _config;
    private static readonly Dictionary<ClimateVariable, (string stdName, string units, TimeStep timeStep)> _standardVariables = new()
    {
        { ClimateVariable.SpecificHumidity, ("huss", "1", TimeStep.Hourly) },
        { ClimateVariable.SurfacePressure, ("ps", "Pa", TimeStep.Hourly) },
        { ClimateVariable.ShortwaveRadiation, ("rsds", "W m-2", TimeStep.Hourly) },
        { ClimateVariable.WindSpeed, ("sfcWind", "m s-1", TimeStep.Hourly) },
        { ClimateVariable.Temperature, ("tas", "degC", TimeStep.Hourly) },
        { ClimateVariable.Precipitation, ("pr", "mm", TimeStep.Daily) }  // Target is accumulated daily precipitation
    };

    public ScriptGenerator(ProcessingConfig config)
    {
        _config = config;
    }

    internal ProcessingCommand GenerateUnitConversionCommand(
        string inputVar,
        string outputVar,
        string inputUnits,
        string targetUnits,
        TimeStep timeStep,
        string inputFile,
        string outputFile)
    {
        var result = UnitConverter.AnalyzeConversion(inputUnits, targetUnits);
        
        if (!result.RequiresConversion && !result.RequiresRenaming)
        {
            return ProcessingCommand.Skip(inputFile, outputFile);
        }

        var commands = new List<string>();

        if (result.RequiresConversion)
        {
            var expression = UnitConverter.GenerateConversionExpression(
                inputVar, 
                outputVar, 
                inputUnits, 
                targetUnits,
                result.RequiresTimeStep ? timeStep : null);
            commands.Add($"cdo -O expr,'{expression}' {inputFile} {outputFile}");
        }

        if (result.RequiresRenaming && !result.RequiresConversion)
        {
            // If we only need to rename units, we can do it in-place with a symlink
            commands.Add($"ln -sf {inputFile} {outputFile}");
            commands.Add($"ncatted -O -a units,{outputVar},m,c,\"{targetUnits}\" {outputFile}");
        }
        else if (result.RequiresRenaming)
        {
            // If we already converted, just update the metadata on the output file
            commands.Add($"ncatted -O -a units,{outputVar},m,c,\"{targetUnits}\" {outputFile}");
        }

        return ProcessingCommand.Process(string.Join("\n", commands));
    }

    internal ProcessingCommand GenerateTimeAggregationCommand(
        ClimateVariable variable,
        string inputFile,
        string outputFile)
    {
        // Only aggregate if input and output timesteps differ
        if (_config.InputTimeStep == _config.OutputTimeStep)
        {
            return ProcessingCommand.Skip(inputFile, outputFile);
        }

        // Calculate the number of timesteps to aggregate
        int stepsToAggregate = _config.OutputTimeStep.Hours / _config.InputTimeStep.Hours;
        
        var aggregationMethod = _config.GetAggregationMethod(variable);
        var @operator = aggregationMethod.ToCdoOperator();

        return ProcessingCommand.Process(
            $"cdo -O -{@operator},{stepsToAggregate} {inputFile} {outputFile}");
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

    private string GenerateVPDCalculation()
    {
        var sb = new StringBuilder();
        sb.AppendLine("# Calculate VPD using equation file");
        
        // Create equation file with selected method
        sb.AppendLine("cat > \"$TMPDIR/vpd_equations.txt\" << 'EOF'");
        sb.AppendLine(GetVPDEquations(_config.VPDMethod));
        sb.AppendLine("EOF");
        sb.AppendLine();
        
        // Calculate VPD using the equation file
        sb.AppendLine("cdo -O exprf,\"$TMPDIR/vpd_equations.txt\" \"$TMPDIR/merged.nc\" \"$TMPDIR/vpd.nc\"");
        sb.AppendLine();
        
        // Clean up equation file
        sb.AppendLine("rm \"$TMPDIR/vpd_equations.txt\"");
        
        return sb.ToString();
    }

    private string GeneratePBSHeader()
    {
        var lines = new List<string>
        {
            "#!/bin/bash",
            $"#PBS -N {_config.JobName}",
            $"#PBS -P {_config.Project}",
            $"#PBS -q {_config.Queue}",
            $"#PBS -l walltime={_config.Walltime}",
            $"#PBS -l ncpus={_config.Ncpus}",
            $"#PBS -l mem={_config.Memory}GB",
            $"#PBS -j oe"
        };

        // Add storage directives if required
        var storageDirectives = _config.GetRequiredStorageDirectives();
        if (storageDirectives.Any())
        {
            lines.Add(PBSStorageHelper.FormatStorageDirectives(storageDirectives));
        }

        lines.Add("");  // Add blank line after header
        return string.Join("\n", lines);
    }

    public string GenerateProcessingScript(IClimateDataset dataset)
    {
        var sb = new StringBuilder();
        
        // Add PBS header
        sb.AppendLine(GeneratePBSHeader());

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
            string targetUnits = _config.GetTargetUnits(variable);

            sb.AppendLine($"# Processing {variable}");
            
            // Merge files and extract variable
            sb.AppendLine($"cdo -O mergetime {string.Join(" ", dataset.GetInputFiles())} merged.nc");
            sb.AppendLine($"cdo -O select,name={varInfo.Name} merged.nc temp.nc");
            sb.AppendLine("rm -f merged.nc");
            sb.AppendLine("START_DATE=$(cdo -s showdate temp.nc | head -n 1)");
            sb.AppendLine("END_DATE=$(cdo -s showdate temp.nc | tail -n 1)");
            sb.AppendLine();

            // Convert units if needed
            var conversionCmd = GenerateUnitConversionCommand(
                varInfo.Name, 
                varInfo.Name,
                varInfo.Units, 
                targetUnits,
                _config.InputTimeStep,
                "temp.nc", 
                "converted.nc");

            if (conversionCmd.RequiresProcessing)
            {
                sb.AppendLine(conversionCmd.Command);
                sb.AppendLine("rm -f temp.nc");
            }
            else
            {
                sb.AppendLine(conversionCmd.Command);
            }

            // Perform time aggregation if needed
            var aggregationCmd = GenerateTimeAggregationCommand(
                variable,
                "converted.nc",
                "aggregated.nc");

            if (aggregationCmd.RequiresProcessing)
            {
                sb.AppendLine(aggregationCmd.Command);
                sb.AppendLine("rm -f converted.nc");
            }
            else
            {
                sb.AppendLine(aggregationCmd.Command);
            }

            // Generate output filename with full date range
            sb.AppendLine($"FILENAME_PATTERN=\"{dataset.GetOutputFilePattern(variable)}\"");
            sb.AppendLine("OUTPUT_FILE=$(echo \"$FILENAME_PATTERN\" | sed \"s/[0-9]\\{8\\}-[0-9]\\{8\\}/${START_DATE}-${END_DATE}/\")");
            sb.AppendLine();

            // Reorder dimensions, improve chunking, and enable compression
            string chunkSpec = $"--cnk_dmn lat/{_config.ChunkSizeSpatial},lon/{_config.ChunkSizeSpatial},time/{_config.ChunkSizeTime}";
            string compressionSpec = _config.CompressOutput ? $"-L{_config.CompressionLevel}" : "";
            var outputPath = Path.Combine(_config.OutputDirectory, "$OUTPUT_FILE");

            sb.AppendLine($"ncpdq -O -a lat,lon,time {chunkSpec} {compressionSpec} aggregated.nc \"{outputPath}\"");
            sb.AppendLine("rm -f aggregated.nc");
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
