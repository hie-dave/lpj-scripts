using Xunit;
using ClimateProcessing.Services;
using ClimateProcessing.Models;
using ClimateProcessing.Units;
using System.Text.RegularExpressions;
using Microsoft.VisualStudio.TestPlatform.ObjectModel.Engine.ClientProtocol;
using ClimateProcessing.Tests.Mocks;
using System.Reflection;
using ClimateProcessing.Configuration;
using Xunit.Abstractions;

namespace ClimateProcessing.Tests.Services;

public class ScriptGeneratorTests : IDisposable
{
    private const string outputDirectoryPrefix = "script_generator_tests_output";
    private readonly ITestOutputHelper outputHelper;
    private readonly string outputDirectory;
    private static readonly string[] cdoTemporalAggregationOperators = [
        "daymin",
        "daymax",
        "daysum",
        "daymean",
        "dayrange",
        "dayavg",
        "daystd",
        "daystd1",
        "dayvar",
        "dayvar1",
        "timselmin",
        "timselmax",
        "timselsum",
        "timselmean",
        "timselrange",
        "timselavg",
        "timselstd",
        "timselstd1",
        "timselvar",
        "timselvar1"
    ];

    private static readonly string[] cdoArithmeticOperators = [
        "addc",
        "subc",
        "mulc",
        "divc",
        "minc",
        "maxc",
        "expr",
    ];

    private readonly NarClim2Config _config;

    private readonly ScriptGenerator _generator;

    public ScriptGeneratorTests(ITestOutputHelper outputHelper)
    {
        this.outputHelper = outputHelper;
        outputDirectory = CreateOutputDirectory();

        _config = new NarClim2Config()
        {
            Project = "test",
            Queue = "normal",
            Walltime = "01:00:00",
            Ncpus = 1,
            Memory = 4,
            OutputDirectory = outputDirectory,
            InputTimeStepHours = 1,
            OutputTimeStepHours = 24
        };

        _generator = new ScriptGenerator(_config);
    }

    /// <summary>
    /// Teardown method - delete the temporary output directory.
    /// </summary>
    public void Dispose()
    {
        try
        {
            Directory.Delete(outputDirectory, true);
        }
        catch (IOException error)
        {
            // Log an error but don't throw an exception.
            Console.Error.WriteLine($"Warning: could not delete temporary output directory used by {GetType().Name}: {error}");
        }
    }

    private static string CreateOutputDirectory()
    {
        return Directory.CreateTempSubdirectory(outputDirectoryPrefix).FullName;
    }

    [Theory]
    [InlineData("K", "K", "temp", false, false)]  // No conversion needed
    [InlineData("K", "degC", "temp", true, true)]  // Conversion needed
    [InlineData("W/m2", "W m-2", "rad", false, true)]  // Only renaming needed
    public void GenerateUnitConversionOperators_GeneratesCorrectOperators(
        string inputUnits,
        string targetUnits,
        string outputVar,
        bool expectsConversion,
        bool expectsRenaming)
    {
        var operators = _generator.GenerateUnitConversionOperators(
            outputVar,
            inputUnits,
            targetUnits,
            TimeStep.Hourly).ToList();

        if (expectsConversion)
        {
            Assert.Contains(operators, op => op.StartsWith("-subc"));
        }
        if (expectsRenaming)
        {
            Assert.Contains(operators, op => op.StartsWith("-setattribute"));
        }
        if (!expectsConversion && !expectsRenaming)
        {
            Assert.Empty(operators);
        }
    }

    [Theory]
    [InlineData(1, 24, true)]   // Hourly to daily
    [InlineData(24, 24, false)] // Daily to daily
    [InlineData(3, 24, true)]   // 3-hourly to daily
    public void GenerateTimeAggregationOperators_GeneratesCorrectOperators(
        int inputHours,
        int outputHours,
        bool requiresAggregation)
    {
        var config = new NarClim2Config
        {
            Project = _config.Project,
            Queue = _config.Queue,
            Walltime = _config.Walltime,
            Ncpus = _config.Ncpus,
            Memory = _config.Memory,
            OutputDirectory = _config.OutputDirectory,
            InputTimeStepHours = inputHours,
            OutputTimeStepHours = outputHours
        };
        var generator = new ScriptGenerator(config);

        string @operator = generator.GenerateTimeAggregationOperator(
            ClimateVariable.Temperature);

        if (requiresAggregation)
        {
            if (outputHours == 24)
                Assert.StartsWith("-daymean", @operator);
            else
                Assert.StartsWith("-timselmean", @operator);
        }
        else
        {
            Assert.Empty(@operator);
        }
    }

    [Theory]
    [InlineData("tas", "tas", false)]      // Same name
    [InlineData("temp", "tas", true)]      // Different name
    [InlineData("pr", "prec", true)]       // Different name
    public void GenerateVariableRenameCommand_DetectsProcessingNeeds(
        string inputVar,
        string outputVar,
        bool requiresProcessing)
    {
        string op = _generator.GenerateRenameOperator(inputVar, outputVar);
        Assert.Equal(requiresProcessing, !string.IsNullOrEmpty(op));
    }

    [Theory]
    [InlineData("kg m-2 s-1", true)]  // Standard notation
    [InlineData("kg/m2/s", true)]     // Division notation
    [InlineData("kg.m-2.s-1", true)]  // Dot notation
    [InlineData("kg m^-2 s^-1", true)] // Caret notation
    [InlineData("KG M-2 S-1", true)]   // Case insensitive
    [InlineData("kg  m-2  s-1", true)] // Extra spaces
    [InlineData("kgm-2s-1", true)]     // No separators
    [InlineData("W", false)]           // No per-area units
    [InlineData("kg s-1", false)]      // Time only, no area
    [InlineData("", false)]            // Empty string
    [InlineData("m2", false)]          // Area but not per-area
    [InlineData("kg/s/m2", true)]      // Different order
    public void HasPerAreaUnits_DetectsUnitsCorrectly(string units, bool expectedResult)
    {
        // Now we can test HasPerAreaUnits directly since it's internal
        var result = ScriptGenerator.HasPerAreaUnits(units);
        Assert.Equal(expectedResult, result);
    }

    [Theory]
    [InlineData(ClimateVariable.Temperature)]
    [InlineData(ClimateVariable.SpecificHumidity)]
    public void GetRemapAlgorithm_NonAreaSensitiveVariables_AlwaysReturnsBilinear(ClimateVariable variable)
    {
        var info = new VariableInfo(Enum.GetName(variable)!, "any_units");
        var result = _generator.GetInterpolationAlgorithm(info, variable);

        Assert.Equal(InterpolationAlgorithm.Bilinear, result);
    }

    [Theory]
    [InlineData(ClimateVariable.Precipitation)]
    [InlineData(ClimateVariable.ShortwaveRadiation)]
    public void GetRemapAlgorithm_AreaSensitiveVariables_DependsOnUnits(ClimateVariable variable)
    {
        // Should use conservative remapping when not per-area
        var nonAreaInfo = new VariableInfo(Enum.GetName(variable)!, "W");
        var nonAreaResult = _generator.GetInterpolationAlgorithm(nonAreaInfo, variable);
        Assert.Equal(InterpolationAlgorithm.Conservative, nonAreaResult);

        // Should use bilinear remapping when per-area
        var areaInfo = new VariableInfo(Enum.GetName(variable)!, "W m-2");
        var areaResult = _generator.GetInterpolationAlgorithm(areaInfo, variable);
        Assert.Equal(InterpolationAlgorithm.Bilinear, areaResult);
    }

    [Theory]
    [InlineData("input dir with spaces")]
    [InlineData("/path/with/special/$chars")]
    [InlineData("/normal/path")]
    public async Task GenerateVariableMergeScript_QuotesVariablesSafely(
        string inputDir)
    {
        NarClim2Config config = new()
        {
            Project = "test",
            Queue = "normal",
            Walltime = "01:00:00",
            Ncpus = 1,
            Memory = 4,
            OutputDirectory = outputDirectory,
            InputDirectory = inputDir,
            InputTimeStepHours = 3,
            OutputTimeStepHours = 24,
        };
        ScriptGenerator generator = new(config);
        StaticMockDataset dataset = new(inputDir);

        string scriptPath = await generator.GenerateVariableMergeScript(
            dataset,
            ClimateVariable.Temperature);
        string scriptContent = await File.ReadAllTextAsync(scriptPath);

        // No unquoted variable references
        ValidateScript(scriptContent);
    }

    private void ValidateScript(string scriptContent)
    {
        // Every variable reference uses braces.
        Assert.DoesNotMatch(@"[^\\]\$[A-Za-z]", scriptContent);

        // No double-escaped braces from string interpolation
        Assert.DoesNotMatch(@"[^\\]\$\{\{", scriptContent);
        Assert.DoesNotMatch(@"[^\\]\$\{\{?[^\}]+\}\}", scriptContent);

        // We can't assume that all variable referenes are quoted, because
        // sometimes they don't need to be or shouldn't be.
    }

    [Theory]
    [InlineData(ClimateVariable.Temperature, 1, 24, "tas", "K", "-daymean,24", "-subc,273.15", "-setattribute,'tas@units=degC")]
    [InlineData(ClimateVariable.Precipitation, 1, 24, "pr", "kg m-2 s-1", "-daysum,24", "-mulc,3600", "-setattribute,'pr@units=mm'")]
    [InlineData(ClimateVariable.Precipitation, 1, 8, "pr", "kg m-2 s-1", "-timselsum,8", "-mulc,3600", "-setattribute,'pr@units=mm'")]
    [InlineData(ClimateVariable.SpecificHumidity, 1, 1, "huss", "1")] // No unit conversion or aggregation
    [InlineData(ClimateVariable.SpecificHumidity, 1, 1, "huss", "kg/kg", null, null, "-setattribute,'huss@units=1'")] // Unit rename, but no unit conversion or aggregation
    [InlineData(ClimateVariable.ShortwaveRadiation, 1, 3, "rsds", "W m-2", "-timselmean,3")] // Aggregation but no unit conversion (intensive variable)
    [InlineData(ClimateVariable.Precipitation, 1, 12, "pr", "mm", "-timselsum,12")] // Aggregation but no unit conversion (extensive variable)
    [InlineData(ClimateVariable.SurfacePressure, 1, 1, "ps", "kPa", null, "-mulc,1000", "-setattribute,'ps@units=Pa'")] // Unit conversion but no aggregation
    public async Task GenerateVariableMergeScript_GeneratesValidCDORemapCommand(
        ClimateVariable variable,
        int inputTimestepHours,
        int outputTimestepHours,
        string varName,
        string inputUnits,
        string? temporalAggregationOperator = null,
        string? unitConversionOperator = null,
        string? unitRenameOperator = null)
    {
        NarClim2Config config = new()
        {
            Project = "test",
            Queue = "normal",
            Walltime = "01:00:00",
            Ncpus = 1,
            Memory = 4,
            OutputDirectory = outputDirectory,
            InputDirectory = "/input",
            InputTimeStepHours = inputTimestepHours,
            OutputTimeStepHours = outputTimestepHours,
            Version = ModelVersion.Dave
        };
        ScriptGenerator generator = new(config);
        StaticMockDataset dataset = new("/input", varName, inputUnits);

        string scriptPath = await generator.GenerateVariableMergeScript(
            dataset,
            variable);
        string scriptContent = await File.ReadAllTextAsync(scriptPath);

        // E.g.
        // cdo -L -O -v -z zip1 -daymean,24 -subc,273.15 -setattribute,'tas@units=degC' -unpack  \"${FILE}\" \"${REMAP_DIR}/$(basename \"${FILE}\")\"
        // "    cdo -L -O -v -z zip1 -daysum,24 -mulc,3600 -setattribute,'pr@units=mm' -unpack  \"${FILE}\" \"${REMAP_DIR}/$(basename \"${FILE}\")\""
        string line = scriptContent.Split("\n").First(l => l.Contains("cdo -"));

        // CDO command should have proper structure

        // Should use -L for thread-safety.
        Assert.Contains("-L", line);

        // Should use -v for progress tracking.
        Assert.Contains("-v", line);

        // Should use -O to overwrite any existing output file.
        Assert.Contains("-O", line);

        // Should use -z zip1 for efficiency.
        Assert.Contains("-z zip1", line);

        // Should unpack data.
        Assert.Contains("-unpack", line);

        // Should apply operators (or not).
        int previousOperatorIndex = -1;
        if (temporalAggregationOperator is null)
            foreach (string @operator in cdoTemporalAggregationOperators)
                Assert.DoesNotContain(@operator, line);
        else
        {
            Assert.Contains(temporalAggregationOperator, line);
            int operatorIndex = line.IndexOf(temporalAggregationOperator);
            Assert.True(operatorIndex > previousOperatorIndex, $"Operator '{temporalAggregationOperator}' appears out of order in CDO command.");
            previousOperatorIndex = operatorIndex;
        }

        if (unitConversionOperator is null)
            foreach (string @operator in cdoArithmeticOperators)
                Assert.DoesNotContain(@operator, line);
        else
        {
            Assert.Contains(unitConversionOperator, line);
            int operatorIndex = line.IndexOf(unitConversionOperator);
            Assert.True(operatorIndex > previousOperatorIndex, $"Operator '{unitConversionOperator}' appears out of order in CDO command.");
            previousOperatorIndex = operatorIndex;
        }

        if (unitRenameOperator is null)
            Assert.DoesNotContain("setattribute", line);
        else
        {
            Assert.Contains(unitRenameOperator, line);
            int operatorIndex = line.IndexOf(unitRenameOperator);
            Assert.True(operatorIndex > previousOperatorIndex, $"Operator '{unitRenameOperator}' appears out of order in CDO command.");
            previousOperatorIndex = operatorIndex;
        }

        // TODO: assert that no additional arguments are present.

        // Rest of script should be valid.
        ValidateScript(scriptContent);
    }

    [Theory]
    [InlineData(ClimateVariable.Precipitation, "precip", "pr", "mm")]
    public async Task GenerateVariableMergeScript_GeneratesValidRenameCommand(
        ClimateVariable variable,
        string inputName,
        string expectedOutputName,
        string inputUnits)
    {
        ScriptGenerator generator = new(_config);
        StaticMockDataset dataset = new("/input", inputName, inputUnits);

        string scriptPath = await generator.GenerateVariableMergeScript(
            dataset,
            variable);
        string scriptContent = await File.ReadAllTextAsync(scriptPath);

        Assert.Contains($"-chname,'{inputName}','{expectedOutputName}'", scriptContent);
    }

    [Theory]
    [InlineData(VPDMethod.Magnus)]
    [InlineData(VPDMethod.Buck1981)]
    [InlineData(VPDMethod.AlduchovEskridge1996)]
    [InlineData(VPDMethod.AllenFAO1998)]
    [InlineData(VPDMethod.Sonntag1990)]
    public async Task GenerateVPDScript_GeneratesCorrectEquations(VPDMethod method)
    {
        NarClim2Config config = new()
        {
            Project = "test",
            Queue = "normal",
            Walltime = "01:00:00",
            Ncpus = 1,
            Memory = 4,
            OutputDirectory = outputDirectory,
            InputDirectory = "/input",
            VPDMethod = method
        };
        ScriptGenerator generator = new(config);

        StringWriter writer = new();
        await generator.WriteVPDEquationsAsync(writer, method);
        string equationContent = writer.ToString();

        // Remove comment lines.
        string sanitised = Regex.Replace(equationContent, @"#.*\n", "");

        // Ensure all lines end with a semicolon.
        Assert.DoesNotMatch(@"[^;\r]\n", sanitised);

        Assert.Matches(@"_e=.*;\n", sanitised);
        Assert.Matches(@"_esat=.*;\n", sanitised);
        Assert.Matches(@"vpd=.*;\n", sanitised);
    }

    [Theory]
    [InlineData(VPDMethod.Magnus)]
    [InlineData(VPDMethod.Buck1981)]
    [InlineData(VPDMethod.AlduchovEskridge1996)]
    [InlineData(VPDMethod.AllenFAO1998)]
    [InlineData(VPDMethod.Sonntag1990)]
    public async Task GenerateVPDScript_GeneratesValidScript(VPDMethod method)
    {
        NarClim2Config config = new()
        {
            Project = "test",
            Queue = "normal",
            Walltime = "01:00:00",
            Ncpus = 1,
            Memory = 4,
            OutputDirectory = outputDirectory,
            InputDirectory = "/input",
            VPDMethod = method
        };
        ScriptGenerator generator = new(config);
        StaticMockDataset dataset = new("/input");

        string scriptPath = await generator.GenerateVPDScript(dataset);

        Assert.True(File.Exists(scriptPath));
        string scriptContent = await File.ReadAllTextAsync(scriptPath);
        ValidateScript(scriptContent);
    }

    [Theory]
    [InlineData(true)]   // With VPD calculation
    [InlineData(false)]  // Without VPD calculation
    public async Task GenerateScriptsAsync_HandlesVPDDependenciesCorrectly(
        bool requiresVPD)
    {
        NarClim2Config config = new()
        {
            Project = "test",
            Queue = "normal",
            Walltime = "01:00:00",
            Ncpus = 1,
            Memory = 4,
            OutputDirectory = outputDirectory,
            InputDirectory = "/input",
            Version = requiresVPD ? ModelVersion.Dave : ModelVersion.Trunk,
            InputTimeStepHours = 1,
            OutputTimeStepHours = 1
        };
        ScriptGenerator generator = new(config);
        DynamicMockDataset dataset = new(config.InputDirectory, config.OutputDirectory);

        string scriptPath = await generator.GenerateScriptsAsync(dataset);

        Assert.True(File.Exists(scriptPath));
        string scriptContent = await File.ReadAllTextAsync(scriptPath);

        // Basic script validation
        ValidateScript(scriptContent);

        if (requiresVPD)
        {
            // Should have all required VPD variables
            Assert.Contains("huss", scriptContent);
            Assert.Contains("ps", scriptContent);
            Assert.Contains("tas", scriptContent);
            Assert.Contains("vpd", scriptContent);

            // Should have proper dependencies
            Assert.Contains("afterok:", scriptContent);
        }
        else
        {
            // Should not have VPD-related content
            Assert.DoesNotContain("vpd", scriptContent);
        }

        // Should always have cleanup job
        Assert.Contains("cleanup_", scriptContent);
    }

    [Theory]
    [InlineData(ClimateVariable.Temperature, "K", ModelVersion.Trunk)]
    [InlineData(ClimateVariable.Temperature, "degC", ModelVersion.Dave)]
    [InlineData(ClimateVariable.MaxTemperature, "K", ModelVersion.Trunk)]
    [InlineData(ClimateVariable.MinTemperature, "K", ModelVersion.Trunk)]
    [InlineData(ClimateVariable.Precipitation, "mm", ModelVersion.Trunk)]
    [InlineData(ClimateVariable.Precipitation, "mm", ModelVersion.Dave)]
    [InlineData(ClimateVariable.ShortwaveRadiation, "W m-2", ModelVersion.Trunk)]
    [InlineData(ClimateVariable.ShortwaveRadiation, "W m-2", ModelVersion.Dave)]
    [InlineData(ClimateVariable.SpecificHumidity, "1", ModelVersion.Trunk)]
    [InlineData(ClimateVariable.SpecificHumidity, "1", ModelVersion.Dave)]
    [InlineData(ClimateVariable.SurfacePressure, "Pa", ModelVersion.Trunk)]
    [InlineData(ClimateVariable.SurfacePressure, "Pa", ModelVersion.Dave)]
    [InlineData(ClimateVariable.WindSpeed, "m s-1", ModelVersion.Trunk)]
    [InlineData(ClimateVariable.WindSpeed, "m s-1", ModelVersion.Dave)]
    public void TestGetTargetUnits_ValidVariable(
        ClimateVariable variable,
        string expectedUnits,
        ModelVersion version)
    {
        _config.Version = version;
        ScriptGenerator generator = new ScriptGenerator(_config);
        string actualUnits = generator.GetTargetUnits(variable);
        Assert.Equal(expectedUnits, actualUnits);
    }

    [Theory]
    [InlineData(ClimateVariable.Temperature, AggregationMethod.Mean)]
    [InlineData(ClimateVariable.MaxTemperature, AggregationMethod.Maximum)]
    [InlineData(ClimateVariable.MinTemperature, AggregationMethod.Minimum)]
    [InlineData(ClimateVariable.Precipitation, AggregationMethod.Sum)]
    [InlineData(ClimateVariable.ShortwaveRadiation, AggregationMethod.Mean)]
    [InlineData(ClimateVariable.SpecificHumidity, AggregationMethod.Mean)]
    [InlineData(ClimateVariable.SurfacePressure, AggregationMethod.Mean)]
    [InlineData(ClimateVariable.WindSpeed, AggregationMethod.Mean)]
    public void TestGetTargetAggregationMethod_ValidVariable(
        ClimateVariable variable,
        AggregationMethod expectedMethod)
    {
        _config.Version = ModelVersion.Trunk;
        ScriptGenerator generator = new ScriptGenerator(_config);
        AggregationMethod actualMethod = generator.GetAggregationMethod(variable);
        Assert.Equal(expectedMethod, actualMethod);

        // Test DAVE version too.
        if (variable != ClimateVariable.MinTemperature && variable != ClimateVariable.MaxTemperature)
        {
            _config.Version = ModelVersion.Dave;
            generator = new ScriptGenerator(_config);
            actualMethod = generator.GetAggregationMethod(variable);
            Assert.Equal(expectedMethod, actualMethod);
        }
    }

    [Theory]
    [InlineData(ClimateVariable.MaxTemperature, ModelVersion.Dave)]
    [InlineData(ClimateVariable.MinTemperature, ModelVersion.Dave)]
    public void GetTargetConfig_ThrowsForInvalidVariable(
        ClimateVariable variable,
        ModelVersion version)
    {
        _config.Version = version;
        ScriptGenerator generator = new ScriptGenerator(_config);
        Assert.ThrowsAny<ArgumentException>(() => generator.GetTargetConfig(variable));
    }

    [Theory]
    [InlineData(true, 0)]
    [InlineData(false, 0)]
    [InlineData(true, 9)]
    [InlineData(false, 9)]
    public async Task GenerateVariableRechunkScript_ConditionallyDisablesCompression(
        bool compressionEnabled,
        int compressionLevel
    )
    {
        StaticMockDataset dataset = new("/input");
        _config.CompressOutput = compressionEnabled;
        _config.CompressionLevel = compressionLevel;
        ScriptGenerator generator = new ScriptGenerator(_config);
        string script = await generator.GenerateVariableRechunkScript(dataset, ClimateVariable.Temperature);
        string[] scriptLines = await File.ReadAllLinesAsync(script);
        IEnumerable<string> ncpdqLines = scriptLines.Where(l => l.Contains("ncpdq"));

        if (compressionEnabled)
            Assert.All(ncpdqLines, l => Assert.Contains($"-L{compressionLevel}", l));
        else
            Assert.All(ncpdqLines, l => Assert.DoesNotContain("-L", l));
    }

    [Fact]
    public async Task TestGenerateWrapperScript()
    {
        string script = Path.GetTempFileName();
        await File.WriteAllLinesAsync(script, ["#!/usr/bin/bash", "echo x"]);
        string wrapper = await ScriptGenerator.GenerateWrapperScript(outputDirectory, [script]);
        string output = await File.ReadAllTextAsync(wrapper);

        // The wrapper script should call each subscript passed into it.
        Assert.Contains($"\n{script}\n", output);
    }

    /// <summary>
    /// Do a full integration test of the script generator under controlled
    /// conditions and do a full string comparison against the expected output.
    /// </summary>
    /// <remarks>
    /// This will be brittle, but thorough. If this fails, the other tests can
    /// be used to figure out why. If no other tests fail, then we need more
    /// tests!
    /// </remarks>
    [Fact]
    public async Task GenerateScriptsAsync_IntegrationTest()
    {
        const string inputDirectory = "/input";
        DynamicMockDataset dataset = new(inputDirectory, outputDirectory);
        NarClim2Config config = new()
        {
            Project = "test",
            Queue = "megamem",
            Walltime = "06:30:00",
            Ncpus = 2,
            Memory = 64,
            InputDirectory = inputDirectory,
            OutputDirectory = outputDirectory,
            InputTimeStepHours = 1,
            OutputTimeStepHours = 3,
            Version = ModelVersion.Dave,
            ChunkSizeTime = 8192,
            ChunkSizeSpatial = 24,
            GridFile = "/home/giraffe/grid.nc",
            CompressionLevel = 8,
            CompressOutput = true,
            DryRun = true,
            Email = "test@example.com",
            JobFS = 128,
            VPDMethod = VPDMethod.AlduchovEskridge1996
        };
        ScriptGenerator generator = new(config);

        // Act.
        await generator.GenerateScriptsAsync(dataset);

        // Assert.
        AssertEmptyDirectory(Path.Combine(outputDirectory, "logs"));
        AssertEmptyDirectory(Path.Combine(outputDirectory, "streams"));
        AssertEmptyDirectory(Path.Combine(outputDirectory, "output", dataset.GetOutputDirectory()));
        AssertEmptyDirectory(Path.Combine(outputDirectory, "tmp", dataset.GetOutputDirectory()));

        string scriptsDirectory = Path.Combine(outputDirectory, "scripts");
        Assert.True(Directory.Exists(scriptsDirectory));
        Assert.NotEmpty(Directory.EnumerateFileSystemEntries(scriptsDirectory));

        string[] expectedScriptNames = [
            "calc_vpd_DynamicMockDataset",
            "cleanup_DynamicMockDataset",
            "mergetime_huss_DynamicMockDataset",
            "mergetime_pr_DynamicMockDataset",
            "mergetime_ps_DynamicMockDataset",
            "mergetime_rsds_DynamicMockDataset",
            "mergetime_sfcWind_DynamicMockDataset",
            "mergetime_tas_DynamicMockDataset",
            "rechunk_huss_DynamicMockDataset",
            "rechunk_pr_DynamicMockDataset",
            "rechunk_ps_DynamicMockDataset",
            "rechunk_rsds_DynamicMockDataset",
            "rechunk_sfcWind_DynamicMockDataset",
            "rechunk_tas_DynamicMockDataset",
            "rechunk_vpd_DynamicMockDataset",
            "submit_DynamicMockDataset"
        ];

        Assert.Equal(expectedScriptNames.Count(), Directory.EnumerateFileSystemEntries(scriptsDirectory).Count());

        // Name of the directory containing this test's data files.
        const string resourcePrefix = "GenerateScriptsAsync_IntegrationTest";

        foreach (string scriptName in expectedScriptNames)
        {
            string actualScriptPath = Path.Combine(scriptsDirectory, scriptName);
            Assert.True(File.Exists(actualScriptPath), $"Script {actualScriptPath} does not exist.");
            string actualScript = await File.ReadAllTextAsync(actualScriptPath);

            // Read expected script from resource in assembly.
            string expectedScript = await ReadResource($"{resourcePrefix}.{scriptName}");
            expectedScript = expectedScript.Replace("@#OUTPUT_DIRECTORY#@", outputDirectory);

            // No custom error messages in xunit, apparently.
            if (expectedScript != actualScript)
                outputHelper.WriteLine($"Script {scriptName} is invalid");

            Assert.Equal(expectedScript, actualScript);
        }
    }

    private async Task<string> ReadResource(string resourceName)
    {
        Assembly assembly = typeof(ScriptGeneratorTests).Assembly;
        string resource = $"ClimateProcessing.Tests.Data.{resourceName}";
        using Stream? stream = assembly.GetManifestResourceStream(resource);
        if (stream is null)
            throw new ArgumentException($"Resource {resource} not found.");

        using StreamReader reader = new StreamReader(stream);
        return await reader.ReadToEndAsync();
    }

    private static void AssertEmptyDirectory(string directory)
    {
        Assert.True(Directory.Exists(directory));
        Assert.Empty(Directory.EnumerateFileSystemEntries(directory));
    }
}
