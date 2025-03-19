using Xunit;
using ClimateProcessing.Services;
using ClimateProcessing.Models;
using ClimateProcessing.Units;
using System.Text.RegularExpressions;

namespace ClimateProcessing.Tests.Services;

public class ScriptGeneratorTests
{
    private readonly NarClim2Config _config = new()
    {
        Project = "test",
        Queue = "normal",
        Walltime = "01:00:00",
        Ncpus = 1,
        Memory = 4,
        OutputDirectory = "/output",
        InputTimeStepHours = 1,
        OutputTimeStepHours = 24
    };

    private readonly ScriptGenerator _generator;

    public ScriptGeneratorTests()
    {
        _generator = new ScriptGenerator(_config);
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
    [InlineData(1, 8, 100, "normal", "01:00:00", "test_job", "p123", "test@example.com")]
    [InlineData(2, 192, 1024, "hugemem", "48:00:00", "asdf_name", "x987", null)]
    public async Task WritesPBSHeader(int ncpu, int memory, int jobfs, string queue, string walltime, string jobName, string project, string? email)
    {
        // Arrange
        StringWriter writer = new();
        NarClim2Config config = new()
        {
            Project = project,
            Queue = queue,
            Walltime = walltime,
            Ncpus = ncpu,
            Memory = memory,
            JobFS = jobfs,
            Email = email ?? string.Empty,
            OutputDirectory = Directory.CreateTempSubdirectory("output").FullName
        };
        ScriptGenerator generator = new(config);

        // Act
        // Don't generate a lightweight script header for this test.
        await generator.WritePBSHeader(writer, jobName, false);
        string result = writer.ToString();

        // Assert
        Assert.Contains($"#PBS -N {jobName}", result);
        Assert.Contains($"#PBS -P {project}", result);
        Assert.Contains($"#PBS -q {queue}", result);
        Assert.Contains($"#PBS -l walltime={walltime}", result);
        Assert.Contains($"#PBS -l ncpus={ncpu}", result);
        Assert.Contains($"#PBS -l mem={memory}GB", result);
        Assert.Contains($"#PBS -l jobfs={jobfs}GB", result);
        Assert.Contains("#PBS -j oe", result);

        if (email != null)
        {
            Assert.Contains($"#PBS -M {email}", result);
            Assert.Contains("#PBS -m abe", result);
        }
        else
        {
            Assert.DoesNotContain("#PBS -M", result);
            Assert.DoesNotContain("#PBS -m", result);
        }

        // Verify proper line endings.
        Assert.DoesNotContain(result, "\r");

        // Verify proper line endings and no empty lines with whitespace.
        string[] lines = result.Split("\n");

        // Get index of last line starting with #PBS.
        string lastPBSLine = lines.Last(line => line.StartsWith("#PBS "));
        int lastPBSLineIndex = Array.LastIndexOf(lines, lastPBSLine);

        // Assert that all lines before the last #PBS line are not empty.
        Assert.All(lines[0..lastPBSLineIndex], Assert.NotEmpty);
    }

    [Theory]
    [InlineData("/data", "/home/asdf/x.y", "/x y z", 0)]  // None of these require a storage directive.
    [InlineData("/g/data/x", "/home/asdf/x.y", "/x y z", 0)]  // /g/data/x isn't really a valid path - it doesn't reference anything inside the gdata file tree
    [InlineData("/scratch/x", "/home/asdf/x.y", "/x y z", 0)]  // /scratch/x isn't really a valid path - it doesn't reference anything inside the scratch file tree
    [InlineData("/g/data/asdf/x.y", "/data", "/data", 1)]    // Single reference to gdata.
    [InlineData("/data", "/g/data/asdf/x.y", "/data", 1)]    // Single reference to gdata.
    [InlineData("/data", "/data", "/g/data/asdf/x.y", 1)]    // Single reference to gdata.
    [InlineData("/scratch/xyz/lkj", "/data", "/data", 1)]    // Single reference to scratch.
    [InlineData("/data", "/scratch/xyz/lkj", "/data", 1)]    // Single reference to scratch.
    [InlineData("/data", "/data", "/scratch/xyz/lkj", 1)]    // Single reference to scratch.
    [InlineData("/g/data/asdf/x.y", "/g/data/asdf/afds", "/data", 1)]    // Same gdata path used twice
    [InlineData("/scratch/xyz/lkj", "/scratch/xyz/fds", "/data", 1)]    // Same scratch path used twice
    [InlineData("/g/data/asdf/fdsa", "/g/data/asdf/asdf", "/g/data/asdf/bql", 1)]    // Same gdata path used thrice
    [InlineData("/scratch/ff/dhwh", "/scratch/ff/wwjd", "/scratch/ff/wwdd", 1)]    // Same scratch paths used thrice
    [InlineData("/g/data/asdf/fdsa", "/g/data/x/y", "/data", 2)] // Two different gdata paths.
    [InlineData("/scratch/ff/ttfn", "/scratch/gg/test", "/home/test/pi.file", 2)] // Two different scratch paths.
    [InlineData("/g/data/asdf/fdsa", "/g/data/k/asdf", "/g/data/llm/bql", 3)] // Three different gdata paths.
    [InlineData("/scratch/ff/wtfn", "/scratch/gg/another_test", "/scratch/hh/wow", 3)] // Three different scratch paths.
    [InlineData("/g/data/jk/hpig", "/scratch/lkjihg/f", "home/glarble", 2)] // Mixture of gdata and scratch paths.
    public async Task WritePBSHeader_GeneratesCorrectStorageDirectives(
        string inputDir,
        string gridlist,
        string outputDir,
        int expectedDirectives)
    {
        // Arrange
        StringWriter writer = new();
        NarClim2Config config = new()
        {
            Project = "p123",
            Queue = "normal",
            Walltime = "01:00:00",
            Ncpus = 2,
            Memory = 8,
            InputDirectory = inputDir,
            OutputDirectory = outputDir,
            GridFile = gridlist
        };
        ScriptGenerator generator = new(config);

        // Act
        await generator.WritePBSHeader(writer, "test_job", true);
        string result = writer.ToString();

        // Assert
        int actualDirectives = Regex.Matches(result, @"#PBS -l storage=").Count;
        Assert.Equal(expectedDirectives > 0 ? 1 : 0, actualDirectives);

        if (expectedDirectives == 0)
            return;

        string line = result.Split("\n").First(line => line.StartsWith("#PBS -l storage="));
        string[] parts = line.Split('+');
        Assert.Equal(expectedDirectives, parts.Length);
    }

    [Theory]
    [InlineData("input dir with spaces")]
    [InlineData("/path/with/special/$chars")]
    [InlineData("/normal/path")]
    public async Task GenerateVariableMergeScript_QuotesVariablesSafely(
        string inputDir)
    {
        // Arrange
        NarClim2Config config = new()
        {
            Project = "test",
            Queue = "normal",
            Walltime = "01:00:00",
            Ncpus = 1,
            Memory = 4,
            OutputDirectory = Directory.CreateTempSubdirectory("script_generator_tests_output").FullName,
            InputDirectory = inputDir,
            InputTimeStepHours = 3,
            OutputTimeStepHours = 24,
        };
        ScriptGenerator generator = new(config);
        MockDataset dataset = new(inputDir);

        // Act
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

    // TODO: This test could be significantly expanded upon.
    [Theory]
    [InlineData(ClimateVariable.Temperature, "tas", "K")]
    [InlineData(ClimateVariable.Precipitation, "pr", "kg m-2 s-1")]
    public async Task GenerateVariableMergeScript_GeneratesValidCDOCommand(
        ClimateVariable variable,
        string varName,
        string inputUnits)
    {
        // Arrange
        NarClim2Config config = new()
        {
            Project = "test",
            Queue = "normal",
            Walltime = "01:00:00",
            Ncpus = 1,
            Memory = 4,
            OutputDirectory = Directory.CreateTempSubdirectory().FullName,
            InputDirectory = "/input",
            InputTimeStepHours = 1,
            OutputTimeStepHours = 24
        };
        ScriptGenerator generator = new(config);
        MockDataset dataset = new("/input", varName, inputUnits);

        // Act
        string scriptPath = await generator.GenerateVariableMergeScript(
            dataset,
            variable);
        string scriptContent = await File.ReadAllTextAsync(scriptPath);

        // E.g.
        // cdo -L -O -v -z zip1 -daymean,24 -subc,273.15 -setattribute,'tas@units=degC' -unpack  \"${FILE}\" \"${REMAP_DIR}/$(basename \"${FILE}\")\"
        string line = scriptContent.Split("\n").First(l => l.Contains("cdo -"));

        // Assert
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

        // Rest of script should be valid.
        ValidateScript(scriptContent);
    }

    // fixme - broken
    [Theory]
    [InlineData(VPDMethod.Magnus)]
    [InlineData(VPDMethod.Buck1981)]
    [InlineData(VPDMethod.AlduchovEskridge1996)]
    [InlineData(VPDMethod.AllenFAO1998)]
    [InlineData(VPDMethod.Sonntag1990)]
    public async Task GenerateVPDScript_GeneratesCorrectEquations(VPDMethod method)
    {
        // Arrange
        NarClim2Config config = new()
        {
            Project = "test",
            Queue = "normal",
            Walltime = "01:00:00",
            Ncpus = 1,
            Memory = 4,
            OutputDirectory = "/output",
            InputDirectory = "/input",
            VPDMethod = method
        };
        ScriptGenerator generator = new(config);
        MockDataset dataset = new("/input");

        // Act
        string scriptPath = await generator.GenerateVPDScript(dataset);
        string scriptContent = await File.ReadAllTextAsync(scriptPath);

        // Assert
        // Variables should be properly quoted
        Assert.Contains("\"${HUSS_FILE}\"", scriptContent);
        Assert.Contains("\"${PS_FILE}\"", scriptContent);
        Assert.Contains("\"${TAS_FILE}\"", scriptContent);
        Assert.Contains("\"${OUT_FILE}\"", scriptContent);
        Assert.Contains("\"${EQN_FILE}\"", scriptContent);

        // Should contain VPD equation components
        Assert.Contains("_esat=", scriptContent);
        Assert.Contains("_e=", scriptContent);
        Assert.Contains("vpd=", scriptContent);

        // Should use CDO exprf operator
        Assert.Contains("cdo", scriptContent);
        Assert.Contains("exprf", scriptContent);
        Assert.Contains("-merge", scriptContent);

        // Should have proper input file handling
        Assert.Contains("${HUSS_FILE}", scriptContent);
        Assert.Contains("${PS_FILE}", scriptContent);
        Assert.Contains("${TAS_FILE}", scriptContent);
    }

    // fixme - broken
    [Theory]
    [InlineData(true)]   // With VPD calculation
    [InlineData(false)]  // Without VPD calculation
    public async Task GenerateScriptsAsync_HandlesVPDDependenciesCorrectly(
        bool requiresVPD)
    {
        // Arrange
        NarClim2Config config = new()
        {
            Project = "test",
            Queue = "normal",
            Walltime = "01:00:00",
            Ncpus = 1,
            Memory = 4,
            OutputDirectory = "/output",
            InputDirectory = "/input",
            Version = requiresVPD ? ModelVersion.Dave : ModelVersion.Trunk
        };
        ScriptGenerator generator = new(config);
        MockDataset dataset = new("/input");

        // Act
        string scriptPath = await generator.GenerateScriptsAsync(dataset);
        string scriptContent = await File.ReadAllTextAsync(scriptPath);

        // Assert
        // Job dependencies
        if (requiresVPD)
        {
            Assert.Contains("VPD_DEPS=", scriptContent);
            Assert.Contains("calc_vpd_", scriptContent);
            Assert.Contains("rechunk_vpd_", scriptContent);
        }
        else
        {
            Assert.DoesNotContain("VPD_DEPS=", scriptContent);
            Assert.DoesNotContain("calc_vpd_", scriptContent);
            Assert.DoesNotContain("rechunk_vpd_", scriptContent);
        }

        // Job dependency syntax
        Assert.Contains("-W depend=afterok:", scriptContent);
        
        // Should always have cleanup job
        Assert.Contains("cleanup_", scriptContent);
    }

    private class MockDataset : IClimateDataset
    {
        private readonly string basePath;
        private readonly string varName;
        private readonly string varUnits;

        public MockDataset(
            string basePath,
            string varName = "tas",
            string varUnits = "K")
        {
            this.basePath = basePath;
            this.varName = varName;
            this.varUnits = varUnits;
        }

        public string DatasetName => "mock";

        public string GetInputFilesDirectory(ClimateVariable variable)
            => basePath;

        public IEnumerable<string> GetInputFiles(ClimateVariable variable)
            => new[] { Path.Combine(basePath, "input.nc") };

        public string GetOutputDirectory() => "mock";

        public VariableInfo GetVariableInfo(ClimateVariable variable)
            => new(varName, varUnits);

        public string GenerateOutputFileName(ClimateVariable variable)
            => "output.nc";
    }
}
