using CommandLine;
using ClimateProcessing.Models;
using ClimateProcessing.Services;

// Parse command line arguments
var result = Parser.Default.ParseArguments<ProcessingConfig>(args);

result.WithParsed(config =>
{
    try
    {
        config.Validate();

        // Create dataset instance based on type
        IClimateDataset dataset = config.DatasetType.ToLower() switch
        {
            "narclim2" => new NarClim2Dataset(config.InputDirectory),
            _ => throw new ArgumentException($"Unsupported dataset type: {config.DatasetType}")
        };

        var scriptGenerator = new ScriptGenerator(config);

        // Generate scripts.
        string submissionScript = scriptGenerator.GenerateScripts(dataset);

        Console.WriteLine($"Processing scripts have been generated in:");
        Console.WriteLine($"{Path.GetDirectoryName(submissionScript)}");
        Console.WriteLine("\nTo submit the job to PBS, run:");
        Console.WriteLine($"{submissionScript}");
    }
    catch (Exception ex)
    {
        Console.Error.WriteLine($"Error: {ex.Message}");
        Environment.Exit(1);
    }
});
