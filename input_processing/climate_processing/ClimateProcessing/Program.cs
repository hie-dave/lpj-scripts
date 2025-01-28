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

        // Generate processing script
        string processingScript = scriptGenerator.GenerateProcessingScript(dataset);
        string processingScriptPath = Path.Combine(config.OutputDirectory, $"process_{dataset.DatasetName}.sh");
        File.WriteAllText(processingScriptPath, processingScript);
        File.SetUnixFileMode(processingScriptPath, UnixFileMode.UserRead | UnixFileMode.UserWrite | UnixFileMode.UserExecute);

        // Generate submission script
        string submissionScript = scriptGenerator.GenerateSubmissionScript(processingScriptPath);
        string submissionScriptPath = Path.Combine(config.OutputDirectory, "submit_job.sh");
        File.WriteAllText(submissionScriptPath, submissionScript);
        File.SetUnixFileMode(submissionScriptPath, UnixFileMode.UserRead | UnixFileMode.UserWrite | UnixFileMode.UserExecute);

        Console.WriteLine($"Processing scripts have been generated in {config.OutputDirectory}:");
        Console.WriteLine($"- Processing script: {Path.GetFileName(processingScriptPath)}");
        Console.WriteLine($"- Submission script: {Path.GetFileName(submissionScriptPath)}");
        Console.WriteLine("\nTo submit the job to PBS, run:");
        Console.WriteLine($"cd {config.OutputDirectory} && ./submit_job.sh");
    }
    catch (Exception ex)
    {
        Console.Error.WriteLine($"Error: {ex.Message}");
        Environment.Exit(1);
    }
});
