using CommandLine;
using ClimateProcessing.Models;
using ClimateProcessing.Services;

// Parse command line arguments
var result = Parser.Default.ParseArguments<ProcessingConfig>(args);

await result.WithParsedAsync(async config =>
{
    try
    {
        config.Validate();

        var processor = new DatasetCombinationProcessor(config);
        await processor.ProcessAllCombinations();

        Console.WriteLine("\nAll dataset combinations have been processed.");
        Console.WriteLine("To submit the jobs to PBS, run each submission script in the output directories.");
    }
    catch (Exception ex)
    {
        Console.Error.WriteLine($"Error: {ex.Message}");
        Environment.Exit(1);
    }
});
