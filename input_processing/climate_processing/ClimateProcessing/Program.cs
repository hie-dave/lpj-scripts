using CommandLine;
using ClimateProcessing.Models;
using ClimateProcessing.Services;
using Microsoft.Extensions.DependencyInjection;

// Parse command line arguments
var result = Parser.Default.ParseArguments<ProcessingConfig>(args);

await result.WithParsedAsync(async config =>
{
    try
    {
        config.Validate();

        // Set up dependency injection
        var services = new ServiceCollection();

        // Register config
        services.AddSingleton(config);

        // Register script generators
        services.AddSingleton<IScriptGeneratorFactory, ScriptGeneratorFactory>();
        services.AddSingleton<ScriptGenerator>();
        services.AddSingleton<IScriptGenerator<NarClim2Dataset>, NarClim2ScriptGenerator>();

        var serviceProvider = services.BuildServiceProvider();

        // Create dataset instance based on type
        IEnumerable<IClimateDataset> datasets = config.DatasetType.ToLower() switch
        {
            "narclim2" => NarClim2Dataset.CreateAll(config),
            _ => throw new ArgumentException($"Unsupported dataset type: {config.DatasetType}")
        };

        // Generate scripts using factory
        var factory = serviceProvider.GetRequiredService<IScriptGeneratorFactory>();
        List<string> scripts = new List<string>();
        foreach (IClimateDataset dataset in datasets)
        {
            string submissionScript = await factory.GenerateScriptsAsync(dataset);
            scripts.Add(submissionScript);
        }

        string wrapper = await ScriptGenerator.GenerateWrapperScript(config.OutputDirectory, scripts);
        Console.WriteLine($"Processing scripts have been generated in:");
        Console.WriteLine($"{Path.GetDirectoryName(config.OutputDirectory)}");
        Console.WriteLine("\nTo submit the jobs to PBS, run:");
        Console.WriteLine($"{wrapper}");
    }
    catch (Exception ex)
    {
        Console.Error.WriteLine($"Error: {ex.Message}");
        Environment.Exit(1);
    }
});
