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
        IClimateDataset dataset = config.DatasetType.ToLower() switch
        {
            "narclim2" => NarClim2Dataset.Create(config),
            _ => throw new ArgumentException($"Unsupported dataset type: {config.DatasetType}")
        };

        // Generate scripts using factory
        var factory = serviceProvider.GetRequiredService<IScriptGeneratorFactory>();
        string submissionScript = await factory.GenerateScriptsAsync(dataset);

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
