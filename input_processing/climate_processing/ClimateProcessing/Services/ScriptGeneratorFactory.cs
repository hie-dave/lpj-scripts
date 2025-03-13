using Microsoft.Extensions.DependencyInjection;
using ClimateProcessing.Models;

namespace ClimateProcessing.Services;

/// <summary>
/// Factory for creating dataset-specific script generators with fallback behavior.
/// </summary>
public class ScriptGeneratorFactory : IScriptGeneratorFactory
{
    private readonly IServiceProvider _serviceProvider;

    public ScriptGeneratorFactory(IServiceProvider serviceProvider)
    {
        _serviceProvider = serviceProvider;
    }

    public async Task<string> GenerateScriptsAsync(IClimateDataset dataset)
    {
        if (dataset is NarClim2Dataset narclim2)
        {
            var generator = _serviceProvider.GetService<IScriptGenerator<NarClim2Dataset>>();
            if (generator != null)
                return await generator.GenerateScriptsAsync(narclim2);
        }

        // Fall back to default generator
        var defaultGenerator = _serviceProvider.GetRequiredService<ScriptGenerator>();
        return await defaultGenerator.GenerateScriptsAsync(dataset);
    }
}
