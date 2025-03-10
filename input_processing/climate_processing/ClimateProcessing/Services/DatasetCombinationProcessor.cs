using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using ClimateProcessing.Models;

namespace ClimateProcessing.Services;

public class DatasetCombinationProcessor
{
    private readonly ProcessingConfig _config;
    private readonly ScriptGenerator _scriptGenerator;

    public DatasetCombinationProcessor(ProcessingConfig config)
    {
        _config = config;
        _scriptGenerator = new ScriptGenerator(config);
    }

    public IEnumerable<IClimateDataset> GenerateDatasetCombinations()
    {
        return _config.DatasetType.ToLower() switch
        {
            "narclim2" => GenerateNarClim2Combinations(),
            _ => throw new ArgumentException($"Unsupported dataset type: {_config.DatasetType}")
        };
    }

    private IEnumerable<IClimateDataset> GenerateNarClim2Combinations()
    {
        var domains = Enum.GetValues<NarClim2Domain>();
        var gcms = Enum.GetValues<NarClim2GCM>();
        var experiments = Enum.GetValues<NarClim2Experiment>();
        var rcms = Enum.GetValues<NarClim2RCM>();

        return from domain in domains
               from gcm in gcms
               from experiment in experiments
               from rcm in rcms
               select new NarClim2Dataset(
                   inputPath: _config.InputDirectory,
                   domain: domain,
                   gcm: gcm,
                   experiment: experiment,
                   rcm: rcm,
                   frequency: NarClim2Constants.ParseFrequency(_config.InputTimeStep.Hours)
               );
    }

    public async Task ProcessAllCombinationsAsync()
    {
        var combinations = GenerateDatasetCombinations();
        var tasks = new List<Task<string>>();

        foreach (var dataset in combinations)
        {
            // Create a new task for processing this combination
            tasks.Add(Task.Run(async () =>
            {
                try
                {
                    string submissionScript = await _scriptGenerator.GenerateScriptsAsync(dataset);
                    Console.WriteLine($"Generated scripts for combination: {dataset.DatasetName}");
                    Console.WriteLine($"Location: {submissionScript}\n");
                    return submissionScript;
                }
                catch (Exception ex)
                {
                    Console.Error.WriteLine($"Error processing {dataset.DatasetName}: {ex.Message}");
                    throw;
                }
            }));
        }

        await Task.WhenAll(tasks);
    }

    public async Task ProcessAllCombinations()
    {
        await ProcessAllCombinationsAsync();
    }
}
