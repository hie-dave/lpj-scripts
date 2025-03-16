using CommandLine;
using ClimateProcessing.Services;

namespace ClimateProcessing.Models;

/// <summary>
/// Contains narclim2-specific CLI options.
/// </summary>
[Verb("narclim2", HelpText = "Process NARCliM2 data.")]
public class NarClim2Config : ProcessingConfig
{
    [Option("domain", HelpText = "Domains to process. Valid values: AUS-18, NARCliM2-0-SEAus-04. Default: process all domains.")]
    public IEnumerable<string>? Domains { get; set; }

    [Option("gcm", HelpText = "Global Climate Models to process. Valid values: ACCESS-ESM1-5, EC-Earth3-Veg, MPI-ESM1-2-HR, NorESM2-MM, UKESM1-0-LL. Default: process all GCMs.")]
    public IEnumerable<string>? GCMs { get; set; }

    [Option("experiment", HelpText = "Experiments to process. Valid values: historical, ssp126, ssp370. Default: process all experiments.")]
    public IEnumerable<string>? Experiments { get; set; }

    [Option("rcm", HelpText = "RCMs to process. Valid values: WRF412R3, WRF412R5. Default: process all RCMs.")]
    public IEnumerable<string>? RCMs { get; set; }

    /// <inheritdoc /> 
    public override IEnumerable<NarClim2Dataset> CreateDatasets()
    {
        IEnumerable<NarClim2Domain> domains = GetDomains();
        IEnumerable<NarClim2GCM> gcms = GetGCMs();
        IEnumerable<NarClim2Experiment> experiments = GetExperiments();
        IEnumerable<NarClim2RCM> rcms = GetRCMs();
        NarClim2Frequency frequency = NarClim2Constants.ParseFrequency(InputTimeStep.Hours);

        if (!domains.Any())
            throw new InvalidOperationException("No domains specified. Use the --domain option to specify at least one domain.");
        if (!gcms.Any())
            throw new InvalidOperationException("No GCMs specified. Use the --gcm option to specify at least one GCM.");
        if (!experiments.Any())
            throw new InvalidOperationException("No experiments specified. Use the --experiment option to specify at least one experiment.");
        if (!rcms.Any())
            throw new InvalidOperationException("No RCMs specified. Use the --rcm option to specify at least one RCM.");

        return from domain in domains
               from gcm in gcms
               from experiment in experiments
               from rcm in rcms
               select new NarClim2Dataset(
                   inputPath: InputDirectory,
                   domain: domain,
                   gcm: gcm,
                   experiment: experiment,
                   rcm: rcm,
                   frequency: frequency);
    }

    /// <inheritdoc />
    public override ScriptGenerator CreateScriptGenerator()
    {
        return new NarClim2ScriptGenerator(this);
    }

    /// <summary>
    /// Get the list of domains to process.
    /// </summary>
    /// <returns>The list of domains to process.</returns>
    private IEnumerable<NarClim2Domain> GetDomains()
    {
        if (Domains == null || !Domains.Any())
            return Enum.GetValues<NarClim2Domain>();
        return Domains.Select(NarClim2Constants.DomainNames.FromString);
    }

    /// <summary>
    /// Get the list of GCMs to process.
    /// </summary>
    /// <returns>The list of GCMs to process.</returns>
    private IEnumerable<NarClim2GCM> GetGCMs()
    {
        if (GCMs == null || !GCMs.Any())
            return Enum.GetValues<NarClim2GCM>();
        return GCMs.Select(NarClim2Constants.GCMNames.FromString);
    }

    /// <summary>
    /// Get the list of experiments to process.
    /// </summary>
    /// <returns>The list of experiments to process.</returns>
    private IEnumerable<NarClim2Experiment> GetExperiments()
    {
        if (Experiments == null || !Experiments.Any())
            return Enum.GetValues<NarClim2Experiment>();
        return Experiments.Select(NarClim2Constants.ExperimentNames.FromString);
    }

    /// <summary>
    /// Get the list of RCMs to process.
    /// </summary>
    /// <returns>The list of RCMs to process.</returns>
    private IEnumerable<NarClim2RCM> GetRCMs()
    {
        if (RCMs == null || !RCMs.Any())
            return Enum.GetValues<NarClim2RCM>();
        return RCMs.Select(NarClim2Constants.RCMNames.FromString);
    }
}
