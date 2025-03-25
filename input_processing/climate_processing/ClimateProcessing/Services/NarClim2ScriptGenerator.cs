using ClimateProcessing.Models;

namespace ClimateProcessing.Services;

/// <summary>
/// Script generator specifically for NarClim2 datasets.
/// </summary>
public class NarClim2ScriptGenerator : ScriptGenerator, IScriptGenerator<NarClim2Dataset>
{
    /// <summary>
    /// Name of the directory in which the files with corrected rlon values are
    /// stored, relative to ${WORK_DIR}.
    /// </summary>
    private const string rlonDir = "corrected_rlon";

    /// <summary>
    /// Creates a new instance of the <see cref="NarClim2ScriptGenerator"/> class.
    /// </summary>
    /// <param name="config">The processing configuration.</param>
    public NarClim2ScriptGenerator(NarClim2Config config) : base(config)
    {
    }

    /// <inheritdoc/>
    public async Task<string> GenerateScriptsAsync(NarClim2Dataset dataset)
    {
        return await base.GenerateScriptsAsync(dataset);
    }

    /// <summary>
    /// Write commands to the script which run before the merge step, which
    /// fix the incorrect rlon values in the input files.
    /// </summary>
    /// <param name="writer"></param>
    /// <param name="dataset"></param>
    /// <param name="variable"></param>
    /// <returns></returns>
    /// <exception cref="ArgumentException"></exception>
    protected override async Task WritePreMerge(TextWriter writer, IClimateDataset dataset, ClimateVariable variable)
    {
        if (dataset is not NarClim2Dataset narclim2)
            // Should never happen.
            throw new ArgumentException($"Expected dataset to be of type {typeof(NarClim2Dataset)}, but got {dataset.GetType().Name}");

        // Some narclim2 files have incorrect rlon values. We can use the
        // setvar.py script to correct them.
        string valuesFile = GetRlonValuesFile(narclim2);

        await writer.WriteLineAsync("# Correct rlon values, which are incorrect in some narclim2 files.");
        await writer.WriteLineAsync($"RLON_VALUES_FILE=\"{valuesFile}\"");
        await writer.WriteLineAsync($"CORRECTED_RLON_DIR=\"${{WORK_DIR}}/{rlonDir}\"");
        await writer.WriteLineAsync($"mkdir -p \"${{CORRECTED_RLON_DIR}}\"");
        await writer.WriteLineAsync("log \"Correcting rlon values...\"");
        await writer.WriteLineAsync($"for FILE in \"${{IN_DIR}}\"/*.nc; do");
        await writer.WriteLineAsync($"    log \"Correcting rlon values in file $(basename \"${{FILE}}\")...\"");
        await writer.WriteLineAsync($"    setvar.py --in-file \"${{FILE}}\" --out-file \"${{CORRECTED_RLON_DIR}}/$(basename \"${{FILE}}\")\" --values-file \"${{RLON_VALUES_FILE}}\" --var rlon");
        await writer.WriteLineAsync("done");
        // The output of this step of the procesing is the new input directory.
        await writer.WriteLineAsync($"{inDirVariable}=\"${{CORRECTED_RLON_DIR}}\"");
        await writer.WriteLineAsync("log \"Successfully corrected all rlon values.\"");
        await writer.WriteLineAsync();
    }

    /// <summary>
    /// Get the absolute path to the plaintext file containing the correct rlon
    /// values.
    /// </summary>
    /// <param name="narclim2">The dataset.</param>
    /// <returns>The absolute path to the file.</returns>
    protected string GetRlonValuesFile(NarClim2Dataset narclim2)
    {
        string directory = narclim2.BasePath;
        string fileName = NarClim2Constants.Files.GetRlonValuesFile(narclim2.Domain);

        // Will fail if output directory is root directory, which should
        // probably never happen. Maybe it's possible if we're running in
        // docker, depending on how the container is set up. In that case, we
        // should add this file name as a CLI option.
        return Path.Combine(directory, "..", fileName);
    }
}
