using ClimateProcessing.Models;
using ClimateProcessing.Constants;
using ClimateProcessing.Configuration;

namespace ClimateProcessing.Services;

/// <summary>
/// Handles generation of PBS-specific script content and headers.
/// </summary>
public class PBSWriter
{
    private readonly PBSConfig config;
    private readonly IPathManager pathManager;

    /// <summary>
    /// Creates a new PBS script generator.
    /// </summary>
    /// <param name="config">The processing configuration.</param>
    public PBSWriter(PBSConfig config, IPathManager pathManager)
    {
        this.config = config;
        this.pathManager = pathManager;
    }

    /// <summary>
    /// Write the PBS header for a job.
    /// </summary>
    /// <param name="writer">The text writer to which the header will be written.</param>
    /// <param name="jobName">The job name.</param>
    /// <param name="storageDirectives">The storage directives required by the job.</param>
    public async Task WritePBSHeader(
        TextWriter writer,
        string jobName,
        IEnumerable<PBSStorageDirective> storageDirectives)
    {
        string logFile = pathManager.GetLogFilePath(jobName);
        string streamFile = pathManager.GetStreamFilePath(jobName);

        await writer.WriteLineAsync("#!/usr/bin/env bash");
        await writer.WriteLineAsync($"#PBS -N {jobName}");
        await writer.WriteLineAsync($"#PBS -o {logFile}");
        await writer.WriteLineAsync($"#PBS -P {config.Project}");
        await writer.WriteLineAsync($"#PBS -q {config.Queue}");
        await writer.WriteLineAsync($"#PBS -l walltime={config.Walltime}");
        await writer.WriteLineAsync($"#PBS -l ncpus={config.Ncpus}");
        await writer.WriteLineAsync($"#PBS -l mem={config.Memory}GB");
        await writer.WriteLineAsync($"#PBS -l jobfs={config.JobFS}GB");
        await writer.WriteLineAsync($"#PBS -j oe");

        if (!string.IsNullOrEmpty(config.Email))
        {
            await writer.WriteLineAsync($"#PBS -M {config.Email}");
            await writer.WriteLineAsync($"#PBS -m abe");
        }

        // Add storage directives if required
        if (storageDirectives.Any())
            await writer.WriteLineAsync(PBSStorageHelper.FormatStorageDirectives(storageDirectives));

        // Add blank line after header
        await writer.WriteLineAsync("");

        await WriteAutoGenerateHeader(writer);

        // Error handling
        await writer.WriteLineAsync("# Exit immediately if any command fails.");
        await writer.WriteLineAsync("set -euo pipefail");
        await writer.WriteLineAsync();

        // Load required modules
        await writer.WriteLineAsync("# Load required modules.");
        await writer.WriteLineAsync("module purge");
        await writer.WriteLineAsync("module load pbs netcdf cdo nco python3/3.12.1");
        await writer.WriteLineAsync();

        // Create temporary directory and cd into it
        await writer.WriteLineAsync("# Create temporary directory and cd into it.");
        await writer.WriteLineAsync("WORK_DIR=\"$(mktemp -d -p \"${PBS_JOBFS}\")\"");
        await writer.WriteLineAsync("cd \"${WORK_DIR}\"");
        await writer.WriteLineAsync();

        // Cleanup on exit
        await writer.WriteLineAsync("# Delete the temporary directory on exit.");
        await writer.WriteLineAsync("trap 'cd \"${PBS_JOBFS}\"; rm -rf \"${WORK_DIR}\"' EXIT");
        await writer.WriteLineAsync();

        // Set up logging that streams all output
        await writer.WriteLineAsync("# Stream all output to a log file without buffering.");
        await writer.WriteLineAsync($"STREAM_FILE=\"{streamFile}\"");
        await writer.WriteLineAsync("rm -f \"${STREAM_FILE}\"");
        await writer.WriteLineAsync("exec 1> >(tee -a \"${STREAM_FILE}\") 2>&1");
        await writer.WriteLineAsync();

        // Add logging function
        await writer.WriteLineAsync("# Print a log message.");
        await writer.WriteLineAsync("log() {");
        await writer.WriteLineAsync("    echo \"[$(date)] $*\"");
        await writer.WriteLineAsync("}");

        // Add blank line after header
        await writer.WriteLineAsync("");
    }

    /// <summary>
    /// Write a comment to a script which indicates that it was automatically generated.
    /// </summary>
    /// <param name="writer">The text writer to which the comment will be written.</param>
    private static async Task WriteAutoGenerateHeader(TextWriter writer)
    {
        await writer.WriteLineAsync("# This script was automatically generated. Do not modify.");
        await writer.WriteLineAsync();
    }
}
