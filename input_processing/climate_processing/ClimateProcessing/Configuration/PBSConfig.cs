using ClimateProcessing.Constants;
using ClimateProcessing.Models;

namespace ClimateProcessing.Configuration;

/// <summary>
/// Configuration options for PBS jobs.
/// </summary>
public class PBSConfig
{
    /// <summary>
    /// The queue to submit the job to.
    /// </summary>
    public string Queue { get; init; }

    /// <summary>
    /// The number of CPUs to use for the job.
    /// </summary>
    public int Ncpus { get; init; }

    /// <summary>
    /// The amount of memory to use for the job, in MB.
    /// </summary>
    public int Memory { get; init; }

    /// <summary>
    /// The amount of job file space to use for the job, in MB.
    /// </summary>
    public int JobFS { get; init; }

    /// <summary>
    /// The project to use for the job.
    /// </summary>
    public string Project { get; init; }

    /// <summary>
    /// The walltime to use for the job.
    /// </summary>
    public PBSWalltime Walltime { get; init; }

    /// <summary>
    /// The email address to use for the job.
    /// </summary>
    public string? Email { get; init; }

    /// <summary>
    /// Create a new instance of the PBSConfig class.
    /// </summary>
    /// <param name="queue">The queue to submit the job to.</param>
    /// <param name="ncpus">The number of CPUs to allocate to the job.</param>
    /// <param name="memory">The amount of memory to allocate to the job, in GiB.</param>
    /// <param name="jobFS">The amount of JobFS space to allocate to the job, in GiB.</param>
    /// <param name="project">The project to which the job will be debited.</param>
    /// <param name="walltime">Maximum walltime the job is allowed to use.</param>
    /// <param name="email">The optional email address to which job notifications will be sent.</param>
    public PBSConfig(
        string queue,
        int ncpus,
        int memory,
        int jobFS,
        string project,
        PBSWalltime walltime,
        string? email = null)
    {
        Queue = queue;
        Ncpus = ncpus;
        Memory = memory;
        JobFS = jobFS;
        Project = project;
        Walltime = walltime;
        Email = email;
    }

    /// <summary>
    /// Creates a lightweight PBS configuration with default values.
    /// </summary>
    /// <param name="jobfs">The amount of JobFS space to allocate to the job, in GiB.</param>
    /// <param name="project">The project to which the job will be debited.</param>
    /// <param name="email">The optional email address to which job notifications will be sent.</param>
    /// <returns>A new PBS configuration instance.</returns>
    public static PBSConfig LightWeight(int jobfs, string project, string? email) => new(
        PBSConstants.QueueNormal,
        PBSConstants.LightweightNcpus,
        PBSConstants.LightweightMemory,
        jobfs,
        project,
        PBSWalltime.MaxValue,
        email
    );
}
