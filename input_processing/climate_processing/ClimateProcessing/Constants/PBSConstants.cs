namespace ClimateProcessing.Constants;

public static class PBSConstants
{
    /// <summary>
    /// Default queue for lightweight jobs.
    /// </summary>
    public const string QueueNormal = "normal";

    /// <summary>
    /// Number of CPUs for lightweight jobs.
    /// </summary>
    public const int LightweightNcpus = 1;

    /// <summary>
    /// Memory limit in GiB for lightweight jobs.
    /// </summary>
    public const int LightweightMemory = 4;
}

