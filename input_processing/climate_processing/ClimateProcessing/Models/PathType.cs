namespace ClimateProcessing.Models;

/// <summary>
/// Represents the type of a path.
/// </summary>
public enum PathType
{
    /// <summary>
    /// Final output files that should be preserved
    /// </summary>
    Output,

    /// <summary>
    /// Temporary/intermediate working files
    /// </summary>
    Working,

    /// <summary>
    /// Generated script files
    /// </summary>
    Script,

    /// <summary>
    /// Log files
    /// </summary>
    Log,

    /// <summary>
    /// Stream output files (stdout/stderr)
    /// </summary>
    Stream
}
