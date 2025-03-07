using System.Text.RegularExpressions;

namespace ClimateProcessing.Models;

public class PBSStorageDirective
{
    /// <summary>
    /// Full path to a resource which requires this storage directive.
    /// </summary>
    public string Path { get; }

    /// <summary>
    /// The kind of storage, e.g. "gdata".
    /// </summary>
    public string Kind { get; }

    /// <summary>
    /// Context required to construct a storage directive. Typically, this will
    /// be a project code.
    /// </summary>
    public string StorageIdentifier { get; }

    /// <summary>
    /// Construct a new PBS storage directive.
    /// </summary>
    /// <param name="path">Full path to a resource which requires this storage directive.</param>
    /// <param name="kind">The kind of storage, e.g. "gdata".</param>
    /// <param name="storageIdentifier">Context required to construct a storage directive. Typically, this will be a project code.</param>
    public PBSStorageDirective(string path, string kind, string storageIdentifier)
    {
        Path = path;
        Kind = kind;
        StorageIdentifier = storageIdentifier;
    }

    public override string ToString() => $"{Kind}/{StorageIdentifier}";
}

public static class PBSStorageHelper
{
    private static readonly Regex gdataRegex = new(@"^/g/data/([^/]+)/.*$");
    private static readonly Regex scratchRegex = new(@"^/scratch/([^/]+)/.*$");

    /// <summary>
    /// Gets the required PBS storage directives for a given path.
    /// </summary>
    /// <param name="path">Absolute path to check</param>
    /// <returns>PBS storage directive if path requires one, null otherwise</returns>
    /// <exception cref="ArgumentException">Thrown when a gdata path doesn't match expected format</exception>
    public static PBSStorageDirective? GetStorageDirective(string path)
    {
        // Check for gdata paths
        var gdataMatch = gdataRegex.Match(path);
        if (gdataMatch.Success)
        {
            var projectId = gdataMatch.Groups[1].Value;
            return new PBSStorageDirective(path, "gdata", projectId);
        }

        // Check for scratch paths
        var scratchMatch = scratchRegex.Match(path);
        if (scratchMatch.Success)
        {
            var projectId = scratchMatch.Groups[1].Value;
            return new PBSStorageDirective(path, "scratch", projectId);
        }

        return null;
    }

    /// <summary>
    /// Gets all required PBS storage directives for a set of paths.
    /// </summary>
    /// <param name="paths">Collection of paths to check</param>
    /// <returns>Collection of unique storage directives required</returns>
    public static IEnumerable<PBSStorageDirective> GetStorageDirectives(IEnumerable<string> paths)
    {
        return paths
            .Select(GetStorageDirective)
            .Where(directive => directive != null)
            .GroupBy(directive => directive!.StorageIdentifier)
            .Select(group => group.First()!);
    }

    /// <summary>
    /// Formats storage directives for inclusion in a PBS script.
    /// </summary>
    /// <param name="directives">Collection of storage directives</param>
    /// <returns>Formatted PBS storage directive line</returns>
    public static string FormatStorageDirectives(IEnumerable<PBSStorageDirective> directives)
    {
        var storageList = string.Join('+', directives.Select(d => d.ToString()));
        return $"#PBS -l storage={storageList}";
    }
}
