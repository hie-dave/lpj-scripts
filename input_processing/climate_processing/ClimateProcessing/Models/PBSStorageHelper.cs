using System.Text.RegularExpressions;

namespace ClimateProcessing.Models;

public class PBSStorageDirective
{
    public string Path { get; }
    public string StorageIdentifier { get; }

    public PBSStorageDirective(string path, string storageIdentifier)
    {
        Path = path;
        StorageIdentifier = storageIdentifier;
    }

    public override string ToString() => $"gdata/{StorageIdentifier}";
}

public static class PBSStorageHelper
{
    private static readonly Regex GDataPathRegex = new(@"^/g/data/([^/]+)/.*$");
    private static readonly Regex ScratchPathRegex = new(@"^/scratch/([^/]+)/.*$");

    /// <summary>
    /// Gets the required PBS storage directives for a given path.
    /// </summary>
    /// <param name="path">Absolute path to check</param>
    /// <returns>PBS storage directive if path requires one, null otherwise</returns>
    /// <exception cref="ArgumentException">Thrown when a gdata path doesn't match expected format</exception>
    public static PBSStorageDirective? GetStorageDirective(string path)
    {
        // Check for gdata paths
        var gdataMatch = GDataPathRegex.Match(path);
        if (gdataMatch.Success)
        {
            var projectId = gdataMatch.Groups[1].Value;
            return new PBSStorageDirective(path, projectId);
        }

        // Check for scratch paths
        var scratchMatch = ScratchPathRegex.Match(path);
        if (scratchMatch.Success)
        {
            var projectId = scratchMatch.Groups[1].Value;
            return new PBSStorageDirective(path, projectId);
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
