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
