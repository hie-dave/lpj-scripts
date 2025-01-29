namespace ClimateProcessing.Extensions;

public static class StringExtensions
{
    /// <summary>
    /// Replaces the first occurrence of a specified string with another string.
    /// </summary>
    /// <param name="text">The string to search in.</param>
    /// <param name="search">The string to search for.</param>
    /// <param name="replace">The string to replace with.</param>
    /// <returns>A new string with the first occurrence replaced, or the original string if the search string is not found.</returns>
    /// <exception cref="ArgumentNullException">Thrown when text or search is null.</exception>
    public static string ReplaceFirst(this string text, string search, string replace)
    {
        if (text == null)
            throw new ArgumentNullException(nameof(text));
        if (search == null)
            throw new ArgumentNullException(nameof(search));

        if (string.IsNullOrEmpty(search))
            return text;

        int pos = text.IndexOf(search);
        if (pos < 0)
            return text;

        return text[..pos] + (replace ?? string.Empty) + text[(pos + search.Length)..];
    }
}
