namespace ClimateProcessing.Models;

/// <summary>
/// Interpolation algorithm that can be used when remapping a variable.
/// </summary>
/// <remarks>
/// This is used by <see cref="ScriptGenerator.GetRemapOperator"/> to generate the
/// remap operator used in the script.
///
/// TODO: support other algorithms? E.g. bicubic, NN, weighted NN, etc.
/// </remarks>
public enum InterpolationAlgorithm
{
    /// <summary>
    /// Bilinear interpolation - for variables that don't need conservation.
    /// </summary>
    Bilinear,

    /// <summary>
    /// Conservative interpolation - for variables that need conservation of
    /// quantities.
    /// </summary>
    Conservative
}
