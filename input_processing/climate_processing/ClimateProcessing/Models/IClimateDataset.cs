namespace ClimateProcessing.Models;

public enum ClimateVariable
{
    SpecificHumidity,  // huss, unitless
    Precipitation,     // pr, mm
    SurfacePressure,   // ps, Pa
    ShortwaveRadiation,// rsds, W m-2
    WindSpeed,         // sfcWind, m s-1
    Temperature        // tas, degC
}

public record VariableInfo(string Name, string Units);

public interface IClimateDataset
{
    string DatasetName { get; }

    /// <summary>
    /// Get the directory containing the input files for a given variable.
    /// </summary>
    /// <param name="variable">The variable.</param>
    string GetInputFilesDirectory(ClimateVariable variable);

    /// <summary>
    /// Get the input files for a given variable, sorted by date ascending.
    /// </summary>
    /// <param name="variable">The variable.</param>
    IEnumerable<string> GetInputFiles(ClimateVariable variable);
    VariableInfo GetVariableInfo(ClimateVariable variable);
    string GenerateOutputFileName(ClimateVariable variable);
    Dictionary<string, string> GetMetadata();
}
