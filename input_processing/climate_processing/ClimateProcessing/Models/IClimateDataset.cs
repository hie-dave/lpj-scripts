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
    IEnumerable<string> GetInputFiles();
    VariableInfo GetVariableInfo(ClimateVariable variable);
    string GetOutputFilePattern(ClimateVariable variable);
    Dictionary<string, string> GetMetadata();
}
