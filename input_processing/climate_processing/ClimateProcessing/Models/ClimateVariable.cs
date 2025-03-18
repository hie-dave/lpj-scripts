namespace ClimateProcessing.Models;

public enum ClimateVariable
{
    SpecificHumidity,  // huss, unitless
    Precipitation,     // pr, mm
    SurfacePressure,   // ps, Pa
    ShortwaveRadiation,// rsds, W m-2
    WindSpeed,         // sfcWind, m s-1
    Temperature,       // tas, degC
    MaxTemperature,    // tasmax, degC
    MinTemperature     // tasmin, degC
}
