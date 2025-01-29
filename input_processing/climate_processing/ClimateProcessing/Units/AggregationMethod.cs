namespace ClimateProcessing.Units;

public enum AggregationMethod
{
    Mean,       // For most instantaneous measurements (temperature, humidity, etc.)
    Sum,        // For accumulation variables (precipitation)
    Maximum,    // For extremes that might be needed
    Minimum     // For extremes that might be needed
}

public static class AggregationMethodExtensions
{
    public static string ToCdoOperator(this AggregationMethod method, TimeStep timeStep) => method switch
    {
        AggregationMethod.Mean => timeStep == TimeStep.Daily ? "daymean" : "timemean",
        AggregationMethod.Sum => timeStep == TimeStep.Daily ? "daysum" : "timesum",
        AggregationMethod.Maximum => timeStep == TimeStep.Daily ? "daymax" : "timemax",
        AggregationMethod.Minimum => timeStep == TimeStep.Daily ? "daymin" : "timemin",
        _ => throw new ArgumentException($"Unknown aggregation method: {method}")
    };
}
