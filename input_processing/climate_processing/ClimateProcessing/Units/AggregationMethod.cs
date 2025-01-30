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
        AggregationMethod.Mean => timeStep == TimeStep.Daily ? "daymean" : "timmean",
        AggregationMethod.Sum => timeStep == TimeStep.Daily ? "daysum" : "timsum",
        AggregationMethod.Maximum => timeStep == TimeStep.Daily ? "daymax" : "timmax",
        AggregationMethod.Minimum => timeStep == TimeStep.Daily ? "daymin" : "timmin",
        _ => throw new ArgumentException($"Unknown aggregation method: {method}")
    };
}
