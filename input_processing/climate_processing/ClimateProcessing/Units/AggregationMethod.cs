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
    public static string ToCdoOperator(this AggregationMethod method) => method switch
    {
        AggregationMethod.Mean => "timemean",
        AggregationMethod.Sum => "timesum",
        AggregationMethod.Maximum => "timemax",
        AggregationMethod.Minimum => "timemin",
        _ => throw new ArgumentException($"Unknown aggregation method: {method}")
    };
}
