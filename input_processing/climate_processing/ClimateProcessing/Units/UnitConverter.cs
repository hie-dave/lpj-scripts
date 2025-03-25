
using System.Diagnostics.CodeAnalysis;

namespace ClimateProcessing.Units;

public static class UnitConverter
{
    private static readonly Dictionary<string, HashSet<string>> UnitSynonyms = new()
    {
        ["1"] = new() { "kg/kg", "kg kg-1", "1" },
        ["W/m2"] = new() { "W/m2", "W/m^2", "W m-2" },
        ["degC"] = new() { "degC", "°C", "Celsius" },
        ["K"] = new() { "K", "Kelvin" },
        ["mm"] = new() { "mm", "kg m-2" }  // 1mm of water = 1 kg/m2
    };

    private record struct ConversionDefinition(
        string Expression,
        bool RequiresTimeStep
    );

    private static readonly Dictionary<(string From, string To), ConversionDefinition> ConversionExpressions = new()
    {
        [("K", "degC")] = new("-subc,273.15", false),
        [("kg m-2 s-1", "mm")] = new("-mulc,{timestep}", true),  // Multiply by seconds in period to get accumulation
        [("kPa", "Pa")] = new("-mulc,1000", false)
    };

    public record ConversionResult(
        bool RequiresConversion,
        bool RequiresRenaming,
        bool RequiresTimeStep,
        string? ConversionExpression = null
    );

    public static ConversionResult AnalyseConversion(string inputUnits, string targetUnits)
    {
        // Check if units are exactly the same (including notation).
        if (inputUnits == targetUnits)
            return new ConversionResult(false, false, false);

        // Normalise both units to their canonical form
        string normalisedInput = NormaliseUnits(inputUnits);
        string normalisedTarget = NormaliseUnits(targetUnits);

        // Check if units are equivalent (different notation but same meaning).
        if (AreUnitsEquivalent(normalisedInput, normalisedTarget))
            return new ConversionResult(false, true, false);

        // Check if we have a conversion expression.
        var conversionKey = (normalisedInput, normalisedTarget);
        if (ConversionExpressions.TryGetValue(conversionKey, out var conversion))
            return new ConversionResult(true, true, conversion.RequiresTimeStep, conversion.Expression);

        throw new ArgumentException($"Unsupported unit conversion from {inputUnits} to {targetUnits}");
    }

    /// <summary>
    /// Normalise a unit string to a canonical form.
    /// </summary>
    /// <param name="units">The unit string to normalise.</param>
    /// <returns>The normalised unit string.</returns>
    internal static string NormaliseUnits(string units)
    {
        // TODO: support spaces:  umol / m2 == umol/m2
        // TODO: support periods: kg.m-2 == kg m-2
        // TODO: support carets:  m2/m2 == m^2/m^2
        // TODO: support slashes: W/m2 == W m-2

        // If no exact match, return as is.
        return units;
    }

    internal static bool AreUnitsEquivalent(string units1, string units2)
    {
        // First check if they're the same after normalisation.
        if (units1 == units2)
            return true;

        // Then check if they belong to the same synonym group.
        foreach (var synonyms in UnitSynonyms.Values)
            if (synonyms.Contains(units1) && synonyms.Contains(units2))
                return true;

        return false;
    }

    public static string GenerateConversionExpression(
        string inputUnits,
        string targetUnits,
        TimeStep timeStep)
    {
        ConversionResult result = AnalyseConversion(inputUnits, targetUnits);

        // Should never happen - this function is only called if a conversion is
        // required.
        if (!result.RequiresConversion)
            return string.Empty;

        // Conversion expression cannot be null here. AnalyseConversion can only
        // return RequiresConversion=true when the conversion expression is
        // non-null.
        string expression = result.ConversionExpression!;
        if (result.RequiresTimeStep)
            expression = expression.Replace("{timestep}", timeStep.GetSecondsInPeriod().ToString());

        return expression;
    }
}
