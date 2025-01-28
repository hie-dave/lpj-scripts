using System.Text.RegularExpressions;

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
        [("K", "degC")] = new("${output}=${input}-273.15", false),
        [("kg m-2 s-1", "mm")] = new("${output}=${input}*${timestep}", true)  // Multiply by seconds in period to get accumulation
    };

    public record ConversionResult(
        bool RequiresConversion,
        bool RequiresRenaming,
        bool RequiresTimeStep,
        string? ConversionExpression = null
    );

    public static ConversionResult AnalyzeConversion(string inputUnits, string targetUnits)
    {
        // Normalize both units to their canonical form
        var normalizedInput = NormalizeUnits(inputUnits);
        var normalizedTarget = NormalizeUnits(targetUnits);

        // Check if units are exactly the same (including notation)
        if (inputUnits == targetUnits)
        {
            return new ConversionResult(false, false, false);
        }

        // Check if units are equivalent (different notation but same meaning)
        if (AreUnitsEquivalent(normalizedInput, normalizedTarget))
        {
            return new ConversionResult(false, true, false);
        }

        // Check if we have a conversion expression
        var conversionKey = (normalizedInput, normalizedTarget);
        if (ConversionExpressions.TryGetValue(conversionKey, out var conversion))
        {
            return new ConversionResult(true, true, conversion.RequiresTimeStep, conversion.Expression);
        }

        throw new ArgumentException($"Unsupported unit conversion from {inputUnits} to {targetUnits}");
    }

    private static string NormalizeUnits(string units)
    {
        // First try to find an exact match in our synonyms
        foreach (var (canonical, synonyms) in UnitSynonyms)
        {
            if (synonyms.Contains(units))
            {
                return canonical;
            }
        }

        // If no exact match, return as is
        return units;
    }

    private static bool AreUnitsEquivalent(string units1, string units2)
    {
        // First check if they're the same after normalization
        if (units1 == units2) return true;

        // Then check if they belong to the same synonym group
        foreach (var synonyms in UnitSynonyms.Values)
        {
            if (synonyms.Contains(units1) && synonyms.Contains(units2))
            {
                return true;
            }
        }

        return false;
    }

    public static string GenerateConversionExpression(
        string inputVar,
        string outputVar,
        string inputUnits,
        string targetUnits,
        TimeStep? timeStep = null)
    {
        var result = AnalyzeConversion(inputUnits, targetUnits);
        
        if (!result.RequiresConversion)
        {
            return $"{outputVar}={inputVar}";
        }

        if (result.ConversionExpression == null)
        {
            throw new InvalidOperationException($"No conversion expression available for {inputUnits} to {targetUnits}");
        }

        if (result.RequiresTimeStep && !timeStep.HasValue)
        {
            throw new ArgumentException($"TimeStep is required for conversion from {inputUnits} to {targetUnits}");
        }

        var expression = result.ConversionExpression
            .Replace("${input}", inputVar)
            .Replace("${output}", outputVar);

        if (result.RequiresTimeStep && timeStep.HasValue)
        {
            expression = expression.Replace("${timestep}", timeStep.Value.GetSecondsInPeriod().ToString());
        }

        return expression;
    }
}
