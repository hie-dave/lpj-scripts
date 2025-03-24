using System;

namespace ClimateProcessing.Models;

/// <summary>
/// Represents a PBS walltime value, which must be between 0 and 48 hours.
/// </summary>
public readonly struct PBSWalltime : IEquatable<PBSWalltime>
{
    /// <summary>
    /// The maximum allowed walltime value (48 hours).
    /// </summary>
    public static readonly PBSWalltime MaxValue = new(48, 0, 0);
    
    /// <summary>
    /// The hours component of the walltime.
    /// </summary>
    public int Hours { get; }

    /// <summary>
    /// The minutes component of the walltime (0-59).
    /// </summary>
    public int Minutes { get; }

    /// <summary>
    /// The seconds component of the walltime (0-59).
    /// </summary>
    public int Seconds { get; }

    /// <summary>
    /// Creates a new PBS walltime value.
    /// </summary>
    /// <param name="hours">The hours component (0-48).</param>
    /// <param name="minutes">The minutes component (0-59).</param>
    /// <param name="seconds">The seconds component (0-59).</param>
    /// <exception cref="ArgumentException">Thrown when the components are invalid or exceed PBS limits.</exception>
    public PBSWalltime(int hours, int minutes, int seconds)
    {
        if (hours < 0 || minutes < 0 || minutes >= 60 || seconds < 0 || seconds >= 60)
            throw new ArgumentException("Invalid time components");
            
        if (hours > 48 || (hours == 48 && (minutes > 0 || seconds > 0)))
            throw new ArgumentException("PBS walltime cannot exceed 48 hours");

        Hours = hours;
        Minutes = minutes;
        Seconds = seconds;
    }

    /// <summary>
    /// Parses a PBS walltime string in the format "HH:MM:SS".
    /// </summary>
    /// <param name="walltime">The walltime string to parse.</param>
    /// <returns>A new PBSWalltime instance.</returns>
    /// <exception cref="FormatException">Thrown when the string format is invalid.</exception>
    public static PBSWalltime Parse(string walltime)
    {
        string[] parts = walltime.Split(':');
        if (parts.Length != 3)
            throw new FormatException("Walltime must be in format 'HH:MM:SS'");

        // Verify each part is exactly 2 digits
        if (parts[0].Length != 2 || parts[1].Length != 2 || parts[2].Length != 2)
            throw new FormatException("Each time component must be exactly 2 digits");

        if (!int.TryParse(parts[0], out int hours) ||
            !int.TryParse(parts[1], out int minutes) ||
            !int.TryParse(parts[2], out int seconds))
            throw new FormatException("Invalid time components");

        return new PBSWalltime(hours, minutes, seconds);
    }

    /// <summary>
    /// Returns a string representation of the walltime in "HH:MM:SS" format.
    /// </summary>
    public override string ToString() => $"{Hours:00}:{Minutes:00}:{Seconds:00}";

    /// <summary>
    /// Converts the PBSWalltime to a TimeSpan.
    /// </summary>
    public static implicit operator TimeSpan(PBSWalltime walltime) =>
        new(walltime.Hours, walltime.Minutes, walltime.Seconds);
        
    /// <summary>
    /// Attempts to convert a TimeSpan to a PBSWalltime.
    /// </summary>
    /// <exception cref="ArgumentException">Thrown when the TimeSpan exceeds PBS limits.</exception>
    public static explicit operator PBSWalltime(TimeSpan timespan)
    {
        int totalMinutes = (int)timespan.TotalMinutes;
        int hours = totalMinutes / 60;
        int minutes = totalMinutes % 60;
        return new PBSWalltime(hours, minutes, timespan.Seconds);
    }

    public bool Equals(PBSWalltime other) =>
        Hours == other.Hours && Minutes == other.Minutes && Seconds == other.Seconds;

    public override bool Equals(object? obj) =>
        obj is PBSWalltime other && Equals(other);

    public override int GetHashCode() =>
        HashCode.Combine(Hours, Minutes, Seconds);

    public static bool operator ==(PBSWalltime left, PBSWalltime right) =>
        left.Equals(right);

    public static bool operator !=(PBSWalltime left, PBSWalltime right) =>
        !left.Equals(right);
}
