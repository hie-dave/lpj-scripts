namespace ClimateProcessing.Units;

public record struct TimeStep
{
    private readonly int hours;

    public int Hours
    {
        get => hours;
        init
        {
            if (value <= 0 || value > 24)
                throw new ArgumentException("TimeStep must be between 1 and 24 hours", nameof(value));
            if (24 % value != 0)
                throw new ArgumentException("TimeStep must divide evenly into 24 hours", nameof(value));
            hours = value;
        }
    }

    public TimeStep(int hours) => Hours = hours;

    public static TimeStep Hourly => new(1);
    public static TimeStep ThreeHourly => new(3);
    public static TimeStep Daily => new(24);

    public int GetSecondsInPeriod() => Hours * 3600; // hours * seconds per hour
}
