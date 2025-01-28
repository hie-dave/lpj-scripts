namespace ClimateProcessing.Services;

public record ProcessingCommand(
    string Command,
    bool RequiresProcessing
)
{
    public static ProcessingCommand Skip(string inputFile, string outputFile) =>
        new($"ln -sf {inputFile} {outputFile}", false);

    public static ProcessingCommand Process(string command) =>
        new(command, true);
}
