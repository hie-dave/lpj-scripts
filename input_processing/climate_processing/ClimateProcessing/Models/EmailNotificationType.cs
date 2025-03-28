namespace ClimateProcessing.Models;

/// <summary>
/// Types of email notifications that can be sent for PBS jobs.
/// </summary>
/// <remarks>
/// PBS also supports sending mail for subjobs, but we don't use subjobs, so
/// this option is not supported here.
/// </remarks>
[Flags]
public enum EmailNotificationType
{
    /// <summary>
    /// No mail is sent.
    /// </summary>
    None = 1,

    /// <summary>
    /// Mail is sent when the job terminates.
    /// </summary>
    After = 2,

    /// <summary>
    /// Mail is sent when the job begins execution.
    /// </summary>
    Before = 4,

    /// <summary>
    /// Mail is sent when the job is aborted by the batch system.
    /// </summary>
    Aborted = 8
}
