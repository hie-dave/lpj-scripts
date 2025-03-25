using ClimateProcessing.Models;
using ClimateProcessing.Services;

namespace ClimateProcessing.Tests.Mocks;

/// <summary>
/// Concrete implementation of ProcessingConfig for testing
/// </summary>
internal class TestProcessingConfig : ProcessingConfig
{
    public override IEnumerable<NarClim2Dataset> CreateDatasets()
    {
        throw new NotImplementedException();
    }

    public override ScriptGenerator CreateScriptGenerator()
    {
        throw new NotImplementedException();
    }
}
