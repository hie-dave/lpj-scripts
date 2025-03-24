# TODO

## Refactor ScriptGenerator

### Extract an IFileWriter interface

- ScriptGenerator is passed an IFileWriter instance
- The writer has methods like WriteLineAsync(), and replaces all of the manual
  TextWriter usage in ScriptGenerator
- Makes the unit tests cleaner by avoiding file IO
- Makes ScriptGenerator more flexible by removing filesystem interactions

### ClimateVariableManager

Handle climate variable configurations

- Current static dictionaries (daveVariables, trunkVariables)
- Methods like GetVariables, GetTargetUnits, GetAggregationMethod

### VPDCalculator

Handle VPD-specific calculations and script generation

- Methods like WriteVPDEquationsAsync, GenerateVPDScript
- VPD configuration and methods could go here

### CDOCommandGenerator

Handle CDO command generation

- Methods like GenerateRenameOperator, GenerateUnitConversionOperators
- CDO-specific constants and configurations

## Additional Unit Tests

- All of the above
