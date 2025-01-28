namespace ClimateProcessing.Models;

public enum VPDMethod
{
    // Default method from Magnus equation
    Magnus,
    
    // From Buck (1981): New equations for computing vapor pressure and enhancement factor
    Buck1981,
    
    // From Alduchov and Eskridge (1996): Improved Magnus form approximation of saturation vapor pressure
    AlduchovEskridge1996,
    
    // From Allen et al. (1998): FAO Irrigation and drainage paper 56
    AllenFAO1998,
    
    // From Sonntag (1990): Important new values of the physical constants of 1986, vapor pressure formulations
    Sonntag1990
}
