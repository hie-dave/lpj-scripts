#
# Helper script providing utility methods for normalising ozflux site names.
#
from ozflux_logging import *

_sites_without_codes = [
    "Arcturus",
    "Longreach",
    "Loxton",
    "RedDirtMelonFarm"
]

_synonyms = {
    "FletcherviewTropicalRangeland": "FletcherView",
    "Wallaby": "WallabyCreek",
}

site_codes = {
    "AU-Ade": "AdelaideRiver",
    "AU-APL": "AlpinePeatland",
    "AU-ASM": "AliceSpringsMulga",
    "AU-Boy": "Boyagin",
    "AU-Cm1": "ColeamballyMaizeWheat",
    "AU-Cm2": "ColeamballyRice",
    "AU-Col": "Collie",
    "AU-Cow": "CowBay",
    "AU-Cpr": "Calperum",
    "AU-Ctr": "CapeTribulation",
    "AU-Cum": "CumberlandPlain",
    "Au-CuP": "CumberlandPlain",
    "AU-DaS": "DalyUncleared",
    "AU-DaP": "DalyPasture",
    "AU-Drg": "Dargo",
    "AU-Dry": "DryRiver",
    "Au-Fle": "FletcherView",
    "AU-Fog": "FoggDam",
    "AU-Gin": "Gingin",
    "AU-GWW": "GreatWesternWoodlands",
    "AU-How": "HowardSprings",
    "AU-Lit": "Litchfield",
    "AU-Nim": "Nimmo",
    "NZ-Oxf": "Oxford",
    "AU-Otw": "Otway",
    "AU-Rgf": "Ridgefield",
    "AU-Rig": "RiggsCreek",
    "AU-Rob": "RobsonCreek",
    "AU-Sam": "Samford",
    "AU-Sil": "SilverPlains",
    "AU-Stp": "SturtPlains",
    "AU-TTE": "TiTreeEast",
    "AU-Tum": "Tumbarumba",
    "AU-Vir": "VirginiaPark",
    "AU-Wac": "WallabyCreek",
    "AU-Whr": "Whroo",
    "AU-Wom": "WombatStateForest",
    "AU-Wrr": "Warra",
    "AU-Ync": "Yanco",
    "AU-Yar": "Yarramundi"
}

_lookup_table = {}

def _pascal(word: str) -> str:
    """
    Captialise the first letter of a word.
    """
    return word[0].upper() + word[1:]

def normalise_site_name(site_name: str) -> str:
    """
    Normalise a site name.
    """
    if site_name in _lookup_table:
        normalised = _lookup_table[site_name]
        # Most lookups will be cached, so use debug verbosity.
        log_debug(f"Using cached normalised site name {site_name}: {normalised}")
        return normalised

    if site_name in site_codes:
        log_diagnostic(f"Identified site code {site_name}: {site_codes[site_name]}")
        _lookup_table[site_name] = site_codes[site_name]
        return site_codes[site_name]

    # Remove whitespace, convert to pascal case.
    normalised = "".join(_pascal(word) for word in site_name.split())
    log_debug(f"Normalised site name {site_name} to {normalised}")

    # Check if normalised site name is a value in the site codes dictionary.
    if normalised in site_codes.values():
        log_diagnostic(f"Identified site code {site_name}: {normalised}")
        _lookup_table[site_name] = normalised
        return normalised

    # Handle sites without codes.
    if normalised in _sites_without_codes:
        log_diagnostic(f"Identified site name {site_name}: {normalised}")
        _lookup_table[site_name] = normalised
        return normalised

    # Handle synonyms/common alternative names.
    if normalised in _synonyms:
        log_diagnostic(f"Identified site name {site_name}: {normalised}")
        _lookup_table[site_name] = normalised
        return normalised

    # Handle site code with fluxnet country prefix.
    code = f"AU-{normalised}"
    if code in site_codes:
        result = site_codes[code]
        log_diagnostic(f"Identified site code {site_name}: {result}")
        _lookup_table[site_name] = result
        return result

    # Log a warning and return original site code.
    log_warning(f"Unknown site name: {site_name}")
    return site_name
