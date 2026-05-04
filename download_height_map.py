#!/usr/bin/env python3
from owslib.wcs import WebCoverageService

OGC_URL = "https://ows.csiro.easi-eo.solutions/"

# Connect to WCS (try 2.0.1 first; if issues, try version="1.0.0")
wcs = WebCoverageService(OGC_URL, version="1.0.0")

# List coverages
for cov_id, meta in list(wcs.contents.items())[:20]:
    print(cov_id, "|", meta.title)

# Choose a coverage that represents vegetation height
coverage_id = "oztreemap_chm_best_pick"

# Request coverage; supply bbox, CRS, and resolution or size
resp = wcs.getCoverage(
    identifier=coverage_id,
    bbox=(89.995, 60.005, 180.005, 5.005),
    crs="EPSG:4326",
    resx=0.01,    # or size=(width,height)
    resy=0.01,
    format="GeoTIFF",
)
with open("veg_height.tif", "wb") as f:
    f.write(resp.read())
