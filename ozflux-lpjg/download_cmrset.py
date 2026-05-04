# https://github.com/ternaustralia/cmrset-examples/blob/main/tds-examples/python/download-aet.py
#!/usr/bin/env python3

# Dataset status (up/down), and important notes for running this script:
# https://github.com/ternaustralia/cmrset-examples/tree/main/tds-examples

from argparse import ArgumentParser
from typing import Optional

import os
import logging
import tempfile
import requests
import defusedxml.ElementTree as ET
from enum import Enum, auto
from operator import itemgetter
from datetime import date, datetime
from dateutil.relativedelta import relativedelta
from sys import argv

class Options:
    """
    Options for the download_cmrset script.

    Attributes:
        api_key (str): TERN API key.
        path_out (str): Output path for the downloaded files.
        update_method (str): Update method for the downloaded files.
        product_code (str): Product code for the CMRSET service ("CMRSET_LANDSAT_V2_2" recommended).
        start (str): Start date for the downloaded files ("YYYY-MM-DD").
        end (str): End date for the downloaded files ("YYYY-MM-DD").
        bands (list[str]): Bands to download ("ETa", "pixel_qa").
        tiles (list[int]): Tiles to download (0-11).
        dryrun (bool): Whether to run in dry run mode.
    """
    def __init__(self, api_key: str, path_out: str, update_method: str,
                 product_code: str, start: str, end: str, bands: list[str],
                 tiles: list[int], dryrun: bool):
        self.api_key = api_key
        self.path_out = path_out
        self.update_method = update_method
        self.product_code = product_code
        self.start = start
        self.end = end
        self.dryrun = dryrun

        self.bands = bands
        if not isinstance(bands, list):
            self.bands = [band.strip() for band in bands.split(",")]

        self.tiles = tiles
        if not isinstance(tiles, list):
            self.tiles = [int(tile.strip()) for tile in tiles.split(",")]

def get_api_key(api_key: Optional[str], api_key_env: Optional[str],
                api_key_file: Optional[str]) -> str:
    """
    Get the API key from the provided arguments.
    """
    if api_key is not None:
        return api_key
    if api_key_env is not None:
        return os.getenv(api_key_env)
    if api_key_file is not None:
        with open(api_key_file, "r") as f:
            return f.read().strip()
    raise ValueError("No API key provided")

def parse_args(args: list[str]) -> Options:
    """
    Parse command line arguments.
    """
    parser = ArgumentParser(description="Download CMRSET data.")
    parser.add_argument("--path-out", "-o", required = True, help = "Output path")
    parser.add_argument("--update-method", "-u", required = True, help = "Update method")
    parser.add_argument("--product-code", "-p", required = True, help = "Product code")
    parser.add_argument("--start", "-s", required = True, help = "Start date")
    parser.add_argument("--end", "-E", required = True, help = "End date")
    parser.add_argument("--bands", "-b", required = True, help = "Bands")
    parser.add_argument("--tiles", "-t", required = True, help = "Tiles")
    parser.add_argument("--dryrun", "-d", action = "store_true", help = "Dry run")
    # Multiple methods of providing API key: directly as CLI arg, env var, or
    # file.
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--api-key", "-k", required = True, help = "API key")
    group.add_argument("--api-key-env", "-e", required = True, help = "API key environment variable")
    group.add_argument("--api-key-file", "-f", required = True, help = "API key file")

    p = parser.parse_args(args)

    api_key = get_api_key(p.api_key, p.api_key_env, p.api_key_file)
    return Options(api_key, p.path_out, p.update_method, p.product_code,
                   p.start, p.end, p.bands, p.tiles, p.dryrun)

# Lookup for available products.
ProductCodes = {
    "CMRSET_LANDSAT_V2_2": {
        "Url": "https://data.tern.org.au/landscapes/aet/v2_2",
        "Start": "2000-02-01",
        "End": None
    },
    "CMRSET_LANDSAT_V2_1": {
        "Url": "https://data.tern.org.au/landscapes/aet/v2_1",
        "Start": "2012-02-01",
        "End": "2021-02-01" # Discontinued
    }
}

# Lookup for tile indicies.
TileLookup = {
    0: "0000000000-0000000000",
    1: "0000000000-0000043776",
    2: "0000000000-0000087552",
    3: "0000000000-0000131328",
    4: "0000043776-0000000000",
    5: "0000043776-0000043776",
    6: "0000043776-0000087552",
    7: "0000043776-0000131328",
    8: "0000087552-0000000000",
    9: "0000087552-0000043776",
    10: "0000087552-0000087552",
    11: "0000087552-0000131328"
}

class UpdateMethod(Enum):
    """
    An enum for the various processing methods.
    """
    # Update missing files within the local archive.
    UPDATE_MISSING = auto()

    # Update missing/outdated files within the local archive.
    UPDATE_NEW = auto()

    # Update all files within the local archive.
    UPDATE_ALL = auto()

def get_months(start, end):
    """
    Get a monthly array of dates between start and end.
    """
    # Ensure dates are 1st of the month.
    start = start.replace(day=1)
    end = end.replace(day=1)

    date = start
    array = []
    while(date <= end):
        array.append(date)
        date = date + relativedelta(months=1)

    return array

def get_vrt_relative_paths(product_code, dates, bands):
    """
    Get the relative paths for the VRT files based upon product, date range and bands.
    """
    date_hash = {}
    for date in dates:
        date_hash[date] = {}
        for band in bands:
            date_str = date.strftime("%Y_%m_%d")
            date_hash[date][band] = f"/{date.year}/{date_str}/{product_code}_{date_str}_{band}.vrt"
    return date_hash

def get_vrt_sources(file):
    """
    Get all the sources referenced within a VRT file.
    """
    # Go to beginning of file.
    file.seek(0)
    xml_doc = file.read()

    # Use wildcard for different source types.
    nodes = ET.fromstring(xml_doc).findall(".//VRTRasterBand/*/SourceFilename")

    # Sort the node values after searching the nodes.
    files = sorted(list(map(lambda node: node.text, nodes)))

    return files

def download_file(url, out_file, dryrun=False):
    """
    Download a file via a session.
    """

    logging.info("Downloading: {url}".format(url=url))
    if dryrun == False:
        # Session accessed from global scope.
        response = Session.get(url, stream = True)
        # Trigger exception for unacceptable status codes.
        response.raise_for_status()

        is_str = isinstance(out_file, str)
        if is_str:
            os.makedirs(os.path.dirname(out_file), exist_ok=True)
            out_file = open(out_file,"wb")

        # Read large files in 10 MiB chunks.
        for chunk in response.iter_content(chunk_size=1024 * 1024 * 10):
            out_file.write(chunk)

        if is_str:
            out_file.close()

        return response

def confirm_download(url, out_file, update_method):
    """
    Determines whether a download should take place based upon the UpdateMethod.
    """
    def update_missing():
        """
        Update missing files within the local archive.
        """
        result = not os.path.exists(out_file)
        if not result: logging.info("Skipping existing file: {url}".format(url=url))
        return result

    def undate_new():
        """
        Update missing/outdated files within the local archive.
        """
        # Return true if the file is missing.
        if not os.path.exists(out_file):
            return True

        dt = datetime.utcfromtimestamp((os.path.getctime(out_file)))

        # File creation date in GMT, for headers.
        date_str = dt.strftime("%a, %d %b %Y %H:%M:%S GMT")
        headers = {'If-Modified-Since': date_str}
        response = Session.head(url, headers=headers, allow_redirects = True)

        # Check if response is a HTTP 304 Not Modified status code.
        result = response.status_code != 304
        if not result:
            logging.info("Skipping up to date file: {url}".format(url=url))
        return result

    def update_all():
        """
        Update all files within the local archive.
        """
        return True

    confirmation = {
        UpdateMethod.UPDATE_MISSING : update_missing,
        UpdateMethod.UPDATE_NEW : undate_new,
        UpdateMethod.UPDATE_ALL : update_all,
    }

    return confirmation[update_method]()

def download_images(base_url, base_dir, relative_paths,
                    tile_ids=list(range(0, 12)),
                    update_method=UpdateMethod.UPDATE_MISSING,
                    dryrun=False):
    """
    Downloads the image tiles for the specified paths.
    """
    for date in relative_paths:
        date_str = date.strftime("%Y-%m-%d")
        nband = len(relative_paths[date])
        logging.info(f"Processing {nband} band(s) for {date_str}...")
        for band in relative_paths[date]:
            logging.info(f"Processing {band} for {date_str}...")

            # Download VRT file that contains references to the files it
            # mosaics.
            try:
                vrt_url = f"{base_url}{relative_paths[date][band]}"

                # Create temporary file for VRT.
                # TODO: replace with "with" construct?
                vrt_file = tempfile.TemporaryFile()

                # Download VRT contents to the temporary file.
                download_file(vrt_url, vrt_file)

                # Read the source files referenced within the VRT.
                files = get_vrt_sources(vrt_file)
            except Exception as error:
                logging.error(error)
                continue
            finally:
                # Delete the temporary file.
                vrt_file.close()

            # Filter the tiles to those specified.
            filtered_files = list(filter(lambda file:
                any(tile_id in file for tile_id in tile_ids), files))

            # Download all the tiles for the filtered space/time/variable
            # parameters.
            n = len(filtered_files)
            logging.info(f"Processing {n} tile(s) for {band} on {date_str}...")
            for file in filtered_files:
                tile_url = f"{base_url}/{date.year}/{date_str}/{file}"
                out_file = f"{base_dir}/{date.year}/{date_str}/{file}"

                # Test whether a file should be downloaded, and do so if True.
                if confirm_download(tile_url, out_file, update_method):
                    try:
                        download_file(tile_url, out_file, dryrun = dryrun)
                    except Exception as error:
                        logging.error(error)
                        continue

def main(opts: Options):
    """
    Main function.
    """
    logging.info("DRYRUN: {dryrun}".format(dryrun=opts.dryrun))

    # A session which contains common settings which will be used for all web requests made.
    # In particular, an X-API-Key auth from a base64 encoded key.
    Session = requests.Session()
    Session.headers.update({"X-API-Key": opts.api_key})

    # Constrain start/end to within dataset temporal bounds.
    product_start = date.fromisoformat(ProductCodes[opts.product_code]["Start"])
    product_end = date.fromisoformat(ProductCodes[opts.product_code]["End"])
    if product_end is None:
        product_end = date.today()

    start = max(date.fromisoformat(opts.start), product_start)
    end = min(date.fromisoformat(opts.end), product_end)

    logging.info("Start: {start}".format(start=start))
    logging.info("End: {end}".format(end=end))

    # Generate the list of dates to download.
    dates = get_months(start, end)
    logging.info("Processing data for the following dates:")
    for d in dates: logging.info(d.strftime("%b %Y"))

    # Get the relative paths for each the VRT files for each date.
    vrt_relative_paths = get_vrt_relative_paths(opts.product_code, dates,
                                                opts.bands)

    # Download all the tiles referenced in each VRT file.
    download_images(ProductCodes[opts.product_code]["Url"],
                    f"{opts.path_out}/{opts.product_code}",
                    vrt_relative_paths, itemgetter(*opts.tiles)(TileLookup),
                    update_method = UpdateMethod[opts.update_method],
                    dryrun = opts.dryrun)

    logging.info("Processing complete")

# Run the script
if __name__=="__main__":
    # TODO: make this configurable?
    logging.getLogger().setLevel(logging.INFO)
    opts = parse_args(argv[1:])
    main(opts)
