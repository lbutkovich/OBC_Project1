#!/usr/bin/env python3
"""
download_fna_files.py
Created 7/6/25

This script uses the MG-RAST API to download a subset of .fna files from the MG-RAST database for 16S rRNA sequencing from David et al. (2014). 
- Paper: https://www.nature.com/articles/nature12820 
- Accession ID: mgp6248
- ID### (plant) and DD### (animal) 16S rRNA sequencing samples for metagenomic analysis.
- MG-RAST API documentation: https://help.mg-rast.org/api.html

# Takes about 40 min. to download 236 target .fna 16S rRNA files from MG-RAST for Accession ID mgp6248.
"""

from __future__ import annotations
from pathlib import Path
from typing import List, Dict, Optional
import requests, sys
import re
# Record time to download each file
import time
start = time.time()


"""
Functions
"""
def get_json(url: str, auth: Optional[str] = None) -> Dict:
    """Fetches MG-RAST API JSON via HTTP GET.
    
    Args:
        url (str): The URL to fetch.
        auth (Optional[str]): Optional authentication key for the request.
    Raises:
        requests.HTTPError: If the HTTP request fails.
    Returns:
        Dict: Parsed JSON response from the API.
    """
    hdr = {"auth": auth} if auth else {}
    r = requests.get(url, headers=hdr, timeout=90)
    r.raise_for_status()
    return r.json()

def list_metagenomes(project_id: str, auth: Optional[str] = None) -> List[dict]:
    """Returns metagenome information [{'id': 'mgm1234.3', 'name': 'DD17'}, …] for every metagenome
    in the project_id.

    Args:
        project_id (str): The ID of the project to list metagenomes for.
        auth (Optional[str]): Optional authentication key for the request.
    Returns:
        List[dict]: A list of dictionaries containing metagenome information.
    """
    url  = f"{BASE_URL}/project/{project_id}?verbosity=full"
    info = get_json(url, auth)
    return [{"id": m["metagenome_id"], "name": m["name"]}
            for m in info["metagenomes"]]

def list_files(mgm_id: str, auth: Optional[str] = None) -> List[Dict]:
    """Generates a catalogue of downloadable artefacts for a metagenome.
    
    Args:
        mgm_id (str): The ID of the metagenome to list files for.
        auth (Optional[str]): Optional authentication key for the request.
    Returns:
        List[Dict]: A list of dictionaries containing file information."""
    return get_json(f"{BASE_URL}/download/{mgm_id}", auth)["data"]

def pick_fna(files: List[Dict], wanted_id: str) -> Optional[Dict]:
    """Finds the .fna file.
    
    Returns the entry that matches wanted_id and ends with '.fna'.
    
    Args:
        files (List[Dict]): List of file dictionaries to search.
        wanted_id (str): The file ID to match.
    Returns:
        Optional[Dict]: The matching file entry or None if not found.
    """
    return next(
        (f for f in files if f["file_id"] == wanted_id and f["file_name"].endswith(".fna")),
        None
    )

def stream_download(url: str, dest: Path, auth: Optional[str] = None) -> None:
    """Streams large files in smaller chunks to keep RAM low.
    
    Args:
        url (str): The URL to download the file from.
        dest (Path): The destination path to save the downloaded file.
        auth (Optional[str]): Optional authentication key for the request.
    Raises:
        requests.HTTPError: If the HTTP request fails.
    Returns:
        None
    """
    hdr = {"auth": auth} if auth else {}
    with requests.get(url, headers=hdr, stream=True, timeout=300) as r:
        r.raise_for_status()
        dest.parent.mkdir(parents=True, exist_ok=True)
        with dest.open("wb") as fh:
            for chunk in r.iter_content(chunk_size=1 << 20):
                fh.write(chunk)


"""
Values
"""
# local folder for saving data
OUTDIR = "downloads"

# MG-RAST-specific values:
# paste your 25-char web-key here (or keep None)
AUTH_KEY = None
# "raw" = original reads; "assembled" = contigs
STAGE = "raw"
FILE_IDS = {
    "raw"      : "050.1",
    "assembled": "299.1",
}
# Base URL for the MG-RAST API           
BASE_URL = "https://api.mg-rast.org/1"

# Project-specific values:
# Project ID for the David et al. (2014) dataset
PROJECT_ID = "mgp6248"
# Acceptable dataset names:  DD1 … DD999  or  ID1 … ID999. Case-insensitive.
NAME_OK = re.compile(r'^(?:DD|ID)\d{1,3}$', re.I)


"""
Download Target .fna files from MG-RAST for Accession ID
"""
if STAGE not in FILE_IDS:
    sys.exit(f"STAGE must be 'raw' or 'assembled' (got: {STAGE!r})")

try:
    metagenomes = list_metagenomes(PROJECT_ID, AUTH_KEY)
    # keep only those dataset names that match the NAME_OK pattern (ie: DD### or ID###)
    metagenomes = [m for m in metagenomes if NAME_OK.match(m["name"])]
except requests.HTTPError as e:
    sys.exit(f"Could not list project metagenomes → {e}")

if not metagenomes:
    sys.exit(f"Project {PROJECT_ID} contains no metagenomes!")

# List of metagenomes to download
print(f"Found {len(metagenomes)} target metagenomes in project {PROJECT_ID}.")
wanted_fileid = FILE_IDS[STAGE]
out_root = Path(OUTDIR).expanduser()

for mgm in metagenomes:
    start_time = time.time()  # Start time for this metagenome download
    mgm_name = mgm["name"]
    mgm_id = mgm["id"]  # e.g. 'mgm1234567.3'

    # 1) discover files available for this metagenome
    files = list_files(mgm_id, AUTH_KEY)

    # 2) pick our target .fna
    entry = pick_fna(files, wanted_fileid)
    if entry is None:
        print(f"[skip] {mgm_name}: no .fna with id {wanted_fileid}")
        continue

    # 3) decide local filename
    mgm_name_ext = mgm_name + ".fna"
    dest = out_root / PROJECT_ID / mgm_name_ext
    if dest.exists():
        print(f"[skip] {dest} (already downloaded)")
        continue

    # 4) download
    print(f"[dl  ] {mgm_name} → {dest}")
    try:
        stream_download(entry["url"], dest, AUTH_KEY)
        print(f"[done] {mgm_name} ({time.time() - start_time:.2f}s)")
    except requests.HTTPError as e:
        print(f"[fail] {mgm_name}: {e}")

print("Finished downloading .fna files.")