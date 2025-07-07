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
    """HTTP GET → JSON (raises on HTTP errors)."""
    hdr = {"auth": auth} if auth else {}
    r = requests.get(url, headers=hdr, timeout=90)
    r.raise_for_status()
    return r.json()

def list_metagenomes(project_id: str, auth: Optional[str] = None) -> List[dict]:
    """
    Return [{'id': 'mgm1234.3', 'name': 'DD17'}, …] for every metagenome
    belonging to *project_id*.
    """
    url  = f"{BASE_URL}/project/{project_id}?verbosity=full"
    info = get_json(url, auth)
    return [{"id": m["metagenome_id"], "name": m["name"]}          # 'name' key
            for m in info["metagenomes"]]

def list_files(mgm_id: str, auth: Optional[str] = None) -> List[Dict]:
    """Catalogue of downloadable artefacts for one metagenome."""
    return get_json(f"{BASE_URL}/download/{mgm_id}", auth)["data"]

def pick_fna(files: List[Dict], wanted_id: str) -> Optional[Dict]:
    """Return the entry that matches *wanted_id* and ends with '.fna'."""
    return next(
        (f for f in files if f["file_id"] == wanted_id and f["file_name"].endswith(".fna")),
        None
    )

def stream_download(url: str, dest: Path, auth: Optional[str] = None) -> None:
    """Stream large files in 1 MiB chunks → disk (keeps RAM low)."""
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
OUTDIR = "downloads"     # local folder for saving data
AUTH_KEY = None            # paste your 25-char web-key here (or keep None)
STAGE = "raw"           # "raw"  → original reads (050.1)
                        # "assembled" → contigs (299.1)
BASE_URL = "https://api.mg-rast.org/1"

FILE_IDS = {
    "raw"      : "050.1",   # *.upload.fna
    "assembled": "299.1",   # *.assembled_contigs.fna
}

PROJECT_ID = "mgp6248"
# Acceptable dataset names:  DD1 … DD999  or  ID1 … ID999
NAME_OK = re.compile(r'^(?:DD|ID)\d{1,3}$', re.I)   # case-insensitive

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
