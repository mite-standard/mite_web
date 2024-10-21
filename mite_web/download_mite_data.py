"""Downloads and unpacks latest version of mite_data.

This script runs as part of the startup
"""

import json
import os
import shutil
import sys
from pathlib import Path

import requests

location = Path(__file__).parent.joinpath("mite_web/data")
record_url = (
    "https://zenodo.org/api/records/13294303"  # Always resolves to the latest version
)
record_path = location.joinpath("record.zip")
version_path = location.joinpath("version.json")
record_path_unzip = location.joinpath("record")

# Create location
if location.exists():
    print(
        f"Location '{location}' exists already - data is presumably already downloaded. SKIP"
    )
    sys.exit(0)
else:
    print(f"Location '{location}' does not exist - started download.")
    location.mkdir(parents=True)

# Download metadata
response = requests.get(record_url)
if response.status_code != 200:
    raise Exception(
        f"Error fetching 'mite_data' record metadata: {response.status_code}"
    )

record_metadata = response.json()

version = record_metadata["metadata"]["version"]
files_url = record_metadata["files"][0]["links"]["self"]

# Download data
response = requests.get(files_url)
if response.status_code == 200:
    with open(record_path, "wb") as f:
        f.write(response.content)

    print(f"Downloaded mite_data record version {version} successfully!")

    with open(version_path, "w") as f:
        f.write(json.dumps({"version": f"{version}"}))
else:
    raise Exception(f"Error downloading 'mite_data' record: {response.status_code}")

# Unpack data
shutil.unpack_archive(filename=record_path, extract_dir=record_path_unzip, format="zip")
if not record_path_unzip.exists():
    raise NotADirectoryError(f"Could not find the directory {record_path_unzip}.")

# Retrieve the data folders and move them for easier access
matching_dirs = list(record_path_unzip.glob("mite-standard-mite_data-*"))
if matching_dirs:
    subdir = matching_dirs[0]
    shutil.move(
        src=record_path_unzip.joinpath(subdir).joinpath("mite_data/data").resolve(),
        dst=location.resolve(),
    )
    shutil.move(
        src=record_path_unzip.joinpath(subdir).joinpath("mite_data/img").resolve(),
        dst=location.resolve(),
    )
    os.remove(record_path)
    shutil.rmtree(record_path_unzip)
    print("Successfully prepared files")
    sys.exit(0)
else:
    raise RuntimeError(
        f"Could not determine data storage location in downloaded directory."
    )
