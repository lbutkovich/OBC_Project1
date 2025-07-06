"""
Tutorial: Illumina Overview for Using QIIME and QIIME 2
"""

import urllib.request
import tarfile
import os

# Download the tutorial file
url = "ftp://ftp.microbio.me/qiime/tutorial_files/moving_pictures_tutorial-1.9.0.tgz"
filename = "moving_pictures_tutorial-1.9.0.tgz"

try:
    print(f"Downloading {filename}...")
    urllib.request.urlretrieve(url, filename)
    print(f"Successfully downloaded {filename}")
except Exception as e:
    print(f"Error downloading file: {e}")

# Extract the tar.gz file
if os.path.exists(filename):
    try:
        print(f"Extracting {filename}...")
        with tarfile.open(filename, 'r:gz') as tar:
            tar.extractall()
        print(f"Successfully extracted {filename}")
    except Exception as e:
        print(f"Error extracting file: {e}")
else:
    print(f"File {filename} not found for extraction")



