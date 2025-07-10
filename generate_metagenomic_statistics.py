#!/usr/bin/env python3
"""
generate_metagenomic_statistics.py
Created 7/9/25


***For Windows users, you can run the script through the Ubuntu terminal, due to QIIME 2 (first line run once per terminal session):
conda activate qiime2-amplicon-2025.4
python3 generate_metagenomic_statistics.py
"""


import os
from os.path import join as pjoin
import pandas as pd
import subprocess
import requests
from pathlib import Path
import tarfile
import shutil

# QIIME 2 imports:
from qiime2 import Metadata
import qiime2.plugins.metadata.actions as metadata_actions


"""
Functions
"""
def classify_taxonomy(merged_seqs_qza, classifier_qza, taxonomy_out_path, threads=2, memory_limit="4G"):
    """Classify the taxonomy of sequences using a pre-trained classifier.
    
    Args:
        merged_seqs_qza (str): Path to the merged feature sequences artifact.
        classifier_qza (str): Path to the pre-trained classifier artifact.
        taxonomy_out_path (str): Output path for the classified taxonomy artifact.
        threads (int): Number of threads to use for classification.
        memory_limit (str): Memory limit per job (e.g., "4G" for 4GB).
    Returns:
        taxonomy_tabulated_path (str): Path to the tabulated taxonomy visualization artifact.
    """
    try:
        print(f"Starting taxonomy classification with {threads} threads and {memory_limit} memory limit per job...")
        # Run QIIME 2 feature-classifier with reduced resources to avoid memory issues
        subprocess.run([
            "qiime", "feature-classifier", "classify-sklearn",
            "--i-classifier", classifier_qza,
            "--i-reads", merged_seqs_qza,
            "--o-classification", taxonomy_out_path,
            "--p-n-jobs", str(threads),
            "--p-reads-per-batch", "5000",  # Process in smaller batches
            "--verbose"
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error in taxonomy classification: {e}")
        print("Try reducing the number of threads or increasing available memory.")
        raise

    # Run QIIME 2 metadata tabulate
    taxonomy_tabulated_path = taxonomy_out_path.replace(".qza", "_tabulated.qzv")
    try:
        subprocess.run([
            "qiime", "metadata", "tabulate",
            "--m-input-file", taxonomy_out_path,
            "--o-visualization", taxonomy_tabulated_path
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error in metadata tabulation: {e}")
        raise
    return taxonomy_tabulated_path



"""
Values
"""
# input folder holds the metadata file
INPUT_FOLDER = "input"
# "metadata_16S.xlsx" for full metadata, "metadata_test.xlsx" for simpler metadata testing
METADATA_FILENAME_XLSX = "metadata_test.xlsx"
# Run descriptor (test or full)
RUN_DESCRIPTOR = "test"
SAMPLE_ID_COL_NAME = "SampleID"

# Create an output sub-folder for this script
OUTPUT_FOLDER = "output"
OUTPUT_SCRIPT_FOLDER = "generate_metagenomic_statistics"
OUTPUT_SCRIPT_FOLDER_DIR = pjoin(OUTPUT_FOLDER, OUTPUT_SCRIPT_FOLDER)
if not os.path.exists(OUTPUT_SCRIPT_FOLDER_DIR):
    os.makedirs(OUTPUT_SCRIPT_FOLDER_DIR)

# Output paths from script cluster_OTUs_from_fna.py
CLUSTER_OTUS_OUTPUT_FOLDER = "cluster_OTUs_from_fna"
CLUSTER_OTUS_OUTPUT_FOLDER_DIR = pjoin(OUTPUT_FOLDER, CLUSTER_OTUS_OUTPUT_FOLDER)
# Use this merged seqs filename for pre-OTU clustering option:
# MERGED_SEQS_FEATURE_QZA_FILENAME = "merged_feature_sequences_" + RUN_DESCRIPTOR + ".qza"
# To use OTU-clustered data, use this filename:
MERGED_SEQS_FEATURE_QZA_FILENAME = "merged-rep-seqs-cr-0.97_" + RUN_DESCRIPTOR + ".qza"
MERGED_SEQS_FEATURE_QZA_PATH = pjoin(CLUSTER_OTUS_OUTPUT_FOLDER_DIR, MERGED_SEQS_FEATURE_QZA_FILENAME)

# Set to True when you want to run all steps, even if the output files already exist
RERUN_EXISTING_DATA = False

# Manually download pre-trained taxonomy classifier silva-138-99-nb-classifier.qza for 16S rRNA gene sequences, and store in input folder
# https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza
CLASSIFIER_FILENAME = "silva-138-99-nb-classifier.qza"
CLASSIFIER_PATH = pjoin(INPUT_FOLDER, CLASSIFIER_FILENAME)

"""
Assign Taxonomy
"""
# Perform naive Bayes classification of sequences using a pre-trained classifier
# Use classifier
taxonomy_out_path = pjoin(OUTPUT_SCRIPT_FOLDER_DIR, "taxonomy_" + RUN_DESCRIPTOR + ".qza")
if not os.path.exists(taxonomy_out_path) or RERUN_EXISTING_DATA:
    if not os.path.exists(CLASSIFIER_PATH):
        raise FileNotFoundError(f"Classifier file not found: {CLASSIFIER_PATH}. Please download it and place it in the input folder.")
    
    # Classify taxonomy
    # Using only 2 threads to reduce memory usage
    taxonomy_tabulated_path = classify_taxonomy(MERGED_SEQS_FEATURE_QZA_PATH, CLASSIFIER_PATH, taxonomy_out_path, threads=2)
    print(f"Taxonomy classification completed. Output saved to: {taxonomy_tabulated_path}")
else:
    print(f"Taxonomy classification already exists at: {taxonomy_out_path}. Skipping classification step.")
    taxonomy_tabulated_path = taxonomy_out_path.replace(".qza", "_tabulated.qzv")


"""
Infer Phylogenetic Tree
"""
# QIIME 2 q2-phylogeny documentation: https://docs.qiime2.org/2024.10/tutorials/phylogeny/


"""
"""