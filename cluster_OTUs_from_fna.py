#!/usr/bin/env python3
"""
cluster_OTUs_from_fna.py
Created 7/8/25

This script clusters OTUs (Operational Taxonomic Units) from 16S rRNA sequences using the VSEARCH tool. This script was built for metagenomic data analysis from David et al. (2014) (https://pmc.ncbi.nlm.nih.gov/articles/PMC3957428/). This script was built using (1) the QIIME 2 "Moving Pictures" tutorial as a template (https://amplicon-docs.qiime2.org/en/latest/tutorials/moving-pictures.html) and (2) the older OTU clustering tutorial for QIIME2 (https://docs.qiime2.org/2024.10/tutorials/otu-clustering/). Use the QIIME 2 viewer to visualize .qza/.qzv file outputs (https://view.qiime2.org/).

Inputs:
- .fna files downloaded from MG-RAST for 16S rRNA sequencing. Note, this script utilizes processed .fna files only, not .fastq files
- metadata file containing sample information

Outputs:
- .qzv file visualizing the metadata
- .qza files for sequences (from .fna) per sample
- .qza files for feature tables and data for sequence frequencies per sample
- .qza files containing (de novo) clustered OTUs per sample

***For Windows users, you can run the script through the Ubuntu terminal, due to QIIME 2 (first line run once per terminal session):
conda activate qiime2-amplicon-2025.4
python3 cluster_OTUs_from_fna.py
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
def import_fna_as_artifact(fna_file, sample_id, qza_output_dir, semantic_type="SampleData[Sequences]"):
    """Import a .fna file as a QIIME 2 Artifact.

    Source: https://docs.qiime2.org/2024.10/tutorials/otu-clustering/
    
    Args:
        fna_file (str): Path to the .fna file.
        sample_id (str): Sample ID for the artifact.
        qza_output_dir (str): Directory for the QIIME 2 artifact output.

    Returns:
        None
    """
    artifact_path = pjoin(qza_output_dir, f"{sample_id}.qza")
    subprocess.run([
        "qiime", "tools", "import",
        "--type", semantic_type,
        "--input-path", fna_file,
        "--output-path", artifact_path,
    ], check=True)
    return

def generate_feature_frequency_table(qza_input_dir, sample_id, output_dir):
    """Run QIIME 2's dereplicate-sequences command to generate a FeatureTable[Frequency] artifact and corresponding FeatureData[Sequence] artifact.

    Source: https://docs.qiime2.org/2024.10/tutorials/otu-clustering/

    Args:
        qza_input_dir (str): Directory containing the input QIIME 2 artifacts.
        sample_id (str): Sample ID for the artifacts.
        output_dir (str): Directory for the output artifacts.

    Returns:
        None
    """
    feature_table_path = pjoin(output_dir, f"{sample_id}_feature_table.qza")
    feature_data_path = pjoin(output_dir, f"{sample_id}_feature_data.qza")

    subprocess.run([
        "qiime", "vsearch", "dereplicate-sequences",
        "--i-sequences", pjoin(qza_input_dir, f"{sample_id}.qza"),
        "--o-dereplicated-table", feature_table_path,
        "--o-dereplicated-sequences", feature_data_path
    ], check=True)
    return

def create_sequences_histogram(qza, output_image):
    """Create a histogram of sequence lengths in a QIIME 2 artifact.

    Args:
        qza (str): Path to a sequences QIIME 2 artifact .qza file.
        output_image (str): Path to save the histogram image.

    Returns:
        None
    """
    subprocess.run([
        "qiime", "feature-table", "tabulate-seqs",
        "--i-data", qza,
        "--o-visualization", output_image
    ], check=True)
    return

def cluster_otus_de_novo(feature_input_dir, sample_id, otu_output_dir,
                 p_ident="0.97", strand="both"):
    """Use QIIME 2's vsearch plugin to perform de novo clustering of OTUs.

    Source: https://docs.qiime2.org/2024.10/tutorials/otu-clustering/

    Args:
        sample_id (str): Sample ID for the artifacts.
        output_dir (str): Directory for the clustered artifacts.
        p_ident (str): Percent identity for clustering (default is "0.97" for 97% identity).
        strand (str): Orientation handling for matching (‘plus’, ‘minus’, or ‘both’).

    Returns:
        None    
    """
    feature_table_path = pjoin(feature_input_dir, f"{sample_id}_feature_table.qza")
    feature_data_path = pjoin(feature_input_dir, f"{sample_id}_feature_data.qza")
    subprocess.run([
        "qiime", "vsearch", "cluster-features-de-novo",
        "--i-table", feature_table_path,
        "--i-sequences", feature_data_path,
        "--p-perc-identity", p_ident,
        "--o-clustered-table", pjoin(otu_output_dir, f"{sample_id}-table-cr-{p_ident}.qza"),
        "--o-clustered-sequences", pjoin(otu_output_dir, f"{sample_id}-rep-seqs-cr-{p_ident}.qza")
    ], check=True)
    return


"""
Values
"""
# input folder holds the metadata file
INPUT_FOLDER = "input"
# "metadata_16S.xlsx" for full metadata, "metadata_test.xlsx" for simpler metadata testing
METADATA_FILENAME_XLSX = "metadata_test.xlsx"
SAMPLE_ID_COL_NAME = "SampleID"

# downloads folder holds the MG-RAST files downloaded by download_fna_files.py
DOWNLOADS_FOLDER = "downloads"
# MG-RAST Accession ID (mgp6248) downloads sub-folder for David et al. (2014) 16S rRNA sequencing data
MGRAST_ACCESSION_FOLDER = "mgp6248"

# Create an overall output folder if it does not exist yet
OUTPUT_FOLDER = "output"
if not os.path.exists(OUTPUT_FOLDER):
    os.makedirs(OUTPUT_FOLDER)
# Create an output sub-folder for this script
OUTPUT_SCRIPT_FOLDER = "cluster_OTUs_from_fna"
OUTPUT_SCRIPT_FOLDER_DIR = pjoin(OUTPUT_FOLDER, OUTPUT_SCRIPT_FOLDER)
if not os.path.exists(OUTPUT_SCRIPT_FOLDER_DIR):
    os.makedirs(OUTPUT_SCRIPT_FOLDER_DIR)

# GreenGenes 97% OTU reference sequences URL for download
URL_GREENGENES_97 = "https://data.qiime2.org/classifiers/greengenes/gg_13_8_otus.tar.gz"
TARGET_GREENGENES_97 = URL_GREENGENES_97.split("/")[-1]  # Extract the filename from the URL
FASTA_FILENAME_GREENGENES_97 = "97_otus.fasta"
# location of the desired 97% OTU reference sequences fasta file in the target downloaded folder
DOWNLOADED_FASTA_DIR_GREENGENES_97 = pjoin("gg_13_8_otus", "rep_set")
V4_PRIMER_F = "GTGCCAGCMGCCGCGGTAA"
V4_PRIMER_R = "GGACTACHVGGGTWTCTAAT"

# Create a name for the GreenGenes 97% OTU reference sequences file
REFERENCE_GREENGENES_97_QZA_FILENAME = "greengenes_97_otus.qza"

# Set to True when you want to run all steps, even if the output files already exist
RERUN_EXISTING_DATA = False


"""
Import Metadata
"""
# Import as pandas. Rewrite as .tsv for QIIME 2 compatibility
sample_metadata_df = pd.read_excel(pjoin(INPUT_FOLDER, METADATA_FILENAME_XLSX))
metadata_filename_tsv = METADATA_FILENAME_XLSX.replace(".xlsx", ".tsv")
# Save as .tsv file for QIIME 2 compatibility
sample_metadata_df.to_csv(pjoin(OUTPUT_SCRIPT_FOLDER_DIR, metadata_filename_tsv), sep="\t", index=False)

# Use QIIME 2's Metadata class to load the metadata from the .tsv file
sample_metadata = Metadata.load(pjoin(OUTPUT_SCRIPT_FOLDER_DIR, metadata_filename_tsv))
# Use QIIME 2's metadata plugin "tabulate" to visualize the metadata
sample_metadata_viz, = metadata_actions.tabulate(sample_metadata)
# Save the metadata visualization to a .qzv file
sample_metadata_viz.save(pjoin(OUTPUT_SCRIPT_FOLDER_DIR, METADATA_FILENAME_XLSX.replace(".xlsx", ".qzv")))

print("Imported metadata.")


"""
Import Data into QIIME 2 Artifacts
"""
# QIIME 2 Data Importing Documentation: https://docs.qiime2.org/2024.10/tutorials/importing/
# Notes: 
# - QIIME 2 currently supports importing the QIIME 1 seqs.fna file format, which consists of a single FASTA file with exactly two lines per record: header and sequence.
# - Each sequence must span exactly one line and cannot be split across multiple lines. The ID in each header must follow the format <sample-id>_<seq-id>.
# - <sample-id> is the identifier of the sample the sequence belongs to, and <seq-id> is an identifier for the sequence within its sample.

# Use subprocess to import .fna files as a QIIME 2 Artifacts

# if the sequence qza folder does not exist, create it
qza_output_path = pjoin(OUTPUT_SCRIPT_FOLDER_DIR, "sequence qza files")
if not os.path.exists(qza_output_path):
    os.makedirs(qza_output_path)

# Iterate through each sample ID in the metadata
for index, row in sample_metadata_df.iterrows():
    sample_id = row[SAMPLE_ID_COL_NAME]
    
    # Construct the path to the .fna file
    fna_file_path = pjoin(DOWNLOADS_FOLDER, MGRAST_ACCESSION_FOLDER, f"{sample_id}.fna")
    
    if not RERUN_EXISTING_DATA:
        # If artifact file already exists, skip importing
        artifact_path = pjoin(OUTPUT_SCRIPT_FOLDER_DIR, "sequence qza files", f"{sample_id}.qza")
        if os.path.exists(artifact_path):
            print(f"Skipping {sample_id}, artifact already exists: {artifact_path}")
            continue

    # Check if the .fna file exists
    if os.path.exists(fna_file_path):
        print(f"Processing sample: {sample_id}")
        
        # Import the .fna file as a QIIME 2 artifact
        try:
            import_fna_as_artifact(fna_file_path, sample_id, qza_output_path)
            print(f"Successfully imported {sample_id} as artifact.")
        except subprocess.CalledProcessError as e:
            print(f"Error importing {sample_id}: {e}")
    else:
        print(f"Warning: .fna file not found for sample {sample_id}: {fna_file_path}")

print("Completed importing all .fna files from the metadata as QIIME 2 artifacts.")


"""
"Dereplicate Sequences" to Generate a Feature Frequency Table for OTU Clustering per Sample
"""
# From QIIME2 OTU Clustering Tutorial: https://docs.qiime2.org/2024.10/tutorials/otu-clustering/
# Notes:
# - After importing the .fna data, you can dereplicate it with the dereplicate-sequences command.
# - This step is still "per sample" (frequency refers to the number of times a sequence appears in one sample, not across samples).
# - Outputs:
# - FeatureTable[Frequency]: contains the frequency of each unique amplicon in the sample.
# - FeatureData[Sequence]: contains mapping of each feature identifier to a corresponding sequence.

# If the feature table output folder does not exist, create it
feature_table_output_path = pjoin(OUTPUT_SCRIPT_FOLDER_DIR, "feature tables")
if not os.path.exists(feature_table_output_path):
    os.makedirs(feature_table_output_path)

# Iterate through each sample ID in the metadata
for index, row in sample_metadata_df.iterrows():
    sample_id = row[SAMPLE_ID_COL_NAME]
    
    if not RERUN_EXISTING_DATA:
        # If artifact file already exists, skip generating feature frequency table
        sample_feature_table_path = pjoin(feature_table_output_path, f"{sample_id}_feature_table.qza")
        if os.path.exists(sample_feature_table_path):
            print(f"Skipping {sample_id}, feature table already exists: {sample_feature_table_path}")
            continue
    
    # Generate the feature frequency table and feature data for the sample
    try:
        generate_feature_frequency_table(qza_output_path, sample_id, feature_table_output_path)
        print(f"Successfully generated feature frequency table and feature data for {sample_id}.")
    except subprocess.CalledProcessError as e:
        print(f"Error generating feature frequency table for {sample_id}: {e}")


"""
OTU Clustering de novo
"""
# Info on picking OTU clustering strategy (de novo vs. closed reference vs. open reference: https://qiime.org/tutorials/otu_picking.html)
# If it does not exist, create the output folder for OTU clustering results
otu_clustering_output_path = pjoin(OUTPUT_SCRIPT_FOLDER_DIR, "OTU_clustering")
if not os.path.exists(otu_clustering_output_path):
    os.makedirs(otu_clustering_output_path)

# Loop over each sample ID in the metadata
# Perform de novo clustering using the vsearch plugin in QIIME 2.
for index, row in sample_metadata_df.iterrows():
    sample_id = row[SAMPLE_ID_COL_NAME]

    if not RERUN_EXISTING_DATA:
        # If artifact file already exists, skip clustering
        otu_table_path = pjoin(otu_clustering_output_path, f"{sample_id}-table-cr-0.97.qza")
        if os.path.exists(otu_table_path):
            print(f"Skipping {sample_id}, OTU clustering already exists: {otu_table_path}")
            continue

    # Cluster the OTUs for the sample
    try:
        cluster_otus_de_novo(feature_table_output_path, sample_id, otu_clustering_output_path)
        print(f"Successfully clustered OTUs for {sample_id}.")
    except subprocess.CalledProcessError as e:
        print(f"Error clustering OTUs for {sample_id}: {e}")