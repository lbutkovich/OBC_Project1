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
import glob

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

def merge_feature_table(feature_tables_dir, output_dir, run_descriptor, batch_size=20):
    """Merge all per-sample feature tables into a single feature table.
    Process files in batches to avoid memory issues.

    Args:
        feature_tables_dir (str): Directory containing per-sample feature tables.
        output_dir (str): Directory to save the merged feature table.
        run_descriptor (str): Descriptor for the run (e.g., "test" or "full").
        batch_size (int): Number of files to process in each batch.

    Returns:
        merged_feature_table_path (str): Path to the merged feature table artifact.
    """
    merged_table_path = pjoin(output_dir, "merged_feature_table_" + run_descriptor + ".qza")
    temp_dir = pjoin(output_dir, "temp_merge")
    
    # Create temp directory if it doesn't exist
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    
    # Get all feature table files in the directory
    feature_table_files = glob.glob(pjoin(feature_tables_dir, "*_feature_table.qza"))
    
    if not feature_table_files:
        raise ValueError(f"No feature table files found in {feature_tables_dir}")
    
    # Process in batches
    batches = [feature_table_files[i:i + batch_size] for i in range(0, len(feature_table_files), batch_size)]
    print(f"Processing {len(feature_table_files)} files in {len(batches)} batches of {batch_size}")
    
    # First round of merges (create intermediate merged files)
    intermediate_files = []
    for i, batch in enumerate(batches):
        intermediate_file = pjoin(temp_dir, f"intermediate_merge_{i}.qza")
        intermediate_files.append(intermediate_file)
        
        # Build the command with feature table files for this batch
        command = ["qiime", "feature-table", "merge"]
        for file in batch:
            command.extend(["--i-tables", file])
        command.extend(["--o-merged-table", intermediate_file])
        
        print(f"Processing batch {i+1} of {len(batches)}...")
        # Run the command for this batch
        subprocess.run(command, check=True)
    
    # If there's only one intermediate file, rename it as the final result
    if len(intermediate_files) == 1:
        shutil.move(intermediate_files[0], merged_table_path)
    else:
        # Second round: merge the intermediate files
        command = ["qiime", "feature-table", "merge"]
        for file in intermediate_files:
            command.extend(["--i-tables", file])
        command.extend(["--o-merged-table", merged_table_path])
        
        print(f"Merging {len(intermediate_files)} intermediate files...")
        subprocess.run(command, check=True)
        
        # Clean up intermediate files
        for file in intermediate_files:
            if os.path.exists(file):
                os.remove(file)
    
    # Remove temp directory if it's empty
    if os.path.exists(temp_dir) and not os.listdir(temp_dir):
        os.rmdir(temp_dir)
    
    return merged_table_path

def merge_feature_sequences(feature_sequences_dir, output_dir, run_descriptor, batch_size=20):
    """Merge all per-sample feature sequences into a single feature sequences artifact.
    Process files in batches to avoid memory issues.

    Args:
        feature_sequences_dir (str): Directory containing per-sample feature sequences.
        output_dir (str): Directory to save the merged feature sequences.
        run_descriptor (str): Descriptor for the run (e.g., "test" or "full").
        batch_size (int): Number of files to process in each batch.

    Returns:
        merged_sequences_path (str): Path to the merged feature sequences artifact.
    """
    merged_sequences_path = pjoin(output_dir, "merged_feature_sequences_" + run_descriptor + ".qza")
    temp_dir = pjoin(output_dir, "temp_merge_seqs")
    
    # Create temp directory if it doesn't exist
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    # Get all feature sequence files in the directory
    feature_sequence_files = glob.glob(pjoin(feature_sequences_dir, "*_feature_data.qza"))
    
    if not feature_sequence_files:
        raise ValueError(f"No feature sequence files found in {feature_sequences_dir}")
    
    # Process in batches
    batches = [feature_sequence_files[i:i + batch_size] for i in range(0, len(feature_sequence_files), batch_size)]
    print(f"Processing {len(feature_sequence_files)} sequence files in {len(batches)} batches of {batch_size}")
    
    # First round of merges (create intermediate merged files)
    intermediate_files = []
    for i, batch in enumerate(batches):
        intermediate_file = pjoin(temp_dir, f"intermediate_merge_seq_{i}.qza")
        intermediate_files.append(intermediate_file)
        
        # Build the command with feature sequence files for this batch
        command = ["qiime", "feature-table", "merge-seqs"]
        for file in batch:
            command.extend(["--i-data", file])
        command.extend(["--o-merged-data", intermediate_file])
        
        print(f"Processing sequence batch {i+1} of {len(batches)}...")
        # Run the command for this batch
        subprocess.run(command, check=True)
    
    # If there's only one intermediate file, rename it as the final result
    if len(intermediate_files) == 1:
        shutil.move(intermediate_files[0], merged_sequences_path)
    else:
        # Second round: merge the intermediate files
        command = ["qiime", "feature-table", "merge-seqs"]
        for file in intermediate_files:
            command.extend(["--i-data", file])
        command.extend(["--o-merged-data", merged_sequences_path])
        
        print(f"Merging {len(intermediate_files)} intermediate sequence files...")
        subprocess.run(command, check=True)
        
        # Clean up intermediate files
        for file in intermediate_files:
            if os.path.exists(file):
                os.remove(file)
    
    # Remove temp directory if it's empty
    if os.path.exists(temp_dir) and not os.listdir(temp_dir):
        os.rmdir(temp_dir)
    
    return merged_sequences_path

def cluster_otus_de_novo_per_sample(feature_input_dir, sample_id, otu_output_dir,
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

def cluster_otus_de_novo_merged(merged_feature_table_path, merged_feature_sequences_path, otu_output_dir, run_descriptor, p_ident="0.97", strand="both"):
    """Use QIIME 2's vsearch plugin to perform de novo clustering of OTUs on a merged feature table.

    Args:
        merged_feature_table_path (str): Path to the merged feature table.
        merged_feature_sequences_path (str): Path to the merged feature sequences.
        otu_output_dir (str): Directory for the clustered artifacts.
        run_descriptor (str): Descriptor for the run (e.g., "test" or "full").
        p_ident (str): Percent identity for clustering (default is "0.97" for 97% identity).
        strand (str): Orientation handling for matching ('plus', 'minus', or 'both').

    Returns:
        None
    """
    subprocess.run([
        "qiime", "vsearch", "cluster-features-de-novo",
        "--i-table", merged_feature_table_path,
        "--i-sequences", merged_feature_sequences_path,
        "--p-perc-identity", p_ident,
        "--o-clustered-table", pjoin(otu_output_dir, f"merged-table-cr-{p_ident}_" + run_descriptor + ".qza"),
        "--o-clustered-sequences", pjoin(otu_output_dir, f"merged-rep-seqs-cr-{p_ident}_" + run_descriptor + ".qza")
    ], check=True)
    return


"""
Values
"""
# input folder holds the metadata file
INPUT_FOLDER = "input"
# "metadata_16S.xlsx" for full metadata, "metadata_test.xlsx" for simpler metadata testing
METADATA_FILENAME_XLSX = "metadata_16S.xlsx"
# Run descriptor (test or full)
RUN_DESCRIPTOR = "full"
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
qza_output_path = pjoin(OUTPUT_SCRIPT_FOLDER_DIR, "sequence_qza_files")
if not os.path.exists(qza_output_path):
    os.makedirs(qza_output_path)

# Iterate through each sample ID in the metadata
for index, row in sample_metadata_df.iterrows():
    sample_id = row[SAMPLE_ID_COL_NAME]
    
    # Construct the path to the .fna file
    fna_file_path = pjoin(DOWNLOADS_FOLDER, MGRAST_ACCESSION_FOLDER, f"{sample_id}.fna")
    
    if not RERUN_EXISTING_DATA:
        # If artifact file already exists, skip importing
        artifact_path = pjoin(OUTPUT_SCRIPT_FOLDER_DIR, "sequence_qza_files", f"{sample_id}.qza")
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
feature_table_output_path = pjoin(OUTPUT_SCRIPT_FOLDER_DIR, "sample_feature_tables")
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
Merge Per-sample Files to Generate (i) a Merged Feature Table and (ii) Merged .qza Sequence File
"""
# Merge all per-sample feature tables into a single feature table.
# This is necessary for de novo OTU clustering on the merged feature table.
# Store in OUTPUT_SCRIPT_FOLDER_DIR (not sub-folder for sample_feature_tables)
merged_feature_table_path = pjoin(OUTPUT_SCRIPT_FOLDER_DIR, "merged_feature_table_" + RUN_DESCRIPTOR + ".qza")
if not os.path.exists(merged_feature_table_path) or RERUN_EXISTING_DATA:
    merged_feature_table_path = merge_feature_table(feature_table_output_path, OUTPUT_SCRIPT_FOLDER_DIR, RUN_DESCRIPTOR)
else:
    print(f"Skipping feature table merge, file already exists: {merged_feature_table_path}")

# Merge all per-sample feature sequences into a single feature sequence file
merged_feature_sequences_path = pjoin(OUTPUT_SCRIPT_FOLDER_DIR, "merged_feature_sequences_" + RUN_DESCRIPTOR + ".qza")
if not os.path.exists(merged_feature_sequences_path) or RERUN_EXISTING_DATA:
    merged_feature_sequences_path = merge_feature_sequences(feature_table_output_path, OUTPUT_SCRIPT_FOLDER_DIR, RUN_DESCRIPTOR)
else:
    print(f"Skipping feature sequences merge, file already exists: {merged_feature_sequences_path}")


# """
# OTU Clustering de novo, per sample basis <-- old version
# """
# # Info on picking OTU clustering strategy (de novo vs. closed reference vs. open reference: https://qiime.org/tutorials/otu_picking.html)
# # If it does not exist, create the output folder for OTU clustering results
# otu_clustering_output_path = pjoin(OUTPUT_SCRIPT_FOLDER_DIR, "sample_OTU_clustering")
# if not os.path.exists(otu_clustering_output_path):
#     os.makedirs(otu_clustering_output_path)

# # Loop over each sample ID in the metadata
# # Perform de novo clustering using the vsearch plugin in QIIME 2.
# for index, row in sample_metadata_df.iterrows():
#     sample_id = row[SAMPLE_ID_COL_NAME]

#     if not RERUN_EXISTING_DATA:
#         # If artifact file already exists, skip clustering
#         otu_table_path = pjoin(otu_clustering_output_path, f"{sample_id}-table-cr-0.97.qza")
#         if os.path.exists(otu_table_path):
#             print(f"Skipping {sample_id}, OTU clustering already exists: {otu_table_path}")
#             continue

#     # Cluster the OTUs for the sample
#     try:
#         cluster_otus_de_novo_per_sample(feature_table_output_path, sample_id, otu_clustering_output_path)
#         print(f"Successfully clustered OTUs for {sample_id}.")
#     except subprocess.CalledProcessError as e:
#         print(f"Error clustering OTUs for {sample_id}: {e}")


"""
OTU Clustering de novo on merged feature table
"""
# Use QIIME 2's vsearch plugin to perform de novo clustering of OTUs on a merged feature table.
# Store in OUTPUT_SCRIPT_FOLDER_DIR
merged_table_cr_path = pjoin(OUTPUT_SCRIPT_FOLDER_DIR, f"merged-table-cr-0.97_{RUN_DESCRIPTOR}.qza")
merged_rep_seqs_cr_path = pjoin(OUTPUT_SCRIPT_FOLDER_DIR, f"merged-rep-seqs-cr-0.97_{RUN_DESCRIPTOR}.qza")

if not os.path.exists(merged_table_cr_path) or not os.path.exists(merged_rep_seqs_cr_path) or RERUN_EXISTING_DATA:
    cluster_otus_de_novo_merged(merged_feature_table_path, merged_feature_sequences_path, OUTPUT_SCRIPT_FOLDER_DIR, RUN_DESCRIPTOR)
    print(f"Successfully clustered OTUs for merged feature table.")
else:
    print(f"Skipping OTU clustering for merged feature table, output files already exist: {merged_table_cr_path} and {merged_rep_seqs_cr_path}")
