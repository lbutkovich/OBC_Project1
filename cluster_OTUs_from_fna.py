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
- .qza file containing clustered OTUs

For Windows users, you can run the script through the Ubuntu terminal, due to QIIME 2 (first line run once per terminal session):
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

def download_greengenes_97_reference(url, target_folder, target_fasta, target_fasta_dir, output_dir):
    """Download the GreenGenes 97% OTU reference .fasta.

    Download the target_folder (zipped) from url. Extract the contents and move the target_fasta (found at target_fasta_dir) to output_dir. Delete the other extracted files.
    
    Args:
        url (str): The URL to download the reference sequences from.
        target_folder (str): The target folder for the downloaded reference sequences.
        target_fasta (str): The target fasta file name for the downloaded reference sequences.
        target_fasta_dir (str): The directory (starting at target_folder) for the target fasta file.
        output_dir (str): The directory to save the downloaded file.

    Returns:
        None
    """
    target_path = Path(target_folder)
    with requests.get(url, stream=True, timeout=30) as r:
        r.raise_for_status() # If the server is down, returns 4xx/5xx error
        with target_path.open("wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)

    # Extract the downloaded tar.gz file
    with tarfile.open(target_path, "r:gz") as tar:
        tar.extractall(path=output_dir)

    base_folder = target_folder.split(".tar.gz")[0] # 'gg_13_8_otus'
    path_fasta = pjoin(output_dir, target_fasta_dir, target_fasta)

    # Use only path_fasta; raise an error if the file is not found
    if os.path.exists(path_fasta):
        extracted_fasta_path = path_fasta
        print("✔ Found greengenes_97_otu FASTA in downloaded folder", extracted_fasta_path)
    else:
        raise FileNotFoundError(
            f"Could not locate {target_fasta} after extraction (looked at {path_fasta})."
        )

    if os.path.exists(path_fasta):
        extracted_fasta_path = path_fasta
        print("✔ Found FASTA at path_a:", extracted_fasta_path)
    else:
        raise FileNotFoundError(f"Could not locate {target_fasta} after extraction.")

    # Move fasta to OUTPUT_SCRIPT_FOLDER_DIR
    shutil.move(extracted_fasta_path, pjoin(output_dir, target_fasta))

    # Clean-up: remove extracted folder and archive
    shutil.rmtree(pjoin(output_dir, base_folder), ignore_errors=True)
    try:
        target_path.unlink()
    except Exception:
        pass

    print(f"Downloaded GreenGenes 97% OTU reference sequences to {pjoin(output_dir, target_fasta)}.")
    return

def trim_reference_sequences_to_V4_region(reference_qza, primer_f, primer_r):
    """Trim the reference sequences to the V4 region using primers and feature-classifier.
    
    Args:
        reference_qza (str): Path to the reference sequences QIIME 2 artifact.
        primer_f (str): Forward primer sequence for the V4 region.
        primer_r (str): Reverse primer sequence for the V4 region.
        
    Returns:
        None
    """
    subprocess.run([
        "qiime", "feature-classifier", "extract-reads",
        "--i-sequences", reference_qza,
        "--p-f-primer", primer_f,
        "--p-r-primer", primer_r,
        "--o-reads", reference_qza.split(".qza")[0] + "_trimmed_V4.qza"
    ], check=True)
    return

def cluster_otus(feature_input_dir, sample_id, reference_qza, otu_output_dir,
                 p_ident="0.97", strand="both"):
    """Use QIIME 2's vsearch plugin to perform closed-reference clustering of OTUs.

    Source: https://docs.qiime2.org/2024.10/tutorials/otu-clustering/

    Args:
        sample_id (str): Sample ID for the artifacts.
        reference_qza (str): Path to the reference sequences QIIME 2 artifact.
        output_dir (str): Directory for the clustered artifacts.
        p_ident (str): Percent identity for clustering (default is "0.97" for 97% identity).
        strand (str): Orientation handling for matching (‘plus’, ‘minus’, or ‘both’).

    Returns:
        None    
    """
    feature_table_path = pjoin(feature_input_dir, f"{sample_id}_feature_table.qza")
    feature_data_path = pjoin(feature_input_dir, f"{sample_id}_feature_data.qza")
    # Example:
    #  qiime vsearch cluster-features-closed-reference \
    #   --i-table pjoin(feature_table_path) \
    #   --i-sequences pjoin(feature_data_path) \
    #   --i-reference-sequences reference_qza \
    #   --p-perc-identity p_ident \
    #   --o-clustered-table f"{sample_id}-table-cr-{p_ident}.qza" \
    #   --o-clustered-sequences f"{sample_id}-rep-seqs-cr-{p_ident}.qza" \
    #   --o-unmatched-sequences f"{sample_id}-unmatched-cr-{p_ident}.qza"
    subprocess.run([
        "qiime", "vsearch", "cluster-features-closed-reference",
        "--i-table", feature_table_path,
        "--i-sequences", feature_data_path,
        "--i-reference-sequences", reference_qza,
        "--p-perc-identity", p_ident,
        "--p-strand", strand,                     # <-- added
        "--o-clustered-table", pjoin(otu_output_dir, f"{sample_id}-table-cr-{p_ident}.qza"),
        "--o-clustered-sequences", pjoin(otu_output_dir, f"{sample_id}-rep-seqs-cr-{p_ident}.qza"),
        "--o-unmatched-sequences", pjoin(otu_output_dir, f"{sample_id}-unmatched-cr-{p_ident}.qza")
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
v4_primer_f = "GTGCCAGCMGCCGCGGTAA"
v4_primer_r = "GGACTACHVGGGTWTCTAAT"

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
OTU Clustering: Download Reference GreenGenes Sequences
"""
# From QIIME2 OTU Clustering Tutorial: https://docs.qiime2.org/2024.10/tutorials/otu-clustering/
# Notes:
# - OTU clustering in QIIME 2 is currently applied to a FeatureTable[Frequency] artifact and a FeatureData[Sequence] artifact. 
# - These artifacts can come from a variety of analysis pipelines, including qiime vsearch dereplicate-sequences (illustrated above), qiime dada2 denoise-*, qiime deblur denoise-*, or one of the clustering processes illustrated below (for example, to recluster data at a lower percent identity).
# - The sequences in the FeatureData[Sequence] artifact are clustered against one another (in de novo clustering) or a reference database (in closed-reference clustering), and then features in the FeatureTable are collapsed, resulting in new features that are clusters of the input features.

# For this study:
# - we generated the FeatureTable[Frequency] and FeatureData[Sequence] artifacts in the previous step (dereplicate-sequences)
# - we will cluster against a reference database (GreenGenes)

# Download the GreenGenes 97% OTU reference sequences file if REFERENCE_GREENGENES_97_QZA_FILENAME does not exist in OUTPUT_SCRIPT_FOLDER_DIR yet
reference_sequences_path = pjoin(OUTPUT_SCRIPT_FOLDER_DIR, FASTA_FILENAME_GREENGENES_97)
reference_sequences_qza_path = pjoin(OUTPUT_SCRIPT_FOLDER_DIR, REFERENCE_GREENGENES_97_QZA_FILENAME)
if not os.path.exists(reference_sequences_qza_path):
    print("Downloading GreenGenes 97% OTU reference sequences...")
    download_greengenes_97_reference(URL_GREENGENES_97, TARGET_GREENGENES_97, FASTA_FILENAME_GREENGENES_97, DOWNLOADED_FASTA_DIR_GREENGENES_97, OUTPUT_SCRIPT_FOLDER_DIR)
    # Import the 97% fasta to a QIIME 2 Artifact (store in OUTPUT_SCRIPT_FOLDER_DIR)
    import_fna_as_artifact(
        reference_sequences_path,
        REFERENCE_GREENGENES_97_QZA_FILENAME.split(".qza")[0], 
        OUTPUT_SCRIPT_FOLDER_DIR,
        semantic_type="FeatureData[Sequence]"
    )
    print(f"Downloaded and imported GreenGenes 97% OTU reference sequences to {reference_sequences_qza_path}.")
else:
    print(f"GreenGenes 97% OTU reference sequences already exist at {reference_sequences_qza_path}. Skipping download.")

# Trim the reference sequences to the V4 region using primers
reference_qza_v4_filename = REFERENCE_GREENGENES_97_QZA_FILENAME.split(".qza")[0] + "_trimmed_V4.qza"
if not os.path.exists(pjoin(OUTPUT_SCRIPT_FOLDER_DIR, reference_qza_v4_filename)):
    print("Trimming Greengenes reference sequences to V4 region...")
    trim_reference_sequences_to_V4_region(reference_sequences_qza_path, v4_primer_f, v4_primer_r)
    print(f"Trimmed reference sequences saved to {pjoin(OUTPUT_SCRIPT_FOLDER_DIR, reference_qza_v4_filename)}.")

# If it does not exist, create the output folder for OTU clustering results
otu_clustering_output_path = pjoin(OUTPUT_SCRIPT_FOLDER_DIR, "OTU_clustering")
if not os.path.exists(otu_clustering_output_path):
    os.makedirs(otu_clustering_output_path)


"""
OTU Clustering
"""
# Loop over each sample ID in the metadata
# Perform closed-reference clustering using the vsearch plugin in QIIME 2.
reference_qza_v4_path = pjoin(OUTPUT_SCRIPT_FOLDER_DIR, reference_qza_v4_filename)
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
        cluster_otus(
            feature_table_output_path,
            sample_id,
            reference_qza_v4_path,
            otu_clustering_output_path,
            p_ident="0.97",
            strand="both"          # ensure both orientations are considered
        )
        print(f"Successfully clustered OTUs for {sample_id}.")
    except subprocess.CalledProcessError as e:
        print(f"Error clustering OTUs for {sample_id}: {e}")
