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

def build_phylogenetic_tree(merged_seqs_qza, tree_out_path, threads=1):
    """Build a phylogenetic tree from the merged feature sequences.
    
    Args:
        merged_seqs_qza (str): Path to the merged feature sequences artifact.
        tree_out_path (str): Output path for the phylogenetic tree artifact.
        threads (int): Number of threads to use for alignment (default: 1).
    Returns:
        dict: Dictionary with paths to all generated tree files
    """
    try:
        # For large datasets, use parttree algorithm which is more memory-efficient
        # See: https://docs.qiime2.org/2024.10/plugins/available/phylogeny/align-to-tree-mafft-fasttree/
        print(f"Building phylogenetic tree using {threads} threads with parttree algorithm...")
        
        # Define the expected output files
        aligned_path = tree_out_path.replace(".qza", "_aligned.qza")
        masked_path = tree_out_path.replace(".qza", "_masked.qza")
        rooted_path = tree_out_path.replace(".qza", "_rooted.qza")
        
        # Run the QIIME 2 command
        subprocess.run([
            "qiime", "phylogeny", "align-to-tree-mafft-fasttree",
            "--i-sequences", merged_seqs_qza,
            "--o-alignment", aligned_path,
            "--o-masked-alignment", masked_path,
            "--o-tree", tree_out_path,
            "--o-rooted-tree", rooted_path,
            "--p-n-threads", str(threads),
            "--p-parttree"  # Use parttree algorithm for large datasets
        ], check=True)
        
        # Verify the output files exist
        output_files = {
            "tree": tree_out_path,
            "rooted_tree": rooted_path,
            "alignment": aligned_path,
            "masked_alignment": masked_path
        }
        
        for file_type, file_path in output_files.items():
            if os.path.exists(file_path):
                print(f"✓ {file_type} created successfully: {file_path}")
            else:
                print(f"⚠ Warning: {file_type} was not created at: {file_path}")
        
        print("Phylogenetic tree building completed successfully.")
        return output_files
        
    except subprocess.CalledProcessError as e:
        print(f"Error in phylogenetic tree building: {e}")
        print("If memory error persists, consider using fasttree directly or reducing your dataset size.")
        raise

def run_core_diversity_analysis(feature_table_qza, rooted_tree_qza, metadata_path, taxonomy_qza, output_dir, sampling_depth=1000, threads=1):
    """Run core diversity analysis using phylogenetic metrics.
    
    Args:
        feature_table_qza (str): Path to the feature table artifact.
        rooted_tree_qza (str): Path to the rooted tree artifact.
        metadata_path (str): Path to the metadata file (TSV format).
        taxonomy_qza (str): Path to the taxonomy artifact.
        output_dir (str): Output directory for diversity analysis results.
        sampling_depth (int): Sampling depth for rarefaction.
        threads (int): Number of threads to use.
    Returns:
        None
    """
    try:
        # Convert Excel metadata to TSV if needed
        if metadata_path.endswith('.xlsx'):
            metadata_df = pd.read_excel(metadata_path)
            tsv_path = metadata_path.replace('.xlsx', '.tsv')
            metadata_df.to_csv(tsv_path, sep='\t', index=False)
            metadata_path = tsv_path
        
        # If the output directory exists and we need to rerun, delete it first
        # Make sure the directory is completely removed before proceeding
        if os.path.exists(output_dir):
            print(f"Removing existing diversity analysis directory: {output_dir}")
            shutil.rmtree(output_dir, ignore_errors=True)
            # Double-check it's gone
            if os.path.exists(output_dir):
                # Try with system command if shutil.rmtree fails
                os.system(f'rm -rf "{output_dir}"')
            # Final check
            if os.path.exists(output_dir):
                raise RuntimeError(f"Could not remove directory: {output_dir}")
                
        # Create a fresh output directory
        os.makedirs(output_dir)
            
        print(f"Running core diversity analysis with sampling depth {sampling_depth}...")
        
        # Run QIIME 2 core-metrics-phylogenetic with added verbose flag
        # and explicit naming of all output files
        alpha_metrics_dir = os.path.join(output_dir, "alpha_metrics")
        beta_metrics_dir = os.path.join(output_dir, "beta_metrics")
        os.makedirs(alpha_metrics_dir, exist_ok=True)
        os.makedirs(beta_metrics_dir, exist_ok=True)
        
        # Create paths for all required output files
        rarefied_table = os.path.join(output_dir, "rarefied_table.qza")
        faith_pd = os.path.join(alpha_metrics_dir, "faith_pd_vector.qza")
        observed_features = os.path.join(alpha_metrics_dir, "observed_features_vector.qza")
        shannon = os.path.join(alpha_metrics_dir, "shannon_vector.qza")
        evenness = os.path.join(alpha_metrics_dir, "evenness_vector.qza")
        
        unweighted_unifrac = os.path.join(beta_metrics_dir, "unweighted_unifrac_distance_matrix.qza")
        weighted_unifrac = os.path.join(beta_metrics_dir, "weighted_unifrac_distance_matrix.qza")
        jaccard = os.path.join(beta_metrics_dir, "jaccard_distance_matrix.qza")
        bray_curtis = os.path.join(beta_metrics_dir, "bray_curtis_distance_matrix.qza")
        
        unweighted_unifrac_pcoa = os.path.join(beta_metrics_dir, "unweighted_unifrac_pcoa_results.qza")
        weighted_unifrac_pcoa = os.path.join(beta_metrics_dir, "weighted_unifrac_pcoa_results.qza")
        jaccard_pcoa = os.path.join(beta_metrics_dir, "jaccard_pcoa_results.qza")
        bray_curtis_pcoa = os.path.join(beta_metrics_dir, "bray_curtis_pcoa_results.qza")
        
        unweighted_unifrac_emperor = os.path.join(beta_metrics_dir, "unweighted_unifrac_emperor.qzv")
        weighted_unifrac_emperor = os.path.join(beta_metrics_dir, "weighted_unifrac_emperor.qzv")
        jaccard_emperor = os.path.join(beta_metrics_dir, "jaccard_emperor.qzv")
        bray_curtis_emperor = os.path.join(beta_metrics_dir, "bray_curtis_emperor.qzv")
        
        # Run with explicit output paths instead of output_dir
        subprocess.run([
            "qiime", "diversity", "core-metrics-phylogenetic",
            "--i-phylogeny", rooted_tree_qza,
            "--i-table", feature_table_qza,
            "--p-sampling-depth", str(sampling_depth),
            "--m-metadata-file", metadata_path,
            "--p-n-jobs-or-threads", str(threads),
            "--o-rarefied-table", rarefied_table,
            "--o-faith-pd-vector", faith_pd,
            "--o-observed-features-vector", observed_features,
            "--o-shannon-vector", shannon,
            "--o-evenness-vector", evenness,
            "--o-unweighted-unifrac-distance-matrix", unweighted_unifrac,
            "--o-weighted-unifrac-distance-matrix", weighted_unifrac,
            "--o-jaccard-distance-matrix", jaccard,
            "--o-bray-curtis-distance-matrix", bray_curtis,
            "--o-unweighted-unifrac-pcoa-results", unweighted_unifrac_pcoa,
            "--o-weighted-unifrac-pcoa-results", weighted_unifrac_pcoa,
            "--o-jaccard-pcoa-results", jaccard_pcoa,
            "--o-bray-curtis-pcoa-results", bray_curtis_pcoa,
            "--o-unweighted-unifrac-emperor", unweighted_unifrac_emperor,
            "--o-weighted-unifrac-emperor", weighted_unifrac_emperor,
            "--o-jaccard-emperor", jaccard_emperor,
            "--o-bray-curtis-emperor", bray_curtis_emperor,
            "--verbose"
        ], check=True)
        
        # Generate alpha diversity significance tests
        for metric in ["faith_pd", "evenness", "shannon"]:
            # Use the correct path from the alpha_metrics subdirectory
            metric_path = os.path.join(alpha_metrics_dir, f"{metric}_vector.qza")
            viz_path = os.path.join(output_dir, f"{metric}_significance.qzv")
            subprocess.run([
                "qiime", "diversity", "alpha-group-significance",
                "--i-alpha-diversity", metric_path,
                "--m-metadata-file", metadata_path,
                "--o-visualization", viz_path
            ], check=True)
        
        # Generate beta diversity significance tests for categorical variables
        # You can customize which metadata columns to test
        metadata_df = pd.read_csv(metadata_path, sep='\t')
        categorical_columns = metadata_df.select_dtypes(include=['object', 'category']).columns.tolist()
        # Remove SampleID column if present
        if SAMPLE_ID_COL_NAME in categorical_columns:
            categorical_columns.remove(SAMPLE_ID_COL_NAME)
            
        for metric in ["unweighted_unifrac", "weighted_unifrac", "jaccard", "bray_curtis"]:
            # Use the correct paths from the beta_metrics subdirectory
            distance_path = os.path.join(beta_metrics_dir, f"{metric}_distance_matrix.qza")
            
            # Generate PCoA plots
            pcoa_path = os.path.join(beta_metrics_dir, f"{metric}_pcoa_results.qza")
            viz_path = os.path.join(beta_metrics_dir, f"{metric}_emperor.qzv")
            subprocess.run([
                "qiime", "emperor", "plot",
                "--i-pcoa", pcoa_path,
                "--m-metadata-file", metadata_path,
                "--o-visualization", viz_path
            ], check=True)
            
            # Test each categorical variable
            for column in categorical_columns[:5]:  # Limit to first 5 columns to avoid excessive computation
                viz_path = os.path.join(beta_metrics_dir, f"{metric}_{column}_significance.qzv")
                subprocess.run([
                    "qiime", "diversity", "beta-group-significance",
                    "--i-distance-matrix", distance_path,
                    "--m-metadata-file", metadata_path,
                    "--m-metadata-column", column,
                    "--p-pairwise",
                    "--o-visualization", viz_path
                ], check=True)
        
        print(f"Core diversity analysis completed. Results saved to: {output_dir}")
        
    except subprocess.CalledProcessError as e:
        print(f"Error in core diversity analysis: {e}")
        raise
    except Exception as e:
        print(f"Error in core diversity analysis: {e}")
        raise

def visualize_taxonomy_barplot(feature_table_qza, taxonomy_qza, metadata_path, output_path):
    """Generate a barplot visualization of taxonomy composition.
    
    Args:
        feature_table_qza (str): Path to the feature table artifact.
        taxonomy_qza (str): Path to the taxonomy artifact.
        metadata_path (str): Path to the metadata file.
        output_path (str): Output path for the barplot visualization.
    Returns:
        None
    """
    try:
        # Convert Excel metadata to TSV if needed
        if metadata_path.endswith('.xlsx'):
            metadata_df = pd.read_excel(metadata_path)
            tsv_path = metadata_path.replace('.xlsx', '.tsv')
            metadata_df.to_csv(tsv_path, sep='\t', index=False)
            metadata_path = tsv_path
            
        # Run QIIME 2 taxa barplot
        subprocess.run([
            "qiime", "taxa", "barplot",
            "--i-table", feature_table_qza,
            "--i-taxonomy", taxonomy_qza,
            "--m-metadata-file", metadata_path,
            "--o-visualization", output_path
        ], check=True)
        
        print(f"Taxonomy barplot created at: {output_path}")
    
    except subprocess.CalledProcessError as e:
        print(f"Error in taxonomy barplot creation: {e}")
        raise
    except Exception as e:
        print(f"Error in taxonomy barplot creation: {e}")
        raise


"""
Values
"""
# input folder holds the metadata file
INPUT_FOLDER = "input"
# "metadata_16S.xlsx" for full metadata, "metadata_test.xlsx" for simpler metadata testing
METADATA_FILENAME_XLSX = "metadata_test.xlsx"
METADATA_PATH = pjoin(INPUT_FOLDER, METADATA_FILENAME_XLSX)
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
Taxonomic Composition Visualization
"""
# QIIME 2 taxa barplot documentation: https://docs.qiime2.org/2024.10/plugins/available/taxa/barplot/

# Generate taxonomy barplot
feature_table_qza_path = pjoin(CLUSTER_OTUS_OUTPUT_FOLDER_DIR, "merged-table-cr-0.97_" + RUN_DESCRIPTOR + ".qza")
taxonomy_barplot_path = pjoin(OUTPUT_SCRIPT_FOLDER_DIR, "taxonomy_barplot_" + RUN_DESCRIPTOR + ".qzv")

if not os.path.exists(taxonomy_barplot_path) or RERUN_EXISTING_DATA:
    print("Generating taxonomy barplot visualization...")
    visualize_taxonomy_barplot(
        feature_table_qza=feature_table_qza_path,
        taxonomy_qza=taxonomy_out_path,
        metadata_path=pjoin(INPUT_FOLDER, METADATA_FILENAME_XLSX),
        output_path=taxonomy_barplot_path
    )
else:
    print(f"Taxonomy barplot already exists at: {taxonomy_barplot_path}. Skipping barplot generation step.")


"""
Infer Phylogenetic Tree
"""
# QIIME 2 q2-phylogeny documentation: https://docs.qiime2.org/2024.10/tutorials/phylogeny/
tree_out_path = pjoin(OUTPUT_SCRIPT_FOLDER_DIR, "phylogenetic_tree_" + RUN_DESCRIPTOR + ".qza")
rooted_tree_out_path = tree_out_path.replace(".qza", "_rooted.qza")

if not os.path.exists(rooted_tree_out_path) or RERUN_EXISTING_DATA:
    # Build phylogenetic tree from the merged feature sequences
    # Use a single thread for memory efficiency, but you can increase this if you have sufficient RAM
    tree_files = build_phylogenetic_tree(MERGED_SEQS_FEATURE_QZA_PATH, tree_out_path, threads=1)
    
    # Update the rooted_tree_out_path with the actual path from the result
    if "rooted_tree" in tree_files:
        rooted_tree_out_path = tree_files["rooted_tree"]
        
    print(f"Phylogenetic tree building completed. Output saved to: {rooted_tree_out_path}")
else:
    print(f"Phylogenetic tree already exists at: {rooted_tree_out_path}. Skipping tree building step.")


"""
Core Diversity Analysis
"""
# QIIME 2 diversity core-metrics-phylogenetic documentation: https://docs.qiime2.org/2024.10/plugins/available/diversity/core-metrics-phylogenetic/
# Perform diversity analysis with every run

# Need both feature table and rooted tree for diversity analysis
feature_table_qza_path = pjoin(CLUSTER_OTUS_OUTPUT_FOLDER_DIR, "merged-table-cr-0.97_" + RUN_DESCRIPTOR + ".qza")
if not os.path.exists(feature_table_qza_path):
    print(f"Error: Feature table not found at {feature_table_qza_path}")
    print("Make sure to run the OTU clustering script first to generate the feature table.")
else:
    # Set up paths for core diversity analysis
    diversity_out_dir = pjoin(OUTPUT_SCRIPT_FOLDER_DIR, "core_diversity_" + RUN_DESCRIPTOR)
    
    # Determine an appropriate sampling depth
    # A higher value gives better resolution but may exclude samples with fewer sequences
    # A lower value includes more samples but with less resolution
    # This would ideally be determined by looking at the feature table summary
    sampling_depth = 1000  # Default conservative value
    
    # Run core diversity analysis
    run_core_diversity_analysis(
        feature_table_qza=feature_table_qza_path,
        rooted_tree_qza=rooted_tree_out_path,
        metadata_path=METADATA_PATH,
        taxonomy_qza=taxonomy_out_path,
        output_dir=diversity_out_dir,
        sampling_depth=sampling_depth,
        threads=1  # Use a single thread for memory efficiency
    )
    


# """
# Differential-abundance & compositional stats
# """
# # Consider ANCOM II method (qiime composition ancom) with add-pseudocounts option for zero counts
# # Or, consider Songbird method (qiime songbird multinomial)


# """
# Functional Inference
# """
# # Consider PICRUSt2 or Tax4Fun2 to predict KEGG pathways or MetaCyc functions from 16S OTUs.


# """
# Time-series Analysis
# """
# # Consider qiime longitudinal plugin for time-series analysis of microbiome data.

# """
# Print instructions for viewing results
# """
# print("\n" + "="*80)
# print("RESULTS SUMMARY")
# print("="*80)
# print("Analysis completed. Here's how to view the results:")

# # Add instructions for viewing the phylogenetic tree
# print("\n1. To view the phylogenetic tree:")
# tree_viz_path = rooted_tree_out_path.replace(".qza", ".qzv")
# print(f"   qiime tools view {tree_viz_path}")

# # Continue with other instructions
# print("\n2. To view the taxonomy classification:")
# print(f"   qiime tools view {taxonomy_out_path.replace('.qza', '_tabulated.qzv')}")
# print("\n3. To view the taxonomy barplot (if generated):")
# taxonomy_barplot_path = pjoin(OUTPUT_SCRIPT_FOLDER_DIR, "taxonomy_barplot_" + RUN_DESCRIPTOR + ".qzv")
# print(f"   qiime tools view {taxonomy_barplot_path}")
# print("\n4. To explore the diversity analysis results:")
# diversity_out_dir = pjoin(OUTPUT_SCRIPT_FOLDER_DIR, "core_diversity_" + RUN_DESCRIPTOR)
# print(f"   cd {diversity_out_dir}")
# print("   qiime tools view faith_pd_significance.qzv  # Alpha diversity")
# print("   qiime tools view weighted_unifrac_emperor.qzv  # Beta diversity PCoA plot")

# print("\nAlternatively, you can use QIIME 2 View online:")
# print("1. Go to https://view.qiime2.org")
# print("2. Drag and drop any of these files to view them in your browser:")
# print(f"   - {tree_viz_path}  # Phylogenetic tree visualization")
# print(f"   - {taxonomy_out_path.replace('.qza', '_tabulated.qzv')}  # Taxonomy classification")
# taxonomy_barplot_path = pjoin(OUTPUT_SCRIPT_FOLDER_DIR, "taxonomy_barplot_" + RUN_DESCRIPTOR + ".qzv")
# print(f"   - {taxonomy_barplot_path}  # Taxonomy barplot (if generated)")
# print("   - [Diversity analysis visualizations from the core_diversity folder]")
# print("\nFor .qza files (like the tree data itself):")
# print(f"1. You can open {rooted_tree_out_path} in QIIME 2 View as well")
# print("2. Or export it to Newick format for use in other tree viewers:")
# print(f"   qiime tools export --input-path {rooted_tree_out_path} --output-path exported_tree")
# print("="*80)
