# Meta-genomic and -transcriptomic Analyses for the Effect of Plant- vs. Animal-based Diets on the Human Gut Microbiome
*A Bioinformatics Portfolio Project*
- Program: Open Bootcamp Collective - Bioinformatics
- Team Members: Lazarina Butkovich, LaShanda Williams, Karl Lundquist, Hitesh Davuluri

## Summary 
- For this project, we replicated some metagenomic and metatranscriptomic analyses of [David et al. (2014)](https://doi.org/10.1038/nature12820).
- Key deliverables:
    - Operational Taxonomic Unit (OTU) Clustering
    - Taxonomic Assignment
    - Alpha and Beta Diversity Metrics
- Main tools:
    - [QIIME2](https://qiime2.org/)
    - [MG-RAST API](https://help.mg-rast.org/api.html)
- The Open Bootcamp Collective [Group Presentation Slides](https://docs.google.com/presentation/d/18d-O28BX4lNezXJkaLLBPzvQ4kP8MUm4zsZS8UC2yn4/edit?usp=sharing) are provided.

## Background
- Previous research in animal models show rapid microbiome shifts with diet changes, but prior to [David et al. (2014)](https://doi.org/10.1038/nature12820), human studies examinded longer timescales.
- In their study, David et al. collected data from 11 healthy adults (aged 21-33 years old) in a crossover study with two diet arms:
    - (1) A plant-based diet with grains, legumes, fruits, vegetables
    - (2) An animal-based diet with meats, eggs, cheeses
- Each diet was consumed for 5 days, separated by baseline and washout periods.
- Fecal samples were collected daily, alongside additional data ()
- For the metagenomic analysis of bacteria in fecal samples, the V4 (high variability) region of 16S rRNA was PCR-amplified. Note that amplicon or targeted metagenomics is distinct from shotgun metagenomics.

## Script Overview for Amplicon Metagenomic Analysis
1. download_fna_files.py
    - Downloads sequence files using [MG-RAST API](https://help.mg-rast.org/api.html)
2. cluster_OTUs_from_fna.py
    - Imports data and metadata for QIIME2 usage
    - Generates "feature frequency" tables per sample, then merges the tables
    - Performs "de novo" OTU clustering over all samples
3. generate_metagenomic_statistics.py
    - Classifies taxonomy of OTUs
    - Builds phylogenetic tree
    - Generates diverstiy statistics (in progress)

- Requirements
    - requests>=2.28.0
    - pandas>=1.3.0
    - openpyxl>=3.1.0
    - numpy>=1.20.0
    - qiime2-amplicon>=2025.4

## Script Overview for Metatranscriptomic Analysis
(in progress)