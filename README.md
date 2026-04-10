# Automated N-Glycoproteomics Data Analysis Pipeline

## Overview

This repository contains a suite of custom Python and MATLAB bioinformatics pipelines designed to process, filter, normalize, and visualize large-scale, multi-dimensional liquid chromatography-mass spectrometry (LC-MS) datasets.Specifically, these tools were developed to decode the N-glycoproteome of murine keratinocytes and exosomes, handling complex tasks such as triplicate variance masking, N-glycosylation site mapping, Total Ion Current (TIC) normalization, and 5-dimensional signature correlation.

## Scripts Overview
This pipeline includes the following Python scripts:

1. glyco_pipeline.py: Main script for running the glycoproteomics analysis, which includes essential preprocessing and analysis steps.

2. glyco_pipeline_batch.py: Script for processing multiple datasets in a batch mode, allowing for high-throughput analysis.

3. site_correlation.py: This script focuses on identifying and correlating glycosylation sites based on their glycosylation profiles.

4. MultiDimensionCorrel.py: A comprehensive tool that allows for multidimensional correlation analysis of glycoproteomics data to uncover hidden relationships at a site-to-site basis.

5. calculate_fastgsite.py: This script calculates glycosylation sites along the protein sequence for downstream analysis.

6. count20.m: This script counts the occurances of the 20 amino acids centered around the N-glycosylation site, enabling downstream motif analysis.

## Key Features & Pipeline Modules

1. Batch Data Processing & Normalization (glyco_pipeline_batch.py)
   
   A robust, end-to-end data processing script designed to handle raw output from Proteome Discoverer / Byonic
   
   Dynamic Bookending: Automatically detects and extracts variable sample abundance columns regardless of naming conventions
   
   Smart Triplicate Masking: Calculates Relative Standard Deviation (RSD) per sample and independently masks failed triplicates (RSD > 30%) with NaN rather than dropping entire protein rows, preserving data integrity
   
   Automated N-Glycosylation Site Mapping: Parses peptide sequences and absolute protein starting positions to extract accurate glycosylation sites using Regex.Human-in-the-loop Database Labeling: Pauses the pipeline to detect and export newly discovered proteins for manual verification/Gene Ontology (GO) filtering, then automatically merges them into a master library
   
   OriginLab Integration: Flattens the final normalized data into a mathematically exact Cartesian matrix, ready for instant Top-10 protein/glycan bubble plot visualization
   
2. Multi-Dimensional Glycosylation Site Correlation (site_correlation.py)
   
   A unique analytical tool that calculates the biological similarity between different glycosylation sites based on their complete glycan composition profiles
   
   5D Dot-Product Correlation: Calculates multi-dimensional correlation flattened across 5 distinct composition dimensions (Hex, HexNAc, NeuAc, NeuGc, Fuc)
   
   Automated Data Pivoting: Auto-aligns disparate glycan lists and pads missing dimensions with zeros for mathematically safe comparison
   
   Dynamic Cutoff Histograms: Automatically generates histograms of glycan counts per site before running the correlation, allowing the user to set a smart threshold (e.g., sites with $\ge$ 5 glycans) to prevent false positives
   
   Site Existence Mapping: Generates a binary "CountIf" summary matrix highlighting which specific glycosylation sites are conserved across all biological conditions vs. which are condition-specific
   
3. Rapid Protein Site Calculation (calculate_fastgsite.py)
   
   A lightweight utility script that parses modified peptide sequences (where target Asparagine residues are substituted) to rapidly calculate absolute sequence positions along the master protein.Installation & Requirements
   
   These scripts are designed to run locally using standard Python data science libraries (pandas numpy matplotlib openpyxl)
   
## Author

Xuyao Zeng, PhD

Postdoctoral Fellow at Indiana University, Dept. of Chemistry

Last updated on 4/9/2026

Contact: xuyzeng@iu.edu
