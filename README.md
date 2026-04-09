# Glycoproteomics Data Processing Pipeline

## Introduction
Glycoproteomics is a branch of proteomics that focuses on the study of glycoproteins, which are proteins that have carbohydrate groups attached to them. The glycosylation of proteins plays a crucial role in various biological processes, including cell signaling, immune response, and disease progression. This pipeline is designed to provide a comprehensive framework for processing glycoproteomics data, integrating various computational methods and workflows to streamline analysis and interpretation of results.

## Scripts Overview
This pipeline includes the following Python scripts:

1. glyco_pipeline.py: Main script for running the glycoproteomics analysis, which includes essential preprocessing and analysis steps.

2. glyco_pipeline_batch.py: Script for processing multiple datasets in a batch mode, allowing for high-throughput analysis.

3. site_correlation.py: This script focuses on identifying and correlating glycosylation sites based on their glycosylation profiles.

4. MultiDimensionCorrel.py: A comprehensive tool that allows for multidimensional correlation analysis of glycoproteomics data to uncover hidden relationships at a site-to-site basis.

5. calculate_fastgsite.py: This script calculates glycosylation sites along the protein sequence for downstream analysis.

## Workflow Examples
To get started, here are some basic workflows:
Single Sample Analysis: Run `glyco_pipeline.py` with a specified input file.
Batch Processing: Use `glyco_pipeline_batch.py` to process all files in a given directory.

## Data Format Requirements
The pipeline requires input data in .xlsx unless otherwise stated.
Ensure that necessary columns are present (e.g., ProteinID, GlycanStructure). Column labels, unless otherwise stated, are formatted based on the output of ProteomeDiscoverer with Byonic plugin.

## Interactive Features
The pipeline includes several interactive features:
Progress trackers for long-running jobs.
Options to visualize intermediate results directly from the command line.

## Output Interpretation
The output will include:
Summary reports of processed data.
Visualizations of glycoprotein structures and correlations.
Detailed logs for error troubleshooting.

## Citations
Please list my name and attach a link to this repository when using this pipeline.

For any further questions or issues, feel free to raise an issue on the repository or contact the maintainer directly!

Last updated on 4/9/2026
