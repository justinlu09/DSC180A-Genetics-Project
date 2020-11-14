# Profiling of Schizophrenia, Bipolar Disorder, and Major Depressive Disorder in Postmortem Brain Tissues

This project demonstrates how researchers used RNA-sequencing of dissected post-mortem tissues in 3 groups of 24 patients, each with schizophrenia, bipolar disorder, major depressive disorder, along with a control group of 24 individuals. This will contain the instructions on how to use the code to extract from the `/datasets/srp073813` folder which contains `SRRXXXXXXX_1.fastq.gz` and `SRRXXXXXXX_2.fastq.gz` files and use the data files for our study. These zipped `.fastq` files are the raw data files of combinations of adenines (A), thymine (T), guanine (G), and cytosine (C) RNA-sequencing reads of individual samples. At this point, the code runs quality control, adapter removal (if necessary), and alignment of each pair of samples in order to quantify gene expression.

## Running the project
* Use the command `launch-180.sh -i ucsdets/dsc180a-genetics:fa20 -G B04_Genetics` in order to have the necessary software (e.g., `FastQC`, `cutadapt`, `Kallisto`) to run data processing and analysis
* symbolic link between DSMLP and home directory

## Building the project using `run.py`
* Use the command `python run.py process` to quality check and process the data (using `FastQC` and `cutadapt`) - right now, this only does the process on around 10 samples (5 pairs)
* Use the command `python run.py align` to align the pairs (using `Kallisto`)

## Group Contributions
* Justin and Joseph worked together (pair-programming) to work on the data ingestion pipeline over Zoom, adding code to run.py, data-params.json, and etl.py. Joseph wrote the README.md, while Justin revised it. 
