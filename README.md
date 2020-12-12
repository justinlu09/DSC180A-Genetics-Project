# Profiling of Schizophrenia, Bipolar Disorder, and Major Depressive Disorder in Postmortem Brain Tissues

This project demonstrates how researchers used RNA-sequencing of dissected post-mortem tissues in 3 groups of 24 patients, each with schizophrenia, bipolar disorder, major depressive disorder, along with a control group of 24 individuals. This will contain the instructions on how to use the code to extract from the `/datasets/srp073813` folder which contains `SRRXXXXXXX_1.fastq.gz` and `SRRXXXXXXX_2.fastq.gz` files and use the data files for our study. These zipped `.fastq` files are the raw data files of combinations of adenines (A), thymine (T), guanine (G), and cytosine (C) RNA-sequencing reads of individual samples. The code runs quality control, adapter removal (if necessary), and alignment of each pair of samples in order to quantify gene expression. After preprocessing the data and aligning the genes, a gene matrix is created from those outputs to be utilized for DESeq2 for differential gene expression analysis. DESeq2 normalizes the counts and makes comparisons of gene expression between two groups; for our purposes, we made 9 comparisons, one comparison of disorder versus control for each of the 3 disorders in each of the 3 brain regions. 

## Running the project
* Use the command `launch.sh -i ucsdets/dsc180a-genetics:fa20` in order to have the necessary software (e.g., `FastQC`, `cutadapt`, `Kallisto`, `DESeq2`) to run data processing and analysis

## Building the project using `run.py`
* Use the command `python run.py process` to quality check and process the data (using `FastQC` and `cutadapt`)
* Use the command `python run.py align` to align the pairs (using `Kallisto`)
* Use the command `python run.py analysis` to start the analysis of the `Kallisto` outputs by creating a `gene_matrix.csv` and splitting the matrix into separate CSVs for DESeq2 gene expression analysis
* Use the command `python run.py data` to create the data folder (this data folder includes all the necessary data for our project); it is created through a symbolic link to our data directory in the /teams directory (only if you are on UC San Diego DSMLP)
* Use the command `python run.py test` to run the pipeline above on a subset of the data; this also creates a report.html that outlines our pipeline
* Use the command `python run.py all` to run the pipeline above on the whole dataset
