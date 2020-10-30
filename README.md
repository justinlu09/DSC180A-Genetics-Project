# Profiling of Schizophrenia, Bipolar Disorder, and Major Depressive Disorder in Postmortem Brain Tissues

This project demonstrates how researchers used RNA-sequencing of dissected post-mortem tissues in 3 groups of 24 patients each with schizophrenia, bipolar disorder, major depressive disorder, and 24 controls. This will contain the instructions on how to use the code to extract from "srp073813" folder which will contain "SRRXXXXXXX.1" files. A SRA toolkit, fastq-dump, is used to convert these files into data of combinations of adenines (A), thymine (T), guanine (G), and cytosine (C) RNA-sequencing reads of individual samples.

## Running the project
* Use Docker for tools, etc.

## Building the project using `run.py`
* Use the command: `run.py data` to get the "SRRXXXXXXX.1" files into a CSV in the `data/output` folder so we can extract these RNA-seq reads using fastq-dump
  * This downloads data from `data/srp073813` in the directory specified in `config/data-params.json`

