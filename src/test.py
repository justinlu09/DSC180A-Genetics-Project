import os
import shutil
from process import quality_check, check_fastqc, clean_adapters, alignment
from analysis import generate_gene_mat

def test(test_data, fastqc_path, outdir_fastqc, kallisto_path, kallisto_idx, outdir_kallisto, outdir_gene_matrix):
    
    if ('.ipynb_checkpoints' in os.listdir(test_data)):
        shutil.rmtree(os.path.join(test_data, '.ipynb_checkpoints'))
    print('Extracting 10000 lines from 8 fastqc files from the DSLP brain region for the different disorders & controls as test data')
    print('Starting pipeline on test data...')
    print('Running FastQC on test data...')
    quality_check(test_data, fastqc_path, outdir_fastqc)
    
    print('Checking FastQC outputs...')
    failed_checks = check_fastqc(outdir_fastqc)
    
    print('Running Kallisto in paired-end mode on test data...')
    alignment(test_data, kallisto_path, kallisto_idx, outdir_kallisto)
    
    print('Generating gene matrix from test data...')
    generate_gene_mat(outdir_kallisto, outdir_gene_matrix)
    
    print('Gene matrix generated from test data can now be found at ./test/out/gene_matrix_test_out.csv')
    return
    