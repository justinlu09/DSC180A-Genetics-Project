import os
import shutil
from process import quality_check, check_fastqc, clean_adapters, alignment
from analysis import generate_gene_mat

def test(test_data, fastqc_path, kallisto_path, kallisto_idx, outdir):
    test_output_dir = os.path.join(outdir, 'out')
    
    #remove data dir if one already exists
    if (os.path.exists(test_output_dir) and os.path.isdir(test_output_dir)):
        shutil.rmtree(test_output_dir)
    
    os.mkdir(test_output_dir)
    
    if ('.ipynb_checkpoints' in os.listdir(test_data)):
        shutil.rmtree(os.path.join(test_data, '.ipynb_checkpoints'))
    
    print('Starting pipeline on test data...')
    print('Test data created using `zcat` on raw fastq files, extracting 10000 lines of sequences from 8 fastq.gz files')
    print('Running FastQC on test data...')
    outdir_fastqc = os.path.join(test_output_dir, 'fastqc_test_out')
    os.mkdir(outdir_fastqc)
    quality_check(test_data, fastqc_path, outdir_fastqc)
    
    print('Checking FastQC outputs...')
    failed_checks = check_fastqc(outdir_fastqc)
    
    print('Running Kallisto in paired-end mode on test data...')
    outdir_kallisto = os.path.join(test_output_dir, 'kallisto_test_out')
    os.mkdir(outdir_kallisto)
    alignment(test_data, kallisto_path, kallisto_idx, outdir_kallisto)
    
    print('Generating gene matrix from test data...')
    outdir_gene_matrix = os.path.join(test_output_dir, 'gene_matrix_test_out.csv')
    generate_gene_mat(outdir_kallisto, outdir_gene_matrix)
    
    print('FastQC outputs can now be found at ./test/out/fastqc_test_out/')
    print('Kallisto outputs can now be found at ./test/out/kallisto_test_out/')
    print('Gene matrix generated from test data can now be found at ./test/out/gene_matrix_test_out.csv')
    return
    