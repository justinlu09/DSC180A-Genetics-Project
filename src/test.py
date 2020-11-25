import os
import shutil
from process import quality_check, check_fastqc, clean_adapters, alignment
from analysis import generate_gene_mat

def test(test_data, fastqc_path, kallisto_path, kallisto_idx, outdir):
    out = os.path.join(outdir, 'out')
    if (os.path.exists(out) and os.path.isdir(out)):
        shutil.rmtree(out)
    os.mkdir(out)
    
    if ('.ipynb_checkpoints' in os.listdir(test_data)):
        shutil.rmtree(os.path.join(test_data, '.ipynb_checkpoints'))
    print('Extracting 10000 lines from 8 fastqc files from the DSLP brain region for the different disorders & controls as test data')
    print('Starting pipeline on test data...')
    print('Running FastQC on test data...')
    fastqc_outdir = os.path.join(out, 'fastqc_test_out')
    os.mkdir(fastqc_outdir)
    quality_check(test_data, fastqc_path, fastqc_outdir)
    
    print('Checking FastQC outputs...')
    failed_checks = check_fastqc(fastqc_outdir)
    
    print('Running Kallisto in paired-end mode on test data...')
    kallisto_outdir = os.path.join(out, 'kallisto_test_out')
    os.mkdir(kallisto_outdir)
    alignment(test_data, kallisto_path, kallisto_idx, kallisto_outdir)
    
    print('Generating gene matrix from test data...')
    gene_matrix_outdir = os.path.join(out, 'gene_matrix_test_out.csv')
    generate_gene_mat(kallisto_outdir, gene_matrix_outdir)
    
    print('Gene matrix generated from test data can now be found at ./test/out/gene_matrix_test_out.csv')
    return
    