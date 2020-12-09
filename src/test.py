import os
import shutil
from process import quality_check, check_fastqc, clean_adapters, alignment
from analysis import generate_gene_mat, split_for_comparison

def test(test_data, fastqc_path, run_name, last_html, kallisto_path, kallisto_idx, num_boots, num_threads, gene_naming_table, sra_run_table, chromosomes_needed, outdir):
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
    quality_check(test_data, fastqc_path, outdir_fastqc, run_name)
    
    print('Checking FastQC outputs...')
    failed_checks = check_fastqc(outdir_fastqc, last_html)
    
    print('Running Kallisto in paired-end mode on test data...')
    outdir_kallisto = os.path.join(test_output_dir, 'kallisto_test_out')
    os.mkdir(outdir_kallisto)
    alignment(test_data, kallisto_path, kallisto_idx, outdir_kallisto, num_boots, num_threads)
    
    print('Generating gene matrix from test data...')
    outdir_gene_matrix = os.path.join(test_output_dir, 'gene_matrix_test_out.csv')
    generate_gene_mat(gene_naming_table, sra_run_table, outdir_kallisto, outdir_gene_matrix, chromosomes_needed)
    
    print('FastQC outputs can now be found at ./test/out/fastqc_test_out/')
    print('Kallisto outputs can now be found at ./test/out/kallisto_test_out/')
    print('Gene matrix generated from test data can now be found at ./test/out/gene_matrix_test_out.csv')
    
    print('Splitting gene matrix up for comparison using R-based tool, DESeq2...')
    
    deseq_output_dir = os.path.join(outdir, 'tmp')
    
    #remove data dir if one already exists
    if (os.path.exists(deseq_output_dir) and os.path.isdir(deseq_output_dir)):
        shutil.rmtree(deseq_output_dir)
    
    os.mkdir(deseq_output_dir)
    
    split_for_comparison("./test/out/gene_matrix_test_out.csv", sra_run_table, deseq_output_dir)  
    return
    