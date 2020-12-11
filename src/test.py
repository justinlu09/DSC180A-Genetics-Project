import os
import shutil
from process import quality_check, check_fastqc, clean_adapters, alignment
from analysis import generate_gene_mat, split_for_comparison, run_deseq

def test(test_data, fastqc_path, run_name, last_html, kallisto_path, kallisto_idx, num_boots, num_threads, gene_naming_table, sra_run_table, chromosomes_needed, outdir, tmp_out_ancg_bpd, tmp_out_ancg_bpd_coldata, tmp_out_ancg_sz, tmp_out_ancg_sz_coldata, tmp_out_ancg_mdd, tmp_out_ancg_mdd_coldata, tmp_out_dlpfc_bpd, tmp_out_dlpfc_bpd_coldata, tmp_out_dlpfc_sz, tmp_out_dlpfc_sz_coldata, tmp_out_dlpfc_mdd, tmp_out_dlpfc_mdd_coldata, tmp_out_nacc_bpd, tmp_out_nacc_bpd_coldata, tmp_out_nacc_sz, tmp_out_nacc_sz_coldata, tmp_out_nacc_mdd, tmp_out_nacc_mdd_coldata, figures_out):
    
    test_output_dir = os.path.join(outdir, 'out')
    
    #remove data dir if one already exists
    if (os.path.exists(test_output_dir) and os.path.isdir(test_output_dir)):
        shutil.rmtree(test_output_dir)
    
    os.mkdir(test_output_dir)
    
    if ('.ipynb_checkpoints' in os.listdir(test_data)):
        shutil.rmtree(os.path.join(test_data, '.ipynb_checkpoints'))
    
    print('Starting pipeline on test data...')
    print('Test data created using `zcat` on raw fastq files, extracting 20 lines of sequences from 24 fastq.gz files')
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
    
    split_for_comparison(outdir_gene_matrix, sra_run_table, deseq_output_dir)
    print('Data needed for DESeq2 can now be found at ./test/tmp/')
    #run_deseq(outdir_gene_matrix, tmp_out_ancg_bpd, tmp_out_ancg_bpd_coldata, tmp_out_ancg_sz, tmp_out_ancg_sz_coldata, tmp_out_ancg_mdd, tmp_out_ancg_mdd_coldata, tmp_out_dlpfc_bpd, tmp_out_dlpfc_bpd_coldata, tmp_out_dlpfc_sz, tmp_out_dlpfc_sz_coldata, tmp_out_dlpfc_mdd, tmp_out_dlpfc_mdd_coldata, tmp_out_nacc_bpd, tmp_out_nacc_bpd_coldata, tmp_out_nacc_sz, tmp_out_nacc_sz_coldata, tmp_out_nacc_mdd, tmp_out_nacc_mdd_coldata, figures_out)
    return
    
def generate_report(notebook_indir):
    os.system("jupyter nbconvert --to html " + notebook_indir)
    print('Please refer to ./notebooks/report.html for more explanation about our test data pipeline')
    