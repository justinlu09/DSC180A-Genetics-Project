#! /usr/bin/env python
import sys
import os
import json
sys.path.insert(0, 'src')
from etl import quality_check, clean_adapters, alignment, run_test
from analysis import generate_gene_mat



def main(targets):
    if 'process' in targets:
        with open('config/process-params.json') as fh:
            data_cfg = json.load(fh)
        
        fastq_data_b = quality_check(data_cfg.get('data_dir'), data_cfg.get('fastqc_path'), data_cfg.get('fq_output_bc'))
        
        if len(fastq_data_b) > 0:
            cutadapt_data = clean_adapters(fastq_data_b, data_cfg.get('cutadapt_output'))
        
    if 'align' in targets:
        with open('config/align-params.json') as fh:
            align_cfg = json.load(fh)
        kallisto_data = alignment(**align_cfg)
        
    
    if 'analysis' in targets:
        with open('config/analysis-params.json') as fh:
            analysis_cfg = json.load(fh)
        
        gene_matrix_data = generate_gene_mat(**analysis_cfg)
        
        
    if 'test' in targets:
        with open('config/test-params.json') as fh:
            test_cfg = json.load(fh)
        
        test_out = run_test(**test_cfg)
        test_gene_mat_out = generate_gene_mat(test_cfg.get('kallisto_out'), test_cfg.get('gene_matrix_out'))
     
    
#     if 'report' in targets:
#         with open('
        
    
    
    
    return
        

    
    
if __name__ == '__main__':
    targets = sys.argv[1]
    main(targets)