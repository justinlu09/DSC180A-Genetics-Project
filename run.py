#! /usr/bin/env python
import sys
import os
import json
sys.path.insert(0, 'src')
from process import quality_check, check_fastqc, clean_adapters, alignment
from etl import get_data
from analysis import generate_gene_mat
from test import test



def main(targets):
    if 'process' in targets:
        with open('config/process-params.json') as fh:
            process_cfg = json.load(fh)
        
        quality_check(process_cfg.get('data_dir'), process_cfg.get('fastqc_path'), process_cfg.get('fq_output_bc'))
        
        failed_checks = check_fastqc(process_cfg.get('fq_output_bc'))
        
        if len(failed_checks) > 0:
            cutadapt_data = clean_adapters(failed_checks, process_cfg.get('cutadapt_output'))
        
    if 'align' in targets:
        with open('config/align-params.json') as fh:
            align_cfg = json.load(fh)
        kallisto_data = alignment(**align_cfg)
        
    
    if 'analysis' in targets:
        with open('config/analysis-params.json') as fh:
            analysis_cfg = json.load(fh)
        
        gene_matrix_data = generate_gene_mat(**analysis_cfg)
        
    if 'data' in targets:
        with open('config/data-params.json') as fh:
            data_cfg = json.load(fh)
        
        getting_data = get_data(**data_cfg)
        
    if 'test' in targets:
        with open('config/test-params.json') as fh:
            test_cfg = json.load(fh)

        test_out = test(**test_cfg)
     
    
#     if 'report' in targets:
#         with open('
        
    
    
    
    return
        

    
    
if __name__ == '__main__':
    targets = sys.argv[1]
    main(targets)