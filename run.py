#! /usr/bin/env python
import sys
import os
import json
sys.path.insert(0, 'src')
from process import quality_check, check_fastqc, clean_adapters, alignment
from etl import get_data
from analysis import generate_gene_mat, split_for_comparison, run_deseq
from test import test



def main(targets):
    if 'process' in targets:
        with open('config/process-params.json') as fh:
            process_cfg = json.load(fh)
        
        quality_check(process_cfg.get('data_dir'), process_cfg.get('fastqc_path'), process_cfg.get('fq_output_bc'), process_cfg.get('run_name'))
        
        failed_checks = check_fastqc(process_cfg.get('fq_output_bc'), process_cfg.get('last_html'))
        
        if len(failed_checks) > 0:
            cutadapt_data = clean_adapters(failed_checks, process_cfg.get('cutadapt_output'), process_cfg.get('adapter1'), process_cfg.get('adapter2'), process_cfg.get('num_cores'))
        
    if 'align' in targets:
        with open('config/align-params.json') as fh:
            align_cfg = json.load(fh)
        kallisto_data = alignment(**align_cfg)
        
    
    if 'analysis' in targets:
        with open('config/analysis-params.json') as fh:
            analysis_cfg = json.load(fh)
        with open('config/deseq-params.json') as fh:
            deseq_cfg = json.load(fh)
        
        #gene_matrix_data = generate_gene_mat(analysis_cfg.get("gene_naming_table"), analysis_cfg.get("sra_run_table"), analysis_cfg.get("kallisto_out"), analysis_cfg.get("gene_matrix_out"), analysis_cfg.get("chromosomes_needed"))
        #split_for_deseq = split_for_comparison(analysis_cfg.get("gene_matrix_out"), analysis_cfg.get("sra_run_table"), analysis_cfg.get("tmp_out"))
        deseq = run_deseq(**deseq_cfg)
        
    if 'data' in targets:
        with open('config/data-params.json') as fh:
            data_cfg = json.load(fh)
        
        getting_data = get_data(**data_cfg)
        
    if 'test' in targets:
        with open('config/test-params.json') as fh:
            test_cfg = json.load(fh)

        test_out = test(**test_cfg)
    
    
    return
        

    
    
if __name__ == '__main__':
    targets = sys.argv[1]
    main(targets)