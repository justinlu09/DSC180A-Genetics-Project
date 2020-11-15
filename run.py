#! /usr/bin/env python
import sys
import os
import json
sys.path.insert(0, 'src')
#import etl
from etl import quality_check, clean_adapters, alignment



def main(targets):
    if 'process' in targets:
        with open('config/process-params.json') as fh:
            data_cfg = json.load(fh)
        
        fastq_data_b = quality_check(data_cfg.get('data_dir'), data_cfg.get('fastqc_path'), data_cfg.get('fq_output_bc'))
        
        
        #check if failed_checks is empty
        #if it is, then dont do cutadapt, if it isn't do cutadapt on those files
        cutadapt_data = clean_adapters(data_cfg.get('data_dir'), data_cfg.get('cutadapt_output'))
        
        
    if 'align' in targets:
        with open('config/align-params.json') as fh:
            align_cfg = json.load(fh)
        kallisto_data = alignment(**align_cfg)
        
    
#     if 'report' in targets:
#         with open('
        
    
    #if 'analysis' in targets:
    
    return
        

    
    
if __name__ == '__main__':
    targets = sys.argv[1]
    main(targets)