#! /usr/bin/env python
import sys
import os
import json
sys.path.insert(0, 'src')
#import etl
from etl import quality_check, clean_adapters



def main(targets):
    if 'data' in targets:
        with open('config/data-params.json') as fh:
            data_cfg = json.load(fh)
        #fastq_data_b = quality_check(data_cfg.get('data_dir'), data_cfg.get('fastqc_path'), data_cfg.get('fq_output_bc'))
        cutadapt_data = clean_adapters(data_cfg.get('data_dir'), data_cfg.get('cutadapt_output'))
        #fastq_data_a = quality_check(data_cfg.get('cutadapt_output'), data_cfg.get('fastqc_path'), data_cfg.get('fq_output_ac'))
        
        
    #if 'analysis' in targets:
        ##
    
    return
        

    
    
if __name__ == '__main__':
    targets = sys.argv[1]
    main(targets)