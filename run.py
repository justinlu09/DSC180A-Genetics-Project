#! /usr/bin/env python
import sys
import os
import json
sys.path.insert(0, 'src')
# import etl
from etl import get_data



def main(targets):
    if 'data' in targets:
        with open('config/data-params.json') as fh:
            data_cfg = json.load(fh)
        data = get_data(data_cfg.get('data_dir'), data_cfg.get('output'))
    return
        

    
    
if __name__ == '__main__':
    targets = sys.argv[1]
    main(targets)