import os
import pandas as pd


def get_data(data_dir, output):
    data = os.listdir(data_dir)
    unprocessed = []
    for i in data:
        process_file = os.path.join(data_dir, i)
        unprocessed.append(process_file)
        # conversion of raw RNA-Seq .1 file into fastq using SRA Toolkit
        # preprocessing using cutadapt or fastqc (depending on what we choose) for quality control
    
    df = pd.DataFrame(unprocessed, columns = ['.1 file names'])
    out_csv = df.to_csv(output, index = False)
    return 