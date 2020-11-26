import pandas as pd
import os
import shutil
import logging

logging.basicConfig(filename="log.txt", filemode='a',
 format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
 datefmt='%H:%M:%S',
 level=logging.DEBUG)


def generate_gene_mat(gene_naming_table, sra_run_table, kallisto_out, gene_matrix_out):
    
    gene_matrix = pd.DataFrame()
    logging.info("Creating gene matrix...")
    if ('.ipynb_checkpoints' in os.listdir(kallisto_out)):
        shutil.rmtree(os.path.join(kallisto_out, '.ipynb_checkpoints'))
    
    sra_run = pd.read_csv(sra_run_table)
    dropped = sra_run.drop_duplicates(subset = ['BioSample'], keep = 'last')
    runs = dropped['Run'].values
    
    for i in sorted(os.listdir(kallisto_out)):
        SRR_path = os.path.join(kallisto_out, i)
        if i in runs:
            abundance_tsv = os.path.join(SRR_path, 'abundance.tsv')
            abundance = pd.read_csv(abundance_tsv, sep = '\t')['tpm']
            gene_matrix[i] = abundance
        else:
            continue

    gene_matrix['genes'] = pd.read_csv(os.path.join(kallisto_out, i, 'abundance.tsv'), sep = '\t')['target_id']
    
    gm = gene_matrix.copy()
    gm = gm[~gm['genes'].str.contains('NR')]
    
    gene_naming = pd.read_csv("/datasets/srp073813/reference/Gene_Naming.csv")
    chromosomes_needed = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
       'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
       'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']

    df_copy = gene_naming.copy()
    non_xy = df_copy[df_copy['chr'].isin(chromosomes_needed)]
    non_xy = non_xy.dropna(subset=['refseq']).reset_index(drop = True)
    ref = non_xy['refseq'].unique()
    
    non_xy_chromosomes = gm[gm['genes'].isin(ref)].set_index('genes')
    non_xy_chromosomes.index.name = '0'
    
    gene_matrix = non_xy_chromosomes
    
    gene_matrix.to_csv(gene_matrix_out)
    logging.info("Finished creating gene matrix.")
    
    #print(gene_matrix.shape)
    return 
