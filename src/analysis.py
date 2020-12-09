import pandas as pd
import os
import shutil
import logging

logging.basicConfig(filename="log.txt", filemode='a',
 format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
 datefmt='%H:%M:%S',
 level=logging.DEBUG)


def generate_gene_mat(gene_naming_table, sra_run_table, kallisto_out, gene_matrix_out, chromosomes_needed):
    
    gene_matrix = pd.DataFrame()
    logging.info("Creating gene matrix...")
    if ('.ipynb_checkpoints' in os.listdir(kallisto_out)):
        shutil.rmtree(os.path.join(kallisto_out, '.ipynb_checkpoints'))
    
    sra_run = pd.read_csv(sra_run_table)
    dups = sra_run[sra_run.duplicated(subset=['BioSample'], keep = False)]
    g = dups.groupby('BioSample')['Run'].apply(list)
    g = pd.DataFrame(g)
    
    for i in sorted(os.listdir(kallisto_out)):
        SRR_path = os.path.join(kallisto_out, i)
        abundance_tsv = os.path.join(SRR_path, 'abundance.tsv')
        abundance = pd.read_csv(abundance_tsv, sep = '\t')['est_counts']
        gene_matrix[i] = abundance

    
    for j in g['Run'].values:
        if j[0] in os.listdir(kallisto_out):
            if len(j) == 2:
                gene_matrix[j[0]] = gene_matrix[j[0]] + gene_matrix[j[1]]
                gene_matrix.drop(j[1], axis = 1, inplace = True)
            else:
                gene_matrix[j[0]] = gene_matrix[j[0]] + gene_matrix[j[1]] + gene_matrix[j[2]]
                gene_matrix.drop(j[1], axis = 1, inplace = True)
                gene_matrix.drop(j[2], axis = 1, inplace = True)
        else:
            break

    gene_matrix['genes'] = pd.read_csv(os.path.join(kallisto_out, i, 'abundance.tsv'), sep = '\t')['target_id']
    
    gm = gene_matrix.copy()
    gm = gm[~gm['genes'].str.contains('NR')]
    
    df_copy = pd.read_csv(gene_naming_table).copy()
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

def split_for_comparison(gene_matrix_out, sra_run_table, tmp_out):
    sra_run = pd.read_csv(sra_run_table)
    
    #getting gene matrix
    gene_matrix = pd.read_csv(gene_matrix_out)
    sra_run = sra_run[sra_run['Run'].isin(gene_matrix.columns)]
    new = sra_run.drop_duplicates(subset = ['BioSample'], keep = 'first')
    
    # get individual groups
    ancg_control = new.groupby('source_name').get_group('AnCg_Control')['Run'].tolist()
    ancg_bpd = new.groupby('source_name').get_group('AnCg_Bipolar Disorder')['Run'].tolist()
    ancg_sz = new.groupby('source_name').get_group('AnCg_Schizophrenia')['Run'].tolist()
    ancg_mdd = new.groupby('source_name').get_group('AnCg_Major Depression')['Run'].tolist()

    nacc_control = new.groupby('source_name').get_group('nAcc_Control')['Run'].tolist()
    nacc_bpd = new.groupby('source_name').get_group('nAcc_Bipolar Disorder')['Run'].tolist()
    nacc_sz = new.groupby('source_name').get_group('nAcc_Schizophrenia')['Run'].tolist()
    nacc_mdd = new.groupby('source_name').get_group('nAcc_Major Depression')['Run'].tolist()

    dlpfc_control = new.groupby('source_name').get_group('DLPFC_Control')['Run'].tolist()
    dlpfc_bpd = new.groupby('source_name').get_group('DLPFC_Bipolar Disorder')['Run'].tolist()
    dlpfc_sz = new.groupby('source_name').get_group('DLPFC_Schizophrenia')['Run'].tolist()
    dlpfc_mdd = new.groupby('source_name').get_group('DLPFC_Major Depression')['Run'].tolist()
    
    
    
    # AnCg Control vs BPD
    ancg_control_df = gene_matrix[ancg_control]
    ancg_bpd_df = gene_matrix[ancg_bpd]
    ancg_control_bpd = pd.concat([ancg_bpd_df, ancg_control_df], axis=1)
    ancg_control_bpd = ancg_control_bpd.fillna(0)
    ancg_control_bpd.to_csv(tmp_out + '/ancg_control_bpd.csv')    
    
    ancg_control_only = new[(new['source_name'] == 'AnCg_Control')][['Run', 'source_name',  'age_at_death', 'post-mortem_interval', 'Brain_pH']]
    ancg_bpd_only = new[(new['source_name'] == 'AnCg_Bipolar Disorder')][['Run', 'source_name',  'age_at_death', 'post-mortem_interval', 'Brain_pH']]
    ancg_control_bpd_coldata = pd.concat([ancg_bpd_only, ancg_control_only])
    ancg_control_bpd_coldata = ancg_control_bpd_coldata.set_index('Run')
    ancg_control_bpd_coldata.columns = ['disorder', 'age', 'PMI', 'pH']
    ancg_control_bpd_coldata.index.name = None
    ancg_control_bpd_coldata['disorder'] = ancg_control_bpd_coldata['disorder'].str.replace('_', '')
    ancg_control_bpd_coldata['disorder'] = ancg_control_bpd_coldata['disorder'].str.replace(' ', '')
    ancg_control_bpd_coldata = ancg_control_bpd_coldata[['age','PMI','pH','disorder']]
    ancg_control_bpd_coldata.to_csv(tmp_out + '/ancg_control_bpd_coldata.csv')
    
    
    # AnCG Control vs SZ
    ancg_control_df = gene_matrix[ancg_control]
    ancg_sz_df = gene_matrix[ancg_sz]
    ancg_control_sz = pd.concat([ancg_sz_df, ancg_control_df], axis=1)
    ancg_control_sz = ancg_control_sz.fillna(0)
    ancg_control_sz.to_csv(tmp_out + '/ancg_control_sz.csv')
    
    ancg_control_only = new[(new['source_name'] == 'AnCg_Control')][['Run', 'source_name',  'age_at_death', 'post-mortem_interval', 'Brain_pH']]
    ancg_sz_only = new[(new['source_name'] == 'AnCg_Schizophrenia')][['Run', 'source_name',  'age_at_death', 'post-mortem_interval', 'Brain_pH']]
    ancg_control_sz_coldata = pd.concat([ancg_sz_only, ancg_control_only])
    ancg_control_sz_coldata = ancg_control_sz_coldata.set_index('Run')
    ancg_control_sz_coldata.columns = ['disorder', 'age', 'PMI', 'pH']
    ancg_control_sz_coldata.index.name = None
    ancg_control_sz_coldata['disorder'] = ancg_control_sz_coldata['disorder'].str.replace('_', '')
    ancg_control_sz_coldata['disorder'] = ancg_control_sz_coldata['disorder'].str.replace(' ', '')
    ancg_control_sz_coldata = ancg_control_sz_coldata[['age','PMI','pH','disorder']]
    ancg_control_sz_coldata['pH'] = ancg_control_sz_coldata['pH'].fillna(ancg_control_sz_coldata['pH'].mean())
    ancg_control_sz_coldata.to_csv(tmp_out + '/ancg_control_sz_coldata.csv')
    
    
    #AnCg Control vs MDD
    ancg_control_df = gene_matrix[ancg_control]
    ancg_mdd_df = gene_matrix[ancg_mdd]
    ancg_control_mdd = pd.concat([ancg_mdd_df, ancg_control_df], axis=1)
    ancg_control_mdd = ancg_control_mdd.fillna(0)
    ancg_control_mdd.to_csv(tmp_out + '/ancg_control_mdd.csv')
    
    ancg_control_only = new[(new['source_name'] == 'AnCg_Control')][['Run', 'source_name',  'age_at_death', 'post-mortem_interval', 'Brain_pH']]
    ancg_mdd_only = new[(new['source_name'] == 'AnCg_Major Depression')][['Run', 'source_name',  'age_at_death', 'post-mortem_interval', 'Brain_pH']]
    ancg_control_mdd_coldata = pd.concat([ancg_mdd_only, ancg_control_only])
    ancg_control_mdd_coldata = ancg_control_mdd_coldata.set_index('Run')
    ancg_control_mdd_coldata.columns = ['disorder', 'age', 'PMI', 'pH']
    ancg_control_mdd_coldata.index.name = None
    ancg_control_mdd_coldata['disorder'] = ancg_control_mdd_coldata['disorder'].str.replace('_', '')
    ancg_control_mdd_coldata['disorder'] = ancg_control_mdd_coldata['disorder'].str.replace(' ', '')
    ancg_control_mdd_coldata = ancg_control_mdd_coldata[['age','PMI','pH','disorder']]
    ancg_control_mdd_coldata['pH'] = ancg_control_mdd_coldata['pH'].fillna(ancg_control_mdd_coldata['pH'].mean())
    ancg_control_mdd_coldata.to_csv(tmp_out + '/ancg_control_mdd_coldata.csv')
    
    
    # nAcc Control vs SZ
    nacc_control_df = gene_matrix[nacc_control]
    nacc_sz_df = gene_matrix[nacc_sz]
    nacc_control_sz = pd.concat([nacc_sz_df, nacc_control_df], axis=1)
    nacc_control_sz = nacc_control_sz.fillna(0)
    nacc_control_sz.to_csv(tmp_out + '/nacc_control_sz.csv')
    
    nacc_control_only = new[(new['source_name'] == 'nAcc_Control')][['Run', 'source_name',  'age_at_death', 'post-mortem_interval', 'Brain_pH']]
    nacc_sz_only = new[(new['source_name'] == 'nAcc_Schizophrenia')][['Run', 'source_name',  'age_at_death', 'post-mortem_interval', 'Brain_pH']]
    nacc_control_sz_coldata = pd.concat([nacc_sz_only, nacc_control_only])
    nacc_control_sz_coldata = nacc_control_sz_coldata.set_index('Run')
    nacc_control_sz_coldata.columns = ['disorder', 'age', 'PMI', 'pH']
    nacc_control_sz_coldata.index.name = None
    nacc_control_sz_coldata['disorder'] = nacc_control_sz_coldata['disorder'].str.replace('_', '')
    nacc_control_sz_coldata['disorder'] = nacc_control_sz_coldata['disorder'].str.replace(' ', '')
    nacc_control_sz_coldata = nacc_control_sz_coldata[['age','PMI','pH','disorder']]
    nacc_control_sz_coldata['pH'] = nacc_control_sz_coldata['pH'].fillna(nacc_control_sz_coldata['pH'].mean())
    nacc_control_sz_coldata.to_csv(tmp_out + '/nacc_control_sz_coldata.csv')
    
    
    # nAcc Control vs BPD
    nacc_control_df = gene_matrix[nacc_control]
    nacc_bpd_df = gene_matrix[nacc_bpd]
    nacc_control_bpd = pd.concat([nacc_bpd_df, nacc_control_df], axis=1)
    nacc_control_bpd = nacc_control_bpd.fillna(0)
    nacc_control_bpd.to_csv(tmp_out + '/nacc_control_bpd.csv')
    
    nacc_control_only = new[(new['source_name'] == 'nAcc_Control')][['Run', 'source_name',  'age_at_death', 'post-mortem_interval', 'Brain_pH']]
    nacc_bpd_only = new[(new['source_name'] == 'nAcc_Bipolar Disorder')][['Run', 'source_name',  'age_at_death', 'post-mortem_interval', 'Brain_pH']]
    nacc_control_bpd_coldata = pd.concat([nacc_bpd_only, nacc_control_only])
    nacc_control_bpd_coldata = nacc_control_bpd_coldata.set_index('Run')
    nacc_control_bpd_coldata.columns = ['disorder', 'age', 'PMI', 'pH']
    nacc_control_bpd_coldata.index.name = None
    nacc_control_bpd_coldata['disorder'] = nacc_control_bpd_coldata['disorder'].str.replace('_', '')
    nacc_control_bpd_coldata['disorder'] = nacc_control_bpd_coldata['disorder'].str.replace(' ', '')
    nacc_control_bpd_coldata = nacc_control_bpd_coldata[['age','PMI','pH','disorder']]
    nacc_control_bpd_coldata['pH'] = nacc_control_bpd_coldata['pH'].fillna(nacc_control_bpd_coldata['pH'].mean())
    nacc_control_bpd_coldata.to_csv(tmp_out + '/nacc_control_bpd_coldata.csv')
    
    
    # nAcc Control vs MDD
    nacc_control_df = gene_matrix[nacc_control]
    nacc_mdd_df = gene_matrix[nacc_mdd]
    nacc_control_mdd = pd.concat([nacc_mdd_df, nacc_control_df], axis=1)
    nacc_control_mdd = nacc_control_mdd.fillna(0)
    nacc_control_mdd.to_csv(tmp_out + '/nacc_control_mdd.csv')
    
    nacc_control_only = new[(new['source_name'] == 'nAcc_Control')][['Run', 'source_name',  'age_at_death', 'post-mortem_interval', 'Brain_pH']]
    nacc_mdd_only = new[(new['source_name'] == 'nAcc_Major Depression')][['Run', 'source_name',  'age_at_death', 'post-mortem_interval', 'Brain_pH']]
    nacc_control_mdd_coldata = pd.concat([nacc_mdd_only, nacc_control_only])
    nacc_control_mdd_coldata = nacc_control_mdd_coldata.set_index('Run')
    nacc_control_mdd_coldata.columns = ['disorder', 'age', 'PMI', 'pH']
    nacc_control_mdd_coldata.index.name = None
    nacc_control_mdd_coldata['disorder'] = nacc_control_mdd_coldata['disorder'].str.replace('_', '')
    nacc_control_mdd_coldata['disorder'] = nacc_control_mdd_coldata['disorder'].str.replace(' ', '')
    nacc_control_mdd_coldata = nacc_control_mdd_coldata[['age','PMI','pH','disorder']]
    nacc_control_mdd_coldata['pH'] = nacc_control_mdd_coldata['pH'].fillna(nacc_control_mdd_coldata['pH'].mean())
    nacc_control_mdd_coldata.to_csv(tmp_out + '/nacc_control_mdd_coldata.csv')
    
    
    # DLPFC Control vs SZ
    dlpfc_control_df = gene_matrix[dlpfc_control]
    dlpfc_sz_df = gene_matrix[dlpfc_sz]
    dlpfc_control_sz = pd.concat([dlpfc_sz_df, dlpfc_control_df], axis=1)
    dlpfc_control_sz = dlpfc_control_sz.fillna(0)
    dlpfc_control_sz.to_csv(tmp_out + '/dlpfc_control_sz.csv')
    
    dlpfc_control_only = new[(new['source_name'] == 'DLPFC_Control')][['Run', 'source_name',  'age_at_death', 'post-mortem_interval', 'Brain_pH']]
    dlpfc_sz_only = new[(new['source_name'] == 'DLPFC_Schizophrenia')][['Run', 'source_name',  'age_at_death', 'post-mortem_interval', 'Brain_pH']]
    dlpfc_control_sz_coldata = pd.concat([dlpfc_sz_only, dlpfc_control_only])
    dlpfc_control_sz_coldata = dlpfc_control_sz_coldata.set_index('Run')
    dlpfc_control_sz_coldata.columns = ['disorder', 'age', 'PMI', 'pH']
    dlpfc_control_sz_coldata.index.name = None
    dlpfc_control_sz_coldata['disorder'] = dlpfc_control_sz_coldata['disorder'].str.replace('_', '')
    dlpfc_control_sz_coldata['disorder'] = dlpfc_control_sz_coldata['disorder'].str.replace(' ', '')
    dlpfc_control_sz_coldata = dlpfc_control_sz_coldata[['age','PMI','pH','disorder']]
    dlpfc_control_sz_coldata['pH'] = dlpfc_control_sz_coldata['pH'].fillna(dlpfc_control_sz_coldata['pH'].mean())
    dlpfc_control_sz_coldata.to_csv(tmp_out + '/dlpfc_control_sz_coldata.csv')
    
    
    # DLPFC Control vs BPD
    dlpfc_control_df = gene_matrix[dlpfc_control]
    dlpfc_bpd_df = gene_matrix[dlpfc_bpd]
    dlpfc_control_bpd = pd.concat([dlpfc_bpd_df, dlpfc_control_df], axis=1)
    dlpfc_control_bpd = dlpfc_control_bpd.fillna(0)
    dlpfc_control_bpd.to_csv(tmp_out + '/dlpfc_control_bpd.csv')
    
    dlpfc_control_only = new[(new['source_name'] == 'DLPFC_Control')][['Run', 'source_name',  'age_at_death', 'post-mortem_interval', 'Brain_pH']]
    dlpfc_bpd_only = new[(new['source_name'] == 'DLPFC_Bipolar Disorder')][['Run', 'source_name',  'age_at_death', 'post-mortem_interval', 'Brain_pH']]
    dlpfc_control_bpd_coldata = pd.concat([dlpfc_bpd_only, dlpfc_control_only])
    dlpfc_control_bpd_coldata = dlpfc_control_bpd_coldata.set_index('Run')
    dlpfc_control_bpd_coldata.columns = ['disorder', 'age', 'PMI', 'pH']
    dlpfc_control_bpd_coldata.index.name = None
    dlpfc_control_bpd_coldata['disorder'] = dlpfc_control_bpd_coldata['disorder'].str.replace('_', '')
    dlpfc_control_bpd_coldata['disorder'] = dlpfc_control_bpd_coldata['disorder'].str.replace(' ', '')
    dlpfc_control_bpd_coldata = dlpfc_control_bpd_coldata[['age','PMI','pH','disorder']]
    dlpfc_control_bpd_coldata['pH'] = dlpfc_control_bpd_coldata['pH'].fillna(dlpfc_control_bpd_coldata['pH'].mean())
    dlpfc_control_bpd_coldata.to_csv(tmp_out + '/dlpfc_control_bpd_coldata.csv')
    
    
    # DLPFC Control vs MDD
    dlpfc_control_df = gene_matrix[dlpfc_control]
    dlpfc_mdd_df = gene_matrix[dlpfc_mdd]
    dlpfc_control_mdd = pd.concat([dlpfc_mdd_df, dlpfc_control_df], axis=1)
    dlpfc_control_mdd = dlpfc_control_mdd.fillna(0)
    dlpfc_control_mdd.to_csv(tmp_out + '/dlpfc_control_mdd.csv')
    
    dlpfc_control_only = new[(new['source_name'] == 'DLPFC_Control')][['Run', 'source_name',  'age_at_death', 'post-mortem_interval', 'Brain_pH']]
    dlpfc_mdd_only = new[(new['source_name'] == 'DLPFC_Major Depression')][['Run', 'source_name',  'age_at_death', 'post-mortem_interval', 'Brain_pH']]
    dlpfc_control_mdd_coldata = pd.concat([dlpfc_mdd_only, dlpfc_control_only])
    dlpfc_control_mdd_coldata = dlpfc_control_mdd_coldata.set_index('Run')
    dlpfc_control_mdd_coldata.columns = ['disorder', 'age', 'PMI', 'pH']
    dlpfc_control_mdd_coldata.index.name = None
    dlpfc_control_mdd_coldata['disorder'] = dlpfc_control_mdd_coldata['disorder'].str.replace('_', '')
    dlpfc_control_mdd_coldata['disorder'] = dlpfc_control_mdd_coldata['disorder'].str.replace(' ', '')
    dlpfc_control_mdd_coldata = dlpfc_control_mdd_coldata[['age','PMI','pH','disorder']]
    dlpfc_control_mdd_coldata['pH'] = dlpfc_control_mdd_coldata['pH'].fillna(dlpfc_control_mdd_coldata['pH'].mean())
    dlpfc_control_mdd_coldata.to_csv(tmp_out + '/dlpfc_control_mdd_coldata.csv')
    
    
    #call deseq.R script here
    
    return
    
    
    
    
