import os
import logging
import shutil

logging.basicConfig(filename="log.txt", filemode='a',
 format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
 datefmt='%H:%M:%S',
 level=logging.DEBUG)


def get_data(indir_fastqc, indir_kallisto, indir_gene_matrix, outdir_fastqc, outdir_kallisto, outdir_gene_matrix):
    directory = "data"
    parent_dir = "./"
    path = os.path.join(parent_dir, directory)
    
    #remove data dir if one already exists
    if (os.path.exists(path) and os.path.isdir(path)):
        shutil.rmtree(path)
    
    os.mkdir(path)
    
    directory1 = "raw"
    parent_dir = "./data/"
    
    os.mkdir(os.path.join(parent_dir, directory1))
    
    #create symlinks
    os.symlink(indir_fastqc, outdir_fastqc)
    os.symlink(indir_kallisto, outdir_kallisto)
    os.symlink(indir_gene_matrix, outdir_gene_matrix)
    
    return #pd.read_csv(outdir_gene_matrix)