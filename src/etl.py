import os

def quality_check(data_dir, fastqc_path, fq_output_bc):
    data = sorted(os.listdir(data_dir))
    count = 0
    for j in range(0, len(data), 2):
        if (data[j] in os.listdir(fq_output_bc)) and (data[j+1] in os.listdir(fq_output_bc)):
            continue
        
        #pair of files (..._1.fastq.gz, ..._2.fastq.gz)
        first_file = os.path.join(data_dir, data[j])
        second_file = os.path.join(data_dir, data[j+1])
        
        #running fastqc on the pair of files before cutadapt
        os.system(fastqc_path + ' ' + first_file + ' ' + second_file + ' --outdir ' + fq_output_bc)
    
        count += 1
        if count == 4:
            break
    return

def clean_adapters(data_dir, cutadapt_output):
    data = sorted(os.listdir(data_dir))
    for i in range(0, len(data), 2):
        
        first_file = os.path.join(data_dir, data[i])
        second_file = os.path.join(data_dir, data[i+1])
        
        os.system('cutadapt -a AACCGGTT -A AACCGGTT -o ' + cutadapt_output + '/out.1.fastq.gz -p ' + cutadapt_output + '/out.2.fastq.gz ' + first_file + ' ' + second_file + ' --cores=16')
        
    
    return
        
            
def alignment(data_dir, kallisto_path, kallisto_idx, kallisto_output):
    data = sorted(os.listdir(data_dir))
    for i in range(0, len(data), 2):
        first_file = os.path.join(data_dir, data[i])
        second_file = os.path.join(data_dir, data[i+1])
        
        os.system(kallisto_path + ' quant -i ' + kallisto_idx + ' -o ' + kallisto_output + ' -b 0 ' + first_file + ' ' + second_file)
        
        break
    
    return
    
        

# def process_data(data_dir, fastqc_path, fq_output_bc, fq_output_ac, cutadapt_output):
#     data = sorted(os.listdir(data_dir))
#     count = 0
    
#     for j in range(0, len(data), 2):
# #         process_file = os.path.join(data_dir, i)
        
#         #pair of files (..._1.fastq.gz, ..._2.fastq.gz)
#         first_file = os.path.join(data_dir, data[j])
#         second_file = os.path.join(data_dir, data[j+1])
        
#         #running fastqc on the pair of files before cutadapt
#         os.system(fastqc_path + ' ' + first_file + ' ' + second_file + ' --outdir ' + fq_output_bc)
        
#         #running cutadapt on pair of files
#         os.system('cutadapt -a AACCGGTT -A AACCGGTT -o ' + cutadapt_output + '/out.1.fastq.gz -p ' + cutadapt_output + '/out.2.fastq.gz ' + first_file + ' ' + second_file + ' --cores=8')
        
#         #running fastqc on the pair of files AFTER cutadapt to make sure cutadapt did its job
#         os.system(fastqc_path + ' ' + cutadapt_output + '/out.1.fastq.gz ' + cutadapt_output + '/out.2.fastq.gz' + ' --outdir ' + fq_output_ac)
        
#         #running Kallisto for alignment
        
        
#         count += 1
#         if count == 5:
#             break
    
#     return 
