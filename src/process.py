import os
import zipfile
import logging

logging.basicConfig(filename="log.txt", filemode='a',
 format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
 datefmt='%H:%M:%S',
 level=logging.DEBUG)

def quality_check(data_dir, fastqc_path, fq_output_bc):
    data = sorted(os.listdir(data_dir))
    count = 0
    lst = []
    
    logging.info("Starting FastQC...")
    for j in range(0, len(data), 2):
        if ((data[j][:12] + '_fastqc.html') in os.listdir(fq_output_bc)) and ((data[j+1] [:12] + '_fastqc.html') in os.listdir(fq_output_bc)):
            #count += 1
            #break
            continue
        
        #pair of files (..._1.fastq.gz, ..._2.fastq.gz)
        first_file = os.path.join(data_dir, data[j])
        second_file = os.path.join(data_dir, data[j+1])
        
        logging.info("Running FastQC on " + str(data[j]) + " and " + str(data[j+1]))
        #running fastqc on the pair of files before cutadapt
        os.system(fastqc_path + ' ' + first_file + ' ' + second_file + ' --outdir ' + fq_output_bc)
        
#         count += 1
#         if count == 30:
#             break
    return
    

def check_fastqc(fq_output_bc):
    zipped = []
    for i in os.listdir(fq_output_bc):
        if i.endswith('zip'):
            zipped.append(i)
    
    zips = sorted(zipped)
    
    #going through zipped results of FastQC to look for pass/fail flag
    logging.info("Checking results of FastQC...")
    failed_checks = []
    for i in range(0, len(zips), 2):
        fq_out_1 = zips[i]
        fq_out_2 = zips[i+1]
        zip_path_1 = os.path.join(fq_output_bc, fq_out_1)
        zip_path_2 = os.path.join(fq_output_bc, fq_out_2)
        
        with zipfile.ZipFile(zip_path_1, 'r') as zip_ref:
            zip_ref.extractall(os.path.join(fq_output_bc, fq_out_1[:-4]))
            
        with zipfile.ZipFile(zip_path_2, 'r') as zip_ref:
            zip_ref.extractall(os.path.join(fq_output_bc, fq_out_2[:-4]))
        
        #print(zips[i])
        text1 = ''
        with open(os.path.join(fq_output_bc, fq_out_1[:-4], fq_out_1[:-4], 'fastqc_data.txt')) as fh:
            text1 = fh.readlines()[1]
        #print(text1)
        
        #print(zips[i+1])
        text2 = ''
        with open(os.path.join(fq_output_bc, fq_out_2[:-4], fq_out_2[:-4], 'fastqc_data.txt')) as fh:
            text2 = fh.readlines()[1]
        #print(text2)
        
        os.remove(zip_path_1)
        os.remove(zip_path_2)
        
        if ('pass' in text1 and 'pass' in text2):
            logging.info(str(fq_out_1[:-4]) + " and " + str(fq_out_2[:-4]) + " passed checks.")
            continue
        else:
            failed_checks.append(fq_out_1)
            failed_checks.append(fq_out_2)
    
    return failed_checks


def clean_adapters(failed, cutadapt_output):
    data = sorted(os.listdir(failed))
    if len(failed) == 0:
        return
    for i in range(0, len(data), 2):
        first_file = os.path.join(data_dir, data[i])
        second_file = os.path.join(data_dir, data[i+1])
        
        os.system('cutadapt -a AACCGGTT -A AACCGGTT -o ' + cutadapt_output + '/' + data[i] + ' -p ' + cutadapt_output + '/' + data[i+1] + ' ' + first_file + ' ' + second_file + ' --cores=16')
        
    return


def alignment(data_dir, kallisto_path, kallisto_idx, kallisto_output):
    
    data = sorted(os.listdir(data_dir))
    
    #count = 0
    logging.info("Starting alignment with Kallisto...")
    for i in range(0, len(data), 2):
        if 'SRR' in data[i]:
            #count += 1 
            first_file = os.path.join(data_dir, data[i])
            second_file = os.path.join(data_dir, data[i+1])
            logging.info("Running pair-end alignment with Kallisto on " + str(data[i]) + " and " + str(data[i+1]))
            os.system(kallisto_path + ' quant -i ' + kallisto_idx + ' -b 8 -o ' + os.path.join(kallisto_output, data[i][:10]) + ' ' + first_file + ' ' + second_file + ' -t 8') #+ ' >>' + os.path.join(kallisto_output, data[i][:10],'pseudo.bam'))
        else:
            break
    
    return
