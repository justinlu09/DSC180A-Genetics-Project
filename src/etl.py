import os
from bs4 import BeautifulSoup


def quality_check(data_dir, fastqc_path, fq_output_bc):
    data = sorted(os.listdir(data_dir))
    count = 0
    lst = []
    for j in range(0, len(data), 2):
        if ((data[j][:12] + '_fastqc.html') in os.listdir(fq_output_bc)) and ((data[j+1] [:12] + '_fastqc.html') in os.listdir(fq_output_bc)):
            count += 1
            continue
        
        
        if count == 5:
            break
        
        #pair of files (..._1.fastq.gz, ..._2.fastq.gz)
        first_file = os.path.join(data_dir, data[j])
        second_file = os.path.join(data_dir, data[j+1])
        
        #running fastqc on the pair of files before cutadapt
        os.system(fastqc_path + ' ' + first_file + ' ' + second_file + ' --outdir ' + fq_output_bc)
        
        #count += 1        
        
    
    htmls = []
    for i in os.listdir(fq_output_bc):
        if i.endswith('html'):
            htmls.append(i)
    
    failed_checks = []
    for i in range(0, len(htmls), 2):
        fq_out_1 = htmls[i]
        fq_out_2 = htmls[i+1]
        html_path_1 = os.path.join(fq_output_bc, fq_out_1)
        html_path_2 = os.path.join(fq_output_bc, fq_out_2)
        soup1 = BeautifulSoup(open(html_path_1), 'html.parser')
        soup2 = BeautifulSoup(open(html_path_2), 'html.parser')
        check1 = soup1.find_all('div', attrs= {'class':"module"})[8].find('p').text
        check2 = soup2.find_all('div', attrs= {'class':"module"})[8].find('p').text
        if (check1 == 'No overrepresented sequences' or check2 == 'No overrepresented sequences'):
            continue
        else:
            failed_checks.append(html_path_1)
            failed_checks.append(html_path_2)
    
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
    count = 0
    for i in range(0, len(data), 2):
        first_file = os.path.join(data_dir, data[i])
        second_file = os.path.join(data_dir, data[i+1])
        os.system(kallisto_path + ' quant -i ' + kallisto_idx + ' -b 0 -o ' + os.path.join(kallisto_output, data[i][:10]) + ' ' + first_file + ' ' + second_file)
        
        count += 1 
        if count == 1:
            break
    
    return
    
