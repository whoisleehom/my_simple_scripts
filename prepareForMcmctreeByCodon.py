__author__ = 'Yuan GAO'
__organization__ = 'Fujian Agriculture and Forestry University HIST'
__doc__ = (
    """
    作者：Yuan GAO  版本:v1.0.0  2020/1/5    本程序所实现的功能：
    将所有的单拷贝protein提取出来，按照orthofinder的结果，将每个orthogroup的单拷贝prot写到一个文件，
    用clustalo来进行比对，然后利用pal2nal将比对结果映射为密码子，用trimal去gap，提取出密码子的三个碱基构成新的序列。
    """
    )

import os
import argparse
def prep_clustalo(msa_dir,single_copy_ortho_list,ortho_csv,all_prot_file,all_cds_file):
#define the function to prepare the files for clustalo
#定义为clustalo进行准备的函数，返回species_list以便modify_for_mcmctree使用
    
    single_ortho_list = []
    prot_dict = {}
    cds_dict = {}
    species_list = []

    with open (single_copy_ortho_list) as in_list:
        for line in in_list:
            single_ortho_list.append(line.strip('\n'))

    with open (all_prot_file) as in_prot:
        gene_name = ''
        for line in in_prot:
            if '>' in line:
                line1 = line.split()
                gene_name = line1[0].strip('>')
                try:
                    prot_dict[gene_name]
                except KeyError:
                    prot_dict[gene_name] = ''
                else:
                    print('ERROR: duplicated {} in protein file'.format(gene_name))
            else:
                prot_dict[gene_name] += line
    
    with open (all_cds_file) as in_cds:
        gene_name = ''
        for line in in_cds:
            if '>' in line:
                line1 = line.split()
                gene_name = line1[0].strip('>')
                try:
                    cds_dict[gene_name]
                except KeyError:
                    cds_dict[gene_name] = ''
                else:
                    print('ERROR: duplicated {} in cds file'.format(gene_name))
            else:
                cds_dict[gene_name] += line

    with open (ortho_csv) as in_csv:
        lines = in_csv.read().splitlines()
        species_names = lines[0].split('\t')
        for line in species_names:#phylip格式有限制，物种名只能有十个字符，因而要对其进行修饰
            species_list.append(line[:10])
        print(species_list)
        for line in lines[1:]:
            infors = line.split('\t')
            if infors[0] in single_ortho_list:
                flag = 1
                for i in range(1,len(infors)):#检查protein和cds序列文件中的单拷贝序列是否完备

                    try:
                        prot_dict[infors[i]]
                    except KeyError:
                        print('{} is not in protein file.'.format(infors[i]))
                        flag = 0
                        break

                    try:
                        cds_dict[infors[i]]
                    except KeyError:
                        print('{} is not in cds file.'.format(infors[i]))
                        flag = 0
                        break

                if flag == 1:
                    with open ('{0}/{1}.pep'.format(msa_dir,infors[0]),'w') as out_prot:
                        for i in range(1,len(infors)):
                            seq = prot_dict[infors[i]].rstrip('\n')
                            species = species_list[i]
                            out_prot.write('>{0}\n{1}\n'.format(species,seq))

                    with open ('{0}/{1}.cds'.format(msa_dir,infors[0]),'w') as out_cds:
                        for i in range(1,len(infors)):
                            seq = cds_dict[infors[i]].rstrip('\n')
                            species = species_list[i]
                            out_cds.write('>{0}\n{1}\n'.format(species,seq))
    
    return (species_list[1:])
    

def msa_by_clustalo(msa_dir,threads_num):
#执行clustalo命令
    rootpath = os.getcwd()
    os.chdir(msa_dir)
    file_list = os.listdir('./')
    for item in file_list:
        infors = item.split('.')
        if infors[-1] == 'pep':
            clustalo_cmd = 'clustalo -i {0} -o {0}.msa --threads={1}'.format(item,threads_num)#clustalo commond
            print('clustalo:\t' + clustalo_cmd)
            check_clustalo = os.system(clustalo_cmd)
            if check_clustalo == 0:
                print('Alignment of {} finished'.format(item))
            else:
                print('ERROR: clustalo broke when align {}.'.format(item))
    os.chdir(rootpath)

def pal2nal_and_trimal(msa_dir):
    rootpath = os.getcwd()
    os.chdir(msa_dir)
    file_list = os.listdir('./')
    for item in file_list:
        infors = item.split('.')
        if infors[-1] == 'msa':
            cds_name = infors[0] + '.cds'
            pal2nal_cmd = 'pal2nal {0} {1} -output fasta > {2}.codon_msa'.format(item,cds_name,infors[0])
            print('pal2nal:\t' + pal2nal_cmd)
            check_pal2nal = os.system(pal2nal_cmd)
            if check_pal2nal == 0:
                print('Pal2nal for {} finished.'.format(infors[0]))
            else:
                print('ERROR:\tpal2nal for {} broke.'.format(infors[0]))

            trimal_cmd = 'trimal -in {0}.codon_msa -out {0}.codon_msa.trimmed -nogaps'.format(infors[0])
            print('trimal:\t' + trimal_cmd)
            check_trimal = os.system(trimal_cmd)
            if check_trimal == 0:
                print('Trimal for {}.codon_msa finished.'.format(infors[0]))
            else:
                print('ERROR:\ttrimal for {}.codon_msa broke.'.format(infors[0]))
    os.chdir(rootpath)

def extractThreeSitesOfCodon(species_list,msa_dir):
    rootpath = os.getcwd()
    os.chdir(msa_dir)
    file_list = os.listdir('./')
    msa_dict = {}
    codon_1_dict = {}
    codon_2_dict = {}
    codon_3_dict = {}
    for species_name in species_list:
        msa_dict[species_name] = ''
        codon_1_dict[species_name] = ''
        codon_2_dict[species_name] = ''
        codon_3_dict[species_name] = ''
    for item in file_list:
        infors = item.split('.')
        if infors[-1] == 'trimmed':
            with open (item) as in_file:
                
                lines = in_file.read().splitlines()
                species_name = ''
                for line in lines:
                    if '>' in line:
                        inf = line.split()
                        species_name = inf[0].strip('>')
                    else:
                        msa_dict[species_name] += line
                        
    for i in range(0,len(msa_dict[species_name]),3):
        for k,v in msa_dict.items():
            codon_1_dict[k] += v[i]
            codon_2_dict[k] += v[i+1]
            codon_3_dict[k] += v[i+2]
    with open ('merged.wholecodon', 'w') as outfile:
        for k,v in msa_dict.items():
            outfile.write('>{0}\n{1}\n'.format(k,v))
    os.system('trimal -in merged.wholecodon -out merged.wholecodon.phylip -phylip')
    with open ('merged.1.codon', 'w') as outfile1:
        for k,v in codon_1_dict.items():
            outfile1.write('>{0}\n{1}\n'.format(k,v))
    os.system('trimal -in merged.1.codon -out merged.1.codon.phylip -phylip')
    with open ('merged.2.codon', 'w') as outfile2:
        for k,v in codon_2_dict.items():
            outfile2.write('>{0}\n{1}\n'.format(k,v))
    os.system('trimal -in merged.2.codon -out merged.2.codon.phylip -phylip')
    with open ('merged.3.codon', 'w') as outfile3:
        for k,v in codon_3_dict.items():
            outfile3.write('>{0}\n{1}\n'.format(k,v))
    os.system('trimal -in merged.3.codon -out merged.3.codon.phylip -phylip')

    os.chdir(rootpath)

def modify_for_mcmctree(msa_result,species_list,outfile):
#将trimal后的结果（phylip格式）调整到适合mcmctree使用的格式
    with open (msa_result) as in_phylip:
        index_flag = 0
        block_flag = 0
        phylip_dict = {}
        for line in in_phylip:
            infors = line.strip('\n').split()
            num_flag = 1
            for infor in infors:
                try:
                    float(infor)
                except ValueError:
                    num_flag = 1 
                else:
                    num_flag = 0
            if num_flag == 0:#如果行为信息行，如 8  2800
                phylip_dict[str(index_flag)] = {'info_line':line.strip('\n'),'species':[],'seq':[]}
                index_flag += 1
                block_flag = -1
                
            elif line.strip() == '':#如果为空行，block_flag重置为-1
                block_flag = -1
            else:#如果不为空行，则有首行和其他行之分
                first_flag = 1
                block_flag += 1
                for species_name in species_list:
                    if species_name in line:
                        first_flag = 0
                        break
                if first_flag == 0:
                    seq = line.strip('\n').replace(species_name,'')
                    phylip_dict[str(index_flag-1)]['species'].append(species_name)
                    phylip_dict[str(index_flag-1)]['seq'].append(seq)
                        
                else:
                    seq = line.strip('\n')
                    phylip_dict[str(index_flag-1)]['seq'][block_flag] += seq
                            

    with open (outfile,'w') as out_phylip:
        for k,v in phylip_dict.items():
            out_phylip.write(v['info_line'] + '\n')
            for i in range(len(v['species'])):
                out_phylip.write('{0}  {1}\n'.format(v['species'][i],v['seq'][i]))
            out_phylip.write('\n')

if __name__ == '__main__':
    #main_process
    parser = argparse.ArgumentParser(description=__doc__,prog='constructSpeciesTree.py')
    parser.add_argument('-scl',action='store',dest='single_copy_ortho_list',required=True,type=str,
        help='full path to single copy orthogroup list')
    parser.add_argument('-og',action='store',dest='orthogroup_csv',required=True,type=str,
        help='full path to orthogroup.tsv generated by orthofinder')
    parser.add_argument('-prot',action='store',dest='merged_prot_file',required=True,type=str,
        help='full path to merged protein sequence file with only one transcription')
    parser.add_argument('-cds',action='store',dest='merged_cds_file',required=True,type=str,
        help='full path to merged cds sequence file with only one transcription')
    parser.add_argument('-t',action='store',dest='threads_num',type=str,default='1',
        help='threads number, default=1')
    
    args = parser.parse_args()
    root_path = os.getcwd()
    msa_dir = root_path + '/single_ortho_seq'
    #step1: make the directory
    print('*********************************\n\tstep_1:\tmake the directory\n*********************************')
    os.mkdir(msa_dir)
    #step2: check and modify the input files
    print('*********************************\n\tstep_2:\textract single copy sequences\n*********************************')
    species_list = prep_clustalo(msa_dir,args.single_copy_ortho_list,args.orthogroup_csv,args.merged_prot_file,args.merged_cds_file)    
    #step3: multiple sequeces alignment of peptide sequnces by clustalo
    print('*********************************\n\tstep_3:\tmultiple sequeces alignment of peptide sequnces by clustalo\n*********************************')
    msa_by_clustalo(msa_dir,args.threads_num)
    #step4: transfer the peptide MSA results to codon MSA results by pal2nal and trim them
    print('*********************************\n\tstep_4:\ttransfer the peptide MSA results to codon MSA results\n*********************************')
    pal2nal_and_trimal(msa_dir)
    #step5: extract 3 codon sites to construct new sequences
    print('*********************************\n\tstep_4:\textract 3 codon sites to\n*********************************')
    extractThreeSitesOfCodon(species_list,msa_dir)
    #step6: modify for mcmctree
    print('*********************************\n\tstep_4:\tmodify for mcmctree\n*********************************')
    modify_for_mcmctree(msa_dir + '/merged.1.codon.phylip',species_list,'merged.1.codon.mod.phylip')
    modify_for_mcmctree(msa_dir + '/merged.2.codon.phylip',species_list,'merged.2.codon.mod.phylip')
    modify_for_mcmctree(msa_dir + '/merged.3.codon.phylip',species_list,'merged.3.codon.mod.phylip')
    modify_for_mcmctree(msa_dir + '/merged.wholecodon.phylip',species_list,'merged.wholecodon.mod.phylip')
    print('prepareForMcmctreeByCodon.py finished!\n')