__author__ = 'Yuan GAO'
__organization__ = 'Fujian Agriculture and Forestry University HIST'
__doc__ = (
    """
    作者：Yuan GAO  版本:v3.0.0  2021/4/3    本程序所实现的功能：
    对所有的序列进行polish以去除序列末尾的非编码符号，将所有的单拷贝序列提取出来，按照orthofinder的结果，将每个orthogroup的单拷贝序列写到一个文件，
    用clustalo来进行比对，然后根据-nogaps指令决定是否通过trimal去除gap。比对结果利用raxml和astral来进行串联和并联的物种树构建。
    """
    )

import os
import argparse
import multiprocessing as mp

def polish(in_file,out_file):
    if out_file == '':
        out_file = '{}.polished'.format(in_file)
    with open (in_file) as in_fasta:
        with open (out_file,'w') as out_fasta:
            for line in in_fasta:
                if '>' in line:
                    out_fasta.write(line)
                elif line == '\n':
                    continue                
                else:
                    out_fasta.write(line.strip('-.*'))

def prep_clustalo(single_copy_ortho_list,ortho_csv,all_prot_file,working_dir):
#define the function to prepare the files for clustalo
#定义为clustalo进行准备的函数，返回species_list以便modify_for_mcmctree使用

    single_ortho_list = []
    prot_dict = {}
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
                    print('ERROR: duplicated gene name')
            else:
                prot_dict[gene_name] += line

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
                for i in range(1,len(infors)):#检查prot文件中的单拷贝序列是否完备
                    try:
                        prot_dict[infors[i]]
                    except KeyError:
                        print('{} is not in prot file.'.format(infors[i]))
                        flag = 0
                        break
                    else:
                        flag = 1
                if flag == 1:
                    with open ('{0}/{1}.prot'.format(working_dir,infors[0]),'w') as out_fasta:
                        for i in range(1,len(infors)):
                            seq = prot_dict[infors[i]].rstrip('\n')
                            species = species_list[i]
                            out_fasta.write('>{0}\n{1}\n'.format(species,seq))
    return (species_list[1:])

def msa_by_clustalo(prot_file):
#执行clustalo命令
    os.system('clustalo -i {0} -o {0}.msa'.format(prot_file))

def triming_by_trimal(msa_file):
    os.system('trimal -in {0} -out {1}.trimmed -noallgaps -automated1'.format(msa_file,msa_file.replace('.msa','.trimmed')))#执行trimal命令
'''
def check_trimmed(msa_dir):
    file_list = os.listdir(msa_dir)
    for item in file_list:
        infors = item.split('.')
        if infors[-1] == 'trimmed':
            with open (msa_dir + '/'  + item) as infile:
                tag = 1
                in_fasta = infile.read().splitlines()
                len_dict = {}
                cur = ''
                for line in in_fasta:
                    if '>' in line:
                        name = line.strip('>')
                        cur = name
                        len_dict[name] = 0
                    else:
                        len_dict[cur] += len(line)
                eq_tag = 0
                eq = 0
                for k,v in len_dict.items():
                    if eq_tag == 0:
                        eq_tag += 1
                        eq = v 
                    elif eq != v:
                        tag = 0
                        break
                if tag == 1:
                    print('{} check passed.'.format(item))
                    os.system('cp {0}/{1} {0}/{1}.checked'.format(msa_dir,item))
                else:                    
                    print('{} check not passed.'.format(item))
                    print(len_dict)
'''     

def construct_by_concatenation_method(species_list,msa_dir,thread,bootstrap,model):
    #Construct the ML species tree with the concatenation method
    os.mkdir('concatenation_result')
    os.chdir('concatenation_result')
    file_list = os.listdir(msa_dir)
    merge_dict = {}
    for species in species_list:
        merge_dict[species] = ''
    for item in file_list:
        infors = item.split('.')
        if infors[-1] == 'trimmed':
            with open(msa_dir + '/' + item) as infile:
                tag = ''
                lines = infile.read().splitlines()
                for line in lines:
                    if '>' in line:
                        infors = line.strip('>').split()
                        tag = infors[0]
                    else:
                        merge_dict[tag] += line
    with open ('merged_trimmed.msa','w') as outfile:
        for k,v in  merge_dict.items():
            outfile.write('>{0}\n{1}\n'.format(k,v))

    os.system('raxml -T {0} -f a -N {1} -m {2} -x 123456 -p 123456 -s merged_trimmed.msa -n concatenation_out.nwk'.format(thread,bootstrap,model))
    os.chdir('../')

def construct_by_coalescent_method(species_list,msa_dir,thread,bootstrap,model,astral):
    #Construct the ML species tree with the coalescent method
    os.mkdir("coalescence_result")
    os.chdir("coalescence_result")
    file_list = os.listdir(msa_dir)
    for item in file_list:
        infors = item.split('.')
        if infors[-1] == 'checked':
            os.system('raxml -T {0} -f a -N {1} -m {2} -x 123456 -p 123456 -s {3}/{4} -n {5}.nwk'.format(thread,bootstrap,model,msa_dir,item,infors[0]))
            os.system('cat RAxML_bipartitions.{} >>allSingleGenes_tree.nwk'.format(infors[0]))
            os.system('echo ./RAxML_bootstrap.{} >>allSingleGenes_bootstrap.txt'.format(infors[0]))    
    os.system('java -jar {0} -i allSingleGenes_tree.nwk -b allSingleGenes_bootstrap.txt -r {1} -o Astral.coalescent_out.result'.format(astral,bootstrap))
    os.system('tail -n 1 Astral.coalescent_out.result >Astral.coalescence_tree.nwk')
    os.chdir("../")


if __name__ == '__main__':
    #main_process
    parser = argparse.ArgumentParser(description=__doc__,prog='constructSpeciesTree.py')
    parser.add_argument('-scl',action='store',dest='single_copy_ortho_list',required=True,type=str,
        help='full path to single copy orthogroup list')
    parser.add_argument('-og',action='store',dest='orthogroup_csv',required=True,type=str,
        help='full path to orthogroup.tsv generated by orthofinder')
    parser.add_argument('-prot',action='store',dest='merged_prot_file',required=True,type=str,
        help='full path to merged protein sequence file with only one transcription')
    parser.add_argument('-t',action='store',dest='threads_num',type=str,default='1',
        help='threads number, default=1')
    parser.add_argument('-m',action='store',dest='model',type=str,default='PROTGAMMAJTT',
        help='model of amino acid substitution, default=PROTGAMMAJTT')
    parser.add_argument('-a',action='store',dest='astral',required=True,type=str,
        help='path to your ASTRAL Java script')
    parser.add_argument('-b',action='store',dest='bootstrap',type=str,default='100',
        help='bootstrap number, default=100')
    args = parser.parse_args()
    root_path = os.getcwd()

    working_dir = 'single_ortho_prot'
    os.mkdir(working_dir)
    species_list = prep_clustalo(args.single_copy_ortho_list,args.orthogroup_csv,args.merged_prot_file,working_dir)

    
    os.chdir(working_dir)
    pep_list = []
    for seq in os.listdir():
        if '.pep' in seq:
            pep_list.append(seq)
    with mp.Pool(args.threads_num) as p:
        p.map(msa_by_clustalo,pep_list)
    
    msa_list = []
    for seq in os.listdir():
        if '.msa' in seq:
            msa_list.append(seq)
    with mp.Pool(args.threads_num) as p:
        p.map(triming_by_trimal,msa_list)
    os.chdir('..')

    #check_trimmed(root_path + '/single_ortho_prot')
    construct_by_concatenation_method(species_list,root_path + '/single_ortho_prot',args.threads_num,args.bootstrap,args.model)
    construct_by_coalescent_method(species_list,root_path + '/single_ortho_prot',args.threads_num,args.bootstrap,args.model,args.astral)
    print('buildSpeciesTree.py finished!\n')