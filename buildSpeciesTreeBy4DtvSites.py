__author__ = 'Yuan GAO'
__organization__ = 'Fujian Agriculture and Forestry University HIST'
__doc__ = (
    """
    作者：Yuan GAO  版本:v2.0.0  2020/10/29    本程序所实现的功能：
    将所有的单拷贝protein提取出来，按照orthofinder的结果，将每个orthogroup的单拷贝prot写到一个文件，
    用clustalo来进行比对，然后利用pal2nal将比对结果映射为密码子，用trimal去gap，提取出4Dtv碱基构成新的序列继而构建物种树。
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
                    print('ERROR: duplicated gene name in protein file')
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
                    print('ERROR: duplicated gene name in cds file')
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

def extract4DtvSites(msa_dir):
    
    codon_dict ={
    'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
    'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
    'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
    'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
    'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
    'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}

    rootpath = os.getcwd()
    os.chdir(msa_dir)
    file_list = os.listdir('./')
    for item in file_list:
        infors = item.split('.')
        if infors[-1] == 'trimmed':
            with open (item) as in_file:
                msa_dict = {}
                FourDtv_dict = {}
                lines = in_file.read().splitlines()
                species_name = ''
                for line in lines:
                    if '>' in line:
                        species_name = line.strip('>')
                        msa_dict[species_name] = ''
                        FourDtv_dict[species_name] = ''
                    else:
                        msa_dict[species_name] += line
                empty_tag = 0
                for i in range(0,len(msa_dict[species_name]),3):
                    peptide = ''
                    same_tag = 0
                    for k,v in msa_dict.items():
                        codon = v[i:i+3]
                        try:
                            codon_dict[codon]
                        except KeyError:#
                            same_tag = 1
                            break
                        else:
                            if peptide == '':
                                peptide = codon_dict[codon]
                            elif peptide != codon_dict[codon]:#
                                same_tag = 1
                                break
                    if same_tag == 0:
                        for k,v in msa_dict.items():
                            FourDtv_dict[k] += v[i+2]
                            empty_tag = 1

                if empty_tag ==  1:
                    with open (infors[0] + '.4dtv','w') as outfile:
                        for k,v in FourDtv_dict.items():
                                outfile.write('>{0}\n{1}\n'.format(k,v))
                        print ('4Dtv sites of {0} were extracted to {0}.4dtv.'.format(infors[0]))
    os.chdir(rootpath)

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
        if infors[-1] == '4dtv':
            with open(msa_dir + '/' + item) as infile:
                tag = ''
                lines = infile.read().splitlines()
                for line in lines:
                    if '>' in line:
                        infors = line.strip('>').split()
                        tag = infors[0]
                    else:
                        merge_dict[tag] += line
    with open ('4dtv.merged','w') as outfile:
        for k,v in  merge_dict.items():
            outfile.write('>{0}\n{1}\n'.format(k,v))
    cmd = 'raxml -T {0} -f a -N {1} -m {2} -x 123456 -p 123456 -s 4dtv.merged -n concatenation_out.nwk'.format(thread,bootstrap,model)
    cmd_check = os.system(cmd)
    if cmd_check == 0:
        print('RAXML:\t{}'.format(cmd))
    else:
        print('ERROR:\traxml for concatenation_result/4dtv.merged broke.')
    os.chdir('../')

def construct_by_coalescent_method(species_list,msa_dir,thread,bootstrap,model,astral):
    #Construct the ML species tree with the coalescent method
    os.mkdir("coalescence_result")
    os.chdir("coalescence_result")
    file_list = os.listdir(msa_dir)
    for item in file_list:
        infors = item.split('.')
        if infors[-1] == '4dtv':
            os.system('raxml -T {0} -f a -N {1} -m {2} -x 123456 -p 123456 -s {3}/{4} -n {5}.nwk'.format(thread,bootstrap,model,msa_dir,item,infors[0]))
            os.system('cat RAxML_bipartitions.{}.nwk >>allSingleGenes_tree.nwk'.format(infors[0]))
            os.system('echo ./RAxML_bootstrap.{}.nwk >>allSingleGenes_bootstrap.txt'.format(infors[0]))    
    os.system('java -jar {0} -i allSingleGenes_tree.nwk -b allSingleGenes_bootstrap.txt -r {1} -o Astral.coalescent_out.result'.format(astral,bootstrap))
    os.system('tail -n 1 Astral.coalescent_out.result >Astral.coalescence_tree.nwk')
    os.chdir("../")

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
    parser.add_argument('-m',action='store',dest='model',type=str,default='PROTGAMMAJTT',
        help='model of amino acid substitution, default=PROTGAMMAJTT')
    parser.add_argument('-a',action='store',dest='astral',required=True,type=str,
        help='path to your ASTRAL Java script')
    parser.add_argument('-b',action='store',dest='bootstrap',type=str,default='100',
        help='bootstrap number, default=100')
    parser.add_argument('-nogaps',action='store_true',
        help='if you want to remove gaps in msa result, add this argument')
    args = parser.parse_args()
    root_path = os.getcwd()
    msa_dir = root_path + '/single_ortho_seq'
    #step1: make the directory
    print('*********************************\n\tstep_1:\tmake the directory\n*********************************')
    #os.mkdir(msa_dir)
    #step2: check and modify the input files
    print('*********************************\n\tstep_2:\textract single copy sequences\n*********************************')
    species_list = prep_clustalo(msa_dir,args.single_copy_ortho_list,args.orthogroup_csv,args.merged_prot_file,args.merged_cds_file)    
    #step3: multiple sequeces alignment of peptide sequnces by clustalo
    print('*********************************\n\tstep_3:\tmultiple sequeces alignment of peptide sequnces by clustalo\n*********************************')
    #msa_by_clustalo(msa_dir,args.threads_num)
    #step4: transfer the peptide MSA results to codon MSA results by pal2nal and trim them
    print('*********************************\n\tstep_4:\ttransfer the peptide MSA results to codon MSA results\n*********************************')
    #pal2nal_and_trimal(msa_dir)
    #step5: extract 4Dtv sites to construct new sequences
    print('*********************************\n\tstep_4:\textract 4Dtv sites to construct new sequences\n*********************************')
    #extract4DtvSites(msa_dir)
    #step6: construct phylogenetic tree by concatenation method
    print('*********************************\n\tstep_4:\tconstruct phylogenetic tree by concatenation method\n*********************************')
    #construct_by_concatenation_method(species_list,msa_dir,args.threads_num,args.bootstrap,args.model)
    modify_for_mcmctree('concatenation_result/4dtv.merged',species_list,'4dtv.mod.phylip')
    #step7: construct phylogenetic tree by coalescent method
    print('*********************************\n\tstep_4:\tconstruct phylogenetic tree by coalescent method\n*********************************')    
    #construct_by_coalescent_method(species_list,msa_dir,args.threads_num,args.bootstrap,args.model,args.astral)
    
    print('buildSpeciesTreeBy4Dtv.py finished!\n')