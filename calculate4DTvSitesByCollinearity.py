#4Dtv (transversion rate on 4-fold degenerated sites)

import csv
import os
import shutil
import re
import argparse
import math
import multiprocessing as mp

codon_dict={
'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
'GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}###This is the DNA codons.

transversion_dict = {
"A":"TC",
"C":"AG",
"G":"TC",
"T":"AG",
}

def cds2pro(in_cds):
    raw_pep = in_cds.replace('.cds','.raw.pep')
    out_pep = in_cds.replace('.cds','.pep')
    os.system('cds2pro.pl {0} > {1}'.format(in_cds,raw_pep))
    with open (raw_pep) as in_fasta:
        with open (out_pep,'w') as out_fasta:
            for line in in_fasta:
                if '>' in line:
                    out_fasta.write(line)
                elif line == '\n':
                    continue                
                else:
                    out_fasta.write(line.rstrip('-.*'))
    os.remove(raw_pep)
    return out_pep

def extractSeqDict(in_fasta):
    seq_dict = {}
    with open (in_fasta) as in_seq:
        seq_name = ''
        seq_infor = in_seq.read().splitlines()
        for line in seq_infor:
            if '>' in line:
                infors = line.split()
                seq_name = infors[0].strip('>')
                seq_dict[seq_name] = ''
            else:
                seq_dict[seq_name] += line
    return seq_dict

def getAnchorSeqFiles(seq_dict,anchor_file,out_dir,seq_type):
    try:
        os.mkdir(out_dir)
    except FileExistsError:
        pass
    with open (anchor_file) as in_anchor:
        pair_num = 0
        for line in in_anchor:
            if '#' not in line:
                infors = line.split()
                with open ('{0}/pair_{1}.{2}'.format(out_dir,pair_num,seq_type),'w') as out_file:
                    out_file.write('>{0}\n{1}\n>{2}\n{3}\n'.format(infors[0],seq_dict[infors[0]],infors[1],seq_dict[infors[1]]))
                pair_num += 1


def alignmentByClustalo(pep_file):
    msa_file = pep_file.replace('.pep','.msa')
    os.system('clustalo -i {0} -o {1}'.format(pep_file,msa_file))

def pal2nal(cds_file):
    msa_file = cds_file.replace('.cds','.msa')
    out_file = cds_file.replace('.cds','.codon')
    os.system('pal2nal {0} {1} -output fasta -nogap -nomismatch > {2}'.format(msa_file,cds_file,out_file))

def calculate_4Dtv(codon_file,codon_d,transversion_d):
    names = ['','']
    seqs = ['','']
    freq = {'A':0,'T':0,'C':0,'G':0,'Y':0,'R':0}
    codon_4d = 0
    codon_4dt = 0
    v = 0.0
    a = 0.0
    b = 0.0
    d = 0.0
    with open(codon_file) as in_file:
        flag = -1
        for line in in_file:
            if '>' in line:
                flag+=1
                names[flag] = line.strip('>\n')
            else:
                seqs[flag] += line.strip()
    for i in range(0,len(seqs[0]),3):
        codon1 = seqs[0][i:i+3].upper()
        codon2 = seqs[1][i:i+3].upper()
        try:
            pep1 = codon_d[codon1]
            pep2 = codon_d[codon2]
        except KeyError:
            continue
        else:
            if pep1 == pep2:
                freq[codon1[-1]]+=1
                freq[codon2[-1]]+=1
                codon_4d+=1
                if codon2[-1] in transversion_d[codon1[-1]]:
                    codon_4dt+=1

    if codon_4d > 0:
        v = codon_4dt/codon_4d##this is raw 4dtv value
		##correction the raw 4dtv values by HKY substitution model
        freq['Y'] = freq['T'] + freq['C']
        freq['R'] = freq['A'] + freq['G']
        new_freq = {}
        for key,value in freq.items():
            new_freq[key] = 0.5*value/codon_4d
        Y=new_freq['Y']
        R=new_freq['R']
        A=new_freq['A']
        T=new_freq['T']
        C=new_freq['C']
        G=new_freq['G']
        if A*T*C*G != 0:
            try:
                a=-1*math.log(1-v*(T*C*R/Y+A*G*Y/R)/(2*(T*C*R+A*G*Y)))
                if (1-v/(2*Y*R)) > 0.0:
                    b=-1*math.log(1-v/(2*Y*R))
                    d=round(2*a*(T*C/Y+A*G/R)-2*b*(T*C*R/Y+A*G*Y/R-Y*R),4)
                    with open(codon_file.replace('.codon','.4dtv'),'w') as out_file:
                        out_file.write('{0};{1},{2},{3},{4},{5}\n'.format(names[0],names[1],d,round(v,4),codon_4d,codon_4dt))
            except ValueError:
                return

    
            
if __name__ == '__main__':
    #main_process
    parser = argparse.ArgumentParser(description=__doc__,prog='constructSpeciesTree.py')
    parser.add_argument('-a',action='store',dest='anchor_file',required=True,type=str,
        help='path to anchor file produced by MCscan')
    parser.add_argument('-pep',action='store',dest='pep_file',type=str,default='',
        help='path to protein sequence file with only one transcription; if not provided, we will translate CDS file')
    parser.add_argument('-cds',action='store',dest='cds_file',required=True,type=str,
        help='path to CDS sequence file with only one transcription')
    parser.add_argument('-t',action='store',dest='threads_num',type=int,default=1,
        help='threads number, default=1')
    parser.add_argument('-o',action='store',dest='out_dir',type=str,default='output',
        help='name of your output directory, default="output"')
    args = parser.parse_args()

    #step1 if protein sequences file is not provided, translate the given CDS file
    if args.pep_file == '':
        pep_file = cds2pro(args.cds_file)
    else:
        pep_file = args.pep_file
    
    #step2 extract seq files' infor
    cds_dict = extractSeqDict(args.cds_file)
    pep_dict = extractSeqDict(pep_file)

    #step3 extract anchored seqs to a directory
    getAnchorSeqFiles(cds_dict,args.anchor_file,args.out_dir,'cds')
    getAnchorSeqFiles(pep_dict,args.anchor_file,args.out_dir,'pep')

    #step4 MSA and arrange CDS to codon seq according to MSA result
    os.chdir(args.out_dir)
    pep_file_list = []
    cds_file_list = []
    codon_file_list = []
    for file_name in os.listdir():
        if '.pep' in file_name:
            pep_file_list.append(file_name)
        elif '.cds' in file_name:
            cds_file_list.append(file_name)

    real_threads = args.threads_num
    with mp.Pool(real_threads) as p:
        p.map(alignmentByClustalo,pep_file_list)
    
    with mp.Pool(real_threads) as p:
        p.map(pal2nal,cds_file_list)

    #step5 calculate 4Dtv value and write the results to 4dtv.csv
  
    for file_name in os.listdir():
        if '.codon' in file_name:
            codon_file_list.append(file_name)

    with mp.Pool(real_threads) as p:
        p.starmap(calculate_4Dtv,[(codon_file,codon_dict,transversion_dict) for codon_file in codon_file_list])
    
    with open ('../4dtv.csv','w') as out_file:
        out_file.write('gene_pair,corrected_4Dtv,raw_4Dtv,codon_4d,codon_4dt\n')
    for file_name in os.listdir():
        if '.4dtv' in file_name:
            os.system('cat {} >> ../4dtv.csv'.format(file_name))
    
    print('Finished! The calculation result has been stored in 4dtv.csv')
    exit()