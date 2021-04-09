__author__ = 'Yuan GAO'
__organization__ = 'Fujian Agriculture and Forestry University HIST'

import argparse
import os
import csv

def getAnchorByDict(species_name,ortho_tsv,split_symbol,out_dir):
    ortho_tsv = open(args.in_tsv)
    out_pep = open ('{0}/{1}.anchor.pep'.format(out_dir,species_name),'w')
    out_cds = open ('{0}/{1}.anchor.cds'.format(out_dir,species_name),'w')
    ortho_dict = csv.DictReader(ortho_tsv,delimiter='\t')
    seq_dict = {}
    with open (species_name + '.pep') as in_pep:
        flag = ''
        for line in in_pep:
            if '>' in line:
                infors = line.split()
                flag = infors[0].strip('>\n')
                seq_dict[flag] = {'pep':'','cds':''}
            else:
                seq_dict[flag]['pep'] += line.strip()
            
    with open (species_name + '.cds') as in_cds:
        flag = ''
        for line in in_cds:
            if '>' in line:
                infors = line.split()
                flag = infors[0].strip('>\n')
            else:
                seq_dict[flag]['cds'] += line.strip()
    for row in ortho_dict:
        if row[species_name] == '':
            continue
        else:
            infors = row[species_name].split(split_symbol)
            if len(infors) != 1:
                for i in range(len(infors)-1):
                    seq_name = infors[i].strip()
                    for j in range(i+1,len(infors)):
                        next_seq_name = infors[j].strip()
                        pep_seq_1 = seq_dict[seq_name]['pep']
                        pep_seq_2 = seq_dict[next_seq_name]['pep']
                        cds_seq_1 = seq_dict[seq_name]['cds']
                        cds_seq_2 = seq_dict[next_seq_name]['cds']
                        out_pep.write('>{0}\n{1}\n>{2}\n{3}\n'.format(seq_name,pep_seq_1,next_seq_name,pep_seq_2))
                        out_cds.write('>{0}\n{1}\n>{2}\n{3}\n'.format(seq_name,cds_seq_1,next_seq_name,cds_seq_2))
    out_pep.close()
    out_cds.close()
    ortho_tsv.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,prog='prepareKaKsByOrthogroup.py')
    parser.add_argument('-t',action='store',dest='in_tsv',required=True,type=str,
        help='your orthogroup tsv file')
    parser.add_argument('-o',action='store',dest='out_dir',type=str,default='kaks_seq',
        help='directory of output files, default=KAKS_seq')
    parser.add_argument('-s',action='store',dest='split_symbol',type=str,default=',',
        help='symbol to split gene names in one orthogroup, default = <,>')
    args = parser.parse_args()

    file_list = os.listdir()
    if args.out_dir not in file_list:
        os.mkdir(args.out_dir)
    #ortho_tsv = open(args.in_tsv)
    #ortho_dict = csv.DictReader(ortho_tsv,delimiter='\t')
    for file_name in file_list:
        print(file_name)
        if '.' in file_name:
            infors = file_name.split('.')
            if infors[-1] == 'pep':
                if (infors[0] + '.cds') in file_list:
                    ortho_tsv = open(args.in_tsv)
                    ortho_dict = csv.DictReader(ortho_tsv,delimiter='\t')
                    ortho_tsv.close()
                    getAnchorByDict(infors[0],args.in_tsv,args.split_symbol,args.out_dir)
                    print(file_name)
                else:
                    print('{}.cds is missing.'.format(infors[0]))
            