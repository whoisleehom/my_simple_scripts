__author__ = 'Yuan GAO'
__organization__ = 'Fujian Agriculture and Forestry University HIST'
__doc__ = (
    """
    作者：Yuan GAO  版本:v1.0.0  2020/10/25    本程序所实现的功能：
    根据gff来重命名基因名
    """
    )

import os
import argparse
import re

def getNewChain(in_gff_file,chr_sign,key_sign,header,out_chain_file):
    unanchoredGene = []
    rename_chain = {}
    unarrangedGene = {}
    with open (in_gff_file) as in_gff:
        lines = in_gff.read().splitlines()
        for line in lines:
            infors = line.split()
            try:
                if infors[2] == 'mRNA':     
                    gene_infors = infors[-1].split(';')
                    temp_dic = {}
                    for ge in gene_infors:
                        ge_inf = ge.split('=')
                        temp_dic[ge_inf[0]] = ge_inf[1]

                    if chr_sign in infors[0]:  
                        try:
                            unarrangedGene[infors[0]]
                        except KeyError:
                            unarrangedGene[infors[0]] = {}
                            unarrangedGene[infors[0]][temp_dic[key_sign]] = int(infors[3])
                        else:
                            unarrangedGene[infors[0]][temp_dic[key_sign]] = int(infors[3])
                    else:
                        unanchoredGene.append(temp_dic[key_sign])
            except IndexError:
                continue
    print(unarrangedGene)
    with open (out_chain_file,'w') as out_file:
        for k,v in unarrangedGene.items():
            chr_num = k.strip(chr_sign)
            sortedGeneList = sorted(v.items(),key=lambda x:x[1])#字典被按照顺序转换为tuple
            start_num = 1
            for gene in sortedGeneList:
                #rename the gene in the format like Mdom01g00016
                rename_chain[gene[0]] = '{0}{1}g{2}'.format(header,chr_num.zfill(2),str(start_num).zfill(5))
                start_num += 1
        unanchored_num = 1
        for gene in unanchoredGene:
            rename_chain[gene] = '{0}00g{1}'.format(header,str(unanchored_num).zfill(5))
            unanchored_num += 1
        print(rename_chain)
        for k,v in rename_chain.items():
            out_file.write(k + '>' + v + '\n')

    return (rename_chain)


def replaceGeneNames_gff(in_gff_file,out_gff_file,rename_chain_dict):
    with open (in_gff_file) as in_gff:
        with open (out_gff_file,'w') as out_gff:
            for line in in_gff:
                for k,v in rename_chain_dict.items():
                    if k in line:
                        new_line = line.replace(k,v)
                        out_gff.write(new_line)
                        break

def replaceGeneNames_fasta(in_fasta_file,out_fasta_file,rename_chain_dict):
    with open (in_fasta_file) as in_fasta:
        with open (out_fasta_file,'w') as out_fasta:
            for line in in_fasta:
                if '>' in line:
                    for k,v in rename_chain_dict.items():
                        if k in line:
                            new_line = line.replace(k,v)
                            out_fasta.write(new_line)
                            break
                else:
                    out_fasta.write(line)
                        

if __name__ == '__main__':
    #main_process
    parser = argparse.ArgumentParser(description=__doc__,prog='constructSpeciesTree.py')
    parser.add_argument('-in_gff',action='store',dest='in_gff',type=str,
        help='path to input gff file')
    parser.add_argument('-chr_sign',action='store',dest='chr_sign',required=True,type=str,
        help='chromsome sign in gff file, like "chr" or "CHR"')
    parser.add_argument('-key',action='store',dest='key',required=True,type=str,
        help='key of your gene infor in gff file, like "ID" or "gene"' )
    parser.add_argument('-header',action='store',dest='header',required=True,type=str,
        help='header of your new gene name, usually the abbreviation of species Latin name' )
    parser.add_argument('-in_prot',action='store',dest='in_prot',type=str,
        help='path to input protein sequences file')
    parser.add_argument('-in_cds',action='store',dest='in_cds',type=str,
        help='path to input cds file')
    parser.add_argument('-out_chain',action='store',dest='out_chain',default='rename.chain',type=str,
        help='path to output chain file')
    args = parser.parse_args()
    rename_chain_dict = getNewChain(args.in_gff,args.chr_sign,args.key,args.header,args.out_chain)
    if args.in_gff:
        replaceGeneNames_gff(args.in_gff,args.in_gff + '.renamed',rename_chain_dict)
    if args.in_prot:
        replaceGeneNames_fasta(args.in_prot,args.in_prot + '.renamed',rename_chain_dict)
    if args.in_cds:
        replaceGeneNames_fasta(args.in_cds,args.in_cds + '.renamed',rename_chain_dict)