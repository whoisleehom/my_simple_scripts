__author__ = 'Yuan GAO'
__organization__ = 'Fujian Agriculture and Forestry University HIST'
__doc__ = (
    """
    作者：Yuan GAO  版本:v3.0.0  2020/12/23    本程序所实现的功能：
    在蛋白质或cds文件中仅保留一个最长的转录本,（可选）删除序列中多余信息
    用法：python3 removeRebundantProteins.py -i in_fasta [-o out_fasta] [-s symbol]
    """
    )

import sys,argparse


def removeRedundant(in_file,out_file,symbol,remove_bool):
    if out_file == '':
        out_file = '{}.simple'.format(in_file)
    
    gene_dic = {}
    flag = ''
    with open (in_file) as in_fasta:
        for line in in_fasta:
            if '>' in line:
                line1 = line.strip('>\n')
                line2 = line1.split()
                line3 = line2[0].split(symbol)
                flag = line2[0].replace(symbol + line3[-1],'',-1)
                try:
                    gene_dic[flag]
                except KeyError:
                    gene_dic[flag] = [[line.strip('>\n'),'']]
                else:
                    gene_dic[flag].append([line.strip('>\n'),''])
            else:
                try:
                    gene_dic[flag][-1][-1] += line.strip('\n*.')
                except KeyError:
                    print(flag)

    with open (out_file,'w') as out_fasta:
        for k,v in gene_dic.items():
            if remove_bool:
                if len(v) == 1:
                    out_fasta.write('>{0}\n{1}\n'.format(k,v[0][-1]))
                else:
                    trans_max = ''
                    for trans in v:
                        a = len(list(trans[-1]))
                        b = len(list(trans_max))
                        if a > b:
                            trans_max = trans[-1]
                    out_fasta.write('>{0}\n{1}\n'.format(k,trans_max))
            else:
                if len(v) == 1:
                    out_fasta.write('>{0}\n{1}\n'.format(v[0][0],v[0][-1]))
                else:
                    trans_max = []
                    for trans in v:
                        a = len(list(trans[1]))
                        b = len(list(trans_max))
                        if a > b:
                            trans_max = trans
                    out_fasta.write('>{0}\n{1}\n'.format(trans_max[0],trans_max[1]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,prog='removeRedundantProteins.py')
    parser.add_argument('-i',action='store',dest='in_file',required=True,type=str,
        help='full path to original fasta file')
    parser.add_argument('-o',action='store',dest='out_file',type=str,default='',
        help='full path to output file, default=#in_file_name#.simple')
    parser.add_argument('-s',action='store',dest='symbol',type=str,default='.',
        help='symbol to split transcription name, default = <.>')
    parser.add_argument('-r',action='store_true',
        help='add this argument if you want to remove the symbol and infors after from gene name')
    args = parser.parse_args()
    print('dealing with {}'.format(args.in_file))
    removeRedundant(args.in_file,args.out_file,args.symbol,args.r)