#! python3
__author__ = 'Yuan GAO'
__organization__ = 'Fujian Agriculture and Forestry University HIST'
__doc__ = (
    """
    作者：Yuan GAO  版本:v1.0.0  2020/8/28    本程序所实现的功能：
    根据mcscan的共线性分析输出文件与cds和pep文件得到kaks分析的输入文件
    """
    )
import os
import argparse

def getFastaDic(in_file_path):
    out_dic = {}
    with open (in_file_path) as in_file:
        flag = ''
        for line in in_file:
            if '>' in line:
                infors = line.split()
                flag = infors[0].strip('>')
                out_dic[flag] = ''
            else:
                out_dic[flag] += line
    return out_dic

def getAnchorList(anchor_file):
    list_1 = []
    list_2 = []

    with open (anchor_file) as in_anchor:
        for line in in_anchor:
            if '#' not in line:
                infors = line.split()
                list_1.append(infors[0])
                list_2.append(infors[1])
    return [list_1,list_2]

def outputMergedFiles(cds_dict_1,cds_dict_2,pep_dict_1,pep_dict_2,list_1,list_2,out_cds_file,out_pep_file):
    out_cds = open(out_cds_file,'w')
    out_pep = open(out_pep_file,'w')
    for i in range(len(list_1)):
        out_cds.write('>' + list_1[i] + '\n' + cds_dict_1[list_1[i]])
        out_cds.write('>' + list_2[i] + '\n' + cds_dict_2[list_2[i]])
        out_pep.write('>' + list_1[i] + '\n' + pep_dict_1[list_1[i]])
        out_pep.write('>' + list_2[i] + '\n' + pep_dict_2[list_2[i]])
    out_cds.close()
    out_pep.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=__doc__,prog='prepareForKaks.py')
    parser.add_argument('--type',action='store',dest='anchor_type',choices=['self','dual'],required=True,type=str,default='self',
        help='type of your anchor file, default = self')
    parser.add_argument('--in_anchor',action='store',dest='anchor_file',required=True,type=str,
        help='input anchor file')    
    parser.add_argument('--in_cds',action='store',dest='in_cds',required=True,type=str,
        help='input cds file or paired cds files, example: <--in_pep A.cds> when --type self and <--in_pep A.cds,B.cds> when --type dual')
    parser.add_argument('--in_pep',action='store',dest='in_pep',required=True,type=str,
        help='input pep file or paired pep files, example: <--in_pep A.pep> when --type self and <--in_pep A.pep,B.pep> when --type dual')
    parser.add_argument('--out_cds',action='store',dest='out_cds',required=False,type=str,default='anchor.cds',
        help='output cds file')
    parser.add_argument('--out_pep',action='store',dest='out_pep',required=False,type=str,default='anchor.pep',
        help='output pep file')
    args = parser.parse_args()

    if args.anchor_type == 'self':
        in_cds = args.in_cds
        in_pep = args.in_pep
        cds_dict = getFastaDic(args.in_cds)
        pep_dict = getFastaDic(args.in_pep)
        anchor_list = getAnchorList(args.anchor_file)
        outputMergedFiles(cds_dict,cds_dict,pep_dict,pep_dict,anchor_list[0],anchor_list[1],args.out_cds,args.out_pep)
        print('Finished!')

    elif args.anchor_type == 'dual':
        in_cds_info = args.in_cds
        in_pep_info = args.in_pep
        in_cds = in_cds_info.split(',')
        in_pep = in_pep_info.split(',')
        cds_dict_1 = getFastaDic(in_cds[0])
        pep_dict_1 = getFastaDic(in_pep[0])
        cds_dict_2 = getFastaDic(in_cds[1])
        pep_dict_2 = getFastaDic(in_pep[1])
        anchor_list = getAnchorList(args.anchor_file)
        outputMergedFiles(cds_dict_1,cds_dict_2,pep_dict_1,pep_dict_2,anchor_list[0],anchor_list[1],args.out_cds,args.out_pep)
        print('Finished!')