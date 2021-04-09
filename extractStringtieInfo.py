__author__ = 'Yuan GAO'
__organization__ = 'Fujian Agriculture and Forestry University HIST'
__doc__ = (
    """
    作者：Yuan GAO  版本:v1.0.1  2020/8/22    本程序所实现的功能：
    根据序列名列表（基因或转录本）从stringtie生成的注释文件中提取必要信息
    """
    )
import os
import argparse

def build_dict(in_count_file,list_kind):
    trans_dict = {}
    gene_dict = {}
    with open (in_count_file) as in_file:
        for line in in_file:
            if '#' not in line:
                infors = line.split('\t')
                if 'reference_id' in infors[8]:
                    in1 = infors[8].split('; ')
                    for ln1 in in1:
                        if 'reference_id' in ln1:
                            ln2 = ln1.split('"')
                            reference_id = ln2[1]
                        elif 'ref_gene_id' in ln1:
                            ln2 = ln1.split('"')
                            ref_gene_id = ln2[1]
                        elif 'FPKM' in ln1:
                            ln2 = ln1.split('"')
                            FPKM = ln2[1]
                        elif 'TPM' in ln1:
                            ln2 = ln1.split('"')
                            TPM = ln2[1]
                        elif 'cov' in ln1:
                            ln2 = ln1.split('"')
                            cov = ln2[1]

                    try:
                        gene_dict[ref_gene_id]
                    except KeyError:
                        try:
                            gene_dict[ref_gene_id] = {reference_id:{'cov':float(cov),'FPKM':float(FPKM),'TPM':float(TPM)}}
                        except UnboundLocalError:
                            gene_dict[ref_gene_id] = {reference_id:{'cov':float(cov),'FPKM':0.0,'TPM':0.0}}
                    else:
                        try:
                            gene_dict[ref_gene_id][reference_id] = {'cov':float(cov),'FPKM':float(FPKM),'TPM':float(TPM)}
                        except UnboundLocalError:
                            gene_dict[ref_gene_id][reference_id] = {'cov':float(cov),'FPKM':0.0,'TPM':0.0}

                    trans_dict[reference_id] = {'cov':float(cov),'FPKM':float(FPKM),'TPM':float(TPM)}

    if list_kind == 'trans':
        return trans_dict
    elif list_kind == 'gene':
        return gene_dict

def build_dict_by_dir(dir_path,list_kind):
    file_names = os.listdir(dir_path)
    trans_dict = {}
    gene_dict = {}
    for in_count_file in file_names:
        with open (dir_path + '/' + in_count_file) as in_file:
            for line in in_file:
                if '#' not in line:
                    infors = line.split('\t')
                    if 'reference_id' in infors[8]:
                        in1 = infors[8].split('; ')
                        for ln1 in in1:
                            if 'reference_id' in ln1:
                                ln2 = ln1.split('"')
                                reference_id = ln2[1]
                            elif 'ref_gene_id' in ln1:
                                ln2 = ln1.split('"')
                                ref_gene_id = ln2[1]
                            elif 'FPKM' in ln1:
                                ln2 = ln1.split('"')
                                FPKM = ln2[1]
                            elif 'TPM' in ln1:
                                ln2 = ln1.split('"')
                                TPM = ln2[1]
                            elif 'cov' in ln1:
                                ln2 = ln1.split('"')
                                cov = ln2[1]

                        try:
                            gene_dict[ref_gene_id]
                        except KeyError:
                            a = {}
                            for fi in file_names:
                                a[fi] = {'cov':0,'FPKM':0,'TPM':0}
                            gene_dict[ref_gene_id] = {reference_id:a}
                            try:
                                gene_dict[ref_gene_id][reference_id][in_count_file] = {'cov':float(cov),'FPKM':float(FPKM),'TPM':float(TPM)}
                            except UnboundLocalError:
                                gene_dict[ref_gene_id][reference_id][in_count_file] = {{'cov':float(cov),'FPKM':0.0,'TPM':0.0}}
                        else:
                            try:
                                gene_dict[ref_gene_id][reference_id]
                            except KeyError:
                                a = {}
                                for fi in file_names:
                                    a[fi] = {'cov':0,'FPKM':0,'TPM':0}
                                gene_dict[ref_gene_id][reference_id]=a
                                try:
                                    gene_dict[ref_gene_id][reference_id][in_count_file] = {'cov':float(cov),'FPKM':float(FPKM),'TPM':float(TPM)}
                                except UnboundLocalError:
                                    gene_dict[ref_gene_id][reference_id][in_count_file] = {'cov':float(cov),'FPKM':0.0,'TPM':0.0}
                            else:
                                try:
                                    gene_dict[ref_gene_id][reference_id][in_count_file] = {'cov':float(cov),'FPKM':float(FPKM),'TPM':float(TPM)}
                                except UnboundLocalError:
                                    gene_dict[ref_gene_id][reference_id][in_count_file] = {'cov':float(cov),'FPKM':0.0,'TPM':0.0}


                        
                        try:
                            trans_dict[reference_id]
                        except KeyError:
                            a = {}
                            for fi in file_names:
                                a[fi] = {'cov':0,'FPKM':0,'TPM':0}
                            trans_dict[reference_id] = a
                            try:
                                trans_dict[reference_id][in_count_file] = {'cov':float(cov),'FPKM':float(FPKM),'TPM':float(TPM)}
                            except UnboundLocalError:
                                trans_dict[reference_id][in_count_file] = {'cov':float(cov),'FPKM':0.0,'TPM':0.0}
                        else:
                            try:
                                trans_dict[reference_id][in_count_file] = {'cov':float(cov),'FPKM':float(FPKM),'TPM':float(TPM)}
                            except UnboundLocalError:
                                trans_dict[reference_id][in_count_file] = {'cov':float(cov),'FPKM':0.0,'TPM':0.0}



    if list_kind == 'trans':
        return trans_dict
    elif list_kind == 'gene':
        return gene_dict

def output_message(in_dict,in_list,list_kind,out_file_name):
    if list_kind == 'trans':
        with open(out_file_name,'w') as outfile:
            with open (in_list) as inlist:
                outfile.write('trans_name,cov,FPKM,TPM\n')
                for line in inlist:
                    trans_name = line.strip()
                    try:
                        a = in_dict[trans_name]
                    except KeyError:
                        outfile.write(trans_name + ',0,0,0\n')
                    else:
                        outfile.write('{0},{1},{2},{3}\n'.format(trans_name,a['cov'],a['FPKM'],a['TPM']))
    else:
        with open(out_file_name,'w') as outfile:
            with open (in_list) as inlist:
                outfile.write('gene_name,trans_name,cov,FPKM,TPM\n')
                for line in inlist:
                    gene_name = line.strip()
                    try:
                        a = in_dict[gene_name]
                    except KeyError:
                        outfile.write(gene_name + ',NA,0,0,0\n')
                    else:
                        for k,v in a.items():
                            outfile.write('{0},{1},{2},{3},{4}\n'.format(gene_name,k,v['cov'],v['FPKM'],v['TPM']))

def output_message_by_dir(in_dict,in_list,dir_path,list_kind,out_file_name):
    file_names = os.listdir(dir_path)
    outfile = open(out_file_name,'w')
    inlist = open (in_list)
    if list_kind == 'trans':
        
        outfile.write('trans_name')
        for line in file_names:
            outfile.write(',' + line)
        outfile.write('\n')

        for line in inlist:
            trans_name = line.strip()
            outfile.write(trans_name)
            for fi in file_names:
                try:
                    a = in_dict[trans_name][fi]
                except KeyError:
                    outfile.write(',0;0;0\n')
                else:
                    outfile.write(',{0};{1};{2},'.format(a['cov'],a['FPKM'],a['TPM']))
            outfile.write('\n')
        
    else:
        outfile.write('gene_name,trans_name')
        for line in file_names:
            outfile.write(',' + line)
        outfile.write('\n')
        for line in inlist:
            gene_name = line.strip()
            try:
                a = in_dict[gene_name]
            except KeyError:
                outfile.write(gene_name + ',NA')
                for fi in file_names:
                    outfile.write(',0;0;0')
                outfile.write('\n')
            else:
                for k,v in a.items():
                    outfile.write('{0},{1}'.format(gene_name,k))
                    for fi in file_names:
                        b = v[fi]
                        outfile.write(',{0};{1};{2}'.format(b['cov'],b['FPKM'],b['TPM']))
                    outfile.write('\n')
    outfile.close()
    inlist.close()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=__doc__,prog='extractStringtieInfo.py')
    parser.add_argument('-s',action='store',dest='in_count_file',required=False,type=str,default='',
        help='full path to stringtie file')
    parser.add_argument('-d',action='store',dest='dir_path',required=False,type=str,default='',
        help='full path to the directory of your stingtie files')
    parser.add_argument('-l',action='store',dest='in_list',required=True,type=str,
        help='full path to gene or transcription list')
    parser.add_argument('-o',action='store',dest='out_file_name',required=False,type=str,default='',
        help='full path to gene or transcription list, default: <IN_STRINGTIE_FILE.csv>')
    parser.add_argument('-g',action='store',dest='list_kind',required=True,type=str,
        help='option to specify your list, use <-g gene> to specify a gene list and <-g trans> to specify a transcription list')
    args = parser.parse_args()

    if args.in_count_file == '':
        if args.dir_path == '':
            print('error:-s or -d is needed')
            exit()
        else:
            if args.out_file_name == '':
                out_file_name = args.dir_path + '.csv'
            else:
                out_file_name = args.out_file_name
            out_dict = build_dict_by_dir(args.dir_path,args.list_kind)
            output_message_by_dir(out_dict,args.in_list,args.dir_path,args.list_kind,out_file_name)
    else:
        if args.dir_path == '':
            if args.out_file_name == '':
                out_file_name = args.in_count_file + '.csv'
            else:
                out_file_name = args.out_file_name
            out_dict = build_dict(args.in_count_file,args.list_kind)
            output_message(out_dict,args.in_list,args.list_kind,out_file_name)
        else:
            print('error:use only -s or -d')
            exit()