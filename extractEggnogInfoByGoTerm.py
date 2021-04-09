__author__ = 'Yuan GAO'
__organization__ = 'Fujian Agriculture and Forestry University HIST'
__doc__ = (
    """
    作者：Yuan GAO  版本:v1.0.0  2020/8/24    本程序所实现的功能：
    根据GO term、富集结果和eggnog注释文件提取基因信息
    """
    )
import os
import argparse

def get_eggnog_dict(eggnog_file):
    eggnog_dict = {}
    with open(eggnog_file) as in_eggnog:
        for line in in_eggnog:
            if '#' not in line:
                infors = line.split('\t')
                query_name = infors[0]
                go_terms = infors[5]
                annotation = infors[-1].strip('\n')
                eggnog_dict[query_name] = {'GO_terms':go_terms,'annotation':annotation}
    return(eggnog_dict)

def extract_by_enrichment(go_term,enrichment_file,eggnog_dict,output_file):
    enrich_dict = {}
    with open (enrichment_file) as in_enrich:
        in_en = in_enrich.read().splitlines()
        for line in in_en[1:]:
            infors = line.split('\t')
            enrich_dict[infors[0]] = infors[8].split('/')
    with open (output_file,'w') as outfile:
        outfile.write('query_name\tannotation\tGO_terms\n')
        for line in enrich_dict[go_term]:
            outfile.write('{0}\t{1}\t{2}\n'.format(line,eggnog_dict[line]['annotation'],eggnog_dict[line]['GO_terms']))

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=__doc__,prog='extractEggnogInfoByGOTerm.py')
    parser.add_argument('-e',action='store',dest='eggnog_file',required=True,type=str,
        help='full path to eggnog annotation file')
    parser.add_argument('-r',action='store',dest='enrichment_file',required=True,type=str,
        help='full path to GO enrichment file')
    parser.add_argument('-g',action='store',dest='go_term',required=True,type=str,
        help='your destination GO term')
    parser.add_argument('-o',action='store',dest='output_file',required=True,type=str,
        help='full path to output file')
    args = parser.parse_args()

    eggnog_dict = get_eggnog_dict(args.eggnog_file)
    extract_by_enrichment(args.go_term,args.enrichment_file,eggnog_dict,args.output_file)

    print('Finished!')
