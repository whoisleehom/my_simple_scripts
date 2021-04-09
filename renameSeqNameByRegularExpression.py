import re
import argparse

def rename(in_fasta,out_fasta,re_formula):
    with open (in_fasta) as infile:
        with open (out_fasta,'w') as outfile:
            for line in infile:
                if '>' in line:
                    match = re.findall(re_formula,line)
                    try:
                        match[0]
                    except IndexError:
                        print('ERROR: Cannot find your dest name by your re formula in ###{}###'.format(line.strip()))
                        outfile.write(line)
                    else:
                        outfile.write('>' + match[0] + '\n')
                else:
                    outfile.write(line)
    print(args.in_fasta + ' is renamed.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,prog='renameSeqNameByRegularExpression.py')
    parser.add_argument('-i',action='store',dest='in_fasta',required=True,type=str,
        help='path to your fasta file')
    parser.add_argument('-o',action='store',dest='out_fasta',type=str,default='',
        help='name of the renamed fasta, default = #in_file_name#.renamed')
    parser.add_argument('-r',action='store',dest='re',required=True,type=str,
        help='the regular expression formula of your destination')
    args = parser.parse_args()
    out_name = args.in_fasta + '.renamed' if args.out_fasta == '' else args.out_fasta
    rename(args.in_fasta,out_name,args.re)
    