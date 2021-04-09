import sys,argparse
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
                    inf = line.strip('\n')
                    out_fasta.write(inf.rstrip('-.*') + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,prog='polish.py')
    parser.add_argument('-i',action='store',dest='in_file',required=True,type=str,
        help='full path to original fasta file')
    parser.add_argument('-o',action='store',dest='out_file',type=str,default='',
        help='full path to output file, default=#in_file_name#.simple')
    args = parser.parse_args()
    polish(args.in_file,args.out_file)

