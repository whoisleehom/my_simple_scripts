__author__ = 'Yuan GAO'
__organization__ = 'AGIS'

import argparse

def generateUniqSeqByKvalue(kvalue:int,seq_type:str,max_len:int):
    seq_list = []
    nuc_list = []
    kmer_dict = {}
    if seq_type == 'DNA':
        nuc_dict = {'0':'C','1':'G','2':'A','3':'T'}
    if seq_type == 'RNA':
        nuc_dict = {'0':'C','1':'G','2':'A','3':'U'}
    for i in range(4**kvalue):
        list = []
        x = i
        while x > 3:
            list.append(str(x%4))
            x = x // 4
        if x:
            list.append(str(x))
        x  = ''.join(reversed(list)).zfill(kvalue)
        kmer_dict[x] = 0
    sc_num = 0
    out_seq = ''.zfill(kvalue)
    for i in range(1,4**kvalue):
        temp_kmer = out_seq[1-kvalue:]

        flag = 0
        for j in range(4):
            kmer = temp_kmer + str(j)
            if kmer_dict[kmer] == 0:
                if len(out_seq) < max_len:
                    kmer_dict[kmer] = 1
                    out_seq += str(j)
                    sc_num += 1
                    flag = 1
                else:
                    seq_list.append(out_seq)
                    out_seq = kmer
                break

        if flag == 0:
            for k,v in kmer_dict.items():
                if v == 0:
                    if len(out_seq) <= max_len-kvalue:
                        out_seq += k
                        kmer_dict[k] = 1
                    else:
                        seq_list.append(out_seq)
                        out_seq = k
                    break
    seq_list.append(out_seq)

    for seq in seq_list:
        nuc_seq = ''
        for i in range(len(seq)):
            nuc_seq += nuc_dict[seq[i]]
        nuc_list.append(nuc_seq)

    return nuc_list


def getArgs():
    parser = argparse.ArgumentParser(description=__doc__,prog='uniqKmerGenerator.py')
    parser.add_argument('-k',action='store',dest='kvalue',required=True,type=int,
        help='your K value')
    parser.add_argument('-n',action='store',dest='nucl_type',required=True,type=str,
        help='the type of your nucleotide, DNA or RNA')
    parser.add_argument('-l',action='store',dest='max_len',required=True,type=int,
        help='max length of each sequence')
    parser.add_argument('-o',action='store',dest='outfile',required=True,type=str,
        help='head of your fasta file names')
    return parser.parse_args()


if __name__ == '__main__':
    args = getArgs()

    out_seq = generateUniqSeqByKvalue(args.kvalue,args.nucl_type,args.max_len)

    for i in range(len(out_seq)):
        with open ('{0}.{1}.fasta'.format(args.outfile,i+1),'w') as out_file:
            out_file.write('>{0}.{1}\t{2}bp\n{3}\n'.format(args.outfile,i+1,len(out_seq[i]),out_seq[i]))