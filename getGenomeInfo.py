import sys,getopt

def usage():
    print('usage:python3 getGenomeInfo.py -i <in_fasta> -o <outfile> <-h>')
    return

def getGenomeInfo(in_fasta_file,out_file):
    chr_len = {}
    contig_len = {}
    with open (in_fasta_file) as in_file:
        flag = 0
        for line in in_file:
            if '>' in line:
                chrname = line.strip('>')
                contig_len[chrname] = []
                chr_len[chrname] = 0
                flag = 0
                
            else:
                li = line.strip()
                chr_len[chrname] += len(li)
                infors = list(li)
                for infor in infors:
                    if infor not in ['n','N']:
                        if flag == 0:
                            contig_len[chrname].append(1)
                            flag = 1
                        else:
                            contig_len[chrname][-1] += 1
                    else:
                        flag = 0
    with open (out_file,'w') as outfile:
        all_contig = []
        for k,v in contig_len.items():
            len_list = contig_len[k]
            len_list.sort()
            all_contig += len_list
            outfile.write(k + '    contig_num:' + str(len(len_list)) + '    chr_len:' + str(chr_len[k]) + '\n')
        all_contig.sort()
        all_contig.reverse()
        print (all_contig)
        total_len = 0
        current_len = 0
        for con in all_contig:
            total_len += con
        for con in all_contig:
            current_len += con
            if current_len/total_len >= 0.5:
                n50 = con
                break

        outfile.write('n50:' + str(n50))

def main(argv):
    
    try:
        opts, args = getopt.getopt(argv,'hi:o:')
    except getopt.GetoptError:
        usage()
        sys.exit()

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt == '-i':
            in_fasta_name = arg
        elif opt == '-o':
            outfile_name = arg
    try:
        getGenomeInfo(in_fasta_name,outfile_name)
    except UnboundLocalError:
        usage()

    return

if __name__ == '__main__':
    main(sys.argv[1:])
