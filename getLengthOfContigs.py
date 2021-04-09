import sys,getopt

def usage():
    print('usage:python3 getLengthOfContig.py -i <in_fasta> -o <outfile> <-h>')
    return

def getLen(in_fasta_name,outfile_name):

    with open (in_fasta_name) as infile:
        in_fasta = infile.read().splitlines()

    len_dict = {}
    cur = ''
    for line in in_fasta:
        if '>' in line:
            name = line.strip('>')
            cur = name
            len_dict[name] = 0
        else:
            len_dict[cur] += len(line)
    with open (outfile_name,'w') as outfile:
        for k,v in len_dict.items():
            outfile.write(k + ' ' + str(v) + '\n')
    
    return

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
        getLen(in_fasta_name,outfile_name)
    except UnboundLocalError:
        usage()

    return

if __name__ == '__main__':
    main(sys.argv[1:])