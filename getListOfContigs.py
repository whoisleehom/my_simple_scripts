import sys,getopt

def usage():
    print('usage:python3 getLengthOfContig.py -i <in_fasta> -o <outfile> <-h>')
    return

def getLen(in_fasta_name,outfile_name):

    with open (in_fasta_name) as infile:
        with open (outfile_name,'w') as outfile:
            for line in infile:
                if '>' in line:
                    name = line.strip('>')
                    outfile.write(name)
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