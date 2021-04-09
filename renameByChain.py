import sys,getopt

def usage():
    print('usage:python3 renameByChain.py -i <old_file> -o <new_file> -c <chain file> <-h> <-r reverse by chain>')
    return

def rename_by_chain(in_file,out_file,chain,flag):
    chain_dic = {}
    with open (chain) as in_chain:
        if flag == 1:
            for line in in_chain:
                li = line.rstrip('\n')
                infors = li.split('>')
                chain_dic[infors[0]] = infors[1]
        else:
            for line in in_chain:
                li = line.rstrip('\n')
                infors = li.split('>')
                chain_dic[infors[1]] = infors[0]
    
    with open (in_file) as infile:
        oldfile = infile.read().splitlines()
        with open (out_file,'w') as outfile:
            for i in range(len(oldfile)):
                fl = 0
                for k,v in chain_dic.items():
                    if k in oldfile[i]:
                        new_line = oldfile[i].replace(k,chain_dic[k])
                        outfile.write(new_line + '\n')
                        fl = 1
                        break
                if fl == 0:
                    outfile.write(oldfile[i] + '\n')

def main(argv):
    
    try:
        opts, args = getopt.getopt(argv,'rhi:o:c:')
    except getopt.GetoptError:
        usage()
        sys.exit()
    flag = 1
    for opt, arg in opts:
    
        if opt == '-h':
            usage()
            sys.exit()
        elif opt == '-i':
            in_fasta_name = arg
        elif opt == '-o':
            outfile_name = arg
        elif opt == '-c':
            chain_file = arg
        elif opt == '-r':
            flag = 0
    try:
        rename_by_chain(in_fasta_name,outfile_name,chain_file,flag)
    except UnboundLocalError:
        usage()

    return

if __name__ == '__main__':
    main(sys.argv[1:])