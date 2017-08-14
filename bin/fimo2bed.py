#!/usr/bin/env python
#convert fimo to BED

import argparse

def arguments():
    #input file (and output file if needed) arguments are defined
    parser = argparse.ArgumentParser(description='takes file ??')
    parser.add_argument("-i", "--input", dest="filename", action='store', required=True,
                        help='input file should be fimo')
    #parser.add_argument("-o", "--output", dest="new", action='store',required=True, help='bed')

    args = parser.parse_args()
    return args

def main():
    #import command line arguments
    args= arguments()
    with open(args.filename + '_.bed', 'w+') as outfile:
        with open(args.filename) as f:
            for line in f:
                if not line.startswith("#"):
                    s = line.split('\t')
                    chrom = s[1]; start = s[2]; end = s[3]; name = s[0]; qval = s[7]; strand = s[4]
                    reorder = chrom + '\t' + start + '\t' + end + '\t' + name + '\t' + qval + '\t' + strand + '\n'
                    outfile.write(reorder)

if __name__ == '__main__':
    main()