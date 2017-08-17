#!/usr/bin/env python
#cuts out gene specific sequences from fasta file, outputs new fasta file

import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse

def arguments():
    #input file (and output file if needed) arguments are defined
    parser = argparse.ArgumentParser(description='takes file ??')
    parser.add_argument("--gff", dest="gff", action='store', required=True,
                        help='needs a gff file of gene info w/slopped 1kb')
    parser.add_argument("--fasta", dest="fasta", action='store', required=True,
                        help='needs a fasta file')
    parser.add_argument("-o", "--output", dest="output", action='store',required=True,
                        help='name the fasta file to write out to')
    args = parser.parse_args()
    return args

def main():
    args = arguments()
    f = SeqIO.parse(args.fasta, "fasta")
    my_dict = SeqIO.to_dict(f)

    my_records = []
    with open(args.gff) as g:
        for line in g:
            splitline = line.split('\t')
            chrom = splitline[0]
            start = int(splitline[3])
            stop = int(splitline[4])
            lesschrom = chrom[3:]
            descrip = splitline[8]
            match = re.compile('ID=(\w*);.*').match(descrip)
            name = match.group(1)
            if lesschrom == 'M':
                look_seq = my_dict['mitochondrion_genome'].seq
            else:
                look_seq = my_dict[lesschrom].seq
            seqslice = look_seq[start - 1:stop]
            new_record = SeqRecord(seqslice, id=name, description="")
            my_records.append(new_record)

    SeqIO.write(my_records, args.output, "fasta")

if __name__ == '__main__':
    main()