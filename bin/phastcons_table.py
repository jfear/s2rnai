#!/usr/bin/env python
#parses transcription factor, gene name/symbol, chrom, start/end, score, & phastcon score
#from input file and creates data frame with mean phastcon scores. writes to output file specified by user
import pandas as pd
import argparse
import re

def arguments():
    #input file and output file arguments are defined
    parser = argparse.ArgumentParser(description='takes file ??')
    parser.add_argument("-i", "--input", dest="filename", action='store', required=True,
                        help='input file obtained from bedtools intersect w/phastcons')
    parser.add_argument("-o", "--output", dest="new", action='store',required=True, help='file to write output to')

    args = parser.parse_args()
    return args


def main():
    #import command line arguments
    args= arguments()
    values = []
    with open(args.filename) as f:
        for line in f:
            pattern = re.compile(
                r'\w*\t(\w*)\t(\w*)\t(FBgn\w*)\t(\S*)\t.*(chr\w*)\tFlyBase\sgene\t\w*\t\w*.*ID=(\w*);Name=(\w*).*\t(\S*)\t\w*')
            match = pattern.match(line)
            TF = match.group(3)
            qval = float(match.group(4))
            chrom = match.group(5)
            start = match.group(1)
            end = match.group(2)
            symbol = match.group(7)
            FBgn = match.group(6)
            phastcon = float(match.group(8))
            reorder = (TF, FBgn, symbol, chrom, start, end, qval, phastcon)
            values.append(reorder)
    df = pd.DataFrame(values, columns=['TF', 'FBgn', 'Symbol', 'Chrom', 'Start', 'End', 'q Value', 'Phastcons'])
    grp = df.groupby(['TF', 'FBgn', 'Symbol', 'Chrom', 'Start', 'End', 'q Value'])
    meanframe = grp.mean()
    meanframe.to_csv(args.new, sep='\t') #write dataframe to output file, tab separated

if __name__ == '__main__':
    main()