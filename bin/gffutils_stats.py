#!/usr/bin/env python

import pandas as pd
import gffutils
import argparse

def arguments():
    #input file and output file arguments are defined
    parser = argparse.ArgumentParser(description='takes file ??')
    parser.add_argument("-i", "--input", dest="infile", action='store', required=True,
                        help='input gffutils database')
    parser.add_argument("-o", "--output", dest="outfile", action='store',required=True, help='file to write output table to')
    args = parser.parse_args()
    return args

def main():
    #import command line arguments
    args = arguments()
    db = gffutils.FeatureDB(args.infile)

    # number of transcripts per gene
    table1 = []
    for gene in db.features_of_type('gene'):
        name = gene.id
        transcript = db.children(gene, featuretype='transcript')
        transcriptnum = len(list(transcript))
        row = [name, transcriptnum]
        table1.append(row)
    df1 = pd.DataFrame(table1, columns=['gene', '#transcripts'])
    newdf = df1.describe()

    # Get list of transcript IDs from mRNA
    tsIds = []
    for ts in list(db.features_of_type('mRNA')):
        tsIds.append(ts.attributes['transcript_id'][0])
    # make list of ts that are mRNA
    tss = []
    for ts in db.features_of_type('transcript'):
        if ts.id in tsIds:
            tss.append(ts)
    # get some stats on mRNA transcripts
    table = []
    for trns in tss:
        name = trns.id
        exons = list(db.children(trns, featuretype='exon'))
        # length of 1st exon for each transcript
        firstex = exons[0]
        exlen = firstex.end - firstex.start
        # number of introns/exons
        numex = len(exons)
        numin = (numex - 1)
        stuff = [name, exlen, numex, numin]
        table.append(stuff)
    df = pd.DataFrame(table, columns=['transcript', '1st_ex_len', '#ex', '#int'])
    tlist = df['transcript']
    # length of coding sequence
    table2 = []
    cdslen = 0
    for x in tlist:
        if list(db.children(x, featuretype='CDS')):
            for i in db.children(x, featuretype='CDS'):
                len1 = (len(i) - 1)
                cdslen += len1
        else:
            cdslen = 0
        newline = [x, cdslen]
        table2.append(newline)
    # length of utrs using transcript
    table3 = []
    utrlen = 0
    for x in tlist:
        if list(db.children(x, featuretype='3UTR')):
            for utr in db.children(x, featuretype='3UTR'):
                len2 = len(utr) - 1
                utrlen += len2
        else:
            utrlen = 0
        newline = [x, utrlen]
        table3.append(newline)
    table4 = []
    utrlen = 0
    for x in tlist:
        if list(db.children(x, featuretype='5UTR')):
            for utr in db.children(x, featuretype='5UTR'):
                len2 = len(utr) - 1
                utrlen += len2
        else:
            utrlen = 0
        newline = [x, utrlen]
        table4.append(newline)
    utr3frame = pd.DataFrame(table3, columns=['transcript', '3UTR_len'])
    utr5frame = pd.DataFrame(table4, columns=['transcript', '5UTR_len'])
    df3 = pd.DataFrame(table2, columns=['transcript', 'cds_len'])
    bigdf = df.merge(df3, how='left', on='transcript').merge(utr3frame, how='left', on='transcript').merge(utr5frame, how='left',
                                                                                                           on='transcript')
    #concatenate:
    concatframe = pd.concat([bigdf.describe(), newdf], axis=1).T
    concatframe.to_csv(args.outfile)


if __name__ == '__main__':
    main()