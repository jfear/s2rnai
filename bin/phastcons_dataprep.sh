#!/bin/bash

#input comes from:
ROOTDIR=/Users/bergeric/Projects/s2rnai
#output into file:
WORKDIR=$ROOTDIR/output/phastcons_workflow
if [ -e $WORKDIR ]; then
mkdir -p $WORKDIR
fi

#slop
bedtools slop \
	-i $ROOTDIR/data/dmel-all-r6.12.gene_only.chr.gff \
	-g $ROOTDIR/data/ChromInfo.txt \
	-l 1000 \
	-r 0 \
	-s \
	> $WORKDIR/gene_only_slop.txt
#convert fimo to bed file type
python fimo2bed.py --i $ROOTDIR/data/motif_alignments_flyFactor_dm6.2L.txt --o $WORKDIR/motif_alignments_flyFactor_dm6.2L_BED.bed
#intersect with gene slop
bedtools intersect -a $WORKDIR/motif_alignments_flyFactor_dm6.2L_BED.bed -b $WORKDIR/gene_only_slop.txt -wb > $WORKDIR/bedtoolsintersect_out.txt
#phastcons steps (ommitted grep for 2L and wig2bed)
#second intersect, this one is output of first intersect with phastcons table
bedtools intersect -a $WORKDIR/bedtoolsintersect_out.txt -b $ROOTDIR/data/small.bed -wo > $WORKDIR/dm6_phastcons_intersect.txt
#create nice output table w mean phastcons scores
python phastcons_table.py --input $WORKDIR/dm6_phastcons_intersect.txt --output $WORKDIR/outfile.txt
