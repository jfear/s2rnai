#!/bin/bash

#input? 
#output into file:
#mkdir /Users/bergeric/data/big_script_out
#slop
bedtools slop -i /Users/bergeric/data/dmel-all-r6.12.gene_only.chr.gff -g /Users/bergeric/data/ChromInfo.txt -l 1000 -r 0 -s > /Users/bergeric/data/big_script_out/gene_only_slop.txt
python fimo2bed.py --i /Users/bergeric/data/motif_alignments_flyFactor_dm6.2L.txt
bedtools intersect -a /Users/bergeric/data/motif_alignments_flyFactor_dm6.2L.txt_.bed -b /Users/bergeric/data/big_script_out/gene_only_slop.txt -wb > /Users/bergeric/data/big_script_out/bedtoolsintersect_out.txt

#phastcons steps (ommitted grep for 2L and wig2bed)
#second intersect, this one is output of first intersect with phastcons table
bedtools intersect -a /Users/bergeric/data/big_script_out/bedtoolsintersect_out.txt -b /Users/bergeric/data/small.bed -wo > /Users/bergeric/data/big_script_out/dm6_phastcons_intersect.txt
#create nice output table w mean phastcons scores
python phastcons_table.py --input /Users/bergeric/data/dm6_phastcons_intersect.txt --output /Users/bergeric/data/big_script_out/practiceoutfile.txt
