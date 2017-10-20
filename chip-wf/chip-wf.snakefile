# ** snakemake ** #
"""chip-seq workflow in progress"""

import pybedtools 
from pybedtools.featurefuncs import gff2bed
import pandas as pd

configfile: "config.yaml"
fbgn = config["fbgn"]
workdir: "../output/chip/" 

rule all:
    input: expand("{NAME}_liftover.bed", NAME=fbgn)

#I think this will only work if I use headers ?? 
rule cut_table:
#separate needed information from metadata
    input: expand("{fbgn}", fbgn=fbgn)
    output: cut = temp("{fbgn}_cut"),
            metadata = "{fbgn}_metadata"
    run:
        big = pd.read_table(input[0])
        bedtrim = big[['chrom', 'start','end','id']]
        metacolumns = ['id']
        for val in big.columns.values:
            if val not in bedtrim.columns.values:
               metacolumns.append(val)
        meta = big[metacolumns]
        bedtrim.to_csv(output.cut, sep='\t', header=None, index=False)
        meta.to_csv(output.metadata, sep='\t', index=False)

rule liftover:
#liftover to current release
    input: oldFile = temp("{fbgn}_cut"),
           mapchain = "../dm6ToDm3.over.chain"
    output: newFile = "{fbgn}_liftover.bed",
            unMapped = "{fbgn}_liftover_unmapped"
    shell: "liftOver {input.oldFile} {input.mapchain} {output.newFile} {output.unMapped}"

rule gff2bed: 
#get gene file as bed
    input: "../dmel-all-r6.12.gene_only.chr.gff"
    output: "dmel6.12.genes.bed"
    run:
        genes = pybedtools.BedTool(input[0])
	beds = genes.each(gff2bed, name_field='ID').saveas().to_dataframe()
        beds.to_csv(output[0], sep='\t', header=None, index=False)

rule intersect:
#intersect to obtain gene information
    input: gene = "dmel6.12.genes.bed",
           peaks = "{fbgn}_liftover.bed" 
    output: "{fbgn}_geneintersect.bed"
    shell: "bedtools intersect -a {input.gene} -b {input.peaks} -wb > {output}"

rule merge_metadata: 
#Necessary because of liftover
    input: geneintersect = "{fbgn}_geneintersect.bed",
           metadata = "{fbgn}_metadata"
    output: "{fbgn}_merged"
    run: 
        intersect = pd.read_table(input.geneintersect, header=None)[[3,6,7,8,9]]
        meta = pd.read_table(input.metadata) 
        intersect.columns= ['gene','chrom','start','end','id']
        intersect.merge(meta, on='id', how='left').to_csv(output[0], sep='\t', index=False)

#rule merge_phantompeaks: 
    #input: "{fbgn}_geneintersect.bed"
    #output: ""
    #run:
    #    probably will do some dataframe stuff but not sure 
