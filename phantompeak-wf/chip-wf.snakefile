# ** snakemake ** #
"""chip-seq workflow in progress"""

import pybedtools 
from pybedtools.featurefuncs import gff2bed
import pandas as pd

configfile: "config.yaml"
fbgn = config["fbgn"]
workdir: "../output/chip/" 

rule all:
    input: expand("out/{NAME}_nophantom", NAME=fbgn)

rule cut_table:
#separate needed information from metadata
    input: expand("{fbgn}.bed", fbgn=fbgn)
    output: cut = temp("{fbgn}_cut"),
            metadata = temp("{fbgn}_metadata")
    run:
        big = pd.read_table(input[0], header=None) #will this mess up if there is header? 
        big['id'] = ['id='+str(x) for x in range(0, len(big[0]))]
        bedtrim = big[[0, 1, 2,'id']]
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
    output: newFile = temp("{fbgn}_liftover.bed"),
            unMapped = "out/{fbgn}_liftover_unmapped" #do we want this? probably
    shell: "liftOver {input.oldFile} {input.mapchain} {output.newFile} {output.unMapped}"

rule gff2bed: 
#get gene file as bed
    input: "../dmel-all-r6.12.gene_only.chr.gff"
    output: "dmel6.12.genes.bed" #should i make this temp? 
    run:
        genes = pybedtools.BedTool(input[0])
	beds = genes.each(gff2bed, name_field='ID').saveas().to_dataframe()
        beds.to_csv(output[0], sep='\t', header=None, index=False)

rule intersect:
#intersect to obtain gene information
    input: gene = "dmel6.12.genes.bed",
           peaks = temp("{fbgn}_liftover.bed") 
    output: temp("{fbgn}_geneintersect.bed")
    shell: "bedtools intersect -a {input.gene} -b {input.peaks} -wb > {output}"

rule merge_metadata: 
#Necessary because of liftover
    input: geneintersect = temp("{fbgn}_geneintersect.bed"),
           metadata = temp("{fbgn}_metadata")
    output: "out/{fbgn}_merged"
    run: 
        intersect = pd.read_table(input.geneintersect, header=None)[[3,6,7,8,9]]
        meta = pd.read_table(input.metadata) 
        intersect.columns= ['gene','chrom','start','end','id']
        intersect.merge(meta, on='id', how='left').to_csv(output[0], sep='\t', index=False)

rule check_phantompeaks:
#determine the amount of phantom peak overlap in the dataset  
    input: merged = "out/{fbgn}_merged",
           phantompeaks = "gkv637_Supplementary_Data/Supplementary_table_3__List_of_Phantom_Peaks.xlsx"
    output: #withpeaks = "out/{fbgn}_withpeaks", 
            nophantom = "out/{fbgn}_nophantom"
    run:
         phantom = pd.read_excel(input.phantompeaks)[['chr ','start','end','Name']]
         merged = pd.read_table(input.merged)
         bed = merged[['chrom','start','end']]
         intersect = pybedtools.BedTool.from_dataframe(bed).intersect(pybedtools.BedTool.from_dataframe(phantom), wo=True).to_dataframe()
         filtered = intersect[intersect.thickEnd >= 50]
         #filtered.to_csv(output.withpeaks, sep='\t', index=False) Don't think we need this
         outermerge = merged.merge(filtered, how='outer', on=['start','end'], indicator=True)
         no_phantom = outermerge[outermerge._merge == 'left_only'].to_csv(output.nophantom, sep='\t', header=None, index=False)
