import gffutils
from gffutils.pybedtools_integration import to_bedtool, featurefuncs

CNT = 1


def main():
    db = gffutils.FeatureDB(snakemake.input.db)
    gene = to_bedtool(db.features_of_type("gene")).saveas()
    slopped = gene.slop(b=100, genome="dm6")
    merged = slopped.sort().merge()
    complement = merged.complement(genome="dm6").saveas()
    bed = complement.each(interName).saveas(snakemake.output.bed)
    bed.each(interGFF).saveas(snakemake.output.gtf)


def interName(feature):
    global CNT

    feature = featurefuncs.extend_fields(feature, 4)
    feature.name = "intergenic{}".format(CNT)
    CNT += 1
    return feature


def interGFF(feature):
    gff = featurefuncs.bed2gff(feature)
    gff[1] = "bedtools"
    gff[2] = "gene"
    gff.attrs["gene_id"] = gff.name
    return gff


if __name__ == "__main__":
    main()
