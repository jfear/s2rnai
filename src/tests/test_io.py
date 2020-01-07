from s2rnai.io import GffRow

GFF_ROW = "3R\tFlyBase\tgene\t7924323\t7967408\t.\t-\t.\tID=FBgn0000504;Name=dsx;fullname=doublesex;Alias=doublesex,DSX,Dmdsx,Dsx,intersex-62c,Hermaphrodite,ix-62c,Hr,CG11094,Doublesex,dsxF,dsxM,double sex,DSXM,DSXF"
GTF_ROW = '3R\tFlyBase\tgene\t7924323\t7967408\t.\t-\t.\tgene_id "FBgn0000504"; gene_symbol "dsx";'


def test_GffRow_gff():
    gff = GffRow(GFF_ROW)
    assert gff.parse_attributes["Name"] == "dsx"


def test_GffRow_gtf():
    gtf = GffRow(GTF_ROW)
    assert gtf.parse_attributes["gene_symbol"] == "dsx"
