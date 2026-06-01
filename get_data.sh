


# Getting the scRNA data from website

mkdir expression
mkdir genotypes
cd ..
cd expression
mkdir -p h5ad
cd h5ad

wget -O onek1k.h5ad "https://datasets.cellxgene.cziscience.com/a3f5651f-cd1a-4d26-8165-74964b79b4f2.h5ad"


# Get real gencode gene annotation file

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz

# Use gencode annotations to create gene_annotation_nochr.tsv file for CONTENT input

zcat gencode.v46.annotation.gtf.gz | \
grep -v '^#' | \
awk '
BEGIN{
    OFS="\t";
    print "gene_id","gene_name","chrom","tss"
}
$3=="gene"{
    match($0,/gene_id "[^"]+"/)
    gene_id=substr($0,RSTART+9,RLENGTH-10)

    match($0,/gene_name "[^"]+"/)
    gene_name=substr($0,RSTART+11,RLENGTH-12)

    sub(/\..*/,"",gene_id)

    tss=($7=="+")?$4:$5

    chrom=$1
    sub(/^chr/,"",chrom)

    print gene_id,gene_name,chrom,tss
}' > gene_annotation_nochr.tsv
