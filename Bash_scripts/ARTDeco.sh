#! /bin/bash

##### PREREQUISTES #####

# 1. Test gtf2bed runs on the GTF, otherwise fix GTF according to https://github.com/sjroth/ARTDeco

#preparing-files

WD="/path/to/working/directory/ARTDeco"
BAMS="/path/to/bam/files/directory"
CORES=10
FASTA="/path/to/GRCh38.p14.ERCC.fa"
GTF="/path/to/gencode.v45.primary_assembly.ERCC.modified_genes.gtf"
META="${WD}/ARTDeco_meta.txt"
COMP="${WD}/ARTDeco_comparisons.txt"

##### Script start #####

# Generate chromosome sizes file
if [[ ! -e ${FASTA}.fai ]]; then
	echo "Indexing FASTA"
	samtools faidx $FASTA
fi
echo "Generating chrom.sizes"
FASTA_ID=$(basename $FASTA)
CHR_SIZES=$(dirname $FASTA)/${FASTA_ID::-3}.chrom.sizes
cut -f1,2 ${FASTA}.fai > $CHR_SIZES

# Create ARTDeco directory
if [[ ! -e $WD ]]; then
	mkdir $WD
fi

printf "\n"
echo "Running ARTDeco with differential expression"
ARTDeco \
	-home-dir $WD \
	-bam-files-dir $BAMS \
	-gtf-file $GTF \
	-cpu $CORES \
	-chrom-sizes-file $CHR_SIZES \
	-layout PE \
	-stranded True \
	-orientation Reverse \
	-skip-bam-summary \
	-meta-file $META \
	-comparisons-file $COMP

printf "\n"

echo "Done."
