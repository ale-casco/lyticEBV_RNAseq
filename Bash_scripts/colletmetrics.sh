#! /bin/bash

## Set inputs

OUTPUT="path/to/output/directory"
BAMPATH="path/to/genome/bams/directory"
REFS="path/to/refs"

# EBV
GTF_EBV="/path/to/M81_DFLR.chrEBV.inverted.gtf"
FASTA_EBV="/path/to/M81_DFLR.chrEBV.inverted.fa"
CHR_EBV="chrEBV" # name of EBV chr

# Host
GTF="/path/to/gencode.v45.primary_assembly.annotation.gtf"
FASTA="/path/to/GRCh38.p14.genome.fa"
CHRS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM" # restrict to these chr

## Set paths to programs

JAVA="/path/to/jdk-21.0.1/bin/java"
PICARD="/path/to/picard/build/libs/picard.jar"
GTF2GPRED="/path/to/gtfToGenePred"

## Prepare filterd host files

# Create refs dir
REFS=${REFS}/metrics
if [[ ! -e $REFS ]]; then
	mkdir $REFS
fi

# Filter FASTA
if [[ ! -e ${FASTA}.fai ]]; then
	samtools faidx $FASTA
fi
samtools faidx $FASTA $CHRS > ${REFS}/host.fa
FASTA=${REFS}/host.fa

# Filter GTF
CHRS_GREP=$(echo ${CHRS// /|})
grep -E -w "$CHRS_GREP" $GTF > ${REFS}/host.gtf
GTF=${REFS}/host.gtf

## Convert GTF to GenePred using USCS tool

# Set temporary output file names
GTFBASE=$(basename $GTF)
REFFLAT="${REFS}/${GTFBASE::-4}.refflat"
GTFBASE_EBV=$(basename $GTF_EBV)
REFFLAT_EBV="${REFS}/${GTFBASE_EBV::-4}.refflat"

# Generate host RefFlat file
if [[ ! -e $REFFLAT ]]; then
	echo "Generating RefFlat file."
	$GTF2GPRED -genePredExt -geneNameAsName2 -ignoreGroupsWithoutExons $GTF /dev/stdout | awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > $REFFLAT
else
	echo "RefFlat file already exists. Using existing file."
fi

# Generate EBV RefFlat file
if [[ ! -e $REFFLAT_EBV ]]; then
	echo "Generating EBV RefFlat file."
	$GTF2GPRED -genePredExt -geneNameAsName2 -ignoreGroupsWithoutExons $GTF_EBV /dev/stdout | awk 'BEGIN { OFS="\t"} {print $1, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > $REFFLAT_EBV
else
	echo "EBV RefFlat file already exists. Using existing file."
fi
printf "\n"

## Get ribosomal intervals

# Calculate chrom sizes
CHROM_SIZES="${REFS}/chrom_sizes"
if [[ ! -e $CHROM_SIZES ]]; then
	echo "Generating chrom sizes file."
	if [[ ! -e ${FASTA}.fai ]]; then
		echo "Indexing FASTA file"
		samtools faidx $FASTA
	fi
	cut -f1,2 ${FASTA}.fai > $CHROM_SIZES
else
	echo "Chrom sizes file already exists. Using existing file."
fi
printf "\n"

# rRNA interval_list file suitable for Picard CollectRnaSeqMetrics
rRNA="${REFS}/${GTFBASE::-4}.rRNA.interval_list"
if [[ ! -e $rRNA ]]; then
	echo "Generating rRNA interval list file suitable for Picard CollectRnaSeqMetrics."
	perl -lane 'print "\@SQ\tSN:$F[0]\tLN:$F[1]\tAS:hg38"' $CHROM_SIZES | grep -v _ >> $rRNA
	grep 'gene_type "rRNA"' $GTF | awk '$3 == "transcript"' | cut -f1,4,5,7,9 | perl -lane '/transcript_id "([^"]+)"/ or die "no transcript_id on $."; print join "\t", (@F[0,1,2,3], $1)' | sort -k1V -k2n -k3n >> $rRNA
else
	echo "rRNA interval list file already exists. Using existing file."
fi
printf "\n"

## Run Picard CollectRnaSeqMetrics

# Set bed file variables
BED=$(basename $FASTA)
BED="${REFS}/${BED::-3}.bed"
BED_EBV=$(basename $FASTA_EBV)
BED_EBV="${REFS}/${BED_EBV::-3}.bed"

# Prepare BED file from host FASTA index
if [[ ! -e $BED ]]; then
	cat ${FASTA}.fai | awk 'BEGIN { OFS="\t"} {print $1, 0, $2}' > $BED
fi

# Prepare BED file from ebv FASTA index
if [[ ! -e $BED_EBV ]]; then
	if [[ ! -e ${FASTA_EBV}.fai ]]; then
		samtools faidx $FASTA_EBV
	fi
	cat ${FASTA_EBV}.fai | awk 'BEGIN { OFS="\t"} {print $1, 0, $2}' > $BED_EBV
fi

# Set header file variables
BEDBASE=$(basename $BED)
CHR_GREP="${REFS}/${BEDBASE::-4}.header_grep.txt"
BEDBASE_EBV=$(basename $BED_EBV)
CHR_GREP_EBV="${REFS}/${BEDBASE_EBV::-4}.header_grep.txt"

# Prepare host header
if [[ ! -e $CHR_GREP ]]; then
	cat $BED | awk '{print $1}' | awk '$1="SN:"$1' > $CHR_GREP
fi

# Prepare ebv header
if [[ ! -e $CHR_GREP_EBV ]]; then
	cat $BED_EBV | awk '{print $1}' | awk '$1="SN:"$1' > $CHR_GREP_EBV
fi

# Set variables for metrics output
BAMS=$( find $BAMPATH -maxdepth 1 -type f -name "*.Aligned.sortedByCoord.out.bam" -exec ls {} + )
BAMVEC=$( ls $BAMS | sed 's!.*/!!' | rev | cut -c 31- | rev | sed 's!.*/!!' )
CHR_HEADER="${REFS}/${BEDBASE::-4}.header.sam"
CHR_HEADER_EBV="${REFS}/${BEDBASE_EBV::-4}.header.sam"

# Create output directories
if [[ ! -e $OUTPUT ]]; then
	mkdir ${OUTPUT}
fi
if [[ ! -e ${OUTPUT}/host ]]; then
	mkdir ${OUTPUT}/host

fi
if [[ ! -e ${OUTPUT}/ebv ]]; then
	mkdir ${OUTPUT}/ebv
fi

# Remove temporary directory if exists
if [[ -e ${OUTPUT}/tmp_dir ]]; then
	rm -r ${OUTPUT}/tmp_dir
fi

# Run Picard CollectRnaSeqMetrics
for samp in $BAMVEC; do
	if [[ ! -e ${OUTPUT}/host/${samp}.host.metrics ]]; then
		echo "Extracting host chromosomes from $samp"
		mkdir ${OUTPUT}/tmp_dir
		samtools view -H ${BAMPATH}/${samp}.Aligned.sortedByCoord.out.bam | grep -w -f $CHR_GREP > $CHR_HEADER
		samtools view -b -h -L $BED ${BAMPATH}/${samp}.Aligned.sortedByCoord.out.bam > ${OUTPUT}/tmp_dir/${samp}.filtered.bam
		samtools reheader $CHR_HEADER ${OUTPUT}/tmp_dir/${samp}.filtered.bam > ${OUTPUT}/tmp_dir/${samp}.filtered.reheader.bam
		echo "Generating metrics for host $samp"
		$JAVA -jar $PICARD CollectRnaSeqMetrics \
			I=${OUTPUT}/tmp_dir/${samp}.filtered.reheader.bam \
			O=${OUTPUT}/host/${samp}.host.metrics \
			REF_FLAT=$REFFLAT \
			STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
			RIBOSOMAL_INTERVALS=$rRNA \
			&> /dev/null
		rm -r ${OUTPUT}/tmp_dir
	fi
	
	if [[ ! -e ${OUTPUT}/ebv/${samp}.ebv.metrics ]]; then
		echo "Extracting ebv chromosome from $samp"
		mkdir ${OUTPUT}/tmp_dir
		samtools view -H ${BAMPATH}/${samp}.Aligned.sortedByCoord.out.bam | grep -w -f $CHR_GREP_EBV > $CHR_HEADER_EBV
		samtools view -b -h -L $BED_EBV ${BAMPATH}/${samp}.Aligned.sortedByCoord.out.bam > ${OUTPUT}/tmp_dir/${samp}.filtered.bam
		cat $CHR_HEADER_EBV <(samtools view ${OUTPUT}/tmp_dir/${samp}.filtered.bam) | samtools view -bo ${OUTPUT}/tmp_dir/${samp}.filtered.reheader.bam
	
		echo "Generating metrics for ebv $samp"
		$JAVA -jar $PICARD CollectRnaSeqMetrics \
			I=${OUTPUT}/tmp_dir/${samp}.filtered.reheader.bam \
			O=${OUTPUT}/ebv/${samp}.ebv.metrics \
			REF_FLAT=$REFFLAT_EBV \
			STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
			&> /dev/null
		rm -r ${OUTPUT}/tmp_dir
	fi
	printf "\n"
done

