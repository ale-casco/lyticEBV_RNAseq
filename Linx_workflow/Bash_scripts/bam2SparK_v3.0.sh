#! /bin/bash

##### Setup

# Inputs
INTERVALS="/path/to/intervals.txt"
BAMDIR="/path/to/genome/bams/directory"
SAMPINFO="/path/to/sampleinfo_IGV.txt" # Run DESeq2_ERCCnorm_SizeFactors.R
STRAND="reverse"	# <reverse|forward>
OUTDIR="/path/to/output/directory"
GTF="/path/to/gencode.v45.primary_assembly.M81_DFLR.chrEBV.inverted.ERCC.gtf"
CORES=14
ALL_GENES="yes"	# <no|yes>
KEEP_BDGS="yes"	# <no|yes>

# Tools 
SPARK="/path/to/SparK/SparK.py"
SAMTOOLS=$(which samtools)
BAMCOVERAGE=$(which bamCoverage)
PYTHON=$(which python)
BEDTOOLS=$(which bedtools)

##### Run

# Create output directory
OUTDIR="${OUTDIR}/SparK"
if [[ ! -e $OUTDIR ]]; then
	mkdir $OUTDIR
fi

# Generate temporary sample info file ordered by sample in column 2
if [[ -e ${OUTDIR}/tmp_sampleinfo.txt ]]; then
	rm ${OUTDIR}/tmp_sampleinfo.txt
fi
head -n 1 $SAMPINFO > ${OUTDIR}/tmp_sampleinfo.txt	# create tmp sample info file with header
CONDITIONS=( $(sed '1d' $SAMPINFO | awk '{print $2}' | awk '!seen[$0]++') ) # uniq without sorting
for ((i=0;i<=$(($(echo ${#CONDITIONS[@]})-1));i++)); do
	cat $SAMPINFO | awk -v var="${CONDITIONS[$i]}" '$2 == var' >> ${OUTDIR}/tmp_sampleinfo.txt
done
SAMPINFO_TMP=${OUTDIR}/tmp_sampleinfo.txt

# Prepare arrays from sampleinfo file
UNIQUE_ID=( $(sed '1d' $SAMPINFO_TMP | awk '{print $1}') )
NORM=( $(sed '1d' $SAMPINFO_TMP | awk '{print $4}') )

# Index bam files
printf "\nIndexing bam files...\n"
for ((i=0;i<=$(($(echo ${#UNIQUE_ID[@]})-1));i++)); do
	bam="$BAMDIR/${UNIQUE_ID[$i]}.Aligned.sortedByCoord.out.bam"
	if [[ ! -e ${bam}.bai ]]; then
		$SAMTOOLS index $bam
	fi
done

# Create arrays from interval file
CHR_LS=( $(sed '1d' $INTERVALS | awk '{print $1}') )
GENE_LS=( $(sed '1d' $INTERVALS | awk '{print $2}') )
START_LS=( $(sed '1d' $INTERVALS | awk '{print $3}') )
END_LS=( $(sed '1d' $INTERVALS | awk '{print $4}') )
GROUP_SCALE_LS=( $(sed '1d' $INTERVALS | awk '{print $5}') )

# Loop through interval file
for ((j=0;j<=$(($(echo ${#CHR_LS[@]})-1));j++)); do
	CHR=${CHR_LS[$j]}
	GENE=${GENE_LS[$j]}
	START=${START_LS[$j]}
	END=${END_LS[$j]}
	GROUP_SCALE=${GROUP_SCALE_LS[$j]}
	
	# Start analysis
	printf "\nInitiating analysis for $GENE gene on ${CHR}:${START}-$END\n"
	
	# Create output directory
	SPARK_OUT=${OUTDIR}/$GENE
	if [[ ! -e $SPARK_OUT ]]; then
		mkdir $SPARK_OUT
	fi

	# Create temporary bedGraph directory
	BDG_OUT="${SPARK_OUT}/bedGraphs"
	if [[ ! -e $BDG_OUT ]]; then
		mkdir $BDG_OUT
	fi
	
	# Set strandedness for second read file
	if [[ $STRAND == reverse ]]; then
		STRAND2="forward"
	else
		STRAND2="reverse"
	fi
	
	# Run bamCoverage to get bedragph files covering specific region
	echo "Generating bedGraphs files..."
	for ((i=0;i<=$(($(echo ${#UNIQUE_ID[@]})-1));i++)); do
		bam="$BAMDIR/${UNIQUE_ID[$i]}.Aligned.sortedByCoord.out.bam"
		$BAMCOVERAGE -b $bam \
			-bs 1 \
			-of bedgraph \
			--numberOfProcessors $CORES \
			--region ${CHR}:${START}:${END} \
			--scaleFactor ${NORM[$i]} \
			--filterRNAstrand $STRAND \
			-o ${BDG_OUT}/${UNIQUE_ID[$i]}.${GENE}-R.bdg \
			&>/dev/null

		$BAMCOVERAGE -b $bam \
			-bs 1 \
			-of bedgraph \
			--numberOfProcessors $CORES \
			--region ${CHR}:${START}:${END} \
			--scaleFactor ${NORM[$i]} \
			--filterRNAstrand $STRAND2 \
			-o ${BDG_OUT}/${UNIQUE_ID[$i]}.${GENE}-F.bdg \
			&>/dev/null
	done

	# Take average
	UNIQUE_CONDITIONS=$(echo "${CONDITIONS[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ')
	echo "Generating averaged bedGraph files..."
	for cond in $UNIQUE_CONDITIONS; do
		UNIQUE_tmp=$(sed '1d' $SAMPINFO_TMP | grep "$cond\b" | awk '{print $1}')
		UNIQUE_tmp=$(printf '%s|' $UNIQUE_tmp | sed 's/|$//')
		files_F=$(find ${BDG_OUT} -maxdepth 1 -type f -name "*F.bdg" -exec ls {} + | grep -E $UNIQUE_tmp)
		files_R=$(find ${BDG_OUT} -maxdepth 1 -type f -name "*R.bdg" -exec ls {} + | grep -E $UNIQUE_tmp)
		
		echo $files_F
			
#		$BEDTOOLS unionbedg -i $files_F | awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' > ${BDG_OUT}/${cond}.avg.${GENE}-F.bdg
#		$BEDTOOLS unionbedg -i $files_R | awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' > ${BDG_OUT}/${cond}.avg.${GENE}-R.bdg
	done
	
	# Get guide for SparK -cg parameter
	for ((i=0;i<=$(($(echo ${#CONDITIONS[@]})-1));i++)); do
		GUIDE+=$(printf "$((i+1)) "%.0s)
	done
	
	# Create lists of output file names
	set -- ${CONDITIONS[@]}
	set -- "${@/#/$BDG_OUT/}" && set -- "${@/%/.avg.${GENE}-R.bdg}"
	OUT_R=$(echo "$@")
	OUT_F=$(echo "$OUT_R" | sed "s/R.bdg/F.bdg/g")

	# Set SparK parameters
	if [[ $ALL_GENES != yes ]]; then
		par1=$(echo "-dg $GENE")
	fi
	if [[ $GROUP_SCALE == yes ]]; then
		par2=$(echo "-gs yes")
	else
		par2=""
	fi
	
	# Run SparK to generate forward, reverse, and forward+reverse tracks
	printf "Generating data tracks...\n"
	$PYTHON $SPARK \
		 -pr ${CHR}:${START}-${END} \
		 -cg $GUIDE \
		 -gl ${CONDITIONS[@]} \
		 $par1 \
		 $par2 \
		 -gtf $GTF \
		 -cf $OUT_F \
		 -o ${SPARK_OUT}/${GENE}_RNAseq_pileup_F \
		 -dt all \
		 > /dev/null
	
	$PYTHON $SPARK \
		-pr ${CHR}:${START}-${END} \
		-cg $GUIDE \
		-gl ${CONDITIONS[@]} \
		$par1 \
		$par2 \
		-gtf $GTF \
		-cf $OUT_R \
		-o ${SPARK_OUT}/${GENE}_RNAseq_pileup_R \
		> /dev/null

	$PYTHON $SPARK \
		-pr ${CHR}:${START}-${END} \
		-cg $GUIDE \
		-tg $GUIDE \
		-gl ${CONDITIONS[@]} \
		$par1 \
		$par2 \
		-gtf $GTF \
		-cf $OUT_R \
		-tf $OUT_F \
		-o ${SPARK_OUT}/${GENE}_RNAseq_pileup \
		-pt sine \
		> /dev/null

	if [[ $KEEP_BDGS = no ]]; then
		rm -r $BDG_OUT
	fi
	GUIDE=""
done
printf "\n"
rm $SAMPINFO_TMP

