#! /bin/bash

# Last updated 05/27/23 by Alejandro Casco
#
# NAME: RNAseq_pipeline
# VERSION: 0.1
# AUTHOR: (c) 2023 Alejandro Casco
# DESCRIPTION: Short-read RNA-seq Consensus Pipeline (RCP) based on NASA GeneLab (PMID: 33870146)
# FEATURES:	- Run Data pre-processing and processing of multiple fastq files at once.
# DEPENDENCIES:	- FastQC, MultiQC, Cutadapt, TrimGalore, STAR, RSEM
#
# LICENSE:      GNU GPLv3 (http://www.gnu.de/documents/gpl-3.0.en.html)
#
# NOTICE:       THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW. 
#               EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES 
#               PROVIDE THE PROGRAM “AS IS” WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR 
#               IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY 
#               AND FITNESS FOR A PARTICULAR PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND 
#               PERFORMANCE OF THE PROGRAM IS WITH YOU. SHOULD THE PROGRAM PROVE DEFECTIVE,
#               YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.
#
#               IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING WILL ANY 
#               COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES AND/OR CONVEYS THE PROGRAM AS 
#               PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, 
#               INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE 
#               THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING RENDERED 
#               INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A FAILURE OF THE 
#               PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH HOLDER OR OTHER 
#               PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
#
# USAGE:	Edit variables under 'Variables and Settings'. May need to change path to 'Programs'.
#		bash RNAseq_pipeline
#		If indexes already present, can overwrite or continue with existing

###### Variables and Settings #######

WD="/path/to/working/directory"			# working directory
IDX_DIR="/path/to/index"		# directory containing/to build STAR and RSEM indexes
FASTA="/path/to//GRCh38.p14.ERCC.M81_DFLR.chrEBV.inverted.fa" # reference genome
GTF="/path/to/gencode.v45.primary_assembly.ERCC.M81_DFLR.chrEBV.inverted.gtf" # reference annotation
SAMPINFO_FILE="/path/to/sampleinfo.txt"	# sample info meta file
FASTQ_DIR="/path/to/fastq/directory"		# directory containing  FASTQ files
RDLT="101"	# Read length
CORES="15"		# Threads
LIBRARY="Paired-End"	# Paired-End or Single-End
FR="_R1"		# (Paired-End only): Forward Read Identifier (e.g., _R1)
RR="_R2"		# (Paired-End only): Reverse Read Ientifiery (e.g., _R2)
STRAND="reverse"	# strandedness: none, revere, or forward

############# Programs ##############

STAR=$(which STAR)
RSEM=$(which rsem-calculate-expression)
RSEM_PREP_REF=$(which rsem-prepare-reference)
SAMTOOLS=$(which samtools)

########### script start #############

##### Usage Checks

# Check working directory exists
if [[ ! -d $WD ]]; then
	echo "Error: working directory does not exist."
	ABORT="YES"
fi

# Check if FASTA file exists
if [[ ! -f $FASTA ]]; then
	echo "Error: FASTA file does not exist."
	ABORT="YES"
fi

# Check if GTF file exists
if [[ ! -f $GTF ]]; then
	echo "Error: GTF file does not exist."
	ABORT="YES"
fi

# Check CORES
re='^[0-9]+$'
if ! [[ $CORES =~ $re ]]; then
	echo "Error: threads input is not a number: $CORES"
	ABORT="YES"
fi

# Check read length
if ! [[ $RDLT =~ $re ]]; then
	echo "Error: read length input is not a number: $RDLT"
fi

# Check if FASTQs directory exists and whether gzipped FASTQs exist within it
if [[ ! -d $FASTQ_DIR ]]; then
	echo "Error: fastq directory does not exist."
	ABORT="YES"
else
	if [[ $(find $FASTQ_DIR -maxdepth 1 -type f -name "*fastq.gz" | wc -l) == 0 ]]; then
		echo "Error: no fastq.gz files detected fastq directory."
		ABORT="YES"
	fi
fi

# Check SAMPINFO file exists
if [[ ! -f $SAMPINFO_FILE ]]; then
	echo "Error: sample treatment matrix file does not exist."
	ABORT="YES"
else
	# Check SAMPINFO file format
	SAMPINFO_CHK=($(sed -re 's/[^\t]//g' "$SAMPINFO_FILE" | expand -t1 | wc -l -c -L))
	SAMPINFO_HEADER=$(head -n1 $SAMPINFO_FILE | tr ' ' '\t' | tr ',' '\t')
	if (( ((SAMPINFO_CHK[0] * 3) == SAMPINFO_CHK[1]) && (SAMPINFO_CHK[2] == 2) )); then
		:
	else
		# does not contain 3 columns
		if (( SAMPINFO_CHK[2] == 3 )); then
			echo "Error: incorrect sample treatment matrix file format. File does not contain 3 columns."
			ABORT="YES"
		fi
		# is not tab-delimited
		if (( (SAMPINFO_CHK[0] * 3) != (SAMPINFO_CHK[1]) )); then
			echo "Error: incorrect sample treatment matrix file format. File is no tab-delimited."
			ABORT="YES"
		fi
		# header does not include UNIQUE_ID Condition Replicate
		if [[ "$SAMPINFO_HEADER" != "$(echo -e "UNIQUE_ID\tCondition\tReplicate\r")" ]]; then
			echo "Error: incorrect sample treatment matrix file format. Column headers do not match \"UNIQUE_ID\", \"Condition\", and \"Replicate\"."
			ABORT="YES"
		fi
	fi
fi

# Check CORES, LIBRARY, and STRAND input are correct
re='^[0-9]+$'
if ! [[ $CORES =~ $re ]]; then
	echo "Error: threads input is not a number: $CORES."
	ABORT="YES"
else
	if [[ $CORES > 20 ]]; then
		echo "Error: number of threads ($CORES) is too high. Input <20."
		ABORT="YES"
	fi
fi
if [[ "$LIBRARY" != "Single-End" && "$LIBRARY" != "Paired-End" ]]; then
	echo "Error: incorrect library type. Input \"Single-End\" or \"Paired-End\"."
	ABORT="YES"
fi
if [[ "$STRAND" != "none" && "$STRAND" != "reverse" && $STRAND != "forward" ]]; then
	echo "Error: incorrect strandedness. Input \"none\", \"reverse\", or \"forward\"."
fi
# Check if STAR exists
if  [[ ! -f $STAR ]]; then
	echo "Error: STAR command does not exist."
	ABORT="YES"
fi
# Check if RSEM exists
if  [[ ! -f $RSEM ]]; then
	echo "Error: rsem-calculate-expression command does not exist."
	ABORT="YES"
fi
# Check if SAMtools exists
if [[ ! -f $SAMTOOLS ]]; then
	echo "Error: samtools command does not exist."
	ABORT="YES"
fi

# if ABORT is YES, exit
if [[ $ABORT == "YES" ]]; then
	exit
fi

# Create Index Directories
if [[ -e $IDX_DIR ]]; then
	printf "\n"
	read -rep $'Index directory already exists, overwite? (y/n):\n' INPUT
	if [[ $INPUT == "y" ]]; then
		printf "\n"
		rm -r $IDX_DIR

		mkdir $IDX_DIR
		mkdir ${IDX_DIR}/STAR_index
		mkdir ${IDX_DIR}/RSEM_index

		# Build STAR index
		$STAR \
			--runThreadN $CORES \
			--runMode genomeGenerate \
			--genomeDir ${IDX_DIR}/STAR_index \
			--genomeFastaFiles $FASTA \
			--sjdbGTFfile $GTF \
			--sjdbOverhang $(($RDLT-1)) \
			--outFileNamePrefix ${IDX_DIR}/STAR_index

		# Build RSEM index
		$RSEM_PREP_REF \
			--gtf $GTF \
			$FASTA \
			${IDX_DIR}/RSEM_index/rsem_index

		# Annotate README file storing info
		if [[ -f ${IDX_DIR}/README ]]; then
			rm ${IDX_DIR}/README
		fi
		DATE=$(date)

cat <<EOT >> ${IDX_DIR}/README
Date generated: $DATE
FASTA file: $FASTA
GTF file: $GTF

$( echo "STAR version:" $($STAR --version) )
cite: PMID: 23104886
indexing non-default parameters:
	--sjdbOverhang $(($RDLT-1))

$( echo "RSEM version:" $($RSEM --version) | cut -d " " -f1,2,6 | cut --complement -c15 )
cite: PMID: 21816040
indexing parameters: default
EOT
	fi
fi

# Prepare sampleinfo list
UNIQUE_ID_LIST=$(cat $SAMPINFO_FILE | sed '1d' | awk '{print $1}')

# Create list of fastq files
FASTQS=$( find $FASTQ_DIR -maxdepth 1 -type f -name "*fastq.gz" -exec ls {} + )

##### Check fastq files for single-end or paired-end
for unique_name in $( echo "$UNIQUE_ID_LIST" ); do
	if [[ $LIBRARY == "Single-End" ]]; then
		# check if possible paired-end
		if [[ $(echo $(ls $FASTQS | grep $unique_name | wc -l)) -gt 1 ]]; then
			echo "Error: UNIQUE_ID is not unique to one fastq file. Data may be paired-end." && exit
		fi
	elif [[ $LIBRARY == "Paired-End" ]]; then
		# check if possible single-end
		if [[ $(echo $(ls -1 $FASTQS | sed 's#.*/##' | grep $unique_name | wc -l )) -lt 2 ]]; then
			echo "Error: one more more paired-end samples contain only one fastq file. Data may be single-end." && exit
		fi
		# check if more than 2 fastq files given unique name
		if [[ $(echo $(ls -1 $FASTQS | sed 's#.*/##' | grep $unique_name | wc -l )) -gt 2 ]]; then
			echo "Error: UNIQUE_ID is present in more than two paired-end fastq files." && exit
		fi
		# Check each sample has a forward and reverse read fastq file
		if [[ $(echo $(ls -1 $FASTQS | sed 's#.*/##' | grep $unique_name | grep $FR | wc -l )) != $(echo $(ls -1 $FASTQS | sed 's#.*/##' | grep $unique_name | grep $RR | wc -l )) ]]; then
			echo "Error: one or more paired-end samples do not contain a Forward and Reverse read fastq file." && exit
		fi
	fi
done

########### Start Analyses ###########

# Record start time and date
STARTDATE=$(date)
SECONDS=0

##### STAR
#Create STAR output directory
if [[ ! -e ${WD}/STAR ]]; then
	mkdir ${WD}/STAR
fi

# STAR mapping
for unique_name in $( echo "$UNIQUE_ID_LIST" ); do
	if [[ $LIBRARY == "Single-End" ]]; then
		R1=$( ls -1 $FASTQS | sed 's#.*/##' | grep $unique_name )
		R2=""
		R1=${FASTQ_DIR}/$R1
		STAROPTION1=""
	fi
	if [[ $LIBRARY == "Paired-End" ]]; then
		R1=$( ls -1 $FASTQS | sed 's#.*/##' | grep $unique_name | grep $FR )
		R2=$( ls -1 $FASTQS | sed 's#.*/##' | grep $unique_name | grep $RR )
		R1=${WD}/trimmed_fastqFiles/$R1
		R2=${WD}/trimmed_fastqFiles/$R2
		R1=$( ls -1 $FASTQS | sed 's#.*/##' | grep $unique_name | grep $FR )
		R2=$( ls -1 $FASTQS | sed 's#.*/##' | grep $unique_name | grep $RR )
		R1=${FASTQ_DIR}/$R1
		R2=${FASTQ_DIR}/$R2
		STAROPTION1="--alignMatesGapMax 1000000"
	fi
	
	# if temp STAR directory exists, delete
	if [[ -d ${WD}/STAR/${unique_name}_tmp ]]; then
		rm -r ${WD}/STAR/${unique_name}_tmp
	fi

	# Run STAR alignment
	STAR \
		--twopassMode Basic \
		--limitBAMsortRAM 0 \
		--genomeDir ${IDX_DIR}/STAR_index \
		--outSAMunmapped Within \
		--outFilterType BySJout \
		--outSAMattributes NH HI AS NM MD MC \
		--outFilterMultimapNmax 20 \
		--outFilterMismatchNmax 999 \
		--alignIntronMin 20 \
		--alignIntronMax 1000000 \
		$STAROPTION1 \
		--alignSJoverhangMin 8 \
		--alignSJDBoverhangMin 1 \
		--sjdbScore 1 \
		--readFilesCommand zcat \
		--runThreadN $CORES \
		--outSAMtype BAM SortedByCoordinate \
		--quantMode TranscriptomeSAM \
		--outSAMheaderHD @HD VN:1.4 SO:coordinate \
		--outFileNamePrefix ${WD}/STAR/${unique_name}. \
		--outTmpDir ${WD}/STAR/${unique_name}_tmp \
		--readFilesIn $R1 $R2
done

# Create list of transcriptome BAM files
BAMS=$( find ${WD}/STAR -maxdepth 1 -type f -name "*.Aligned.toTranscriptome.out.bam" -exec ls {} + )

##### RSEM
# Create RSEM output directory
if [[ ! -e ${WD}/RSEM ]]; then
	mkdir ${WD}/RSEM
fi
if [[ $LIBRARY == "Single-End" ]]; then
	RSEMLIB=""
elif [[ $LIBRARY == "Paired-End" ]]; then
	RSEMLIB="--paired-end"
fi

# Run RSEM transcript quantification
for unique_name in $( echo "$UNIQUE_ID_LIST" ); do
	rsem-calculate-expression \
		--num-threads $CORES \
		--alignments \
		--bam $RSEMLIB \
		--seed 12345 \
		--estimate-rspd \
		--no-bam-output \
		--strandedness $STRAND \
		${WD}/STAR/${unique_name}.Aligned.toTranscriptome.out.bam \
		${IDX_DIR}/RSEM_index/rsem_index  \
		${WD}/RSEM/$unique_name
done

ENDDATE=$(date)
duration=$SECONDS

##### Annotate output reference file
# If reference file exists, delete
DATE=$(date +%m)$(date +%d)$(date +%y)
if [[ -f ${WD}/RNAseq_completed_${DATE}.log ]]; then
	rm ${WD}/RNAseq_completed_${DATE}.log
fi

cat <<EOT >> ${WD}/RNAseq_completed_${DATE}.log
START DATE: $STARTDATE
END DATE: $ENDDATE
Duration of run: $(($duration / 3600)) hours, $((($duration / 60) % 60)) minutes and $(($duration % 60)) seconds.

Working directory: $WD
Fastq directory: $FASTQ_DIR
Index directory: $IDX_DIR
Threads: $CORES
$(cat ${IDX_DIR}/README | grep -A1 "FASTA")

$( echo "STAR version:" $($STAR --version) )
cite: PMID: 23104886
$(cat ${IDX_DIR}/README | grep -A1 "indexing non-default")
mapping non-default parameters:
	--twopassMode Basic
	--outSAMunmapped Within
	--outFilterType BySJout
	--outSAMattributes NH HI AS NM MD MC
	--outFilterMultimapNmax 20
	--outFilterMismatchNmax 999
	--outFilterMismatchNoverReadLmax 0.04
	--alignIntronMin 20
	--alignIntronMax 1000000
	--alignSJoverhangMin 8
	--alignSJDBoverhangMin 1
	--sjdbScore 1
	--quantMode TranscriptomeSAM
	--outSAMheaderHD @HD VN1:1.4 SO:coordinate
	$STAROPTION1

$( echo "RSEM version:" $($RSEM --version) | cut -d " " -f1,2,6 | cut --complement -c15 )
cite: PMID: 21816040
indexing parameters: default
quantification non-default parameters:
	--alignments
	--seed 12345
	--estimate-rspd
	--strandedness $STRAND

EOT




