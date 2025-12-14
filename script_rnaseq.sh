##!bin/bash
SAMPLES=("mut_lib1" "mut_lib2" "WT_lib1" "WT_lib2")
N_THREADS=4
ADAPTER_FILE="TruSeq3-PE.fa"
HISAT2_INDEX="DM6"
TRIMMOMATIC_JAR="/usr/local/bin/trimmomatic.jar"

RAW_DIR="data/raw"
QC_DIR="qc"
TRIMMED_DIR="trimmed"
ALIGNED_DIR="alignment"
LOGS_DIR="logs"
MULTI_QC_DIR="multiqc"
QC_TRIMMED_DIR=$QC_DIR/"QC_trimmed"
MULTI_QC_TRIMMED_DIR=$MULTI_QC_DIR/"multi_qc_trimmed"
ANNOTATION_GTF="data/reference/Drosophila_melanogaster.BDGP6.46.110.gtf"
COUNT_DIR="comptage"

#créer les dossiers si pas existants
mkdir -p $RAW_DIR $QC_DIR $TRIMMED_DIR $ALIGNED_DIR $LOGS_DIR $MULTI_QC_DIR $QC_TRIMMED_DIR $MULTI_QC_TRIMMED_DIR \
	 $COUNT_DIR
#=============
#boucle de traitement
#===================
for SAMPLE in "${SAMPLES[@]}"; do
	R1_RAW="$RAW_DIR/${SAMPLE}_R1.fq.gz"
	R2_RAW="$RAW_DIR/${SAMPLE}_R2.fq.gz"
	R1_PAIRED="$TRIMMED_DIR/${SAMPLE}_R1_paired.fq.gz"
	R2_PAIRED="$TRIMMED_DIR/${SAMPLE}_R2_paired.fq.gz"
	R1_UNPAIRED="$TRIMMED_DIR/${SAMPLE}_R1_unpaired.fq.gz" # Non utilisé pour alignement
	R2_UNPAIRED="$TRIMMED_DIR/${SAMPLE}_R2_unpaired.fq.gz" # Non utilisé pour alignement	

	LOG_TRIM="$LOGS_DIR/${SAMPLE}_trimmomatic.log"
	LOG_HISAT="$LOGS_DIR/${SAMPLE}_hisat2.log"
	FEATURE_COUNTS_LOG="$LOGS_DIR/${SAMPLE}_featurecount.log"
	SAM_OUTPUT="$ALIGNED_DIR/${SAMPLE}.sam"
	BAM_OUTPUT="$ALIGNED_DIR/${SAMPLE}.bam"
	BAM_OUTPUT_SORTED="$ALIGNED_DIR/${SAMPLE}_sorted.bam"
	COUNT_OUTPUT="$COUNT_DIR/${SAMPLE}_counts.txt"
# Controle qualité avant trimming

	fastqc $R1_RAW $R2_RAW -o $QC_DIR
	echo "controle qualité terminé"


#nettoyage des données

	java -jar  $TRIMMOMATIC_JAR PE -threads $N_THREADS -phred33 \
		$R1_RAW $R2_RAW \
		$R1_PAIRED $R1_UNPAIRED \
		$R2_PAIRED $R2_UNPAIRED \
		ILLUMINACLIP:$ADAPTER_FILE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
		2> $LOG_TRIM

	echo "Nettoyage terminé"

 controle qualité après trimming
	fastqc $R1_PAIRED $R2_PAIRED -o $QC_TRIMMED_DIR

	# alignement

	hisat2 -p $N_THREADS -x $HISAT2_INDEX \
	-1 $R1_PAIRED \
	-2 $R2_PAIRED \
	-S $SAM_OUTPUT \
	2> $LOG_HISAT

	echo "alignement terminé"

# conversion sam to bam
	samtools view -bS $SAM_OUTPUT -o $BAM_OUTPUT
	
	echo "bam crées"
#lancemement featureCounts
	featureCounts -T $N_THREADS -p \
		-a $ANNOTATION_GTF \
		-o $COUNT_OUTPUT \
		$BAM_OUTPUT
		2> $FEATURE_COUNTS_LOG
	echo "feature counts terminé"
done

multiqc $QC_DIR -o $MULTI_QC_DIR

multiqc $QC_TRIMMED_DIR -o $MULTI_QC_TRIMMED_DIR


