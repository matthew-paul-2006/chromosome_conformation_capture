#! /bin/bash
#SBATCH --verbose
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --time=24:00:00
#SBATCH --job-name=repeat_in-silico
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=mrp420@nyu.edu
#SBATCH --output=/scratch/%u/reports/%x_%j.out
#SBATCH --error=/scratch/%u/reports/%x_%j.err

#------------------------------------------------------------------------------#
#                                INSTRUCTIONS                                  #
#------------------------------------------------------------------------------#

# Will take paired end data conformational data and map it to the C. elegans rDNA. 
# Will then report mapped 

### Argument options:
# EXPID     	Custom ID for output files
# RUNDIR    	Path to directory to run script and save output in
# GENNAME   	Basename of reference genome files preceded by absolute path to
#           	directory containing the files. Must include a FASTA file and a
#           	matching GFF file with the same basename.
#          		An existing Bowtie index with a basename matching the files' is used
#           	if found in the same directory; otherwise a new index is built
# REPGENNAME	Basename of repetitive region of the genome, with FASTA file and 
#				potentially bowtie index.
# FQ1        	Absolute path to input fastq file 1
# FQ2			Absolute path to input fastq file 1
# DIGENDS		Bed file showing the position of every restriction site
# READLEN		How long are the reads for the fastq file

### EXAMPLE:
#sbatch --export \
#EXPID="N2B22015",\
#RUNDIR="/scratch/mrp420/insilico",\
#GENNAME="/home/mrp420/worms/genomes/WS220/WS220_ucsc",\
#REPGENNAME="/home/mrp420/worms/genomes/rDNA/rDNA",\
#FQ1="/scratch/cgsb/ercan/GEO/2015_meyer/N2B22015_X_X_R1_X.fastq",\
#FQ2="/scratch/cgsb/ercan/GEO/2015_meyer/N2B22015_X_X_R2_X.fastq",\
#DIGENDS='GATC',\
#READLEN=100 \
#~/worms/scripts/repeat_in-silico4C.sh
 
#------------------------------------------------------------------------------#
#                                  Functions                                   #
#------------------------------------------------------------------------------#

function elapsed_time() {
    ENDTIME=$(date +%s)

    TIME=$(($ENDTIME - $1))
    if [ $TIME -lt 60 ]
    then
        echo "$TIME sec"
    elif [ $TIME -ge 60 ]  && [ $TIME -lt 3600 ]
    then
        echo "$(($TIME / 60)) min"
    else
        echo "$(($TIME / 60 / 60)) hr"
    fi
}

function check_arg() {
    if [ -z "$1" ]
    then
        echo ">>>>> Please provide values for all required arguments"
        exit 1
    else 
    	echo "Input argument is $1"
    fi
}

#------------------------------------------------------------------------------#
#                                  IO checks                                   #
#------------------------------------------------------------------------------#

# Check arguments
check_arg $EXPID
check_arg $RUNDIR 
check_arg $FQ1
check_arg $FQ2
check_arg $DIGENDS
check_arg $READLEN

# Check input files / dirs
[ -f $FQ ] || { echo "Could not find file: $FQ"; exit 1; }
[ -d $RUNDIR ] || { echo "Could not find directory: $RUNDIR"; exit 1; }
GENDIR=$(dirname "$GENNAME")
GENNAME=$(basename "$GENNAME")
[ -d $GENDIR ] || { echo "Could not find directory: $GENDIR"; exit 1; }
REPGENDIR=$(dirname "$REPGENNAME")
REPGENNAME=$(basename "$REPGENNAME")
[ -d $REPGENDIR ] || { echo "Could not find directory: $REPGENDIR"; exit 1; }


# Search for reference genome fasta file; exit if not found
FA=$(find $GENDIR -iname "${GENNAME}.fa*")

if [ -z "$FA" ]
then
    echo ">>>>> Could not find reference genome FASTA file"
    exit 1
fi

# Search for Bowtie index and build one if not found
# (a file named according to rule "fasta_base_name.1.ebwt")
# The following code will return the full basename
# (the provided $GENNAME might not include it in full)
IX=$(basename $FA)                              # Drop path to file
IX=${IX%.*}                                     # Drop extension
CHECKIX=$(find $GENDIR -iname "${IX}.1.ebwt")   # Search file
#IX=$(basename $IX | cut -d '.' -f 1)

if [ -z "$CHECKIX" ]
then
    echo ">>>>> Building Bowtie index..."
    module purge
    module load bowtie/gnu/1.2.1.1
    # Build index
    cd $GENDIR
    bowtie-build -f $FA $IX
fi


# Search for repeat fasta file; exit if not found
REPFA=$(find $REPGENDIR -iname "${REPGENNAME}.fa*")

if [ -z "$REPFA" ]
then
    echo ">>>>> Could not find repeat FASTA file"
    exit 1
fi

# Search for Bowtie index and build one if not found
# (a file named according to rule "fasta_base_name.1.ebwt")
# The following code will return the full basename
# (the provided $GENNAME might not include it in full)
REPIX=$(basename $REPFA)                              # Drop path to file
REPIX=${REPIX%.*}                                     # Drop extension
CHECKREPIX=$(find $REPGENDIR -iname "${REPIX}.1.ebwt")   # Search file
#REPIX=$(basename $REPIX | cut -d '.' -f 1)

if [ -z "$CHECKREPIX" ]
then
    echo ">>>>> Building Bowtie index..."
    module purge
    module load bowtie/gnu/1.2.1.1
    # Build index
    cd $REPGENDIR
    bowtie-build -f $REPFA $REPIX
fi

#------------------------------------------------------------------------------#
#                                                                              #
#                                Run pipeline                                  #
#                                                                              #
#------------------------------------------------------------------------------#

STARTTIME=$(date +%s)
echo \
"------------------------------------------------------------------------------"
echo ">>>>> Started Repeat In-silico: $EXPID"
echo \
"------------------------------------------------------------------------------"
date

#------------------------------------------------------------------------------#
#                  Align reads to reference genome with Bowtie                 #
#------------------------------------------------------------------------------#
cd $RUNDIR

# Abort if output directory already exists
if [ -d "$EXPID" ]
then
    echo ">>>>> Output directory already exists"
    exit 1
fi

#Define sub directories for storing mapped reads
mkdir ${EXPID}/
cd ${EXPID}/
mkdir R1_mapping
mkdir R2_mapping

echo ">>>>> Map reads with Bowtie:"
module purge
module load bowtie/gnu/1.2.1.1

#Define the amount of times we will need to map, based on the length of the reads
START=1
END=$(((${READLEN}-20)/10))

echo "      Iterative mapping of the first reads to reference..."

bowtie -q -5 1 -v 0 -m 1 -p 8 --seed=123 -S \
	--un R1_mapping/unmapped.fastq \
	$GENDIR/$IX ${FQ1} \
	R1_mapping/${EXPID}_R1_${READLEN}.sam

for (( c=$START; c<=$END; c++ ))
do
	echo -n "$c "
	#Define how much to trim 3' for mapping
	cut3prime=$(($c*5))
	#Create tmp file for unmapped reads
	cp R1_mapping/unmapped.fastq R1_mapping/tmp.fastq
	#Conduct mapping
	bowtie -q -5 1 -3 ${cut3prime} -v 0 -m 1 -p 8 --seed=123 -S \
	--un R1_mapping/unmapped.fastq \
	$GENDIR/$IX R1_mapping/tmp.fastq \
	R1_mapping/${EXPID}_R1_$((${READLEN}-${cut3prime})).sam
	sleep 1
done

echo "      Iterative mapping of the second reads to reference..."

bowtie -q -5 1 -v 0 -m 1 -p 8 --seed=123 -S \
	--un R1_mapping/unmapped.fastq \
	$GENDIR/$IX ${FQ2} \
	R2_mapping/${EXPID}_R2_${READLEN}.sam

for (( c=$START; c<=$END; c++ ))
do
	echo -n "$c "
	#Define how much to trim 3' for mapping
	cut3prime=$(($c*5))
	#Create tmp file for unmapped reads
	cp R2_mapping/unmapped.fastq R2_mapping/tmp.fastq
	#Conduct mapping
	bowtie -q -5 1 -3 ${cut3prime} -v 0 -m 1 -p 8 --seed=123 -S \
	--un R2_mapping/unmapped.fastq \
	$GENDIR/$IX R2_mapping/tmp.fastq \
	R2_mapping/${EXPID}_R2_$((${READLEN}-${cut3prime})).sam
	sleep 1
done


echo "      Iterative mapping of the first reads to repeat..."

cp R1_mapping/unmapped.fastq R1_mapping/tmp.fastq
bowtie -q -5 1 -v 0 -p 8 --seed=123 -S \
	--un R1_mapping/unmapped.fastq \
	$REPGENDIR/$REPIX R1_mapping/tmp.fastq \
	R1_mapping/${EXPID}_REP_R1_${READLEN}.sam

for (( c=$START; c<=$END; c++ ))
do
	echo -n "$c "
	#Define how much to trim 3' for mapping
	cut3prime=$(($c*5))
	#Create tmp file for unmapped reads
	cp R1_mapping/unmapped.fastq R1_mapping/tmp.fastq
	#Conduct mapping
	bowtie -q -5 1 -3 ${cut3prime} -v 0 -p 8 --seed=123 -S \
	--un R1_mapping/unmapped.fastq \
	$GENDIR/$IX R1_mapping/tmp.fastq \
	R1_mapping/${EXPID}_REP_R1_$((${READLEN}-${cut3prime})).sam
	sleep 1
done

echo "      Iterative mapping of the second reads to repeat..."

bowtie -q -5 1 -v 0 -p 8 --seed=123 -S \
	--un R1_mapping/unmapped.fastq \
	$REPGENDIR/$REPIX ${FQ2} \
	R2_mapping/${EXPID}_REP_R2_${READLEN}.sam

for (( c=$START; c<=$END; c++ ))
do
	echo -n "$c "
	#Define how much to trim 3' for mapping
	cut3prime=$(($c*5))
	#Create tmp file for unmapped reads
	cp R2_mapping/unmapped.fastq R2_mapping/tmp.fastq
	#Conduct mapping
	bowtie -q -5 1 -3 ${cut3prime} -v 0 -p 8 --seed=123 -S \
	--un R2_mapping/unmapped.fastq \
	$REPGENDIR/$REPIX R2_mapping/tmp.fastq \
	R2_mapping/${EXPID}_REP_R2_$((${READLEN}-${cut3prime})).sam
	sleep 1
done

#------------------------------------------------------------------------------#
#                             Convert SAM to BED files                         #
#------------------------------------------------------------------------------#

module load samtools/intel/1.3.1
module load bedtools/intel/2.26.0 

START=0

echo "      Converting first read Bam2Bed..."

for (( c=$START; c<=$END; c++ ))
do
	echo -n "$c "
	#Need the trim number for accessing different files
	cut3prime=$(($c*5))
	#Convert bam file
	samtools view -b \
	R1_mapping/${EXPID}_R1_$((${READLEN}-${cut3prime})).sam \
	> R1_mapping/${EXPID}_R1_$((${READLEN}-${cut3prime})).bam
	#Remove the header and give only the mapped reads
	samtools view -h -F 4 -b \
	R1_mapping/${EXPID}_R1_$((${READLEN}-${cut3prime})).bam \
	> R1_mapping/${EXPID}_MAPPED_R1_$((${READLEN}-${cut3prime})).bam
	#Convert to a bed file
	bedtools bamtobed -i \
	R1_mapping/${EXPID}_MAPPED_R1_$((${READLEN}-${cut3prime})).bam \
	> R1_mapping/${EXPID}_MAPPED_R1_$((${READLEN}-${cut3prime})).bed
	sleep 1
done

echo "      Converting second read Bam2Bed..."

for (( c=$START; c<=$END; c++ ))
do
	echo -n "$c "
	#Need the trim number for accessing different files
	cut3prime=$(($c*5))
	#Convert bam file
	samtools view -b \
	R2_mapping/${EXPID}_R2_$((${READLEN}-${cut3prime})).sam \
	> R2_mapping/${EXPID}_R2_$((${READLEN}-${cut3prime})).bam
	#Remove the header and give only the mapped reads
	samtools view -h -F 4 -b \
	R2_mapping/${EXPID}_R2_$((${READLEN}-${cut3prime})).bam \
	> R2_mapping/${EXPID}_MAPPED_R2_$((${READLEN}-${cut3prime})).bam
	#Convert to a bed file
	bedtools bamtobed -i \
	R2_mapping/${EXPID}_MAPPED_R2_$((${READLEN}-${cut3prime})).bam \
	> R2_mapping/${EXPID}_MAPPED_R2_$((${READLEN}-${cut3prime})).bed
	sleep 1
done

echo "      Converting first read repeat Bam2Bed..."

for (( c=$START; c<=$END; c++ ))
do
	echo -n "$c "
	#Need the trim number for accessing different files
	cut3prime=$(($c*5))
	#Convert bam file
	samtools view -b \
	R1_mapping/${EXPID}_REP_R1_$((${READLEN}-${cut3prime})).sam \
	> R1_mapping/${EXPID}_REP_R1_$((${READLEN}-${cut3prime})).bam
	#Remove the header and give only the mapped reads
	samtools view -h -F 4 -b \
	R1_mapping/${EXPID}_REP_R1_$((${READLEN}-${cut3prime})).bam \
	> R1_mapping/${EXPID}_REP_MAPPED_R1_$((${READLEN}-${cut3prime})).bam
	#Convert to a bed file
	bedtools bamtobed -i \
	R1_mapping/${EXPID}_REP_MAPPED_R1_$((${READLEN}-${cut3prime})).bam \
	> R1_mapping/${EXPID}_REP_MAPPED_R1_$((${READLEN}-${cut3prime})).bed
	sleep 1
done

echo "      Converting second read repeat Bam2Bed..."

for (( c=$START; c<=$END; c++ ))
do
	echo -n "$c "
	#Need the trim number for accessing different files
	cut3prime=$(($c*5))
	#Convert bam file
	samtools view -b \
	R2_mapping/${EXPID}_REP_R2_$((${READLEN}-${cut3prime})).sam \
	> R2_mapping/${EXPID}_REP_R2_$((${READLEN}-${cut3prime})).bam
	#Remove the header and give only the mapped reads
	samtools view -h -F 4 -b \
	R2_mapping/${EXPID}_REP_R2_$((${READLEN}-${cut3prime})).bam \
	> R2_mapping/${EXPID}_REP_MAPPED_R2_$((${READLEN}-${cut3prime})).bam
	#Convert to a bed file
	bedtools bamtobed -i \
	R2_mapping/${EXPID}_REP_MAPPED_R2_$((${READLEN}-${cut3prime})).bam \
	> R2_mapping/${EXPID}_REP_MAPPED_R2_$((${READLEN}-${cut3prime})).bed
	sleep 1
done

echo "      Catenate the myriad bed files together..."

cat R1_mapping/${EXPID}_MAPPED_R1_*.bed > ${EXPID}_MAPPED_R1.bed
cat R2_mapping/${EXPID}_MAPPED_R2_*.bed > ${EXPID}_MAPPED_R2.bed
cat R1_mapping/${EXPID}_REP_MAPPED_R1_*.bed > ${EXPID}_REP_MAPPED_R1.bed
cat R2_mapping/${EXPID}_REP_MAPPED_R2_*.bed > ${EXPID}_REP_MAPPED_R2.bed

#------------------------------------------------------------------------------#
#                               Merge files                                    #
#------------------------------------------------------------------------------#

module load perl/intel/5.24.0

echo "      Modify files to allow them to be merged..."
awk '{OFS="\t"; print $4,$1,$2,$6}' ${EXPID}_REP_MAPPED_R1.bed > MERGEREADY_${EXPID}_REP_R1.bed
awk '{OFS="\t"; print $4,$1,$2,$6}' ${EXPID}_REP_MAPPED_R2.bed > MERGEREADY_${EXPID}_REP_R2.bed
awk '{OFS="\t"; print $4,$1,$2,$6}' ${EXPID}_MAPPED_R1.bed > MERGEREADY_${EXPID}_REP_R1.bed
awk '{OFS="\t"; print $4,$1,$2,$6}' ${EXPID}_MAPPED_R2.bed > MERGEREADY_${EXPID}_REP_R2.bed

echo "      Merge files using perl script. Uses the ID for each read to find its respective pair...."
perl /home/mrp420/worms/scripts/bed_pe_merge.pl MERGEREADY_${EXPID}_REP_R1.bed MERGEREADY_${EXPID}_REP_R2.bed ${EXPID}_REP_intra_MERGED.bed
perl /home/mrp420/worms/scripts/bed_pe_merge.pl MERGEREADY_${EXPID}_R1.bed MERGEREADY_${EXPID}_R2.bed ${EXPID}_intra_MERGED.bed
perl /home/mrp420/worms/scripts/bed_pe_merge.pl MERGEREADY_${EXPID}_REP_R1.bed MERGEREADY_${EXPID}_R2.bed ${EXPID}_inter1_MERGED.bed
perl /home/mrp420/worms/scripts/bed_pe_merge.pl MERGEREADY_${EXPID}_R1.bed MERGEREADY_${EXPID}_REP_R2.bed ${EXPID}_inter2_MERGED.bed
cat ${EXPID}_inter1_MERGED.bed ${EXPID}_inter2_MERGED.bed > ${EXPID}_inter_MERGED.bed
 
#------------------------------------------------------------------------------#
#                Filter reads for intra versus inter                           #
#------------------------------------------------------------------------------# 
 
 echo "      Interactions are filtered to either inter-chromosomal or intra-chromosomal ...."
awk '{OFS="\t"} {if ($2 == $5) print $2,$3,$5,$6,$4,$7,$1}' ${EXPID}_intra_MERGED.bed > ${EXPID}_intra.bed

awk '{OFS="\t"} {if ($2 < $4) print $1,$2,$4,$5,$6,$7;
else
print $1,$4,$2,$5,$6,$7;
}' ${EXPID}_intra.bed > ${EXPID}_intra_sorted.bed

awk '{OFS="\t"} {if ($2 != $5) print $2,$3,$5,$6,$4,$7,$1}' ${EXPID}_intra_MERGED.bed > ${EXPID}_inter.bed

#------------------------------------------------------------------------------#
#                                  Clean up                                    #
#------------------------------------------------------------------------------#
echo ">>>>> Delete unnecessary files..."
#rm R1_mapping/*
#rm R2_mapping/*

#------------------------------------------------------------------------------#
ELAPSEDTIME=$(elapsed_time $STARTTIME)
echo
echo "-----"
echo "-----"
echo "Completed pipeline in $ELAPSEDTIME"
echo \
"------------------------------------------------------------------------------"





