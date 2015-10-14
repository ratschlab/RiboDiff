#!/bin/bash

set -e

STAR=$1
genome=$2
fastq1=$3
fastq2=$4
bamDir=$5
logDir=$6

sample=`basename $Fq1 | sed 's/.fastq.gz//'`

cd $bamDir

if [ ! -f $logDir/$sample.Log.final.out ]
then

	$STAR --genomeDir $genome --readFilesIn $fastq1 $fastq2 --readFilesCommand zcat --runThreadN 1 --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 0 --outFilterMismatchNmax 3 --alignMatesGapMax 10000 --alignIntronMax 500000 --outFilterMatchNmin 20 --genomeLoad NoSharedMemory --outFileNamePrefix $sample. --seedSearchStartLmax 30

	cat $sample.Aligned.out.sam | samtools view -Sb - > $sample.bam

	rm $sample.Aligned.out.sam
	rm $sample.Log.progress.out
	rm $sample.SJ.out.tab
	rm $sample.Log.out
	mv $sample.Log.final.out $logDir

	echo $sample: Done!

else
	echo $sample.bam already exists!
fi
