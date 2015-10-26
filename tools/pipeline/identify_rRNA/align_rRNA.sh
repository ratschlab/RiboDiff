#!/bin/bash

set -e

ARGS=4
if [ $# -ne $ARGS ]
then
	echo Usage: `basename $0` \<STAR excutable\> \<rRNA reference\> \<FASTQ\> \<output Dir\>
	exit
fi

STAR=$1
rRNAdb=$2
fastq=$3
outputDir=$4

sample=`basename $fastq | sed 's/.fastq.gz//'`

cd $outputDir

if [ ! -f $outputDir/$sample.rRNAdb.Log.final.out ]
then

	$STAR --genomeDir $rRNAdb --readFilesIn $fastq --readFilesCommand zcat --runThreadN 2 --clip3pAdapterSeq CTGTAGGCAC --clip3pAdapterMMp 0.1 --outFilterMultimapNmax 1000000 --outFilterMultimapScoreRange 0 --outFilterMismatchNmax 2 --alignIntronMax 9998 --outFilterMatchNmin 15 --outFileNamePrefix $sample. --seedSearchStartLmax 15 --outSJfilterOverhangMin 28 6 6 6 --genomeLoad NoSharedMemory

	cat $sample.Aligned.out.sam | samtools view -Sb - > $sample.rRNAdb.bam 

	rm $sample.Aligned.out.sam
	rm $sample.Log.progress.out
	rm $sample.SJ.out.tab
	rm $sample.Log.out

	mv $sample.Log.final.out $sample.rRNAdb.Log.final.out

else
	echo $sample.rRNAdb.bam already exists!
fi
