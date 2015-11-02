#!/bin/bash

set -e

ARGS=1
if [ $# -ne $ARGS ]
then
	echo Usage: `basename $0` \<STAR excutable\>
	exit
fi

STAR=$1

$STAR --runMode genomeGenerate --genomeDir ./rRNA_ref_human/ --genomeFastaFiles ./rRNA_ref_human/biomart_silva_human_rRNA.fasta --runThreadN 1
