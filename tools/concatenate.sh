#!/bin/bash

set -e

cntRiboFile=$1
repRiboCondA=$2
repRiboCondB=$3
cntRnaFile=$4
repRnaCondA=$5
repRnaCondB=$6
outputName=$7

sed 's/\tCondition/\tRbCondition/g' $cntRiboFile > Tmp.txt; mv Tmp.txt $cntRiboFile
sed 's/\tCondition/\tRnaCondition/g' $cntRnaFile > Tmp.txt; mv Tmp.txt $cntRnaFile

seqNumRibo=`echo $(seq -s, 1 $(($repRiboCondA+$repRiboCondB+1))) | sed 's/,$//'`
seqNumRna=`echo $(seq -s, 2 $(($repRnaCondA+$repRnaCondB+1))) | sed 's/,$//'`
paste <(cut -f $seqNumRibo $cntRiboFile) <(cut -f $seqNumRna $cntRnaFile) > $outputName
