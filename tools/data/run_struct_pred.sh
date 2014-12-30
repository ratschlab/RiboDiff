#/bin/bash

set -e

seqPath=/cbio/grlab/projects/RibosomeFootprint/GQ/Seq
seqFiles=`find $seqPath -type f -name 'AllProtGenes.all.fasta'`

for IN_FILE in $seqFiles
do
	OUT_FILE_QUAD=`echo $IN_FILE | sed -e 's/Seq/Struc/' -e 's/fasta/struc/'`

	#echo "Predicting secondary structure"
	#(/cbio/grlab/share/software/ViennaRNA/ViennaRNA-2.1.2-bin/bin/RNAfold --noPS < $IN_FILE) > $OUT_FILE_STRUCT

	#echo "Predicting base pair probabilities"
	#(/cbio/grlab/share/software/ViennaRNA/ViennaRNA-2.1.2-bin/bin/RNAplfold < $IN_FILE) > $OUT_FILE_PROB

	echo "Predicting secondary structure with G quadruplexes"
	(/cbio/grlab/share/software/ViennaRNA/ViennaRNA-2.1.2-bin/bin/RNAfold --noPS -g < $IN_FILE) > $OUT_FILE_QUAD

	#cd /cbio/grlab/projects/GuidoWendel/SecStruct
	#/cbio/grlab/share/software/ViennaRNA/ViennaRNA-2.1.2-bin/bin/RNAplot < $OUT_FILE_QUAD

	echo "$OUT_FILE_QUAD: Done!"
done
