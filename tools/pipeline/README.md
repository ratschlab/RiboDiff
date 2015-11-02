Data preprocessing
----------

###DESCRIPTION

Scripts for preprocessing ribosome footprint and RNA-Seq data. We use STAR (Dobin A. Bioinformatics. 2013) 
to align reads to references.

###Pipeline

* identify_rRNA
  * Build your local rRNA reference:

    ```
	generate_rRNA_ref.sh <Path To STAR Excutable>
    ```

  * Align reads to rRNA reference:

	```
	align_rRNA.sh <STAR Excutable> <rRNA Reference> <FASTQ> <Output Dir>
	```

  * Identify rRNA reads:

	```
	python get_rRNA_reads.py <Input Bam> <Output File>
	```

	The rRNA read IDs will be stored in the output file.

* align_reads
  * Align RNA-Seq reads:

    ```
	align_RNA.sh <STAR Excutable> <Genome Reference> <FASTQ 1> <FASTQ 2> <Output Dir>
	```

  * Align ribosome footprint reads:

	```
	align_RF.sh <STAR Excutable> <Genome Reference> <FASTQ> <Output Dir>
	```

* filter_reads
  * Filter rRNA read alignments from RNA-Seq Bam:

	```
	python filtering_RNA.py <Input Bam> <rRNA IDs File>
	```

  * Filter rRNA read alignments from ribosome footprint Bam:

	```
	python filtering_RF.py <Input Bam> <rRNA IDs File>
	```

	Note: the filtered Bam file will be generated in the same directory as the input Bam file.

* count_reads
  * Count reads for genes given the GTF file containing the gene annotation:

	```
	python count_expression.pysam.py -A <Bam file> -a <GTF File> -o <Output Count File> -v -b
	```

	Note: count_expression.pysam.py requires pysam module installed.
