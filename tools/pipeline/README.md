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
	align_RNA.sh <STAR excutable> <genome reference> <FASTQ 1> <FASTQ 2> <output Dir>
	```
