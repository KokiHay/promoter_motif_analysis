# **Promoter motif analysis**

This is set of tools to perform promoter motif analysis from a genome (fasta file), a structural annotation (gff file), and a motif database, provided by the user. By running this program, the user will obtain a matrix describing which transcription factor binding motifs (TFBMs) are found in each of the promoter regions of every gene in the genome. If your genome has N different genes and the TFBM database has M motifs, your result will come out as a matrix of N rows and M columns.

To use the scripts first clone this repository to any deirectory.
~~~bash
git clone https://github.com/KokiHay/promoter_motif_analysis.git
~~~

The following command will extract promoter sequences from the user provided genome sequence and the structural annotation file.

~~~bash
Rscript ./promoter_motif_analysis/extract_proseq.R [geneome.fa] [annotation.gff] [promoter_length]
~~~

**Arguments**  
genome.fa: genome sequunce in fasta format  
annotation.gff: structrual annotaion file  
promoter/length: integer defining promoter_length
