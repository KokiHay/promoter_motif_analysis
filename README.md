# **Promoter motif analysis**
**Koki Hayashi**

This is a tool to perform promoter motif analysis from a genome (fasta file), a structural annotation (gff file), and a motif database, provided by the user. By running this program, the user will obtain a matrix describing which transcription factor binding motifs (TFBMs) are found in each of the promoter regions of every gene in the genome. If your genome has N different genes and the TFBM database has M motifs, your result will come out as a matrix of N rows and M columns.

~~~bash
Rscript extract_proseq.R [geneome.fa] [annotation.gff]
~~~
