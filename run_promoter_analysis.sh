#!/bin/bash

# extract promoter seq
Rscript extract_proseq.R $1 $2 $3

# build hisat2 index for mapping
hisat2-build -p 8 \
$1 \
./hisat2_index/genome_hisat2_index

# mapping the promoter seqs to the genome
hisat2 -p 8 \
--new-summary \
-x ./hisat2_index/genome_hisat2_index \
-f \
-U proseq.full.includes.rev.comp.cds.fa \
-S proseq_hisat2.sam

# sort and change into bam file
samtools sort -@ 10 -o proseq_hisat2_sorted.bam proseq_hisat2.sam

# create index for visualization eg IGV
samtools index -@ 5 proseq_hisat2_sorted.bam

# remove sam file as it is unnecessary 
rm proseq_hisat2.sam

fasta-get-markov \
-m 0 \
-dna \
-norc \
proseq.full.includes.rev.comp.cds.fa | sed -n 3,7p > markov_0_order.txt

fimo --o out_fimo \
--norc \
--bfile markov_0_order.txt \
PlantTFDB_TF_binding_motifs_from_experiments \
proseq.full.includes.rev.comp.cds.fa

Rscript change_fimo.tsv_to_gene_motif_matrix.R ./out_fimo/fimo.tsv



