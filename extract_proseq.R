rm(list = ls(all = TRUE))             #remove all objects

options(repos="https://cran.ism.ac.jp/") ## use japanes CRAN mirror
required_packages = c("readxl","stringr","seqinr")
# check if the required packages are installed
# if not install the packeges
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) {install.packages(new_packages,quiet = T)}

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
BiocManager::install("Biostrings")
}

# load the packages
library(readxl,quietly = T)
library(stringr,quietly = T)
library(Biostrings,quietly = T)
library(seqinr,quietly = T)

files = commandArgs(trailingOnly = TRUE)  ## the user should provide 2 files
                                          ## the order should be .fa then .gff3

######  enter necessary information  ##############
in_f1 = files[1]       #input genome seq (fasta format)
in_f2 = files[2]      #input gff file (gff3)

#### prepare the genome ####
gen_seq = readDNAStringSet(in_f1, format="fasta")    # read text file as fasta
gen_seq = as.data.frame(gen_seq)                     # convert to dataframe
gen_seq = as.matrix(gen_seq)                         # convert into matrix

# parameters to use #
param1 = rownames(gen_seq) ##caution!##this doesn't necessarilly correspond to chromosome number
param2 = param1  # enter from which chromosome you want the promoter sequece
pro_len = 2000                                                  #define promoter length
out_f = paste("pro_seq_f_cds",Sys.Date(),".fa",sep = "_")                 #forward strand output file name
out_r = paste("pro_seq_r_cds",Sys.Date(),".fa",sep = "_")                 #reverse strand output file name
out_r_revcomp = paste("pro_seq_r_rev_comp_cds",Sys.Date(),sep = "_")  #reverse complement of reverse strand output file name
out_full = "proseq.full.includes.rev.comp.cds.fa" #output file name for total promoter sequence including reverse strand as reverse complement

#### preparation (1) ####
gff_file = as.matrix(read.table(file = in_f2, header = F))   # read gff file
gff_CDS = gff_file[which(gff_file[,3] == "CDS"),] # extract rows with "CDS"
gff_CDS_uniqueid = gff_CDS[!duplicated(gff_CDS[,9]),] # take the first CDS of all genes as genes are made of multiple CDSs

gff_CDS_uniqueid_f = gff_CDS_uniqueid[which(gff_CDS_uniqueid[,7] == "+"), ] ## separate "+" 
gff_CDS_uniqueid_r = gff_CDS_uniqueid[which(gff_CDS_uniqueid[,7] == "-"), ] # and "-"

gff_f = gff_CDS_uniqueid_f[!duplicated(gff_CDS_uniqueid_f[,c(1,4)]),] ## if there are sequences starting from the 
gff_r = gff_CDS_uniqueid_r[!duplicated(gff_CDS_uniqueid_r[,c(1,5)]),] # same location remove them 


######################################################################################
##### extract forward strand promoter sequences ####
######################################################################################
pro_seq_f = apply(gff_f,MARGIN = 1, function(x){
    if ( as.numeric(x[4]) < 2000 ) {                    # when the starting location of the gene is smaller than 2000
      pro_seq = str_sub(gen_seq[rownames(gen_seq) == x[1],], start = 1, end = as.numeric(x[4])-1)  # extract from the first nucleotide of the chromosome
    }else{
      pro_seq = str_sub(gen_seq[rownames(gen_seq) == x[1],], start = as.numeric(x[4]) - pro_len , end = as.numeric(x[4]) - 1)  # extract 2000bp upstram of the gene
    }
  return(pro_seq)
    }    
)
names(pro_seq_f) = gff_f[,9]  # name the promoter sequences
### promoter sequences are stored in pro_seq_f 
#######################################################################################
print("extraction of forward strand completed")

#######################################################################################
#### extract reverse strand promoter sequences #####
#######################################################################################
pro_seq_r = apply(gff_r, MARGIN = 1, function(x){
  pro_seq = str_sub(gen_seq[rownames(gen_seq) == x[1],], start = as.numeric(x[5]) + 1, end = as.numeric(x[5]) + pro_len)  # Even if the promoter length on the 3prime side of the gene is less than 2000, R will just extract untill the last nucleotide
})
names(pro_seq_r) = gff_r[,9] # name the promoter sequences
### promoter sequences are stored in pro_seq_r
#######################################################################################
print("extraction of reverse strand completed")

write.fasta(sequences = as.list(pro_seq_f), name = names(pro_seq_f), file.out = out_f)     # output forward strand promoter sequnces as fasta file
write.fasta(sequences = as.list(pro_seq_r), name = names(pro_seq_r), file.out = out_r)     # output reverse strand promoter sequnces as fasta file

pro_seq_r_dna_object = readDNAStringSet(filepath = out_r, format = "fasta")  #take in reverse strand output file as DNAStringSet
pro_seq_r_reversed = reverseComplement(pro_seq_r_dna_object)           # get reverse complement

pro_seq_f_dna_object = readDNAStringSet(filepath = out_f, format = "fasta")   # read the output file of forward strand
proseq_full = c(pro_seq_f_dna_object, pro_seq_r_reversed)  # bind the two vectors
writeXStringSet(x = proseq_full,filepath = out_full,format = "fasta")  # output as fasta file