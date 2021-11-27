####The sequencing data was obtained from two different runs and processed separately using standard dada2 pipeline. The sequence tables (output from dada2) were stored as two rds files and imported here. 
library(dada2)
library(phangorn)
library(ggplot2)
library(BiocStyle)
library(devtools)
library(gridExtra)
library(knitr)
library(phyloseq)
library(DECIPHER)
library(dplyr)
library(microbiome)
library(vegan)
library(viridis)
library(RColorBrewer)
library(DESeq2)
library(decontam)
library(DAtest)
library(samr)
library(rcompanion)
library(dplyr) 
library(ggplot2)
library(ggpubr)
library(rstatix)
library(ComplexHeatmap)
library(zCompositions)
library(reshape)
library(ZIBR)
library(lme4)
library(nlme)
library(reshape)
library(fantaxtic)
library(FSA)
library(plotly)
library(tidyverse)
library(splinectomeR)
library(tibble)
library(reshape2)
library(tidyr)
setwd("/Users/data_analysis/")
set.seed(555)

####Merging the two sequence tables generated from DADA2

st1 <- readRDS("/Users/seqtab_batch1.rds")
st2 <- readRDS("/Users/seqtab_batch2.rds")

st.all <- mergeSequenceTables(st1, st2)

#Check the merged sequence table
dim(st.all)
table(nchar(getSequences(st.all)))


##Chimera removal
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
#Checking the dimension of the chimera free sequence table
dim(seqtab)
asvtable<- table(nchar(getSequences(seqtab)))
asvtable

##Taxonomic classification using HIT-DB
ref_fasta_hit <- "hitdb_v1.00.fa"
taxtab_hit <- assignTaxonomy(seqtab, refFasta = ref_fasta_hit, multithread = TRUE, verbose = TRUE)
colnames(taxtab_hit) <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

#Importing the metadata file and merging all the separate objects created above as a single phyloseq object
meta_table<-read.csv("metadata.csv", row.names=1, check.names=FALSE)
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
                    sample_data(meta_table), 
                    tax_table(taxtab_hit))
ps

####Decontaminatiion steps using decontam package
#defining the control samples
sample_data(pstotal_non_DNA_control)$is.neg <- sample_data(pstotal_non_DNA_control)$SampleIdentity == "Control Sample"

#determining the true contaminants based on prevalence method
contamdf.prev <- isContaminant(pstotal_non_DNA_control, method="prevalence", neg="is.neg", batch = "Run", threshold = 0.1)
#checking the number of true and false contaminants detected
table(contamdf.prev$contaminant)

#filtering the true contaminant ASVs (detected in above step) from the phyloseq object
pstotal_noncontam <- prune_taxa(!contamdf.prev$contaminant, pstotal_non_DNA_control)
pstotal_non_DNA_control
pstotal_noncontam

###Negative control samples were removed the phyloseq object and taxa not present in true samples were pruned
psnoncontam_tru <- prune_taxa(taxa_sums(pstotal_noncontam)>0, pstotal_noncontam)
psnoncontam_tru

#filtering out unwanted (non-specific) taxa and unclassified phyla at phylum level from the phyloseq object
psnoncontam_tru_filtered <- subset_taxa(psnoncontam_tru,Phylum!= "Euryarchaeota" & Phylum!="NA" & Class!="Chloroplast" & Order!="Mitochondria")
psnoncontam_tru
psnoncontam_tru_filtered

##Further low read samples were removed
#Also Removed taxa not seen more than 1 times in at least 10% of the samples. This protects against an ASV with small mean & trivially large C.V.
pstrial <- psnoncontam_tru_final
pstrial
psnoncontam_tru_final_filtered <-  filter_taxa(pstrial, function(x) sum(x > 1) > (0.10 * length(x)),  TRUE)
psnoncontam_tru_final_filtered
summarize_phyloseq(psnoncontam_tru_final_filtered)
##Almost 93% reads were retained here, while the sparsity reduced by 23.7%, all singletons were removed. Used this phyloseq object later in the analysis.

#rarefaction of the data to the minimum number of reads (n=7578)
psnoncontam_tru_final_rare <- phyloseq::rarefy_even_depth(psnoncontam_tru_final_unrarified, rngseed = 123, replace = FALSE)
psnoncontam_tru_final_rare

###Saved both unrarified and rarified version of phyloseq objects for further analysis
