######Alpha and beta diversity analysis
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

###Alpha diversity
psnoncontam_tru_final_unrarified
ps_total_div_all <- estimate_richness(psnoncontam_tru_final_unrarified)
write.csv(ps_total_div_all, "alpha_div_all.csv")

###Keeping only the measures/indices that I need to plot and importing it here
alpha_time_data <- read.table("alpha_div_plot_time_data.csv", sep=",", check.names = F, header = T)
alpha_time_data

#melting the data frame to use it for box plot and statistical test
alpha_time_data_melt <- melt(alpha_time_data)
alpha_time_data_melt

#Changing the column name variable to alpha, as term variable is used internally in the code, so to avoid conflict..
colnames(alpha_time_data_melt)[3] <- "Alpha"
alpha_time_data_melt

#Rstatix package for doing a pairwise comparison
stat.test.alpha <- alpha_time_data_melt %>%
  group_by(Alpha) %>%
  wilcox_test(value ~ TimePoint) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()
stat.test.alpha

#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/rstatix")

#Adding the xy position for plotting p-values on the box plot
stat.test.alpha <- stat.test.alpha %>% add_xy_position(x = "TimePoint", scales = "free")
stat.test.alpha

#Box plot using the ggboxplot function from ggpubr, using ggplot function gives an error for adding the p-values
bxp.alpha <- ggboxplot(data= alpha_time_data_melt, x = "TimePoint", y= "value", fill = "TimePoint", xlab = "Time Point", ylab= "Alpha diversity measure/index")
bxp.alpha <- facet(p = bxp.alpha,facet.by = "Alpha",scales = "free_y")

bxp.alpha + scale_fill_manual(values=c("#ffb6b9", "#fae3d9", "#bbded6", "#8ac6d3")) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="none") + scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + geom_point(position = position_dodge(0.6))

bxp.alpha + stat_pvalue_manual(stat.test.alpha, label= "p.adj.signif", tip.length=0.03, hide.ns = T, step.increase = 0.3) + scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

#multiple group comparison
kruskal.test(x= alpha_time_data$Chao1, g= alpha_time_data$TimePoint, p.adjust.method = "fdr")
kruskal.test(alpha_time_data$Shannon, alpha_time_data$TimePoint, p.adjust.method = "fdr")
kruskal.test(alpha_time_data$InvSimpson, alpha_time_data$TimePoint, p.adjust.method = "fdr")

#Multiple pairwise comparison
pairwise.wilcox.test(ps_total_div_all$Observed, sample_data(psnoncontam_tru_final_unrarified)$Description, p.adjust.method = "fdr", paired = FALSE)

pairwise.wilcox.test(ps_total_div_all$Shannon, sample_data(psnoncontam_tru_final_unrarified)$Description, p.adjust.method = "fdr", paired = FALSE)

pairwise.wilcox.test(ps_total_div_all$InvSimpson, sample_data(psnoncontam_tru_final_unrarified)$Description, p.adjust.method = "fdr", paired = FALSE)


#######For beta diversity analysis CLR transformed dataset was used using following steps

#Write the phyloseq object as a CSV table with ASVs and counts
psnoncontam_tru_final_rare
write_phyloseq(psnoncontam_tru_final_rare)

#Read the CSV file
abund_table <- read.csv("otu_table_noncontam_tru_final_rare.csv", row.names=1, check.names= FALSE) #Added header "Features" to the first clumn of ASVs
abund_table_t<-t(abund_table)
#View(abund_table_t)

#install.packages("zCompositions")
#library (zCompositions)
abund_table_r <- t(cmultRepl((abund_table_t), method="CZM", output="prop")) #prop calculates the imputed proportion is similar to "counts" and then taking proportions in the previous version; CZM: count zero multiplicative

#Check the before and after files
#View(abund_table)
#View(abund_table_r)

#Apply CLR transformation
abund_clr <- t(apply(abund_table_r, 2, function(x){log(x) - mean (log(x))}))
#View(abund_clr)

#Replacing original counts with CLR transformed values
ps_noncontam_zclr_rare <- psnoncontam_tru_final_rare
otu_table(ps_noncontam_zclr_rare) <- otu_table(abund_clr, taxa_are_rows = F)
ps_noncontam_zclr_rare
ps_noncontam_zclr_rare_pseudo <- ps_noncontam_zclr_rare
otu_table(ps_noncontam_zclr_rare_pseudo) <- otu_table(ps_noncontam_zclr_rare_pseudo) + 3
ps_noncontam_clr_plot_rare_bray <- plot_ordination(ps_noncontam_zclr_rare_pseudo, ordinate(ps_noncontam_zclr_rare_pseudo, "MDS", "bray"), color = "Description") + geom_point(size = 3)
ps_noncontam_clr_plot_rare_bray + theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_color_manual(values=c("#ffb6b9", "#fae3d9", "#bbded6", "#8ac6d3")) + theme(legend.position = "right")

#PERMANOVA

#comparing the samples at different time points
metadata5 <- as(sample_data(ps_noncontam_zclr_rare_pseudo), "data.frame")
adonis(phyloseq::distance(ps_noncontam_zclr_rare_pseudo, method = "bray")~Description, data = metadata5, p.adjust.methods= "fdr", permutations = 99999, strata = metadata5$SubjectID)###added infant IDs as strata to account for repeated measures.

#beta-dispersion test
ps_total_noncontam_dist <- phyloseq::distance(ps_noncontam_zclr_rare_pseudo, method = "bray")
sampledf5 <- data.frame(sample_data(ps_noncontam_zclr_rare_pseudo))
beta5 <- betadisper(ps_total_noncontam_dist, sampledf5$Description, bias.adjust = TRUE)
permutest(beta5, p.adjust.methods= "fdr", permutations = 99999)

#Jaccard distance on presence/absence of ASVs. Used rarified but untranformed phyloseq file here
ps_noncontam_rare_plot_jaccard <- plot_ordination(psnoncontam_tru_final_rare, ordinate(psnoncontam_tru_final_rare, "MDS", "jaccard"), color = "Description") + geom_point(size = 3)
ps_noncontam_rare_plot_jaccard + theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_color_manual(values=c("#ffb6b9", "#fae3d9", "#bbded6", "#8ac6d3")) + theme(legend.position = "right")

#PERMANOVA
#comparing the samples at different time points
metadata5 <- as(sample_data(psnoncontam_tru_final_rare), "data.frame")
adonis(phyloseq::distance(psnoncontam_tru_final_rare, method = "jaccard")~Description, data = metadata5, p.adjust.methods= "fdr", permutations = 99999, strata = metadata5$SubjectID)###added infant IDs as strata to account for repeated measures.

#beta-dispersion test
ps_total_noncontam_dist <- phyloseq::distance(psnoncontam_tru_final_rare, method = "jaccard")
sampledf5 <- data.frame(sample_data(psnoncontam_tru_final_rare))
beta5 <- betadisper(ps_total_noncontam_dist, sampledf5$Description, bias.adjust = TRUE)
permutest(beta5, p.adjust.methods= "fdr", permutations = 99999)

######################################Bacterial abundance over time##################################

##Extracted phylum and genus level OTU tables from phyloseq and added additional information

##Extracted phylum and genus level OTU tables from phyloseq and added additional information


##Aggolomerating at phyla level
phylum_data_total <- tax_glom(psnoncontam_tru_final_rare,taxrank = "Phylum")
phylum_data_total

#Relative abundances
#phylum_data_total.rel <- transform(phylum_data_total, "compositional")
#phylum_data_total.rel
#Using CLR transformation here
#Write the phyloseq object as a CSV table with ASVs and counts
phylum_data_total
write_phyloseq(phylum_data_total)

#Read the CSV file
abund_table <- read.csv("phylum_rare_clr.csv", row.names=1, check.names= FALSE) #Added header "Features" to the first clumn of ASVs
abund_table_t<-t(abund_table)
#View(abund_table_t)

#install.packages("zCompositions")
#library (zCompositions)
abund_table_r <- t(cmultRepl((abund_table_t), method="CZM", output="prop")) #prop calculates the imputed proportion is similar to "counts" and then taking proportions in the previous version; CZM: count zero multiplicative

#Check the before and after files
#View(abund_table)
#View(abund_table_r)

#Apply CLR transformation
abund_clr <- t(apply(abund_table_r, 2, function(x){log(x) - mean (log(x))}))
#View(abund_clr)

#Replacing original counts with CLR transformed values
phylum_data_total_zclr <- phylum_data_total
otu_table(phylum_data_total_zclr) <- otu_table(abund_clr, taxa_are_rows = F)
phylum_data_total_zclr

#replacing ASVs with short string
taxa_names(phylum_data_total_zclr) <- paste0("SV", seq(ntaxa(phylum_data_total_zclr)))
#View(tax_table(phylum_data_total_zclr))

#Extracted the phyloseq file in csv
write_phyloseq(phylum_data_total_zclr, type= 'all')

#Replaced string names with phyla names, transposed and added metadata for time point
##Note: Made some modifications in excel here

#importing the abundance table
phylum_time_data <- read.table("otu_table_phylum_total_data_clr_new.csv", header = TRUE, sep= ",", check.names = F)
phylum_time_data

#class(phylum_time_data)
trans_phylum_time_data <- melt(phylum_time_data)
trans_phylum_time_data

#Changing the column name from 'variable' to 'Phylum'
colnames(trans_phylum_time_data)[3] <- "Phylum"
trans_phylum_time_data
#class(trans_phylum_time_data$variable)

#Rstatix
stat.test.phyla <- trans_phylum_time_data %>%
  group_by(Phylum) %>%
  wilcox_test(value ~ TimePoint) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()
stat.test.phyla

#Adding XY position
stat.test.phyla <- stat.test.phyla %>% add_xy_position(x = "TimePoint")
stat.test.phyla

#Box plot
bxp <- ggboxplot(data= trans_phylum_time_data, x = "TimePoint", y= "value", facet.by = "Phylum", xlab = "Time Points", ylab= "Abundance (CLR)", fill="Phylum")

bxp + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="bottom") + stat_pvalue_manual(stat.test.phyla, label= "p.adj.signif", tip.length=0.03, hide.ns = T) + scale_fill_manual(values=c("#5B7E00", "#C3CFA2", "#3A5200"))

##There are no significant differences after fdr correction, so no brackets seen in the plot below


##Multiple group comparison
kruskal.test(phylum_time_data$Actinobacteria, phylum_time_data$TimePoint)
kruskal.test(phylum_time_data$Firmicutes, phylum_time_data$TimePoint)
kruskal.test(phylum_time_data$Proteobacteria, phylum_time_data$TimePoint)


#Posthoc using wilcoxon sum rank test
pairwise.wilcox.test(phylum_time_data$Actinobacteria, phylum_time_data$TimePoint, p.adjust.method = "fdr", paired = FALSE)

pairwise.wilcox.test(phylum_time_data$Firmicutes, phylum_time_data$TimePoint, p.adjust.method = "fdr", paired = FALSE)

pairwise.wilcox.test(phylum_time_data$Proteobacteria, phylum_time_data$TimePoint, p.adjust.method = "fdr", paired = FALSE)

##Genus level abundance over time
##Extracted phylum and genus level OTU tables from phyloseq and added additional information
psnoncontam_tru_final_rare
##Aggolomerating at genus level
genus_data_total <- tax_glom(psnoncontam_tru_final_rare,taxrank = "Genus")
genus_data_total

#Relative abundances
#genus_data_total.rel <- transform(genus_data_total, "compositional")
#genus_data_total.rel

#Using CLR transformed abundances here
#Write the phyloseq object as a CSV table with ASVs and counts
genus_data_total
write_phyloseq(genus_data_total)

#Read the CSV file
abund_table <- read.csv("genus_rare_clr_new.csv", row.names=1, check.names= FALSE) #Added header "Features" to the first clumn of ASVs
abund_table_t<-t(abund_table)
#View(abund_table_t)

#install.packages("zCompositions")
#library (zCompositions)
abund_table_r <- t(cmultRepl((abund_table_t), method="CZM", output="prop")) #prop calculates the imputed proportion is similar to "counts" and then taking proportions in the previous version; CZM: count zero multiplicative

#Check the before and after files
#View(abund_table)
#View(abund_table_r)

#Apply CLR transformation
abund_clr <- t(apply(abund_table_r, 2, function(x){log(x) - mean (log(x))}))
#View(abund_clr)

#Replacing original counts with CLR transformed values
genus_data_total_zclr <- genus_data_total
otu_table(genus_data_total_zclr) <- otu_table(abund_clr, taxa_are_rows = F)
genus_data_total_zclr

#replacing ASVs with short string
taxa_names(genus_data_total_zclr) <- paste0("SV", seq(ntaxa(genus_data_total_zclr)))
#View(tax_table(genus_data_total_zclr))

#Extracting short strings and their respective taxonomy
genus_list<- tax_table(genus_data_total_zclr)@.Data
write.csv(genus_list, "genus_total_list.csv")

#extracting as csv and replacing string names with phyla names, transposed and added metadata
write_phyloseq(genus_data_total_zclr, type= 'all')


genus_time_data <- read.csv("otu_table_genus_total_data_clr_new.csv", header = TRUE, sep= ",")
genus_time_data
#library(reshape)
#class(genus_time_data)
trans_genus_time_data <- melt.data.frame(genus_time_data)
trans_genus_time_data


##multiple group comparison using Kruskal wallis
kruskal.test(genus_time_data$Staphylococcus, genus_time_data$TimePoint)
kruskal.test(genus_time_data$Enterobacteriaceae..F., genus_time_data$TimePoint)
kruskal.test(genus_time_data$Erwinia, genus_time_data$TimePoint)
kruskal.test(genus_time_data$Clostridium, genus_time_data$TimePoint)
kruskal.test(genus_time_data$Enterococcus, genus_time_data$TimePoint)
kruskal.test(genus_time_data$Bifidobacterium, genus_time_data$TimePoint)
kruskal.test(genus_time_data$Propionibacterium, genus_time_data$TimePoint)
kruskal.test(genus_time_data$Veillonella, genus_time_data$TimePoint)
kruskal.test(genus_time_data$Clostridiales.XI..F., genus_time_data$TimePoint)
kruskal.test(genus_time_data$Peptostreptococcaceae..F., genus_time_data$TimePoint)
kruskal.test(genus_time_data$Pasteurellaceae..F., genus_time_data$TimePoint)
kruskal.test(genus_time_data$Lactobacillus, genus_time_data$TimePoint)
kruskal.test(genus_time_data$Bacillus, genus_time_data$TimePoint)
kruskal.test(genus_time_data$Streptococcus, genus_time_data$TimePoint)
kruskal.test(genus_time_data$Dialister, genus_time_data$TimePoint)
kruskal.test(genus_time_data$Peptoniphilus, genus_time_data$TimePoint)
kruskal.test(genus_time_data$Lachnoclostridium, genus_time_data$TimePoint)
kruskal.test(genus_time_data$Corynebacterium, genus_time_data$TimePoint)
kruskal.test(genus_time_data$Acinetobacter, genus_time_data$TimePoint)
kruskal.test(genus_time_data$Anaerococcus, genus_time_data$TimePoint)
kruskal.test(genus_time_data$Varibaculum, genus_time_data$TimePoint)
kruskal.test(genus_time_data$Comamonas, genus_time_data$TimePoint)
kruskal.test(genus_time_data$Gemella, genus_time_data$TimePoint)
kruskal.test(genus_time_data$Lactobacillales..O., genus_time_data$TimePoint)


#####Pairwise comparison for time points
pairwise.wilcox.test(genus_time_data$Staphylococcus, genus_time_data$TimePoint, p.adjust.method = "fdr", paired = FALSE)
pairwise.wilcox.test(genus_time_data$Enterobacteriaceae..F., genus_time_data$TimePoint, p.adjust.method = "fdr", paired = FALSE)
pairwise.wilcox.test(genus_time_data$Erwinia, genus_time_data$TimePoint, p.adjust.method = "fdr", paired = FALSE)
pairwise.wilcox.test(genus_time_data$Clostridium, genus_time_data$TimePoint, p.adjust.method = "fdr", paired = FALSE)
pairwise.wilcox.test(genus_time_data$Enterococcus, genus_time_data$TimePoint, p.adjust.method = "fdr", paired = FALSE)
pairwise.wilcox.test(genus_time_data$Enterococcus, genus_time_data$TimePoint, p.adjust.method = "fdr", paired = FALSE)
pairwise.wilcox.test(genus_time_data$Bifidobacterium, genus_time_data$TimePoint, p.adjust.method = "fdr", paired = FALSE)
pairwise.wilcox.test(genus_time_data$Propionibacterium, genus_time_data$TimePoint, p.adjust.method = "fdr", paired = FALSE)
pairwise.wilcox.test(genus_time_data$Veillonella, genus_time_data$TimePoint, p.adjust.method = "fdr", paired = FALSE)
pairwise.wilcox.test(genus_time_data$Clostridiales.XI..F., genus_time_data$TimePoint, p.adjust.method = "fdr", paired = FALSE)
pairwise.wilcox.test(genus_time_data$Peptostreptococcaceae..F., genus_time_data$TimePoint, p.adjust.method = "fdr", paired = FALSE)
pairwise.wilcox.test(genus_time_data$Pasteurellaceae..F., genus_time_data$TimePoint, p.adjust.method = "fdr", paired = FALSE)
pairwise.wilcox.test(genus_time_data$Lactobacillus, genus_time_data$TimePoint, p.adjust.method = "fdr", paired = FALSE)
pairwise.wilcox.test(genus_time_data$Bacillus, genus_time_data$TimePoint, p.adjust.method = "fdr", paired = FALSE)
pairwise.wilcox.test(genus_time_data$Streptococcus, genus_time_data$TimePoint, p.adjust.method = "fdr", paired = FALSE)
pairwise.wilcox.test(genus_time_data$Dialister, genus_time_data$TimePoint, p.adjust.method = "fdr", paired = FALSE)
pairwise.wilcox.test(genus_time_data$Peptoniphilus, genus_time_data$TimePoint, p.adjust.method = "fdr", paired = FALSE)
pairwise.wilcox.test(genus_time_data$Lachnoclostridium, genus_time_data$TimePoint, p.adjust.method = "fdr", paired = FALSE)
pairwise.wilcox.test(genus_time_data$Corynebacterium, genus_time_data$TimePoint, p.adjust.method = "fdr", paired = FALSE)
pairwise.wilcox.test(genus_time_data$Acinetobacter, genus_time_data$TimePoint, p.adjust.method = "fdr", paired = FALSE)
pairwise.wilcox.test(genus_time_data$Anaerococcus, genus_time_data$TimePoint, p.adjust.method = "fdr", paired = FALSE)
pairwise.wilcox.test(genus_time_data$Varibaculum, genus_time_data$TimePoint, p.adjust.method = "fdr", paired = FALSE)
pairwise.wilcox.test(genus_time_data$Comamonas, genus_time_data$TimePoint, p.adjust.method = "fdr", paired = FALSE)
pairwise.wilcox.test(genus_time_data$Gemella, genus_time_data$TimePoint, p.adjust.method = "fdr", paired = FALSE)
pairwise.wilcox.test(genus_time_data$Lactobacillales..O., genus_time_data$TimePoint, p.adjust.method = "fdr", paired = FALSE)

#Replotting box plots with genera

#importing the abundance table
genus_time_data <- read.table("otu_table_genus_total_data_clr_new.csv", header = TRUE, sep= ",", check.names = F)
genus_time_data

#Melting the data frame to use it for boxplot and statistics
trans_genus_time_data <- melt(genus_time_data)
trans_genus_time_data

#Changing the column name from 'variable' to 'Genus'
colnames(trans_genus_time_data)[3] <- "Genus"
trans_genus_time_data

#Rstatix
stat.test.genus <- trans_genus_time_data %>%
  group_by(Genus) %>%
  wilcox_test(value ~ TimePoint) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()
stat.test.genus

#Adding XY position
stat.test.genus <- stat.test.genus %>% add_xy_position(x = "TimePoint")
stat.test.genus

#Box plot
bxp.genus <- ggboxplot(data= trans_genus_time_data, x = "TimePoint", y= "value", facet.by = "Genus", xlab = "Time Points", ylab="Abundance (CLR)", fill = "Genus", lwd=0.3, fatten=0.3) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.3)))

#scale_fill_manual(values=c("#FDD349", "#80D7E2", "#D8AE92", "#F6CDC4", "#C36200", "#B8EE8A", "#004B8D", "#C89F00", "#537345", "#FF4739", "#CBA96E", "#FFB55C", "#54565A", "#4288A4", "#9E0000", "#E80076", "#C3CB3C", "#6600CC", "#773F05"))  

bxp.genus

#modifying the plot to use individual y-scale for each facet
bxp.genus <- facet(p = bxp.genus,facet.by = "Genus",scales = "free_y") + theme(legend.position = "bottom") + scale_fill_manual(values = c("#F7D363", "#20899F", "#94D5E0", "#D2AF96", "#F0CEC6", "#B56602",  "#C3EB95", "#044B88", "#935272", "#0DA703", "#653723",  "#C29F00", "#5B724B", "#F35733", "#C6AA76", "#FFB86C", "#53565A", "#EAF0F3", "#B00D33", "#D52F75", "#C3CB3C", "#6600CC", "#773F05", "#DDA0DD"))
bxp.genus

#To add significance for pairwise comparisons
#bxp.genus + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="none") + stat_pvalue_manual(stat.test.genus, label="p.adj.signif", tip.length=0.02, hide.ns = T)


##################Section 2: Gestational age and changes in microbial diversity######################


#Using AGA as a continous gradient to colour
#For CLR and bray curtis

ps_noncontam_zclr_rare
ps_noncontam_zclr_rare_pseudo <- ps_noncontam_zclr_rare
otu_table(ps_noncontam_zclr_rare_pseudo) <- otu_table(ps_noncontam_zclr_rare_pseudo) + 3
ps_noncontam_clr_plot_rare_bray <- plot_ordination(ps_noncontam_zclr_rare_pseudo, ordinate(ps_noncontam_zclr_rare_pseudo, "MDS", "bray"), color = "GestationalAge") + geom_point(size = 3)

ps_noncontam_clr_plot_rare_bray + theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_color_gradient(low= "#FE3700", high = "#FAE4E3") + theme(legend.position = "right")

#comparing the samples at different time points
metadata5 <- as(sample_data(ps_noncontam_zclr_rare_pseudo), "data.frame")
adonis(phyloseq::distance(ps_noncontam_zclr_rare_pseudo, method = "bray")~GestationalAge, data = metadata5, p.adjust.methods= "fdr", permutations = 99999)###do a pairwise adonis here

#beta-dispersion test
#ps_total_noncontam_dist <- phyloseq::distance(ps_noncontam_zclr_rare_pseudo, method = "bray")
#sampledf5 <- data.frame(sample_data(ps_noncontam_zclr_rare_pseudo))
#beta5 <- betadisper(ps_total_noncontam_dist, sampledf5$GestationalAge, bias.adjust = TRUE)
#permutest(beta5, p.adjust.methods= "fdr", permutations = 99999)


#Jaccard distance on presence/absence of ASVs

ps_noncontam_rare_plot_jaccard <- plot_ordination(psnoncontam_tru_final_rare, ordinate(psnoncontam_tru_final_rare, "MDS", "jaccard"), color = "GestationalAge") + geom_point(size = 3)

ps_noncontam_rare_plot_jaccard + theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_color_gradient(low= "#FE3700", high = "#FAE4E3") + theme(legend.position = "right")

#comparing the samples at different time points
metadata5 <- as(sample_data(psnoncontam_tru_final_rare), "data.frame")
adonis(phyloseq::distance(psnoncontam_tru_final_rare, method = "jaccard")~GestationalAge, data = metadata5, p.adjust.methods= "fdr", permutations = 99999)###do a pairwise adonis here

#beta-dispersion test
#ps_total_noncontam_dist <- phyloseq::distance(psnoncontam_tru_final_rare, method = "jaccard")
#sampledf5 <- data.frame(sample_data(psnoncontam_tru_final_rare))
#beta5 <- betadisper(ps_total_noncontam_dist, sampledf5$GestationalAge, bias.adjust = TRUE)
#permutest(beta5, p.adjust.methods= "fdr", permutations = 99999)

###########################Section 3: RCT groups and microbial composition##########################

####Composition of each group at phylum level

psnoncontam_tru_final_rare
#subset control group samples
ps_control_final <- subset_samples(psnoncontam_tru_final_rare,Group%in%c("Control_group"))
ps_control_final
ps_control_final <- prune_taxa(taxa_sums(ps_control_final)>0, ps_control_final)
ps_control_final


####further separating the samples from control group at different time points

#before fortification samples
ps_control_bf <- subset_samples(ps_control_final,Description%in%c("T1"))
ps_control_bf
ps_control_bf <- prune_taxa(taxa_sums(ps_control_bf)>0, ps_control_bf)
ps_control_bf
#during fortification samples (after 7 days)
ps_control_df <- subset_samples(ps_control_final,Description%in%c("T2"))
ps_control_df
ps_control_df <- prune_taxa(taxa_sums(ps_control_df)>0, ps_control_df)
ps_control_df
#at the end of fortification samples (week 33)
ps_control_ef <- subset_samples(ps_control_final,Description%in%c("T3"))
ps_control_ef
ps_control_ef <- prune_taxa(taxa_sums(ps_control_ef)>0, ps_control_ef)
ps_control_ef
#follow up samples (week 35)
ps_control_ff <- subset_samples(ps_control_final,Description%in%c("T4"))
ps_control_ff
ps_control_ff <- prune_taxa(taxa_sums(ps_control_ff)>0, ps_control_ff)
ps_control_ff


####plotting composition of control samples at each time point

ps_control_bf.rel <- microbiome::transform(ps_control_bf, "compositional")
plot_t1_c_phylum <- plot_bar(ps_control_bf.rel, fill = "Phylum") + geom_bar(aes( fill=Phylum), stat="identity", position = "stack") + scale_fill_manual(values=c("#40530A","#637C33", "#C4CEA7")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90), aspect.ratio = 4/4)
plot_t1_c <- plot_t1_c_phylum + labs(x="", y= "Relative abundance")
plot_t1_c

ps_control_df.rel <- microbiome::transform(ps_control_df, "compositional")
plot_t2_c_phylum <- plot_bar(ps_control_df.rel, fill = "Phylum") + geom_bar(aes( fill=Phylum), stat="identity", position = "stack") + scale_fill_manual(values=c("#40530A","#637C33", "#C4CEA7")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90), aspect.ratio = 4/4)
plot_t2_c <- plot_t2_c_phylum + labs(x="", y= "Relative abundance")
plot_t2_c

ps_control_ef.rel <- microbiome::transform(ps_control_ef, "compositional")
plot_t3_c_phylum <- plot_bar(ps_control_ef.rel, fill = "Phylum") + geom_bar(aes( fill=Phylum), stat="identity", position = "stack") + scale_fill_manual(values=c("#40530A","#637C33", "#C4CEA7")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90), aspect.ratio = 4/4)
plot_t3_c <- plot_t3_c_phylum + labs(x="", y= "Relative abundance")
plot_t3_c


ps_control_ff.rel <- microbiome::transform(ps_control_ff, "compositional")
plot_t4_c_phylum <- plot_bar(ps_control_ff.rel, fill = "Phylum") + geom_bar(aes( fill=Phylum), stat="identity", position = "stack") + scale_fill_manual(values=c("#40530A","#637C33", "#C4CEA7")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90), aspect.ratio = 4/4)
plot_t4_c <- plot_t4_c_phylum + labs(x="", y= "Relative abundance")
plot_t4_c


####subset treatment group samples and separating the samples from treatment group at different time points


#subset treatment group samples
ps_treatment_final <- subset_samples(psnoncontam_tru_final_rare,Group%in%c("Treatment_group"))
ps_treatment_final
ps_treatment_final <- prune_taxa(taxa_sums(ps_treatment_final)>0, ps_treatment_final)
ps_treatment_final

#before fortification samples
ps_treatment_bf <- subset_samples(ps_treatment_final,Description%in%c("T1"))
ps_treatment_bf
ps_treatment_bf <- prune_taxa(taxa_sums(ps_treatment_bf)>0, ps_treatment_bf)
ps_treatment_bf
#during fortification samples (after 7 days)
ps_treatment_df <- subset_samples(ps_treatment_final,Description%in%c("T2"))
ps_treatment_df
ps_treatment_df <- prune_taxa(taxa_sums(ps_treatment_df)>0, ps_treatment_df)
ps_treatment_df
#at the end of fortification samples (week 33)
ps_treatment_ef <- subset_samples(ps_treatment_final,Description%in%c("T3"))
ps_treatment_ef
ps_treatment_ef <- prune_taxa(taxa_sums(ps_treatment_ef)>0, ps_treatment_ef)
ps_treatment_ef
#follow up samples (week 35)
ps_treatment_ff <- subset_samples(ps_treatment_final,Description%in%c("T4"))
ps_treatment_ff
ps_treatment_ff <- prune_taxa(taxa_sums(ps_treatment_ff)>0, ps_treatment_ff)
ps_treatment_ff


####plotting composition of treatment samples at each time point

ps_treatment_bf.rel <- microbiome::transform(ps_treatment_bf, "compositional")
plot_t1_t_phylum <- plot_bar(ps_treatment_bf.rel, fill = "Phylum") + geom_bar(aes( fill=Phylum), stat="identity", position = "stack") + scale_fill_manual(values=c("#40530A","#637C33", "#C4CEA7")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90), aspect.ratio = 4/4)
plot_t1_t <- plot_t1_t_phylum + labs(x="", y= "Relative abundance")
plot_t1_t

ps_treatment_df.rel <- microbiome::transform(ps_treatment_df, "compositional")
plot_t2_t_phylum <- plot_bar(ps_treatment_df.rel, fill = "Phylum") + geom_bar(aes( fill=Phylum), stat="identity", position = "stack") + scale_fill_manual(values=c("#40530A","#637C33", "#C4CEA7")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90), aspect.ratio = 4/4)
plot_t2_t <- plot_t2_t_phylum + labs(x="", y= "Relative abundance")
plot_t2_t

ps_treatment_ef.rel <- microbiome::transform(ps_treatment_ef, "compositional")
plot_t3_t_phylum <- plot_bar(ps_treatment_ef.rel, fill = "Phylum") + geom_bar(aes( fill=Phylum), stat="identity", position = "stack") + scale_fill_manual(values=c("#40530A","#637C33", "#C4CEA7")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90), aspect.ratio = 4/4)
plot_t3_t <- plot_t3_t_phylum + labs(x="", y= "Relative abundance")
plot_t3_t

ps_treatment_ff.rel <- microbiome::transform(ps_treatment_ff, "compositional")
plot_t4_t_phylum <- plot_bar(ps_treatment_ff.rel, fill = "Phylum") + geom_bar(aes( fill=Phylum), stat="identity", position = "stack") + scale_fill_manual(values=c("#40530A","#637C33", "#C4CEA7")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90), aspect.ratio = 4/4)
plot_t4_t <- plot_t4_t_phylum + labs(x="", y= "Relative abundance")
plot_t4_t


####Trying out the combined mean abundance bar chart RCT group wise at phylum level
psnoncontam_tru_final_rare


#Seperating out samples based on time points
ps_bf_total <- subset_samples(psnoncontam_tru_final_rare,Description%in%c("T1"))
ps_bf_total <- prune_taxa(taxa_sums(ps_bf_total)>0, ps_bf_total)
ps_bf_total
ps_df_total <- subset_samples(psnoncontam_tru_final_rare,Description%in%c("T2"))
ps_df_total <- prune_taxa(taxa_sums(ps_df_total)>0, ps_df_total)
ps_df_total
ps_ef_total <- subset_samples(psnoncontam_tru_final_rare,Description%in%c("T3"))
ps_ef_total <- prune_taxa(taxa_sums(ps_ef_total)>0, ps_ef_total)
ps_ef_total
ps_ff_total <- subset_samples(psnoncontam_tru_final_rare,Description%in%c("T4"))
ps_ff_total <- prune_taxa(taxa_sums(ps_ff_total)>0, ps_ff_total)
ps_ff_total

#for first time point at phyla level
ps_bf_total
ps_bf_total_merge <- merge_samples(ps_bf_total, "Group")
ps_bf_total_merge
ps_bf_total_merge.rel <- microbiome::transform(ps_bf_total_merge, "compositional")
plot_t1_phylum <- plot_bar(ps_bf_total_merge.rel, fill = "Phylum") + geom_bar(stat="identity", geom_params = list(width=.5)) + theme(panel.background = element_blank()) + scale_fill_manual(values=c("#40530A","#637C33", "#C4CEA7")) + theme(legend.position = "right") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 0, hjust= 0.5), aspect.ratio = 4/4)
plot_t1<- plot_t1_phylum + labs(x="", y= "Relative abundance")
plot_t1


#for second time point at phyla level
ps_df_total
ps_df_total_merge <- merge_samples(ps_df_total, "Group")
ps_df_total_merge
ps_df_total_merge.rel <- microbiome::transform(ps_df_total_merge, "compositional")
plot_t2_phylum <- plot_bar(ps_df_total_merge.rel, fill = "Phylum") + geom_bar(aes( fill=Phylum), stat="identity", position = "stack") + theme(panel.background = element_blank()) + scale_fill_manual(values=c("#40530A","#637C33", "#C4CEA7")) + theme(legend.position = "right") +theme(legend.position = "right") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 0, hjust= 0.5), aspect.ratio = 4/4)
plot_t2<- plot_t2_phylum + labs(x="", y= "Relative abundance")
plot_t2

#for third time point at phyla level
ps_ef_total
ps_ef_total_merge <- merge_samples(ps_ef_total, "Group")
ps_ef_total_merge
ps_ef_total_merge.rel <- microbiome::transform(ps_ef_total_merge, "compositional")
plot_t3_phylum<- plot_bar(ps_ef_total_merge.rel, fill = "Phylum") + geom_bar(aes( fill=Phylum), stat="identity", position = "stack") + theme(panel.background = element_blank()) + scale_fill_manual(values=c("#40530A","#637C33", "#C4CEA7")) + theme(legend.position = "right") + theme(legend.position = "right") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 0, hjust= 0.5), aspect.ratio = 4/4)
plot_t3<- plot_t3_phylum + labs(x="", y= "Relative abundance")
plot_t3

#for fourth time point at phyla level
ps_ff_total
ps_ff_total_merge <- merge_samples(ps_ff_total, "Group")
ps_ff_total_merge
ps_ff_total_merge.rel <- microbiome::transform(ps_ff_total_merge, "compositional")
plot_t4_phylum <- plot_bar(ps_ff_total_merge.rel, fill = "Phylum") + geom_bar(aes( fill=Phylum), stat="identity", position = "stack") + theme(panel.background = element_blank()) +  scale_fill_manual(values=c("#40530A","#637C33", "#C4CEA7")) + theme(legend.position = "right") + theme(legend.position = "right") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 0, hjust= 0.5), aspect.ratio = 4/4)
plot_t4<- plot_t4_phylum + labs(x="", y= "Relative abundance")
plot_t4

###Combining all the phylum level plots for RCT groups

ggarrange(plotlist= list(plot_t1_c, plot_t1_t, plot_t1, plot_t2_c, plot_t2_t, plot_t2, plot_t3_c, plot_t3_t, plot_t3, plot_t4_c, plot_t4_t, plot_t4), ncol=3, nrow = 4, legend = "none")


###########Genus composition at sample level and group level (group mean) for each group (either RCT or MOM groups) at discrete time points. For this analysis it is important that the agglomeration at genus level followed by relative abundance transformation takes place independently. So I have plotted them separately and tried to maintain the order of genus in the per sample plot similar to the group mean plot. This was done comparing the list of genera in individual sample plots with the group mean plots. The order was then set accordinly by sorting the list and adding 'other genus' at the end (details below in the code)

##For Control and treatment samples at 4 time points, group mean and indvidual sample bar plots
####As first step setting up a color palette for all the genera

###I have created a color palette in excel with name of genera and hex codes and imported them as a data frame, so that my plots can pick up the same color for each respective genus every time

##There are other options available for working on color palettes such as the ggsci package, however the number of colors are limited and picking up distinct colours was important to avoid confusion.

##Setting up color palette
color_palette <- read.csv("/Users/shreyas/Documents/Prolacta study docs/sequence data/combined_analysis/color_palette.csv", sep=",", header = T)
color_palette$Genus <- as.character(color_palette$Genus)
color_palette$HEX <- as.character(color_palette$HEX)
color_genus <- setNames(color_palette$HEX, color_palette$Genus)

###As a next step seperating out samples based on time points
ps_bf_total <- subset_samples(psnoncontam_tru_final_rare,Description%in%c("T1"))
ps_bf_total <- prune_taxa(taxa_sums(ps_bf_total)>0, ps_bf_total)
ps_bf_total
ps_df_total <- subset_samples(psnoncontam_tru_final_rare,Description%in%c("T2"))
ps_df_total <- prune_taxa(taxa_sums(ps_df_total)>0, ps_df_total)
ps_df_total
ps_ef_total <- subset_samples(psnoncontam_tru_final_rare,Description%in%c("T3"))
ps_ef_total <- prune_taxa(taxa_sums(ps_ef_total)>0, ps_ef_total)
ps_ef_total
ps_ff_total <- subset_samples(psnoncontam_tru_final_rare,Description%in%c("T4"))
ps_ff_total <- prune_taxa(taxa_sums(ps_ff_total)>0, ps_ff_total)
ps_ff_total

#for first time point split by RCT groups
ps_bf_total
ps_bf_rct_merge <- merge_samples(ps_bf_total, "Group")
ps_bf_rct_merge

#for second time point split by RCT groups
ps_df_total
ps_df_rct_merge <- merge_samples(ps_df_total, "Group")
ps_df_rct_merge

#for third time point split by RCT groups
ps_ef_total
ps_ef_rct_merge <- merge_samples(ps_ef_total, "Group")
ps_ef_rct_merge

#for fourth time point split by RCT groups
ps_ff_total
ps_ff_rct_merge <- merge_samples(ps_ff_total, "Group")
ps_ff_rct_merge

#####Group mean T1, Combined plot###

psbf.rct.gglom <- tax_glom(physeq = ps_bf_rct_merge,taxrank = 'Genus')
#convert to relative abundance and select top 10 taxa
psbf.rct.gglom.prop <- transform_sample_counts(psbf.rct.gglom, function(otu) otu/sum(otu))  

top10 <- names(sort(taxa_sums(psbf.rct.gglom.prop), decreasing=TRUE))[1:10]
ps.top10 <- prune_taxa(top10, psbf.rct.gglom.prop)

otu.tab <- as.data.frame(t(otu_table(ps.top10)))
tax.tab <- as.matrix(tax_table(ps.top10))

#for ordering purpose
temp.levels <- as.character(tax.tab[order(rowSums(otu.tab),decreasing = T),][,'Genus'])

#insert others
otu.tab['Other',] <- as.numeric(1 - colSums(otu.tab))
tax.tab <- rbind(tax.tab,paste0('Other ',rank_names(ps.top10)))
rownames(tax.tab) <- rownames(otu.tab)

#build phyloseq
ps.final <- phyloseq(otu_table(otu.tab,taxa_are_rows = T),tax_table(tax.tab),sample_data(sample_data(ps_bf_rct_merge)))

# combine plot
ps.final.melt <- psmelt(ps.final)
ps.final.melt$Genus <- factor(ps.final.melt$Genus,levels = c(temp.levels,'Other Genus'))
plot_t1_rct_g_group <- ggplot(data = ps.final.melt,mapping = aes(x = Sample, y = Abundance,fill=Genus)) + geom_bar(stat = 'identity') + theme(panel.background = element_blank()) + scale_fill_manual(values = color_genus) + theme(legend.position = "right") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 0, hjust= 0.5), aspect.ratio = 4/4) + labs(x="", y= "Relative abundance")
plot_t1_rct_g_group

###################

####per sample bar plot control
ps.glom <- tax_glom(physeq = ps_control_bf,taxrank = 'Genus')
#convert to relative abundance and select top 10 taxa
ps.glom.prop <- transform_sample_counts(ps.glom, function(otu) otu/sum(otu))  

top10 <- names(sort(taxa_sums(ps.glom.prop), decreasing=TRUE))[1:10]
ps.top10 <- prune_taxa(top10, ps.glom.prop)

otu.tab <- as.data.frame(t(otu_table(ps.top10)))
tax.tab <- as.matrix(tax_table(ps.top10))

#for ordering purpose
#temp.levels <- as.character(tax.tab[order(rowSums(otu.tab),decreasing = T),][,'Genus'])

#insert others
otu.tab['Other',] <- as.numeric(1 - colSums(otu.tab))
tax.tab <- rbind(tax.tab,paste0('Other ',rank_names(ps.top10)))
rownames(tax.tab) <- rownames(otu.tab)

#build phyloseq
ps.final <- phyloseq(otu_table(otu.tab,taxa_are_rows = T),tax_table(tax.tab),sample_data(sample_data(ps_control_bf)))

# per sample genus bar plot for control at T1
ps.final.melt <- psmelt(ps.final)
ps.final.melt$Genus = as.character(ps.final.melt$Genus)

###In order to pick up the same order for genera in individual plots I am picking up the levels from the group mean plot and filter/add the genera that are missing after comparing the list of genera in individual plot with the group mean plot above
rct_genus_t1_levels <- temp.levels[!(temp.levels %in% "Other Genus")]
target_genus = as.character(unique(ps.final.melt$Genus)[!(unique(ps.final.melt$Genus) %in% "Other Genus")])
temp = c(rct_genus_t1_levels[rct_genus_t1_levels %in% target_genus], target_genus[!(target_genus %in% rct_genus_t1_levels)])

ps.final.melt$Genus <- factor(ps.final.melt$Genus,levels = c(temp, "Other Genus"))
levels(ps.final.melt$Genus)
plot_t1_c_genus <- ggplot(data = ps.final.melt,mapping = aes(x = SampleName,y = Abundance,fill=Genus)) + geom_bar(stat = 'identity') + scale_fill_manual(values = color_genus[c(temp, "Other Genus")]) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90), aspect.ratio = 4/4) + labs(x="", y= "Relative abundance")
plot_t1_c_genus

# per sample genus bar plot for treatment at T1

ps.glom <- tax_glom(physeq = ps_treatment_bf,taxrank = 'Genus')
#convert to relative abundance and select top 10 taxa
ps.glom.prop <- transform_sample_counts(ps.glom, function(otu) otu/sum(otu))  

top10 <- names(sort(taxa_sums(ps.glom.prop), decreasing=TRUE))[1:10]
ps.top10 <- prune_taxa(top10, ps.glom.prop)

otu.tab <- as.data.frame(t(otu_table(ps.top10)))
tax.tab <- as.matrix(tax_table(ps.top10))

#for ordering purpose. I used it previously to order the taxa according to their abundance.
#temp.levels <- as.character(tax.tab[order(rowSums(otu.tab),decreasing = T),][,'Genus'])

#insert others
otu.tab['Other',] <- as.numeric(1 - colSums(otu.tab))
tax.tab <- rbind(tax.tab,paste0('Other ',rank_names(ps.top10)))
rownames(tax.tab) <- rownames(otu.tab)

#build phyloseq
ps.final <- phyloseq(otu_table(otu.tab,taxa_are_rows = T),tax_table(tax.tab),sample_data(sample_data(ps_treatment_bf)))

# plot for treatment at T1
ps.final.melt <- psmelt(ps.final)
ps.final.melt$Genus = as.character(ps.final.melt$Genus)

###Again as done for the control plot, in order to pick up the same order for genera in individual plots I am picking up the levels from the group mean plot and filter/add the genera that are missing after comparing the list of genera in individual plot with the group mean plot above

rct_genus_t1_levels <- temp.levels[!(temp.levels %in% "Other Genus")]
target_genus = as.character(unique(ps.final.melt$Genus)[!(unique(ps.final.melt$Genus) %in% "Other Genus")])
temp = c(rct_genus_t1_levels[rct_genus_t1_levels %in% target_genus], target_genus[!(target_genus %in% rct_genus_t1_levels)])

ps.final.melt$Genus <- factor(ps.final.melt$Genus,levels = c(temp, "Other Genus"))
levels(ps.final.melt$Genus)
plot_t1_t_genus <- ggplot(data = ps.final.melt,mapping = aes(x = SampleName,y = Abundance,fill=Genus)) + geom_bar(stat = 'identity') + scale_fill_manual(values = color_genus[c(temp, "Other Genus")]) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90), aspect.ratio = 4/4) + labs(x="", y= "Relative abundance")
plot_t1_t_genus


###For time point T2
#####Group mean T2, Combined plot###

psdf.rct.gglom <- tax_glom(physeq = ps_df_rct_merge,taxrank = 'Genus')
#convert to relative abundance and select top 10 taxa
psdf.rct.gglom.prop <- transform_sample_counts(psdf.rct.gglom, function(otu) otu/sum(otu))  

top10 <- names(sort(taxa_sums(psdf.rct.gglom.prop), decreasing=TRUE))[1:10]
ps.top10 <- prune_taxa(top10, psdf.rct.gglom.prop)

otu.tab <- as.data.frame(t(otu_table(ps.top10)))
tax.tab <- as.matrix(tax_table(ps.top10))

#for ordering purpose
temp.levels <- as.character(tax.tab[order(rowSums(otu.tab),decreasing = T),][,'Genus'])

#insert others
otu.tab['Other',] <- as.numeric(1 - colSums(otu.tab))
tax.tab <- rbind(tax.tab,paste0('Other ',rank_names(ps.top10)))
rownames(tax.tab) <- rownames(otu.tab)

#build phyloseq
ps.final <- phyloseq(otu_table(otu.tab,taxa_are_rows = T),tax_table(tax.tab),sample_data(sample_data(ps_df_rct_merge)))

# combine plot
ps.final.melt <- psmelt(ps.final)
ps.final.melt$Genus <- factor(ps.final.melt$Genus,levels = c(temp.levels,'Other Genus'))
plot_t2_rct_g_group <- ggplot(data = ps.final.melt,mapping = aes(x = Sample, y = Abundance,fill=Genus)) + geom_bar(stat = 'identity') + theme(panel.background = element_blank()) + scale_fill_manual(values = color_genus) + theme(legend.position = "right") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 0, hjust= 0.5), aspect.ratio = 4/4) + labs(x="", y= "Relative abundance")
plot_t2_rct_g_group

###################

####per sample bar plot control
ps.glom <- tax_glom(physeq = ps_control_df,taxrank = 'Genus')
#convert to relative abundance and select top 10 taxa
ps.glom.prop <- transform_sample_counts(ps.glom, function(otu) otu/sum(otu))  

top10 <- names(sort(taxa_sums(ps.glom.prop), decreasing=TRUE))[1:10]
ps.top10 <- prune_taxa(top10, ps.glom.prop)

otu.tab <- as.data.frame(t(otu_table(ps.top10)))
tax.tab <- as.matrix(tax_table(ps.top10))

#for ordering purpose
#temp.levels <- as.character(tax.tab[order(rowSums(otu.tab),decreasing = T),][,'Genus'])

#insert others
otu.tab['Other',] <- as.numeric(1 - colSums(otu.tab))
tax.tab <- rbind(tax.tab,paste0('Other ',rank_names(ps.top10)))
rownames(tax.tab) <- rownames(otu.tab)

#build phyloseq
ps.final <- phyloseq(otu_table(otu.tab,taxa_are_rows = T),tax_table(tax.tab),sample_data(sample_data(ps_control_df)))

# per sample genus bar plot for control at T2
ps.final.melt <- psmelt(ps.final)
ps.final.melt$Genus = as.character(ps.final.melt$Genus)

###In order to pick up the same order for genera in individual plots I am picking up the levels from the group mean plot and filter/add the genera that are missing after comparing the list of genera in individual plot with the group mean plot above
rct_genus_t2_levels <- temp.levels[!(temp.levels %in% "Other Genus")]
target_genus = as.character(unique(ps.final.melt$Genus)[!(unique(ps.final.melt$Genus) %in% "Other Genus")])
temp = c(rct_genus_t2_levels[rct_genus_t2_levels %in% target_genus], target_genus[!(target_genus %in% rct_genus_t2_levels)])

ps.final.melt$Genus <- factor(ps.final.melt$Genus,levels = c(temp, "Other Genus"))
levels(ps.final.melt$Genus)
plot_t2_c_genus <- ggplot(data = ps.final.melt,mapping = aes(x = SampleName,y = Abundance,fill=Genus)) + geom_bar(stat = 'identity') + scale_fill_manual(values = color_genus[c(temp, "Other Genus")]) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90), aspect.ratio = 4/4) + labs(x="", y= "Relative abundance")
plot_t2_c_genus


# per sample genus bar plot for treatment at T2
ps.glom <- tax_glom(physeq = ps_treatment_df,taxrank = 'Genus')
#convert to relative abundance and select top 10 taxa
ps.glom.prop <- transform_sample_counts(ps.glom, function(otu) otu/sum(otu))  

top10 <- names(sort(taxa_sums(ps.glom.prop), decreasing=TRUE))[1:10]
ps.top10 <- prune_taxa(top10, ps.glom.prop)

otu.tab <- as.data.frame(t(otu_table(ps.top10)))
tax.tab <- as.matrix(tax_table(ps.top10))

#for ordering purpose. I used it previously to order the taxa according to their abundance.
#temp.levels <- as.character(tax.tab[order(rowSums(otu.tab),decreasing = T),][,'Genus'])

#insert others
otu.tab['Other',] <- as.numeric(1 - colSums(otu.tab))
tax.tab <- rbind(tax.tab,paste0('Other ',rank_names(ps.top10)))
rownames(tax.tab) <- rownames(otu.tab)

#build phyloseq
ps.final <- phyloseq(otu_table(otu.tab,taxa_are_rows = T),tax_table(tax.tab),sample_data(sample_data(ps_treatment_df)))

# plot for treatment at T2
ps.final.melt <- psmelt(ps.final)
ps.final.melt$Genus = as.character(ps.final.melt$Genus)

###Again as done for the control plot, in order to pick up the same order for genera in individual plots I am picking up the levels from the group mean plot and filter/add the genera that are missing after comparing the list of genera in individual plot with the group mean plot above

rct_genus_t2_levels <- temp.levels[!(temp.levels %in% "Other Genus")]
target_genus = as.character(unique(ps.final.melt$Genus)[!(unique(ps.final.melt$Genus) %in% "Other Genus")])
temp = c(rct_genus_t2_levels[rct_genus_t2_levels %in% target_genus], target_genus[!(target_genus %in% rct_genus_t2_levels)])

ps.final.melt$Genus <- factor(ps.final.melt$Genus,levels = c(temp, "Other Genus"))
levels(ps.final.melt$Genus)
plot_t2_t_genus <- ggplot(data = ps.final.melt,mapping = aes(x = SampleName,y = Abundance,fill=Genus)) + geom_bar(stat = 'identity') + scale_fill_manual(values = color_genus[c(temp, "Other Genus")]) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90), aspect.ratio = 4/4) + labs(x="", y= "Relative abundance")
plot_t2_t_genus


##For T3 by RCT groups
#####Group mean T3, Combined plot###

psef.rct.gglom <- tax_glom(physeq = ps_ef_rct_merge,taxrank = 'Genus')
#convert to relative abundance and select top 10 taxa
psef.rct.gglom.prop <- transform_sample_counts(psef.rct.gglom, function(otu) otu/sum(otu))  

top10 <- names(sort(taxa_sums(psef.rct.gglom.prop), decreasing=TRUE))[1:10]
ps.top10 <- prune_taxa(top10, psef.rct.gglom.prop)

otu.tab <- as.data.frame(t(otu_table(ps.top10)))
tax.tab <- as.matrix(tax_table(ps.top10))

#for ordering purpose
temp.levels <- as.character(tax.tab[order(rowSums(otu.tab),decreasing = T),][,'Genus'])

#insert others
otu.tab['Other',] <- as.numeric(1 - colSums(otu.tab))
tax.tab <- rbind(tax.tab,paste0('Other ',rank_names(ps.top10)))
rownames(tax.tab) <- rownames(otu.tab)

#build phyloseq
ps.final <- phyloseq(otu_table(otu.tab,taxa_are_rows = T),tax_table(tax.tab),sample_data(sample_data(ps_ef_rct_merge)))

# combine plot
ps.final.melt <- psmelt(ps.final)
ps.final.melt$Genus <- factor(ps.final.melt$Genus,levels = c(temp.levels,'Other Genus'))
plot_t3_rct_g_group <- ggplot(data = ps.final.melt,mapping = aes(x = Sample, y = Abundance,fill=Genus)) + geom_bar(stat = 'identity') + theme(panel.background = element_blank()) + scale_fill_manual(values = color_genus) + theme(legend.position = "right") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 0, hjust= 0.5), aspect.ratio = 4/4) + labs(x="", y= "Relative abundance")
plot_t3_rct_g_group

###################

####per sample bar plot control
ps.glom <- tax_glom(physeq = ps_control_ef,taxrank = 'Genus')
#convert to relative abundance and select top 10 taxa
ps.glom.prop <- transform_sample_counts(ps.glom, function(otu) otu/sum(otu))  

top10 <- names(sort(taxa_sums(ps.glom.prop), decreasing=TRUE))[1:10]
ps.top10 <- prune_taxa(top10, ps.glom.prop)

otu.tab <- as.data.frame(t(otu_table(ps.top10)))
tax.tab <- as.matrix(tax_table(ps.top10))

#for ordering purpose if to be plotted separately, this oredering is according to the abundance
#temp.levels <- as.character(tax.tab[order(rowSums(otu.tab),decreasing = T),][,'Genus'])

#insert others
otu.tab['Other',] <- as.numeric(1 - colSums(otu.tab))
tax.tab <- rbind(tax.tab,paste0('Other ',rank_names(ps.top10)))
rownames(tax.tab) <- rownames(otu.tab)

#build phyloseq
ps.final <- phyloseq(otu_table(otu.tab,taxa_are_rows = T),tax_table(tax.tab),sample_data(sample_data(ps_control_ef)))

# per sample genus bar plot for control at T3
ps.final.melt <- psmelt(ps.final)
ps.final.melt$Genus = as.character(ps.final.melt$Genus)

###In order to pick up the same order for genera in individual plots I am picking up the levels from the group mean plot and filter/add the genera that are missing after comparing the list of genera in individual plot with the group mean plot above
rct_genus_t3_levels <- temp.levels[!(temp.levels %in% "Other Genus")]
target_genus = as.character(unique(ps.final.melt$Genus)[!(unique(ps.final.melt$Genus) %in% "Other Genus")])
temp = c(rct_genus_t3_levels[rct_genus_t3_levels %in% target_genus], target_genus[!(target_genus %in% rct_genus_t3_levels)])

ps.final.melt$Genus <- factor(ps.final.melt$Genus,levels = c(temp, "Other Genus"))
levels(ps.final.melt$Genus)
plot_t3_c_genus <- ggplot(data = ps.final.melt,mapping = aes(x = SampleName,y = Abundance,fill=Genus)) + geom_bar(stat = 'identity') + scale_fill_manual(values = color_genus[c(temp, "Other Genus")]) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90), aspect.ratio = 4/4) + labs(x="", y= "Relative abundance")
plot_t3_c_genus


# per sample genus bar plot for treatment at T3

ps.glom <- tax_glom(physeq = ps_treatment_ef,taxrank = 'Genus')
#convert to relative abundance and select top 10 taxa
ps.glom.prop <- transform_sample_counts(ps.glom, function(otu) otu/sum(otu))  

top10 <- names(sort(taxa_sums(ps.glom.prop), decreasing=TRUE))[1:10]
ps.top10 <- prune_taxa(top10, ps.glom.prop)

otu.tab <- as.data.frame(t(otu_table(ps.top10)))
tax.tab <- as.matrix(tax_table(ps.top10))

#for ordering purpose. I used it previously to order the taxa according to their abundance.
#temp.levels <- as.character(tax.tab[order(rowSums(otu.tab),decreasing = T),][,'Genus'])

#insert others
otu.tab['Other',] <- as.numeric(1 - colSums(otu.tab))
tax.tab <- rbind(tax.tab,paste0('Other ',rank_names(ps.top10)))
rownames(tax.tab) <- rownames(otu.tab)

#build phyloseq
ps.final <- phyloseq(otu_table(otu.tab,taxa_are_rows = T),tax_table(tax.tab),sample_data(sample_data(ps_treatment_ef)))

# plot for treatment at T3
ps.final.melt <- psmelt(ps.final)
ps.final.melt$Genus = as.character(ps.final.melt$Genus)

###Again as done for the control plot, in order to pick up the same order for genera in individual plots I am picking up the levels from the group mean plot and filter/add the genera that are missing after comparing the list of genera in individual plot with the group mean plot above

rct_genus_t3_levels <- temp.levels[!(temp.levels %in% "Other Genus")]
target_genus = as.character(unique(ps.final.melt$Genus)[!(unique(ps.final.melt$Genus) %in% "Other Genus")])
temp = c(rct_genus_t3_levels[rct_genus_t3_levels %in% target_genus], target_genus[!(target_genus %in% rct_genus_t3_levels)])

ps.final.melt$Genus <- factor(ps.final.melt$Genus,levels = c(temp, "Other Genus"))
levels(ps.final.melt$Genus)
plot_t3_t_genus <- ggplot(data = ps.final.melt,mapping = aes(x = SampleName,y = Abundance,fill=Genus)) + geom_bar(stat = 'identity') + scale_fill_manual(values = color_genus[c(temp, "Other Genus")]) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90), aspect.ratio = 4/4) + labs(x="", y= "Relative abundance")
plot_t3_t_genus


##For T4 by RCT

#####Group mean T4, Combined plot###

psff.rct.gglom <- tax_glom(physeq = ps_ff_rct_merge,taxrank = 'Genus')
#convert to relative abundance and select top 10 taxa
psff.rct.gglom.prop <- transform_sample_counts(psff.rct.gglom, function(otu) otu/sum(otu))  

top10 <- names(sort(taxa_sums(psff.rct.gglom.prop), decreasing=TRUE))[1:10]
ps.top10 <- prune_taxa(top10, psff.rct.gglom.prop)

otu.tab <- as.data.frame(t(otu_table(ps.top10)))
tax.tab <- as.matrix(tax_table(ps.top10))

#for ordering purpose
temp.levels <- as.character(tax.tab[order(rowSums(otu.tab),decreasing = T),][,'Genus'])

#insert others
otu.tab['Other',] <- as.numeric(1 - colSums(otu.tab))
tax.tab <- rbind(tax.tab,paste0('Other ',rank_names(ps.top10)))
rownames(tax.tab) <- rownames(otu.tab)

#build phyloseq
ps.final <- phyloseq(otu_table(otu.tab,taxa_are_rows = T),tax_table(tax.tab),sample_data(sample_data(ps_ff_rct_merge)))

# combine plot
ps.final.melt <- psmelt(ps.final)
ps.final.melt$Genus <- factor(ps.final.melt$Genus,levels = c(temp.levels,'Other Genus'))
plot_t4_rct_g_group <- ggplot(data = ps.final.melt,mapping = aes(x = Sample, y = Abundance,fill=Genus)) + geom_bar(stat = 'identity') + theme(panel.background = element_blank()) + scale_fill_manual(values = color_genus) + theme(legend.position = "right") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 0, hjust= 0.5), aspect.ratio = 4/4) + labs(x="", y= "Relative abundance")
plot_t4_rct_g_group

###################

####per sample bar plot control
ps.glom <- tax_glom(physeq = ps_control_ff,taxrank = 'Genus')
#convert to relative abundance and select top 10 taxa
ps.glom.prop <- transform_sample_counts(ps.glom, function(otu) otu/sum(otu))  

top10 <- names(sort(taxa_sums(ps.glom.prop), decreasing=TRUE))[1:10]
ps.top10 <- prune_taxa(top10, ps.glom.prop)

otu.tab <- as.data.frame(t(otu_table(ps.top10)))
tax.tab <- as.matrix(tax_table(ps.top10))

#for ordering purpose if to be plotted separately, this oredering is according to the abundance
#temp.levels <- as.character(tax.tab[order(rowSums(otu.tab),decreasing = T),][,'Genus'])

#insert others
otu.tab['Other',] <- as.numeric(1 - colSums(otu.tab))
tax.tab <- rbind(tax.tab,paste0('Other ',rank_names(ps.top10)))
rownames(tax.tab) <- rownames(otu.tab)

#build phyloseq
ps.final <- phyloseq(otu_table(otu.tab,taxa_are_rows = T),tax_table(tax.tab),sample_data(sample_data(ps_control_ff)))

# per sample genus bar plot for control at T1
ps.final.melt <- psmelt(ps.final)
ps.final.melt$Genus = as.character(ps.final.melt$Genus)

###In order to pick up the same order for genera in individual plots I am picking up the levels from the group mean plot and filter/add the genera that are missing after comparing the list of genera in individual plot with the group mean plot above
rct_genus_t4_levels <- temp.levels[!(temp.levels %in% "Other Genus")]
target_genus = as.character(unique(ps.final.melt$Genus)[!(unique(ps.final.melt$Genus) %in% "Other Genus")])
temp = c(rct_genus_t4_levels[rct_genus_t4_levels %in% target_genus], target_genus[!(target_genus %in% rct_genus_t4_levels)])

ps.final.melt$Genus <- factor(ps.final.melt$Genus,levels = c(temp, "Other Genus"))
levels(ps.final.melt$Genus)
plot_t4_c_genus <- ggplot(data = ps.final.melt,mapping = aes(x = SampleName,y = Abundance,fill=Genus)) + geom_bar(stat = 'identity') + scale_fill_manual(values = color_genus[c(temp, "Other Genus")]) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90), aspect.ratio = 4/4) + labs(x="", y= "Relative abundance")
plot_t4_c_genus


# per sample genus bar plot for treatment at T3

ps.glom <- tax_glom(physeq = ps_treatment_ff,taxrank = 'Genus')
#convert to relative abundance and select top 10 taxa
ps.glom.prop <- transform_sample_counts(ps.glom, function(otu) otu/sum(otu))  

top10 <- names(sort(taxa_sums(ps.glom.prop), decreasing=TRUE))[1:10]
ps.top10 <- prune_taxa(top10, ps.glom.prop)

otu.tab <- as.data.frame(t(otu_table(ps.top10)))
tax.tab <- as.matrix(tax_table(ps.top10))

#for ordering purpose. I used it previously to order the taxa according to their abundance.
#temp.levels <- as.character(tax.tab[order(rowSums(otu.tab),decreasing = T),][,'Genus'])

#insert others
otu.tab['Other',] <- as.numeric(1 - colSums(otu.tab))
tax.tab <- rbind(tax.tab,paste0('Other ',rank_names(ps.top10)))
rownames(tax.tab) <- rownames(otu.tab)

#build phyloseq
ps.final <- phyloseq(otu_table(otu.tab,taxa_are_rows = T),tax_table(tax.tab),sample_data(sample_data(ps_treatment_ff)))

# plot for treatment at T3
ps.final.melt <- psmelt(ps.final)
ps.final.melt$Genus = as.character(ps.final.melt$Genus)

###Again as done for the control plot, in order to pick up the same order for genera in individual plots I am picking up the levels from the group mean plot and filter/add the genera that are missing after comparing the list of genera in individual plot with the group mean plot above

rct_genus_t4_levels <- temp.levels[!(temp.levels %in% "Other Genus")]
target_genus = as.character(unique(ps.final.melt$Genus)[!(unique(ps.final.melt$Genus) %in% "Other Genus")])
temp = c(rct_genus_t4_levels[rct_genus_t4_levels %in% target_genus], target_genus[!(target_genus %in% rct_genus_t4_levels)])

ps.final.melt$Genus <- factor(ps.final.melt$Genus,levels = c(temp, "Other Genus"))
levels(ps.final.melt$Genus)
plot_t4_t_genus <- ggplot(data = ps.final.melt,mapping = aes(x = SampleName,y = Abundance,fill=Genus)) + geom_bar(stat = 'identity') + scale_fill_manual(values = color_genus[c(temp, "Other Genus")]) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90), aspect.ratio = 4/4) + labs(x="", y= "Relative abundance")
plot_t4_t_genus


###Combining all the plots now
ggarrange(plotlist = list(plot_t3_c_genus, plot_t3_t_genus, plot_t3_rct_g_group, plot_t2_c_genus, plot_t2_t_genus, plot_t2_rct_g_group, plot_t3_c_genus, plot_t3_t_genus, plot_t3_rct_g_group, plot_t4_c_genus, plot_t4_t_genus, plot_t4_rct_g_group), ncol=3, nrow=4, legend = "none")


####Repeating the same for MOM groups. Plotting genera level composition at each time point at sample and group level
###I had already separated the phyloseq objects for each time point above, so I am using the objects here and now merging them using MOM_Group variable

#for first time point split by MOM groups
ps_bf_total
ps_bf_mom_merge <- merge_samples(ps_bf_total, "MOM_Group")
ps_bf_mom_merge

#for second time point split by MOM groups
ps_df_total
ps_df_mom_merge <- merge_samples(ps_df_total, "MOM_Group")
ps_df_mom_merge

#for third time point split by MOM groups
ps_ef_total
ps_ef_mom_merge <- merge_samples(ps_ef_total, "MOM_Group")
ps_ef_mom_merge

#for fourth time point split by MOM groups
ps_ff_total
ps_ff_mom_merge <- merge_samples(ps_ff_total, "MOM_Group")
ps_ff_mom_merge

####Now separating out based on MOM groups at each time point

#subset LMOM group samples
ps_LMOM_final <- subset_samples(psnoncontam_tru_final_rare,MOM_Group%in%c("LMOM"))
ps_LMOM_final
ps_LMOM_final <- prune_taxa(taxa_sums(ps_LMOM_final)>0, ps_LMOM_final)
ps_LMOM_final

#subset HMOM group samples

ps_HMOM_final <- subset_samples(psnoncontam_tru_final_rare,MOM_Group%in%c("HMOM"))
ps_HMOM_final
ps_HMOM_final <- prune_taxa(taxa_sums(ps_HMOM_final)>0, ps_HMOM_final)
ps_HMOM_final

##Now subsetting each LMOM and HMOM group for all time points

###First for LMOM groups
ps_LMOM_bf <- subset_samples(ps_LMOM_final,Description%in%c("T1"))
ps_LMOM_bf
ps_LMOM_bf <- prune_taxa(taxa_sums(ps_LMOM_bf)>0, ps_LMOM_bf)
ps_LMOM_bf

ps_LMOM_df <- subset_samples(ps_LMOM_final,Description%in%c("T2"))
ps_LMOM_df
ps_LMOM_df <- prune_taxa(taxa_sums(ps_LMOM_df)>0, ps_LMOM_df)
ps_LMOM_df

ps_LMOM_ef <- subset_samples(ps_LMOM_final,Description%in%c("T3"))
ps_LMOM_ef
ps_LMOM_ef <- prune_taxa(taxa_sums(ps_LMOM_ef)>0, ps_LMOM_ef)
ps_LMOM_ef

ps_LMOM_ff <- subset_samples(ps_LMOM_final,Description%in%c("T4"))
ps_LMOM_ff
ps_LMOM_ff <- prune_taxa(taxa_sums(ps_LMOM_ff)>0, ps_LMOM_ff)
ps_LMOM_ff

###for HMOM groups
ps_HMOM_bf <- subset_samples(ps_HMOM_final,Description%in%c("T3"))
ps_HMOM_bf
ps_HMOM_bf <- prune_taxa(taxa_sums(ps_HMOM_bf)>0, ps_HMOM_bf)
ps_HMOM_bf

ps_HMOM_df <- subset_samples(ps_HMOM_final,Description%in%c("T2"))
ps_HMOM_df
ps_HMOM_df <- prune_taxa(taxa_sums(ps_HMOM_df)>0, ps_HMOM_df)
ps_HMOM_df

ps_HMOM_ef <- subset_samples(ps_HMOM_final,Description%in%c("T3"))
ps_HMOM_ef
ps_HMOM_ef <- prune_taxa(taxa_sums(ps_HMOM_ef)>0, ps_HMOM_ef)
ps_HMOM_ef

ps_HMOM_ff <- subset_samples(ps_HMOM_final,Description%in%c("T4"))
ps_HMOM_ff
ps_HMOM_ff <- prune_taxa(taxa_sums(ps_HMOM_ff)>0, ps_HMOM_ff)
ps_HMOM_ff

#####Group mean T1, Combined plot###

psbf.mom.gglom <- tax_glom(physeq = ps_bf_mom_merge,taxrank = 'Genus')
#convert to relative abundance and select top 10 taxa
psbf.mom.gglom.prop <- transform_sample_counts(psbf.mom.gglom, function(otu) otu/sum(otu))  

top10 <- names(sort(taxa_sums(psbf.mom.gglom.prop), decreasing=TRUE))[1:10]
ps.top10 <- prune_taxa(top10, psbf.mom.gglom.prop)

otu.tab <- as.data.frame(t(otu_table(ps.top10)))
tax.tab <- as.matrix(tax_table(ps.top10))

#for ordering purpose
temp.levels <- as.character(tax.tab[order(rowSums(otu.tab),decreasing = T),][,'Genus'])

#insert others
otu.tab['Other',] <- as.numeric(1 - colSums(otu.tab))
tax.tab <- rbind(tax.tab,paste0('Other ',rank_names(ps.top10)))
rownames(tax.tab) <- rownames(otu.tab)

#build phyloseq
ps.final <- phyloseq(otu_table(otu.tab,taxa_are_rows = T),tax_table(tax.tab),sample_data(sample_data(ps_bf_mom_merge)))

# combine plot
ps.final.melt <- psmelt(ps.final)
ps.final.melt$Genus <- factor(ps.final.melt$Genus,levels = c(temp.levels,'Other Genus'))
ps.final.melt$Sample <- factor(ps.final.melt$Sample,levels = c("LMOM", "HMOM"))
plot_t1_mom_g_group <- ggplot(data = ps.final.melt,mapping = aes(x = Sample, y = Abundance,fill=Genus)) + geom_bar(stat = 'identity') + theme(panel.background = element_blank()) + scale_fill_manual(values = color_genus) + theme(legend.position = "right") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 0, hjust= 0.5), aspect.ratio = 4/4) + labs(x="", y= "Relative abundance")
plot_t1_mom_g_group

###################

####per sample bar plot LMOM
ps.glom <- tax_glom(physeq = ps_LMOM_bf,taxrank = 'Genus')
#convert to relative abundance and select top 10 taxa
ps.glom.prop <- transform_sample_counts(ps.glom, function(otu) otu/sum(otu))  

top10 <- names(sort(taxa_sums(ps.glom.prop), decreasing=TRUE))[1:10]
ps.top10 <- prune_taxa(top10, ps.glom.prop)

otu.tab <- as.data.frame(t(otu_table(ps.top10)))
tax.tab <- as.matrix(tax_table(ps.top10))

#for ordering purpose
#temp.levels <- as.character(tax.tab[order(rowSums(otu.tab),decreasing = T),][,'Genus'])

#insert others
otu.tab['Other',] <- as.numeric(1 - colSums(otu.tab))
tax.tab <- rbind(tax.tab,paste0('Other ',rank_names(ps.top10)))
rownames(tax.tab) <- rownames(otu.tab)

#build phyloseq
ps.final <- phyloseq(otu_table(otu.tab,taxa_are_rows = T),tax_table(tax.tab),sample_data(sample_data(ps_LMOM_bf)))

# per sample genus bar plot for LMOM at T1
ps.final.melt <- psmelt(ps.final)
ps.final.melt$Genus = as.character(ps.final.melt$Genus)

###In order to pick up the same order for genera in individual plots I am picking up the levels from the group mean plot and filter/add the genera that are missing after comparing the list of genera in individual plot with the group mean plot above
MOM_genus_t1_levels <- temp.levels[!(temp.levels %in% "Other Genus")]
target_genus = as.character(unique(ps.final.melt$Genus)[!(unique(ps.final.melt$Genus) %in% "Other Genus")])
temp = c(MOM_genus_t1_levels[MOM_genus_t1_levels %in% target_genus], target_genus[!(target_genus %in% MOM_genus_t1_levels)])

ps.final.melt$Genus <- factor(ps.final.melt$Genus,levels = c(temp, "Other Genus"))
levels(ps.final.melt$Genus)
plot_t1_l_genus <- ggplot(data = ps.final.melt,mapping = aes(x = SampleName,y = Abundance,fill=Genus)) + geom_bar(stat = 'identity') + scale_fill_manual(values = color_genus[c(temp, "Other Genus")]) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90), aspect.ratio = 4/2) + labs(x="", y= "Relative abundance")
plot_t1_l_genus

# per sample genus bar plot for HMOM at T1

ps.glom <- tax_glom(physeq = ps_HMOM_bf,taxrank = 'Genus')
#convert to relative abundance and select top 10 taxa
ps.glom.prop <- transform_sample_counts(ps.glom, function(otu) otu/sum(otu))  

top10 <- names(sort(taxa_sums(ps.glom.prop), decreasing=TRUE))[1:10]
ps.top10 <- prune_taxa(top10, ps.glom.prop)

otu.tab <- as.data.frame(t(otu_table(ps.top10)))
tax.tab <- as.matrix(tax_table(ps.top10))

#for ordering purpose. I used it previously to order the taxa according to their abundance.
#temp.levels <- as.character(tax.tab[order(rowSums(otu.tab),decreasing = T),][,'Genus'])

#insert others
otu.tab['Other',] <- as.numeric(1 - colSums(otu.tab))
tax.tab <- rbind(tax.tab,paste0('Other ',rank_names(ps.top10)))
rownames(tax.tab) <- rownames(otu.tab)

#build phyloseq
ps.final <- phyloseq(otu_table(otu.tab,taxa_are_rows = T),tax_table(tax.tab),sample_data(sample_data(ps_HMOM_bf)))

# plot for treatment at T1
ps.final.melt <- psmelt(ps.final)
ps.final.melt$Genus = as.character(ps.final.melt$Genus)

###Again as done for the control plot, in order to pick up the same order for genera in individual plots I am picking up the levels from the group mean plot and filter/add the genera that are missing after comparing the list of genera in individual plot with the group mean plot above

MOM_genus_t1_levels <- temp.levels[!(temp.levels %in% "Other Genus")]
target_genus = as.character(unique(ps.final.melt$Genus)[!(unique(ps.final.melt$Genus) %in% "Other Genus")])
temp = c(MOM_genus_t1_levels[MOM_genus_t1_levels %in% target_genus], target_genus[!(target_genus %in% MOM_genus_t1_levels)])

ps.final.melt$Genus <- factor(ps.final.melt$Genus,levels = c(temp, "Other Genus"))
levels(ps.final.melt$Genus)
plot_t1_h_genus <- ggplot(data = ps.final.melt,mapping = aes(x = SampleName,y = Abundance,fill=Genus)) + geom_bar(stat = 'identity') + scale_fill_manual(values = color_genus[c(temp, "Other Genus")]) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90), aspect.ratio = 4/4) + labs(x="", y= "Relative abundance")
plot_t1_h_genus


####For T2 by MOM groups
#####Group mean T2, Combined plot###

psdf.mom.gglom <- tax_glom(physeq = ps_df_mom_merge,taxrank = 'Genus')
#convert to relative abundance and select top 10 taxa
psdf.mom.gglom.prop <- transform_sample_counts(psdf.mom.gglom, function(otu) otu/sum(otu))  

top10 <- names(sort(taxa_sums(psdf.mom.gglom.prop), decreasing=TRUE))[1:10]
ps.top10 <- prune_taxa(top10, psdf.mom.gglom.prop)

otu.tab <- as.data.frame(t(otu_table(ps.top10)))
tax.tab <- as.matrix(tax_table(ps.top10))

#for ordering purpose
temp.levels <- as.character(tax.tab[order(rowSums(otu.tab),decreasing = T),][,'Genus'])

#insert others
otu.tab['Other',] <- as.numeric(1 - colSums(otu.tab))
tax.tab <- rbind(tax.tab,paste0('Other ',rank_names(ps.top10)))
rownames(tax.tab) <- rownames(otu.tab)

#build phyloseq
ps.final <- phyloseq(otu_table(otu.tab,taxa_are_rows = T),tax_table(tax.tab),sample_data(sample_data(ps_df_mom_merge)))

# combine plot
ps.final.melt <- psmelt(ps.final)
ps.final.melt$Genus <- factor(ps.final.melt$Genus,levels = c(temp.levels,'Other Genus'))
ps.final.melt$Sample <- factor(ps.final.melt$Sample,levels = c("LMOM", "HMOM"))
plot_t2_mom_g_group <- ggplot(data = ps.final.melt,mapping = aes(x = Sample, y = Abundance,fill=Genus)) + geom_bar(stat = 'identity') + theme(panel.background = element_blank()) + scale_fill_manual(values = color_genus) + theme(legend.position = "right") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 0, hjust= 0.5), aspect.ratio = 4/4) + labs(x="", y= "Relative abundance")
plot_t2_mom_g_group

###################

####per sample bar plot LMOM
ps.glom <- tax_glom(physeq = ps_LMOM_df,taxrank = 'Genus')
#convert to relative abundance and select top 10 taxa
ps.glom.prop <- transform_sample_counts(ps.glom, function(otu) otu/sum(otu))  

top10 <- names(sort(taxa_sums(ps.glom.prop), decreasing=TRUE))[1:10]
ps.top10 <- prune_taxa(top10, ps.glom.prop)

otu.tab <- as.data.frame(t(otu_table(ps.top10)))
tax.tab <- as.matrix(tax_table(ps.top10))

#for ordering purpose
#temp.levels <- as.character(tax.tab[order(rowSums(otu.tab),decreasing = T),][,'Genus'])

#insert others
otu.tab['Other',] <- as.numeric(1 - colSums(otu.tab))
tax.tab <- rbind(tax.tab,paste0('Other ',rank_names(ps.top10)))
rownames(tax.tab) <- rownames(otu.tab)

#build phyloseq
ps.final <- phyloseq(otu_table(otu.tab,taxa_are_rows = T),tax_table(tax.tab),sample_data(sample_data(ps_LMOM_df)))

# per sample genus bar plot for LMOM at T2
ps.final.melt <- psmelt(ps.final)
ps.final.melt$Genus = as.character(ps.final.melt$Genus)

###In order to pick up the same order for genera in individual plots I am picking up the levels from the group mean plot and filter/add the genera that are missing after comparing the list of genera in individual plot with the group mean plot above
MOM_genus_t2_levels <- temp.levels[!(temp.levels %in% "Other Genus")]
target_genus = as.character(unique(ps.final.melt$Genus)[!(unique(ps.final.melt$Genus) %in% "Other Genus")])
temp = c(MOM_genus_t2_levels[MOM_genus_t2_levels %in% target_genus], target_genus[!(target_genus %in% MOM_genus_t2_levels)])

ps.final.melt$Genus <- factor(ps.final.melt$Genus,levels = c(temp, "Other Genus"))
levels(ps.final.melt$Genus)
plot_t2_l_genus <- ggplot(data = ps.final.melt,mapping = aes(x = SampleName,y = Abundance,fill=Genus)) + geom_bar(stat = 'identity') + scale_fill_manual(values = color_genus[c(temp, "Other Genus")]) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90), aspect.ratio = 4/2) + labs(x="", y= "Relative abundance")
plot_t2_l_genus


# per sample genus bar plot for HMOM at T2

ps.glom <- tax_glom(physeq = ps_HMOM_df,taxrank = 'Genus')
#convert to relative abundance and select top 10 taxa
ps.glom.prop <- transform_sample_counts(ps.glom, function(otu) otu/sum(otu))  

top10 <- names(sort(taxa_sums(ps.glom.prop), decreasing=TRUE))[1:10]
ps.top10 <- prune_taxa(top10, ps.glom.prop)

otu.tab <- as.data.frame(t(otu_table(ps.top10)))
tax.tab <- as.matrix(tax_table(ps.top10))

#for ordering purpose. I used it previously to order the taxa according to their abundance.
#temp.levels <- as.character(tax.tab[order(rowSums(otu.tab),decreasing = T),][,'Genus'])

#insert others
otu.tab['Other',] <- as.numeric(1 - colSums(otu.tab))
tax.tab <- rbind(tax.tab,paste0('Other ',rank_names(ps.top10)))
rownames(tax.tab) <- rownames(otu.tab)

#build phyloseq
ps.final <- phyloseq(otu_table(otu.tab,taxa_are_rows = T),tax_table(tax.tab),sample_data(sample_data(ps_HMOM_df)))

# plot for treatment at T2
ps.final.melt <- psmelt(ps.final)
ps.final.melt$Genus = as.character(ps.final.melt$Genus)

###Again as done for the control plot, in order to pick up the same order for genera in individual plots I am picking up the levels from the group mean plot and filter/add the genera that are missing after comparing the list of genera in individual plot with the group mean plot above

MOM_genus_t2_levels <- temp.levels[!(temp.levels %in% "Other Genus")]
target_genus = as.character(unique(ps.final.melt$Genus)[!(unique(ps.final.melt$Genus) %in% "Other Genus")])
temp = c(MOM_genus_t2_levels[MOM_genus_t2_levels %in% target_genus], target_genus[!(target_genus %in% MOM_genus_t2_levels)])

ps.final.melt$Genus <- factor(ps.final.melt$Genus,levels = c(temp, "Other Genus"))
levels(ps.final.melt$Genus)
plot_t2_h_genus <- ggplot(data = ps.final.melt,mapping = aes(x = SampleName,y = Abundance,fill=Genus)) + geom_bar(stat = 'identity') + scale_fill_manual(values = color_genus[c(temp, "Other Genus")]) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90), aspect.ratio = 4/5) + labs(x="", y= "Relative abundance")
plot_t2_h_genus


###For T3 by MOM groups
#####Group mean T3, Combined plot###

psef.mom.gglom <- tax_glom(physeq = ps_ef_mom_merge,taxrank = 'Genus')
#convert to relative abundance and select top 10 taxa
psef.mom.gglom.prop <- transform_sample_counts(psef.mom.gglom, function(otu) otu/sum(otu))  

top10 <- names(sort(taxa_sums(psef.mom.gglom.prop), decreasing=TRUE))[1:10]
ps.top10 <- prune_taxa(top10, psef.mom.gglom.prop)

otu.tab <- as.data.frame(t(otu_table(ps.top10)))
tax.tab <- as.matrix(tax_table(ps.top10))

#for ordering purpose
temp.levels <- as.character(tax.tab[order(rowSums(otu.tab),decreasing = T),][,'Genus'])

#insert others
otu.tab['Other',] <- as.numeric(1 - colSums(otu.tab))
tax.tab <- rbind(tax.tab,paste0('Other ',rank_names(ps.top10)))
rownames(tax.tab) <- rownames(otu.tab)

#build phyloseq
ps.final <- phyloseq(otu_table(otu.tab,taxa_are_rows = T),tax_table(tax.tab),sample_data(sample_data(ps_ef_mom_merge)))

# combine plot
ps.final.melt <- psmelt(ps.final)
ps.final.melt$Genus <- factor(ps.final.melt$Genus,levels = c(temp.levels,'Other Genus'))
ps.final.melt$Sample <- factor(ps.final.melt$Sample,levels = c("LMOM", "HMOM"))
plot_t3_mom_g_group <- ggplot(data = ps.final.melt,mapping = aes(x = Sample, y = Abundance,fill=Genus)) + geom_bar(stat = 'identity') + theme(panel.background = element_blank()) + scale_fill_manual(values = color_genus) + theme(legend.position = "right") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 0, hjust= 0.5), aspect.ratio = 4/4) + labs(x="", y= "Relative abundance")
plot_t3_mom_g_group

###################

####per sample bar plot LMOM
ps.glom <- tax_glom(physeq = ps_LMOM_ef,taxrank = 'Genus')
#convert to relative abundance and select top 10 taxa
ps.glom.prop <- transform_sample_counts(ps.glom, function(otu) otu/sum(otu))  

top10 <- names(sort(taxa_sums(ps.glom.prop), decreasing=TRUE))[1:10]
ps.top10 <- prune_taxa(top10, ps.glom.prop)

otu.tab <- as.data.frame(t(otu_table(ps.top10)))
tax.tab <- as.matrix(tax_table(ps.top10))

#for ordering purpose
#temp.levels <- as.character(tax.tab[order(rowSums(otu.tab),decreasing = T),][,'Genus'])

#insert others
otu.tab['Other',] <- as.numeric(1 - colSums(otu.tab))
tax.tab <- rbind(tax.tab,paste0('Other ',rank_names(ps.top10)))
rownames(tax.tab) <- rownames(otu.tab)

#build phyloseq
ps.final <- phyloseq(otu_table(otu.tab,taxa_are_rows = T),tax_table(tax.tab),sample_data(sample_data(ps_LMOM_ef)))

# per sample genus bar plot for LMOM at T3
ps.final.melt <- psmelt(ps.final)
ps.final.melt$Genus = as.character(ps.final.melt$Genus)

###In order to pick up the same order for genera in individual plots I am picking up the levels from the group mean plot and filter/add the genera that are missing after comparing the list of genera in individual plot with the group mean plot above
MOM_genus_t3_levels <- temp.levels[!(temp.levels %in% "Other Genus")]
target_genus = as.character(unique(ps.final.melt$Genus)[!(unique(ps.final.melt$Genus) %in% "Other Genus")])
temp = c(MOM_genus_t3_levels[MOM_genus_t3_levels %in% target_genus], target_genus[!(target_genus %in% MOM_genus_t3_levels)])

ps.final.melt$Genus <- factor(ps.final.melt$Genus,levels = c(temp, "Other Genus"))
levels(ps.final.melt$Genus)
plot_t3_l_genus <- ggplot(data = ps.final.melt,mapping = aes(x = SampleName,y = Abundance,fill=Genus)) + geom_bar(stat = 'identity') + scale_fill_manual(values = color_genus[c(temp, "Other Genus")]) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90), aspect.ratio = 4/3) + labs(x="", y= "Relative abundance")
plot_t3_l_genus


# per sample genus bar plot for HMOM at T3

ps.glom <- tax_glom(physeq = ps_HMOM_ef,taxrank = 'Genus')
#convert to relative abundance and select top 10 taxa
ps.glom.prop <- transform_sample_counts(ps.glom, function(otu) otu/sum(otu))  

top10 <- names(sort(taxa_sums(ps.glom.prop), decreasing=TRUE))[1:10]
ps.top10 <- prune_taxa(top10, ps.glom.prop)

otu.tab <- as.data.frame(t(otu_table(ps.top10)))
tax.tab <- as.matrix(tax_table(ps.top10))

#for ordering purpose. I used it previously to order the taxa according to their abundance.
#temp.levels <- as.character(tax.tab[order(rowSums(otu.tab),decreasing = T),][,'Genus'])

#insert others
otu.tab['Other',] <- as.numeric(1 - colSums(otu.tab))
tax.tab <- rbind(tax.tab,paste0('Other ',rank_names(ps.top10)))
rownames(tax.tab) <- rownames(otu.tab)

#build phyloseq
ps.final <- phyloseq(otu_table(otu.tab,taxa_are_rows = T),tax_table(tax.tab),sample_data(sample_data(ps_HMOM_ef)))

# plot for treatment at T3
ps.final.melt <- psmelt(ps.final)
ps.final.melt$Genus = as.character(ps.final.melt$Genus)

###Again as done for the control plot, in order to pick up the same order for genera in individual plots I am picking up the levels from the group mean plot and filter/add the genera that are missing after comparing the list of genera in individual plot with the group mean plot above

MOM_genus_t3_levels <- temp.levels[!(temp.levels %in% "Other Genus")]
target_genus = as.character(unique(ps.final.melt$Genus)[!(unique(ps.final.melt$Genus) %in% "Other Genus")])
temp = c(MOM_genus_t3_levels[MOM_genus_t3_levels %in% target_genus], target_genus[!(target_genus %in% MOM_genus_t3_levels)])

ps.final.melt$Genus <- factor(ps.final.melt$Genus,levels = c(temp, "Other Genus"))
levels(ps.final.melt$Genus)
plot_t3_h_genus <- ggplot(data = ps.final.melt,mapping = aes(x = SampleName,y = Abundance,fill=Genus)) + geom_bar(stat = 'identity') + scale_fill_manual(values = color_genus[c(temp, "Other Genus")]) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90), aspect.ratio = 4/4) + labs(x="", y= "Relative abundance")
plot_t3_h_genus


#######For T4 by MOM
#####Group mean T4, Combined plot###

psff.mom.gglom <- tax_glom(physeq = ps_ff_mom_merge,taxrank = 'Genus')
#convert to relative abundance and select top 10 taxa
psff.mom.gglom.prop <- transform_sample_counts(psff.mom.gglom, function(otu) otu/sum(otu))  

top10 <- names(sort(taxa_sums(psff.mom.gglom.prop), decreasing=TRUE))[1:10]
ps.top10 <- prune_taxa(top10, psff.mom.gglom.prop)

otu.tab <- as.data.frame(t(otu_table(ps.top10)))
tax.tab <- as.matrix(tax_table(ps.top10))

#for ordering purpose
temp.levels <- as.character(tax.tab[order(rowSums(otu.tab),decreasing = T),][,'Genus'])

#insert others
otu.tab['Other',] <- as.numeric(1 - colSums(otu.tab))
tax.tab <- rbind(tax.tab,paste0('Other ',rank_names(ps.top10)))
rownames(tax.tab) <- rownames(otu.tab)

#build phyloseq
ps.final <- phyloseq(otu_table(otu.tab,taxa_are_rows = T),tax_table(tax.tab),sample_data(sample_data(ps_ff_mom_merge)))

# combine plot
ps.final.melt <- psmelt(ps.final)
ps.final.melt$Genus <- factor(ps.final.melt$Genus,levels = c(temp.levels,'Other Genus'))
ps.final.melt$Sample <- factor(ps.final.melt$Sample,levels = c("LMOM", "HMOM"))
plot_t4_mom_g_group <- ggplot(data = ps.final.melt,mapping = aes(x = Sample, y = Abundance,fill=Genus)) + geom_bar(stat = 'identity') + theme(panel.background = element_blank()) + scale_fill_manual(values = color_genus) + theme(legend.position = "right") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 0, hjust= 0.5), aspect.ratio = 4/4) + labs(x="", y= "Relative abundance")
plot_t4_mom_g_group

###################

####per sample bar plot LMOM
ps.glom <- tax_glom(physeq = ps_LMOM_ff,taxrank = 'Genus')
#convert to relative abundance and select top 10 taxa
ps.glom.prop <- transform_sample_counts(ps.glom, function(otu) otu/sum(otu))  

top10 <- names(sort(taxa_sums(ps.glom.prop), decreasing=TRUE))[1:10]
ps.top10 <- prune_taxa(top10, ps.glom.prop)

otu.tab <- as.data.frame(t(otu_table(ps.top10)))
tax.tab <- as.matrix(tax_table(ps.top10))

#for ordering purpose
#temp.levels <- as.character(tax.tab[order(rowSums(otu.tab),decreasing = T),][,'Genus'])

#insert others
otu.tab['Other',] <- as.numeric(1 - colSums(otu.tab))
tax.tab <- rbind(tax.tab,paste0('Other ',rank_names(ps.top10)))
rownames(tax.tab) <- rownames(otu.tab)

#build phyloseq
ps.final <- phyloseq(otu_table(otu.tab,taxa_are_rows = T),tax_table(tax.tab),sample_data(sample_data(ps_LMOM_ff)))

# per sample genus bar plot for LMOM at T4
ps.final.melt <- psmelt(ps.final)
ps.final.melt$Genus = as.character(ps.final.melt$Genus)

###In order to pick up the same order for genera in individual plots I am picking up the levels from the group mean plot and filter/add the genera that are missing after comparing the list of genera in individual plot with the group mean plot above
MOM_genus_t4_levels <- temp.levels[!(temp.levels %in% "Other Genus")]
target_genus = as.character(unique(ps.final.melt$Genus)[!(unique(ps.final.melt$Genus) %in% "Other Genus")])
temp = c(MOM_genus_t4_levels[MOM_genus_t4_levels %in% target_genus], target_genus[!(target_genus %in% MOM_genus_t4_levels)])

ps.final.melt$Genus <- factor(ps.final.melt$Genus,levels = c(temp, "Other Genus"))
levels(ps.final.melt$Genus)
plot_t4_l_genus <- ggplot(data = ps.final.melt,mapping = aes(x = SampleName,y = Abundance,fill=Genus)) + geom_bar(stat = 'identity') + scale_fill_manual(values = color_genus[c(temp, "Other Genus")]) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90), aspect.ratio = 4/3.5) + labs(x="", y= "Relative abundance")
plot_t4_l_genus


# per sample genus bar plot for HMOM at T4

ps.glom <- tax_glom(physeq = ps_HMOM_ff,taxrank = 'Genus')
#convert to relative abundance and select top 10 taxa
ps.glom.prop <- transform_sample_counts(ps.glom, function(otu) otu/sum(otu))  

top10 <- names(sort(taxa_sums(ps.glom.prop), decreasing=TRUE))[1:10]
ps.top10 <- prune_taxa(top10, ps.glom.prop)

otu.tab <- as.data.frame(t(otu_table(ps.top10)))
tax.tab <- as.matrix(tax_table(ps.top10))

#for ordering purpose. I used it previously to order the taxa according to their abundance.
#temp.levels <- as.character(tax.tab[order(rowSums(otu.tab),decreasing = T),][,'Genus'])

#insert others
otu.tab['Other',] <- as.numeric(1 - colSums(otu.tab))
tax.tab <- rbind(tax.tab,paste0('Other ',rank_names(ps.top10)))
rownames(tax.tab) <- rownames(otu.tab)

#build phyloseq
ps.final <- phyloseq(otu_table(otu.tab,taxa_are_rows = T),tax_table(tax.tab),sample_data(sample_data(ps_HMOM_ff)))

# plot for treatment at T4
ps.final.melt <- psmelt(ps.final)
ps.final.melt$Genus = as.character(ps.final.melt$Genus)

###Again as done for the control plot, in order to pick up the same order for genera in individual plots I am picking up the levels from the group mean plot and filter/add the genera that are missing after comparing the list of genera in individual plot with the group mean plot above

MOM_genus_t4_levels <- temp.levels[!(temp.levels %in% "Other Genus")]
target_genus = as.character(unique(ps.final.melt$Genus)[!(unique(ps.final.melt$Genus) %in% "Other Genus")])
temp = c(MOM_genus_t4_levels[MOM_genus_t4_levels %in% target_genus], target_genus[!(target_genus %in% MOM_genus_t4_levels)])

ps.final.melt$Genus <- factor(ps.final.melt$Genus,levels = c(temp, "Other Genus"))
levels(ps.final.melt$Genus)
plot_t4_h_genus <- ggplot(data = ps.final.melt,mapping = aes(x = SampleName,y = Abundance,fill=Genus)) + geom_bar(stat = 'identity') + scale_fill_manual(values = color_genus[c(temp, "Other Genus")]) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 90), aspect.ratio = 4/4) + labs(x="", y= "Relative abundance")
plot_t4_h_genus


####Combining all plots for MOM groups
ggarrange(plotlist = list(plot_t3_l_genus, plot_t3_h_genus, plot_t3_mom_g_group, plot_t2_l_genus, plot_t2_h_genus, plot_t2_mom_g_group, plot_t3_l_genus, plot_t3_h_genus, plot_t3_mom_g_group, plot_t4_l_genus, plot_t4_h_genus, plot_t4_mom_g_group ), ncol=3, nrow= 4, legend = "none")

#group combined genus level bar plots

#For first time point

psbf.gglom <- tax_glom(physeq = ps_bf_total_merge,taxrank = 'Genus')
#convert to relative abundance and select top 10 taxa
psbf.gglom.prop <- transform_sample_counts(psbf.gglom, function(otu) otu/sum(otu))  

top10 <- names(sort(taxa_sums(psbf.gglom.prop), decreasing=TRUE))[1:10]
ps.top10 <- prune_taxa(top10, psbf.gglom.prop)

otu.tab <- as.data.frame(t(otu_table(ps.top10)))
tax.tab <- as.matrix(tax_table(ps.top10))

#for ordering purpose
temp.levels <- as.character(tax.tab[order(rowSums(otu.tab),decreasing = T),][,'Genus'])

#insert others
otu.tab['Other',] <- as.numeric(1 - colSums(otu.tab))
tax.tab <- rbind(tax.tab,paste0('Other ',rank_names(ps.top10)))
rownames(tax.tab) <- rownames(otu.tab)

#build phyloseq
ps.final <- phyloseq(otu_table(otu.tab,taxa_are_rows = T),tax_table(tax.tab),sample_data(sample_data(ps_bf_total_merge)))

# plot2
ps.final.melt <- psmelt(ps.final)
ps.final.melt$Genus <- factor(ps.final.melt$Genus,levels = c(temp.levels,'Other Genus'))
plot_t1_genus <- ggplot(data = ps.final.melt,mapping = aes(x = Sample, y = Abundance,fill=Genus)) + geom_bar(stat = 'identity') + theme(panel.background = element_blank()) + scale_fill_manual(values=c("#F7D363", "#94D5E0", "#20899F", "#F0CEC6", "#D2AF96", "#C3EB95", "#EAF0F3", "#935272", "#F35733", "#5B724B", "#000000")) + theme(legend.position = "right") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 0, hjust= 0.5))
plot_t1g<- plot_t1_genus + labs(x="", y= "Relative abundance")
plot_t1g
#rct_genus_t1_levels <- levels(ps.final.melt$Genus)


#For second time point
psdf.gglom <- tax_glom(physeq = ps_df_total_merge,taxrank = 'Genus')
#convert to relative abundance and select top 10 taxa
psdf.gglom.prop <- transform_sample_counts(psdf.gglom, function(otu) otu/sum(otu))  

top10 <- names(sort(taxa_sums(psdf.gglom.prop), decreasing=TRUE))[1:10]
ps.top10 <- prune_taxa(top10, psdf.gglom.prop)

otu.tab <- as.data.frame(t(otu_table(ps.top10)))
tax.tab <- as.matrix(tax_table(ps.top10))

#for ordering purpose
temp.levels <- as.character(tax.tab[order(rowSums(otu.tab),decreasing = T),][,'Genus'])

#insert others
otu.tab['Other',] <- as.numeric(1 - colSums(otu.tab))
tax.tab <- rbind(tax.tab,paste0('Other ',rank_names(ps.top10)))
rownames(tax.tab) <- rownames(otu.tab)

#build phyloseq
ps.final <- phyloseq(otu_table(otu.tab,taxa_are_rows = T),tax_table(tax.tab),sample_data(sample_data(ps_df_total_merge)))

# plot2
ps.final.melt <- psmelt(ps.final)
ps.final.melt$Genus <- factor(ps.final.melt$Genus,levels = c(temp.levels,'Other Genus'))
plot_t2_genus <- ggplot(data = ps.final.melt,mapping = aes(x = Sample, y = Abundance,fill=Genus)) + geom_bar(stat = 'identity') + theme(panel.background = element_blank()) + scale_fill_manual(values=c("#20899F", "#F7D363", "#D2AF96", "#94D5E0", "#044B88", "#F0CEC6", "#C3EB95", "#C29F00", "#935272", "#B56602", "#000000")) + theme(legend.position = "right") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 0, hjust= 0.5))
plot_t2g<- plot_t2_genus + labs(x="", y= "Relative abundance")
plot_t2g


##For third time point
psef.gglom <- tax_glom(physeq = ps_ef_total_merge,taxrank = 'Genus')
#convert to relative abundance and select top 10 taxa
psef.gglom.prop <- transform_sample_counts(psef.gglom, function(otu) otu/sum(otu))  

top10 <- names(sort(taxa_sums(psef.gglom.prop), decreasing=TRUE))[1:10]
ps.top10 <- prune_taxa(top10, psef.gglom.prop)

otu.tab <- as.data.frame(t(otu_table(ps.top10)))
tax.tab <- as.matrix(tax_table(ps.top10))

#for ordering purpose
temp.levels <- as.character(tax.tab[order(rowSums(otu.tab),decreasing = T),][,'Genus'])

#insert others
otu.tab['Other',] <- as.numeric(1 - colSums(otu.tab))
tax.tab <- rbind(tax.tab,paste0('Other ',rank_names(ps.top10)))
rownames(tax.tab) <- rownames(otu.tab)

#build phyloseq
ps.final <- phyloseq(otu_table(otu.tab,taxa_are_rows = T),tax_table(tax.tab),sample_data(sample_data(ps_ef_total_merge)))

# plot2
ps.final.melt <- psmelt(ps.final)
ps.final.melt$Genus <- factor(ps.final.melt$Genus,levels = c(temp.levels,'Other Genus'))
plot_t3_genus <- ggplot(data = ps.final.melt,mapping = aes(x = Sample, y = Abundance,fill=Genus)) + geom_bar(stat = 'identity') + theme(panel.background = element_blank()) + scale_fill_manual(values=c("#D2AF96", "#20899F", "#044B88", "#F7D363", "#B56602", "#94D5E0", "#F0CEC6", "#653723", "#C3EB95", "#C29F00", "#000000")) + theme(legend.position = "right") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 0, hjust= 0.5))
plot_t3g<- plot_t3_genus + labs(x="", y= "Relative abundance")
plot_t3g


##For fourth time point
psff.gglom <- tax_glom(physeq = ps_ff_total_merge,taxrank = 'Genus')
#convert to relative abundance and select top 10 taxa
psff.gglom.prop <- transform_sample_counts(psff.gglom, function(otu) otu/sum(otu))  

top10 <- names(sort(taxa_sums(psff.gglom.prop), decreasing=TRUE))[1:10]
ps.top10 <- prune_taxa(top10, psff.gglom.prop)

otu.tab <- as.data.frame(t(otu_table(ps.top10)))
tax.tab <- as.matrix(tax_table(ps.top10))

#for ordering purpose
temp.levels <- as.character(tax.tab[order(rowSums(otu.tab),decreasing = T),][,'Genus'])

#insert others
otu.tab['Other',] <- as.numeric(1 - colSums(otu.tab))
tax.tab <- rbind(tax.tab,paste0('Other ',rank_names(ps.top10)))
rownames(tax.tab) <- rownames(otu.tab)

#build phyloseq
ps.final <- phyloseq(otu_table(otu.tab,taxa_are_rows = T),tax_table(tax.tab),sample_data(sample_data(ps_ff_total_merge)))

# plot2
ps.final.melt <- psmelt(ps.final)
ps.final.melt$Genus <- factor(ps.final.melt$Genus,levels = c(temp.levels,'Other Genus'))
plot_t4_genus <- ggplot(data = ps.final.melt,mapping = aes(x = Sample, y = Abundance,fill=Genus)) + geom_bar(stat = 'identity') + theme(panel.background = element_blank()) + scale_fill_manual(values=c("#D2AF96", "#20899F", "#94D5E0", "#044B88", "#B56602", "#F0CEC6", "#C3EB95", "#C29F00", "#ca624b", "#53565A", "#000000")) + theme(legend.position = "right") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 0, hjust= 0.5))
plot_t4g<- plot_t4_genus + labs(x="", y= "Relative abundance")
plot_t4g

#combining plots


genus_bar_plot <- ggarrange(plot_t1g, plot_t2g, plot_t3g, plot_t4g,
                            common.legend = FALSE, legend = "none",
                            ncol = 2, nrow = 2)
genus_bar_plot

################################Comaprison of alpha diversity between RCT groups#####################


##Normal Repeated measures analysis https://stats.idre.ucla.edu/r/seminars/repeated-measures-analysis-with-r/

###For alpha diversity analysis repeated measures analysis using unrarified dataset
#For repeated measures analysis removed subjects that did not have samples at all time points

psnoncontam_tru_final_unrarified

psnoncontam_tru_final_unrarified_rm <- subset_samples(psnoncontam_tru_final_unrarified, sample_data(psnoncontam_tru_final_unrarified)$SubjectID!= "W28" & sample_data(psnoncontam_tru_final_unrarified)$SubjectID!= "W30" & sample_data(psnoncontam_tru_final_unrarified)$SubjectID!= "W06")


#filtering out taxa not present in the true samples
psnoncontam_tru_final_unrarified_rm <- prune_taxa(taxa_sums(psnoncontam_tru_final_unrarified_rm)>0, psnoncontam_tru_final_unrarified_rm)
psnoncontam_tru_final_unrarified_rm


##First trying out repreated measures analysis on the alpha diversity indices between RCT groups at four time points


###Alpha diversity estimates
psnoncontam_tru_final_unrarified_rm

ps_total_div_all <- estimate_richness(psnoncontam_tru_final_unrarified_rm)
write.csv(ps_total_div_all, "alpha_div_pstotal_revised_new.csv")

##Added other metadata (Group and timepoint info) in the file and imported it as a data frame
pstotal.div <- read.table(file = "alpha_div_pstotal_revised_new_March2021.csv", header = TRUE, sep = ",")
pstotaldiv.new <- as.data.frame(pstotal.div)
pstotaldiv.new


###First for Chao1 ASVs

##Post hoc analysis with Wilcox t-test (For group)
stat.test.observed <- pstotaldiv.new %>%
  group_by(TimePoint) %>%
  wilcox_test(data =., Chao1 ~ Group, paired = FALSE) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "TimePoint") %>%
  p_round(digits = 2)
stat.test.observed

##Post hoc analysis with Wilcox t-test (For TimePoint)
stat.test.time <- pstotaldiv.new %>%
  group_by(Group) %>%
  wilcox_test(data =., Chao1 ~ TimePoint, paired=TRUE) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "TimePoint") %>%
  p_round(digits = 2)
stat.test.time

###Plot the data with the stats
observed_plot <- ggplot(pstotaldiv.new,aes(x=TimePoint, y=Chao1)) +
  geom_boxplot(aes(fill=Group), width=0.6) + scale_fill_manual(values = c("#80D7DD", "#F3DB80")) +
  geom_point(aes(color = Group), position = position_dodge(0.6)) + stat_summary(
    fun = median,
    geom = 'line',
    aes(group = Group, colour = Group),
    position = position_dodge(width = 0.6) #this has to be added
  ) +
  scale_color_manual(values = c("#354360", "#7D6608")) + xlab("Time points") +
  ylab("Chao1 estimate")
observed_plot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))


#For shannon index
##Post hoc analysis with wilcox test (for group)
stat.test.shannon <- pstotaldiv.new %>%
  group_by(TimePoint) %>%
  wilcox_test(data =., Shannon ~ Group, paired=FALSE) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "TimePoint") %>%
  p_round(digits = 2)
stat.test.shannon

##Post hoc analysis with wilcox test (for time)

stat.test.time <- pstotaldiv.new %>%
  group_by(Group) %>%
  wilcox_test(data =., Shannon ~ TimePoint, paired=TRUE) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "TimePoint") %>%
  p_round(digits = 2)
stat.test.time

###Plot the data with the stats
Shannon_plot <- ggplot(pstotaldiv.new,aes(x=TimePoint, y=Shannon)) +
  geom_boxplot(aes(fill=Group), width=0.6) + scale_fill_manual(values = c("#80D7DD", "#F3DB80")) +
  geom_point(aes(color = Group), position = position_dodge(0.6)) + stat_summary(
    fun = median,
    geom = 'line',
    aes(group = Group, colour = Group),
    position = position_dodge(width = 0.6) #this has to be added
  ) +
  scale_color_manual(values = c("#354360", "#7D6608")) + xlab("Time points") +
  ylab("Shannon index")
shannon_plot <- Shannon_plot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

###Splinectome (longitudinal analysis) for shannon index
##replacing the alpha-numeric characters with numbers instead for the timepoints, as splinectome needs continous/numeric scale variable.

pstotaldiv.new$TimePoint <- as.numeric(gsub(pattern = 'T',replacement = '',x = pstotaldiv.new$TimePoint))


shannon.rct.result <- permuspliner(data = pstotaldiv.new, xvar = 'TimePoint',yvar = 'Shannon', perms = 999,category = 'Group', cases = 'SampleID',quiet = T)

p.shannon.rct <- permuspliner.plot.permsplines(data = shannon.rct.result,
                                               xvar = 'timepoint',
                                               yvar = 'Shannon index')
p.shannon.rct <- p.shannon.rct + scale_colour_manual(values = c("#80D7DD", "#F3DB80")) + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

p.shannon.rct

shannon_rct_dist <- permuspliner.plot.permdistance(shannon.rct.result, xlabel = 'timepoint') + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

shannon_rct_dist

shannon.rct.result$pval

#for InvSimpson
##Post hoc analysis with wilcox sum rank t-test for the group factor
stat.test.invsimp <- pstotaldiv.new %>%
  group_by(TimePoint) %>%
  wilcox_test(data =., InvSimpson ~ Group, paired=FALSE) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "TimePoint") %>%
  p_round(digits = 2)
stat.test.invsimp

##Post hoc analysis with wilcox test (for time)

stat.test.time <- pstotaldiv.new %>%
  group_by(Group) %>%
  wilcox_test(data =., InvSimpson ~ TimePoint, paired=TRUE) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "TimePoint") %>%
  p_round(digits = 2)
stat.test.time


###Plot the data with the stats
InvSimp_plot <- ggplot(pstotaldiv.new,aes(x=TimePoint, y=InvSimpson)) +
  geom_boxplot(aes(fill=Group), width=0.6) + scale_fill_manual(values = c("#80D7DD", "#F3DB80")) +
  geom_point(aes(color = Group), position = position_dodge(0.6)) + stat_summary(
    fun = median,
    geom = 'line',
    aes(group = Group, colour = Group),
    position = position_dodge(width = 0.6) #this has to be added
  ) +
  scale_color_manual(values = c("#354360", "#7D6608")) + xlab("Time points") +
  ylab("InvSimpson index")
InvSimp_plot <- InvSimp_plot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

####Longitudinal analysis for Inverse Simpson index
##replacing the alpha-numeric characters with whole numbers instead for the timepoints, as splinectome needs continous/numeric scale variable.

pstotaldiv.new$TimePoint <- as.numeric(gsub(pattern = 'T',replacement = '',x = pstotaldiv.new$TimePoint))


invsimp.rct.result <- permuspliner(data = pstotaldiv.new, xvar = 'TimePoint',yvar = 'InvSimpson', perms = 999,category = 'Group', cases = 'SampleID',quiet = T)

p.invsimp.rct <- permuspliner.plot.permsplines(data = invsimp.rct.result,
                                               xvar = 'timepoint',
                                               yvar = 'Inverse Simpson index')
p.invsimp.rct <- p.invsimp.rct + scale_colour_manual(values = c("#80D7DD", "#F3DB80")) + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

p.invsimp.rct

invsimp_rct_dist <- permuspliner.plot.permdistance(invsimp.rct.result, xlabel = 'timepoint') + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

invsimp_rct_dist

invsimp.rct.result$pval

####Combining these images

plot_alpha_rct_trend <- ggarrange(plotlist = list(shannon_plot, p.shannon.rct, shannon_rct_dist, InvSimp_plot, p.invsimp.rct, invsimp_rct_dist), ncol = 3, nrow = 2, legend = "none")
plot_alpha_rct_trend

##Note: These alpha diversity box plots were not included in manuscript, however splienctome plots for shannon and inverse simpson were added

###########################Beta diversity analysis for RCT group and microbiome######################


##Using relative abundances with clr transformation prior to PCoA (Ref: https://msystems.asm.org/content/msys/4/3/e00036-39.full.pdf and https://www.nature.com/articles/s43467-020-36224-6.pdf)
#Using rarified phyloseq for CLR transformation on count data

###Time point T1

#Write the phyloseq object as a CSV table with ASVs and counts
ps_bf_total
write_phyloseq(ps_bf_total)

#Read the CSV file
abund_table <- read.csv("otu_table_bf_total_rare_revised_new.csv", row.names=1, check.names= FALSE) #Added header "Features" to the first clumn of ASVs
abund_table_t<-t(abund_table)
#View(abund_table_t)

#install.packages("zCompositions")
#library (zCompositions)
abund_table_r <- t(cmultRepl((abund_table_t), method="CZM", output="prop")) #prop calculates the imputed proportion is similar to "counts" and then taking proportions in the previous version; CZM: count zero multiplicative

#Check the before and after files
#View(abund_table)
#View(abund_table_r)

#Apply CLR transformation
abund_clr <- t(apply(abund_table_r, 2, function(x){log(x) - mean (log(x))}))
#View(abund_clr)

#Replacing original counts with CLR transformed values
ps_bf_total_rare_zclr <- ps_bf_total
otu_table(ps_bf_total_rare_zclr) <- otu_table(abund_clr, taxa_are_rows = F)
ps_bf_total_rare_zclr
#write_phyloseq(ps_bf_zclr) Just checking if the replacement went well

#CLR with Bray distance (adding a pseudocount)

ps_bf_total_rare_zclr
ps_bf_total_zclr_pseudo <- ps_bf_total_rare_zclr
otu_table(ps_bf_total_zclr_pseudo) <- otu_table(ps_bf_total_zclr_pseudo) + 5
ps_bftotal_clr_plot_bray <- plot_ordination(ps_bf_total_zclr_pseudo, ordinate(ps_bf_total_zclr_pseudo, "MDS", "bray"), color = "Group") + geom_point(size = 3)

ps_bftotal_clr_plot_bray + scale_color_manual(values=c("#00AFBB", "#E7B800")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + stat_ellipse(type = "norm")

#comparing two groups
metadata3 <- as(sample_data(ps_bf_total_zclr_pseudo), "data.frame")
adonis(phyloseq::distance(ps_bf_total_zclr_pseudo, method = "bray")~Group, data = metadata3, p.adjust.methods= "fdr", permutations = 99999)

#beta-dispersion test
ps_total_bf_dist <- phyloseq::distance(ps_bf_total_zclr_pseudo, method = "bray")
sampledf3 <- data.frame(sample_data(ps_bf_total_zclr_pseudo))
beta3 <- betadisper(ps_total_bf_dist, sampledf3$Group, bias.adjust = TRUE)
permutest(beta3, p.adjust.methods= "fdr", permutations=99999)


#Jaccard (presence/absence)

ps_bftotal_plot_jaccard <- plot_ordination(ps_bf_total, ordinate(ps_bf_total, "MDS", "jaccard"), color = "Group") + geom_point(size = 3)

ps_bftotal_plot_jaccard + scale_color_manual(values=c("#00AFBB", "#E7B800")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + stat_ellipse(type = "norm")

#comparing two groups
metadata3 <- as(sample_data(ps_bf_total), "data.frame")
adonis(phyloseq::distance(ps_bf_total, method = "jaccard")~Group, data = metadata3, p.adjust.methods= "fdr", permutations = 99999)

#beta-dispersion test
ps_total_bf_dist <- phyloseq::distance(ps_bf_total, method = "jaccard")
sampledf3 <- data.frame(sample_data(ps_bf_total))
beta3 <- betadisper(ps_total_bf_dist, sampledf3$Group, bias.adjust = TRUE)
permutest(beta3, p.adjust.methods= "fdr", permutations=99999)

####Second time point

#CLR transformation on relative abundances

#Write the phyloseq object as a CSV table with ASVs and counts
ps_df_total
write_phyloseq(ps_df_total)

#Read the CSV file
abund_table <- read.csv("otu_table_df_total_rare_revised_new.csv", row.names=1, check.names= FALSE) #Added header "Features" to the first clumn of ASVs
abund_table_t<-t(abund_table)
#View(abund_table_t)

#install.packages("zCompositions")
#library (zCompositions)
abund_table_r <- t(cmultRepl((abund_table_t), method="CZM", output="prop")) #prop calculates the imputed proportion is similar to "counts" and then taking proportions in the previous version; CZM: count zero multiplicative

#Check the before and after files
#View(abund_table)
#View(abund_table_r)

#Apply CLR transformation
abund_clr <- t(apply(abund_table_r, 2, function(x){log(x) - mean (log(x))}))
#View(abund_clr)

#Replacing original counts with CLR transformed values
ps_df_total_zclr <- ps_df_total
otu_table(ps_df_total_zclr) <- otu_table(abund_clr, taxa_are_rows = F)
ps_df_total_zclr
#write_phyloseq(ps_bf_zclr) Just checking if the replacement went well

#CLR with Bray distance

ps_df_total_zclr
ps_df_total_zclr_pseudo <- ps_df_total_zclr
otu_table(ps_df_total_zclr_pseudo) <- otu_table(ps_df_total_zclr_pseudo) + 5
ps_dftotal_clr_plot_bray <- plot_ordination(ps_df_total_zclr_pseudo, ordinate(ps_df_total_zclr_pseudo, "MDS", "bray"), color = "Group") + geom_point(size = 3)

ps_dftotal_clr_plot_bray + scale_color_manual(values=c("#00AFBB", "#E7B800")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + stat_ellipse(type = "t")

#comparing two groups
metadata3 <- as(sample_data(ps_df_total_zclr_pseudo), "data.frame")
adonis(phyloseq::distance(ps_df_total_zclr_pseudo, method = "bray")~Group, data = metadata3, p.adjust.methods= "fdr", permutations = 99999)

#beta-dispersion test
ps_total_df_dist <- phyloseq::distance(ps_df_total_zclr_pseudo, method = "bray")
sampledf3 <- data.frame(sample_data(ps_df_total_zclr_pseudo))
beta3 <- betadisper(ps_total_df_dist, sampledf3$Group, bias.adjust = TRUE)
permutest(beta3, p.adjust.methods= "fdr", permutations=99999)

#Jaccard (presence/absence)

ps_dftotal_plot_jaccard <- plot_ordination(ps_df_total, ordinate(ps_df_total, "MDS", "jaccard"), color = "Group") + geom_point(size = 3)

ps_dftotal_plot_jaccard + scale_color_manual(values=c("#00AFBB", "#E7B800")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + stat_ellipse(type = "t")


#comparing two groups
metadata3 <- as(sample_data(ps_df_total), "data.frame")
adonis(phyloseq::distance(ps_df_total, method = "jaccard")~Group, data = metadata3, p.adjust.methods= "fdr", permutations = 99999)

#beta-dispersion test
ps_total_df_dist <- phyloseq::distance(ps_df_total, method = "jaccard")
sampledf3 <- data.frame(sample_data(ps_df_total))
beta3 <- betadisper(ps_total_df_dist, sampledf3$Group, bias.adjust = TRUE)
permutest(beta3, p.adjust.methods= "fdr", permutations=99999)


#####Third time point

#Using rarified phyloseq for CLR transformation on count data

#Write the phyloseq object as a CSV table with ASVs and counts
ps_ef_total
write_phyloseq(ps_ef_total)

#Read the CSV file
abund_table <- read.csv("otu_table_ef_total_rare_revised_new.csv", row.names=1, check.names= FALSE) #Added header "Features" to the first clumn of ASVs
abund_table_t<-t(abund_table)
#View(abund_table_t)

#install.packages("zCompositions")
#library (zCompositions)
abund_table_r <- t(cmultRepl((abund_table_t), method="CZM", output="prop")) #prop calculates the imputed proportion is similar to "counts" and then taking proportions in the previous version; CZM: count zero multiplicative

#Check the before and after files
#View(abund_table)
#View(abund_table_r)

#Apply CLR transformation
abund_clr <- t(apply(abund_table_r, 2, function(x){log(x) - mean (log(x))}))
#View(abund_clr)

#Replacing original counts with CLR transformed values
ps_ef_total_rare_zclr <- ps_ef_total
otu_table(ps_ef_total_rare_zclr) <- otu_table(abund_clr, taxa_are_rows = F)
ps_ef_total_rare_zclr
#write_phyloseq(ps_bf_zclr) Just checking if the replacement went well

#CLR with Bray distance (adding a pseudocount)

ps_ef_total_rare_zclr
ps_ef_total_zclr_pseudo <- ps_ef_total_rare_zclr
otu_table(ps_ef_total_zclr_pseudo) <- otu_table(ps_ef_total_zclr_pseudo) + 5
ps_eftotal_clr_plot_bray <- plot_ordination(ps_ef_total_zclr_pseudo, ordinate(ps_ef_total_zclr_pseudo, "MDS", "bray"), color = "Group") + geom_point(size = 3)

ps_eftotal_clr_plot_bray + scale_color_manual(values=c("#00AFBB", "#E7B800")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + stat_ellipse(type = "norm")

#comparing two groups
metadata3 <- as(sample_data(ps_ef_total_zclr_pseudo), "data.frame")
adonis(phyloseq::distance(ps_ef_total_zclr_pseudo, method = "bray")~Group, data = metadata3, p.adjust.methods= "fdr", permutations = 99999)

#beta-dispersion test
ps_total_ef_dist <- phyloseq::distance(ps_ef_total_zclr_pseudo, method = "bray")
sampledf3 <- data.frame(sample_data(ps_ef_total_zclr_pseudo))
beta3 <- betadisper(ps_total_ef_dist, sampledf3$Group, bias.adjust = TRUE)
permutest(beta3, p.adjust.methods= "fdr", permutations=99999)


#Jaccard (presence/absence)


ps_eftotal_plot_jaccard <- plot_ordination(ps_ef_total, ordinate(ps_ef_total, "MDS", "jaccard"), color = "Group") + geom_point(size = 3)

ps_eftotal_plot_jaccard + scale_color_manual(values=c("#00AFBB", "#E7B800")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + stat_ellipse(type = "norm")

#comparing two groups
metadata3 <- as(sample_data(ps_ef_total), "data.frame")
adonis(phyloseq::distance(ps_ef_total, method = "jaccard")~Group, data = metadata3, p.adjust.methods= "fdr", permutations = 99999)

#beta-dispersion test
ps_total_ef_dist <- phyloseq::distance(ps_ef_total, method = "jaccard")
sampledf3 <- data.frame(sample_data(ps_ef_total))
beta3 <- betadisper(ps_total_ef_dist, sampledf3$Group, bias.adjust = TRUE)
permutest(beta3, p.adjust.methods= "fdr", permutations=99999)


#####fourth time point

#Using rarified phyloseq for CLR transformation on count data

#Write the phyloseq object as a CSV table with ASVs and counts
ps_ff_total
write_phyloseq(ps_ff_total)

#Read the CSV file
abund_table <- read.csv("otu_table_ff_total_rare_revised.csv", row.names=1, check.names= FALSE) #Added header "Features" to the first clumn of ASVs
abund_table_t<-t(abund_table)
#View(abund_table_t)

#install.packages("zCompositions")
#library (zCompositions)
abund_table_r <- t(cmultRepl((abund_table_t), method="CZM", output="prop")) #prop calculates the imputed proportion is similar to "counts" and then taking proportions in the previous version; CZM: count zero multiplicative

#Check the before and after files
#View(abund_table)
#View(abund_table_r)

#Apply CLR transformation
abund_clr <- t(apply(abund_table_r, 2, function(x){log(x) - mean (log(x))}))
#View(abund_clr)

#Replacing original counts with CLR transformed values
ps_ff_total_rare_zclr <- ps_ff_total
otu_table(ps_ff_total_rare_zclr) <- otu_table(abund_clr, taxa_are_rows = F)
ps_ff_total_rare_zclr
#write_phyloseq(ps_bf_zclr) Just checking if the replacement went well

#CLR with Bray distance (adding a pseudocount)

ps_ff_total_rare_zclr
ps_ff_total_zclr_pseudo <- ps_ff_total_rare_zclr
otu_table(ps_ff_total_zclr_pseudo) <- otu_table(ps_ff_total_zclr_pseudo) + 5
ps_fftotal_clr_plot_bray <- plot_ordination(ps_ff_total_zclr_pseudo, ordinate(ps_ff_total_zclr_pseudo, "MDS", "bray"), color = "Group") + geom_point(size = 3)

ps_fftotal_clr_plot_bray + scale_color_manual(values=c("#00AFBB", "#E7B800")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + stat_ellipse(type = "norm")

#comparing two groups
metadata3 <- as(sample_data(ps_ff_total_zclr_pseudo), "data.frame")
adonis(phyloseq::distance(ps_ff_total_zclr_pseudo, method = "bray")~Group, data = metadata3, p.adjust.methods= "fdr", permutations = 99999)

#beta-dispersion test
ps_total_ff_dist <- phyloseq::distance(ps_ff_total_zclr_pseudo, method = "bray")
sampledf3 <- data.frame(sample_data(ps_ff_total_zclr_pseudo))
beta3 <- betadisper(ps_total_ff_dist, sampledf3$Group, bias.adjust = TRUE)
permutest(beta3, p.adjust.methods= "fdr", permutations=99999)

#Jaccard (presence/absence)

ps_fftotal_plot_jaccard <- plot_ordination(ps_ff_total, ordinate(ps_ff_total, "MDS", "jaccard"), color = "Group") + geom_point(size = 3)

ps_fftotal_plot_jaccard + scale_color_manual(values=c("#00AFBB", "#E7B800")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + stat_ellipse(type = "norm")

#comparing two groups
metadata3 <- as(sample_data(ps_ff_total), "data.frame")
adonis(phyloseq::distance(ps_ff_total, method = "jaccard")~Group, data = metadata3, p.adjust.methods= "fdr", permutations = 99999)

#beta-dispersion test
ps_total_ff_dist <- phyloseq::distance(ps_ff_total, method = "jaccard")
sampledf3 <- data.frame(sample_data(ps_ff_total))
beta3 <- betadisper(ps_total_ff_dist, sampledf3$Group, bias.adjust = TRUE)
permutest(beta3, p.adjust.methods= "fdr", permutations=99999)

###No differences in RCT groups seen in beta-diversity analysis in any time points, except marginal difference at T2 for beta dispersion.


################## MOM and Microbiome analysis #######################

#######################beta diversity with multivariate analysis####################################

#For first time point T1

#Bray curtis

#unconstrained ordination
psnewt1_plot <- plot_ordination(ps_bf_total_zclr_pseudo, ordinate(ps_bf_total_zclr_pseudo, "MDS", "bray"), color = "MOM_Group_T1", shape = "Group") + geom_point(size = 3)

psnewt1_plot + scale_color_manual(values=c("#DC354B", "#A9BFCC")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position = "right")  + scale_shape_manual(values=c(1, 16)) + stat_ellipse(type = "norm", linetype =1, aes(group= MOM_Group_T1))


##PERMANOVA using bray curtis

#Multivariable analysis
metadata35 <- as(sample_data(ps_bf_total_zclr_pseudo), "data.frame")
adonis(phyloseq::distance(ps_bf_total_zclr_pseudo, method = "bray") ~ MOM_Group_T1 + Group + MOM_Group_T1*Group, data = metadata35, p.adjust.methods= "fdr", permutations = 99999)


##Jaccard distances for T1

#unconstrained ordination
psnewt1_plot <- plot_ordination(ps_bf_total, ordinate(ps_bf_total, "MDS", "jaccard"), color = "MOM_Group_T1", shape = "Group") + geom_point(size = 3)

psnewt1_plot + scale_color_manual(values=c("#DC354B", "#A9BFCC")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position = "right")  + scale_shape_manual(values=c(1, 16)) + stat_ellipse(type = "norm", linetype =3, aes(group= MOM_Group_1))


#Multivariable analysis
metadata35 <- as(sample_data(ps_bf_total), "data.frame")
adonis(phyloseq::distance(ps_bf_total, method = "jaccard") ~ MOM_Group_T1 + Group + MOM_Group_T1*Group, data = metadata35, p.adjust.methods= "fdr", permutations = 99999)


# for T2

#unconstrained ordination
psnewt2_plot <- plot_ordination(ps_df_total_zclr_pseudo, ordinate(ps_df_total_zclr_pseudo, "MDS", "bray"), color = "MOM_Group_T2", shape = "Group") + geom_point(size = 3)

psnewt2_plot + scale_color_manual(values=c("#DC354B", "#A9BFCC")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position = "right")  + scale_shape_manual(values=c(1, 16)) + stat_ellipse(type = "norm", linetype =3, aes(group= MOM_Group_T2))


##PERMANOVA using bray curtis
#Multivariable analysis
metadata35 <- as(sample_data(ps_df_total_zclr_pseudo), "data.frame")
adonis(phyloseq::distance(ps_df_total_zclr_pseudo, method = "bray") ~ MOM_Group_T2 + Group + MOM_Group_T2*Group, data = metadata35, p.adjust.methods= "fdr", permutations = 99999)


##With Jaccard for T2
#unconstrained ordination
psnewt2_plot <- plot_ordination(ps_df_total, ordinate(ps_df_total_zclr_pseudo, "MDS", "jaccard"), color = "MOM_Group_T2", shape = "Group") + geom_point(size = 3)

psnewt2_plot + scale_color_manual(values=c("#DC354B", "#A9BFCC")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position = "right")  + scale_shape_manual(values=c(1, 16)) + stat_ellipse(type = "norm", linetype =3, aes(group= MOM_Group_T2))

#Multivariable analysis
metadata35 <- as(sample_data(ps_df_total), "data.frame")
adonis(phyloseq::distance(ps_df_total, method = "jaccard") ~ MOM_Group_T2 + Group + MOM_Group_T2*Group, data = metadata35, p.adjust.methods= "fdr", permutations = 99999)

#For T3

#unconstrained ordination
psnewt3_plot <- plot_ordination(ps_ef_total_zclr_pseudo, ordinate(ps_ef_total_zclr_pseudo, "MDS", "bray"), color = "MOM_Group_T3", shape = "Group") + geom_point(size = 3)
psnewt3_plot + scale_color_manual(values=c("#DC354B", "#A9BFCC")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_shape_manual(values=c(1, 16)) + theme(legend.position = "right") + stat_ellipse(type = "norm", linetype =3, aes(group= MOM_Group))

##PERMANOVA using bray curtis

#Multivariable analysis
metadata36 <- as(sample_data(ps_ef_total_zclr_pseudo), "data.frame")
adonis(phyloseq::distance(ps_ef_total_zclr_pseudo, method = "bray") ~ MOM_Group_T3 + Group + MOM_Group_T3*Group, data = metadata36, p.adjust.methods= "fdr", permutations = 99999)

#With Jaccard for T3

#unconstrained ordination
psnewt3_plot <- plot_ordination(ps_ef_total, ordinate(ps_ef_total, "MDS", "jaccard"), color = "MOM_Group_T3", shape = "Group") + geom_point(size = 3)

psnewt3_plot + scale_color_manual(values=c("#DC354B", "#A9BFCC")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_shape_manual(values=c(1, 16)) + theme(legend.position = "right") + stat_ellipse(type = "norm", linetype =3, aes(group= MOM_Group_T3))

#Multivariable analysis
metadata36 <- as(sample_data(ps_ef_total), "data.frame")
adonis(phyloseq::distance(ps_ef_total, method = "jaccard") ~ MOM_Group_T3 + Group + MOM_Group_T3*Group, data = metadata36, p.adjust.methods= "fdr", permutations = 99999)

##For T4

#unconstrained ordination
psnewt4_plot <- plot_ordination(ps_ff_total_zclr_pseudo, ordinate(ps_ff_total_zclr_pseudo, "MDS", "bray"), color = "MOM_Group_T4", shape = "Group") + geom_point(size = 3)

psnewt4_plot + scale_color_manual(values=c("#DC354B", "#A9BFCC")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_shape_manual(values=c(1, 16)) + theme(legend.position = "right") + stat_ellipse(type = "t", linetype =3, aes(group= MOM_Group_T4))

##PERMANOVA using bray curtis

#Multivariable analysis
metadata37 <- as(sample_data(ps_ff_total_zclr_pseudo), "data.frame")
adonis(phyloseq::distance(ps_ff_total_zclr_pseudo, method = "bray") ~ MOM_Group_T4 + Group + MOM_Group_T4*Group, data = metadata37, p.adjust.methods= "fdr", permutations = 99999)


##Using Jaccard distances

#unconstrained ordination
psnewt4_plot <- plot_ordination(ps_ff_total, ordinate(ps_ff_total, "MDS", "jaccard"), color = "MOM_Group_T4", shape = "Group") + geom_point(size = 3)

psnewt4_plot + scale_color_manual(values=c("#DC354B", "#A9BFCC")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_shape_manual(values=c(1, 16)) + theme(legend.position = "right") + stat_ellipse(type = "t", linetype =3, aes(group= MOM_Group_T4))


#Multivariable analysis
metadata37 <- as(sample_data(ps_ff_total), "data.frame")
adonis(phyloseq::distance(ps_ff_total, method = "jaccard") ~ MOM_Group_T4 + Group + MOM_Group_T4*Group, data = metadata37, p.adjust.methods= "fdr", permutations = 99999)





