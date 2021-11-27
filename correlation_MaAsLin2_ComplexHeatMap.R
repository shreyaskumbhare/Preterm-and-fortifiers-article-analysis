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

####Correlation analysis######################

########## T1###############


#ps1 <- readRDS(file = 'ps1')
ps_bf_total
#agglomerating at phylum level
ps1.glom <- tax_glom(ps_bf_total, "Phylum")
abund_table_z <- t(cmultRepl((otu_table(ps1.glom)), method="CZM", output="prop"))
#Apply CLR transformation
abund_clr <- t(apply(abund_table_z, 2, function(x){log(x) - mean (log(x))}))
#Replacing original counts with CLR transformed values
otu_table(ps1.glom) <- otu_table(abund_clr,taxa_are_rows = F)
#separating objects to be used for correlation analysis in the next step
otu.tab.ps1 <- otu_table(ps1.glom)
tax.tab.ps1 <- tax_table(ps1.glom)
metadata_t1 <- as.matrix(read.table(file = "metadata_T1_corr_final.csv", header = T,sep = ",",row.names = 1))
#microbiome package association fucntion
correlations.ps1 <- associate(otu.tab.ps1, metadata_t1, method = "spearman", mode = "table", p.adj.method = "BH")
correlations.ps1 <- merge(x=correlations.ps1,y=data.frame(ASV=rownames(t(otu.tab.ps1)),Phylum=tax.tab.ps1[,'Phylum']),by.x='X1',by.y='ASV')
correlations.ps1$TimePoint <- "T1"



########## T2 ###############

#ps2 <- readRDS(file = 'ps2')
ps_df_total
ps2.glom <- tax_glom(ps_df_total, "Phylum")
#library (zCompositions)
abund_table_z <- t(cmultRepl((otu_table(ps2.glom)), method="CZM", output="prop"))
#Apply CLR transformation
abund_clr <- t(apply(abund_table_z, 2, function(x){log(x) - mean (log(x))}))
#Replacing original counts with CLR transformed values
otu_table(ps2.glom) <- otu_table(abund_clr,taxa_are_rows = F)
#separating objects to be used for correlation analysis in the next step
otu.tab.ps2 <- otu_table(ps2.glom)
tax.tab.ps2 <- tax_table(ps2.glom)
metadata_t2 <- as.matrix(read.table(file = "metadata_T2_corr_final.csv", header = T,sep = ",",row.names = 1))
#microbiome package association fucntion
correlations.ps2 <- associate(otu.tab.ps2, metadata_t2, method = "spearman", mode = "table", p.adj.method = "BH")
correlations.ps2 <- merge(x=correlations.ps2,y=data.frame(ASV=rownames(t(otu.tab.ps2)),Phylum=tax.tab.ps2[,'Phylum']),by.x='X1',by.y='ASV')
correlations.ps2$TimePoint <- "T2"


########## T3 ###############


#ps3 <- readRDS(file = 'ps3')
ps3.glom <- tax_glom(ps_ef_total, "Phylum")
#library (zCompositions)
abund_table_z <- t(cmultRepl((otu_table(ps3.glom)), method="CZM", output="prop"))
#Apply CLR transformation
abund_clr <- t(apply(abund_table_z, 2, function(x){log(x) - mean (log(x))}))
#Replacing original counts with CLR transformed values
otu_table(ps3.glom) <- otu_table(abund_clr,taxa_are_rows = F)
#separating objects to be used for correlation analysis in the next step
otu.tab.ps3 <- otu_table(ps3.glom)
tax.tab.ps3 <- tax_table(ps3.glom)
metadata_t3 <- as.matrix(read.table(file = "metadata_T3_corr_final.csv", header = T,sep = ",",row.names = 1))
#microbiome package association fucntion
correlations.ps3 <- associate(otu.tab.ps3, metadata_t3, method = "spearman", mode = "table", p.adj.method = "BH")
correlations.ps3 <- merge(x=correlations.ps3,y=data.frame(ASV=rownames(t(otu.tab.ps3)),Phylum=tax.tab.ps3[,'Phylum']),by.x='X1',by.y='ASV')
correlations.ps3$TimePoint <- "T3"



########## T4 ###############

#ps4 <- readRDS(file = 'ps4')
ps4.glom <- tax_glom(ps_ff_total, "Phylum")
#library (zCompositions)
abund_table_z <- t(cmultRepl((otu_table(ps4.glom)), method="CZM", output="prop"))
#Apply CLR transformation
abund_clr <- t(apply(abund_table_z, 2, function(x){log(x) - mean (log(x))}))
#Replacing original counts with CLR transformed values
otu_table(ps4.glom) <- otu_table(abund_clr,taxa_are_rows = F)
#separating objects to be used for correlation analysis in the next step
otu.tab.ps4 <- otu_table(ps4.glom)
tax.tab.ps4 <- tax_table(ps4.glom)
metadata_t4 <- as.matrix(read.table(file = "metadata_T4_corr_final.csv", header = T,sep = ",",row.names = 1))
#microbiome package association fucntion
correlations.ps4 <- associate(otu.tab.ps4, metadata_t4, method = "spearman", mode = "table", p.adj.method = "BH")
correlations.ps4 <- merge(x=correlations.ps4,y=data.frame(ASV=rownames(t(otu.tab.ps4)),Phylum=tax.tab.ps4[,'Phylum']),by.x='X1',by.y='ASV')
correlations.ps4$TimePoint <- "T4"

##Adding all correlations data into a single object
correlations_phylum <- rbind(correlations.ps1,correlations.ps2,correlations.ps3,correlations.ps4)



######### There are multiple variables, but we were interested in a few that are known to be associated with microbome and are part of the study hypothesis. The following steps selects out correlations for variables of interest.


correlations_phylum <- correlations_phylum[correlations_phylum$X2 %in% c("MOM.intake","TEV","Nutritional.support", "Weight","WeightGain","GestationalAge","Antibiotics","Calprotectin","F2IsoP"),]



##Reading/saving the object in RDS file
#correlations_phylum <- readRDS("phylum_corr_data.rds")

View(correlations_phylum)

############### plot the phylum level correlations####################


cor_plot_phylum <- ggplot(data = correlations_phylum, mapping = aes(x = X2,y = Phylum,fill=Correlation)) + geom_tile(color='white', size=3) + scale_fill_gradient2(low = 'blue',mid = 'white',high = 'red') + theme_bw()

cor_plot_phylum <- cor_plot_phylum + coord_flip() + xlab("") + facet_wrap(~TimePoint,ncol = 4)

cor_plot_phylum


####Save data file

#save.image("phylum_correlation_final_data.RData")


#At Genus level
#For T1


#Agglomerate at genus level
psbftotal_genus <- tax_glom(ps_bf_total, "Genus")
psbftotal_genus
abund_table_z <- t(cmultRepl((otu_table(psbftotal_genus)), method="CZM", output="prop"))
#Apply CLR transformation
abund_clr <- t(apply(abund_table_z, 2, function(x){log(x) - mean (log(x))}))
#Replacing original counts with CLR transformed values
otu_table(psbftotal_genus) <- otu_table(abund_clr,taxa_are_rows = F)
#separating objects to be used for correlation analysis in the next step
otu.tab.ps1 <- otu_table(psbftotal_genus)
tax.tab.ps1 <- tax_table(psbftotal_genus)
metadata_t1 <- as.matrix(read.table(file = "metadata_T1_corr_final.csv", header = T,sep = ",",row.names = 1))

#microbiome package association fucntion
correlations.ps1 <- associate(otu.tab.ps1, metadata_t1, method = "spearman", mode = "table", p.adj.method = "BH")
correlations.ps1 <- merge(x=correlations.ps1,y=data.frame(ASV=rownames(t(otu.tab.ps1)),Genus=tax.tab.ps1[,'Genus']),by.x='X1',by.y='ASV')
View(correlations.ps1)


#For T2
#At Genus level

#Agglomerate at genus level
ps_df_total
psdftotal_genus <- tax_glom(ps_df_total, "Genus")
psdftotal_genus
abund_table_z <- t(cmultRepl((otu_table(psdftotal_genus)), method="CZM", output="prop"))
#Apply CLR transformation
abund_clr <- t(apply(abund_table_z, 2, function(x){log(x) - mean (log(x))}))
#Replacing original counts with CLR transformed values
otu_table(psdftotal_genus) <- otu_table(abund_clr,taxa_are_rows = F)
#separating objects to be used for correlation analysis in the next step
otu.tab.ps2 <- otu_table(psdftotal_genus)
tax.tab.ps2 <- tax_table(psdftotal_genus)
metadata_t2 <- as.matrix(read.table(file = "metadata_T2_corr_final.csv", header = T,sep = ",",row.names = 1))

#microbiome package association fucntion
correlations.ps2 <- associate(otu.tab.ps2, metadata_t2, method = "spearman", mode = "table", p.adj.method = "BH")
correlations.ps2 <- merge(x=correlations.ps2,y=data.frame(ASV=rownames(t(otu.tab.ps2)),Genus=tax.tab.ps2[,'Genus']),by.x='X1',by.y='ASV')
View(correlations.ps2)


#At T3

#For Genus level

#Agglomerate at genus level
ps_ef_total
pseftotal_genus <- tax_glom(ps_ef_total, "Genus")
pseftotal_genus
abund_table_z <- t(cmultRepl((otu_table(pseftotal_genus)), method="CZM", output="prop"))
#Apply CLR transformation
abund_clr <- t(apply(abund_table_z, 2, function(x){log(x) - mean (log(x))}))
#Replacing original counts with CLR transformed values
otu_table(pseftotal_genus) <- otu_table(abund_clr,taxa_are_rows = F)
#separating objects to be used for correlation analysis in the next step
otu.tab.ps3 <- otu_table(pseftotal_genus)
tax.tab.ps3 <- tax_table(pseftotal_genus)
metadata_t3 <- as.matrix(read.table(file = "metadata_T3_corr_final.csv", header = T,sep = ",",row.names = 1))

#microbiome package association fucntion
correlations.ps3 <- associate(otu.tab.ps3, metadata_t3, method = "spearman", mode = "table", p.adj.method = "fdr")
correlations.ps3 <- merge(x=correlations.ps3,y=data.frame(ASV=rownames(t(otu.tab.ps3)),Genus=tax.tab.ps3[,'Genus']),by.x='X1',by.y='ASV')
View(correlations.ps3)


#For T4

#At genus level


#Agglomerate at genus level
ps_ff_total
psfftotal_genus <- tax_glom(ps_ff_total, "Genus")
psfftotal_genus
abund_table_z <- t(cmultRepl((otu_table(psfftotal_genus)), method="CZM", output="prop"))
#Apply CLR transformation
abund_clr <- t(apply(abund_table_z, 2, function(x){log(x) - mean (log(x))}))
#Replacing original counts with CLR transformed values
otu_table(psfftotal_genus) <- otu_table(abund_clr,taxa_are_rows = F)
#separating objects to be used for correlation analysis in the next step
otu.tab.ps4 <- otu_table(psfftotal_genus)
tax.tab.ps4 <- tax_table(psfftotal_genus)
metadata_t4 <- as.matrix(read.table(file = "metadata_T4_corr_final.csv", header = T,sep = ",",row.names = 1))

#microbiome package association fucntion
correlations.ps4 <- associate(otu.tab.ps4, metadata_t4, method = "spearman", mode = "table", p.adj.method = "BH")
correlations.ps4 <- merge(x=correlations.ps4,y=data.frame(ASV=rownames(t(otu.tab.ps4)),Genus=tax.tab.ps4[,'Genus']),by.x='X1',by.y='ASV')
View(correlations.ps4)




######### filtering unwanted variables and then merging them using a for loop
correlations.ps1.filtered <- correlations.ps1[correlations.ps1$X2 %in% c("MOM.intake","TEV","Weight","WeightGain","GestationalAge","Antibiotics","Calprotectin","F2IsoP"),]
correlations.ps2.filtered <- correlations.ps2[correlations.ps2$X2 %in% c("MOM.intake","TEV","Weight","WeightGain","GestationalAge","Antibiotics","Calprotectin","F2IsoP"),]
correlations.ps3.filtered <- correlations.ps3[correlations.ps3$X2 %in% c("MOM.intake","TEV","Weight","WeightGain","GestationalAge","Antibiotics","Calprotectin","F2IsoP"),]
correlations.ps4.filtered <- correlations.ps4[correlations.ps4$X2 %in% c("MOM.intake","TEV","Weight","WeightGain","GestationalAge","Antibiotics","Calprotectin","F2IsoP"),]





#### Using the Complex heatmap to combine the genera level correlations together with hirarcheal clustering
#temp data frame
#library(ComplexHeatmap)

##we need to store all the correlations for all variables vs. all taxa for all time points together in one data frame. So, first create an empty data frame

data_temp <- data.frame()


###Now the following code is a for loop to store correlations one time point at a time, so we need to run the following chunk 4 times with each time point (use ps1, ps2, ps3, ps4 and t1, t2, t3, t4 respectively) and only then go to the next chunk. Note: please use appropriate combinations ps1 with t1....and so on

##this is a for loop for adding the correlation plots together in one single object, run the loop for all time points. Start with T1 and keep replacing the time points one after the other
for (i in 1:length(rownames(correlations.ps4.filtered))) {
  m <- correlations.ps4.filtered[i,'X2']
  m <- as.character(paste0('t4_',m))
  g <- as.character(correlations.ps4.filtered[i,'Genus'])
  c <- correlations.ps4.filtered[i,'Correlation']
  data_temp[g,m] <- c
}



#Saving this data frame as a separate object with correlation data added with the for loop above from all time points
saveRDS(object = data_temp,file = 'correlation_all.rds')


final_df <- as.matrix(readRDS(file = 'correlation_all.rds'))
final_df <- replace(final_df,is.na(final_df),0)

cols <- colnames(final_df)
cols <- gsub(pattern = 't1_.*',replacement = 'T1',x = cols,perl = T)
cols <- gsub(pattern = 't2_.*',replacement = 'T2',x = cols,perl = T)
cols <- gsub(pattern = 't3_.*',replacement = 'T3',x = cols,perl = T)
cols <- gsub(pattern = 't4_.*',replacement = 'T4',x = cols,perl = T)

colnames(final_df) <- gsub(pattern = 't\\d\\_',replacement = '',x = colnames(final_df),perl = T)
#pvalues
#data_temp_p <- data.frame()
#for (i in 3:length(rownames(correlation_t2))) {
#  m <- correlation_t2[i,'X2']
#  m <- as.character(paste0('t2_',m))
#  g <- as.character(correlation_t2[i,'Genus'])
#  c <- correlation_t2[i,'p.adj']
#  data_temp_p[g,m] <- c
#}

temp_row <- rownames(final_df)
temp_row <- gsub(pattern = "Unannotated ", replacement = '', x=temp_row)
temp_row <- gsub(pattern = ".IncertaeSedis", replacement = '', x=temp_row)
rownames(final_df)<- temp_row

library(circlize)
f1 = colorRamp2(seq(min(final_df), max(final_df), length = 1), c("blue", "#EEEEEE", "red"))
#f <- function(j, i, x, y, width, height, fill) {
#  grid.text(ifelse(data_temp_p[i, j] <= 0.05, sprintf("%.1f", #data_temp_p[i, j]), x, y, gp = gpar(fontsize = 10),""))
#}
genus_cor_plot <- Heatmap(matrix = final_df,column_split = cols,col=f1,row_names_side = 'left',row_dend_side = 'right', heatmap_legend_param =  list(title= "Spearman correlation"), clustering_distance_rows = function(x, y) 1 - cor(x, y), row_km = 2, cluster_column_slices = F)
genus_cor_plot





# save plot to file without using ggsave
#tiff(filename="genus_cor_plot_with_col_clust.tiff", width=25, height=32, units="in", compression="lzw", bg="white", res=300)
#print(genus_cor_plot)
#dev.off()



###Extracting genus level information for preparing a table


genus_rct_data <- read.csv("otu_table_genus_rm_clr_new.csv", sep = ",", header = T)
genus_rct_data

#Subset by time point
genus_data_t1 <- genus_rct_data[ which(genus_rct_data$TimePoint=='T1'), ]
genus_data_t2 <- genus_rct_data[ which(genus_rct_data$TimePoint=='T2'), ]
genus_data_t3 <- genus_rct_data[ which(genus_rct_data$TimePoint=='T3'), ]
genus_data_t4 <- genus_rct_data[ which(genus_rct_data$TimePoint=='T4'), ]

#transpose data frame and add a column variable header 'Genus'
#Melting the data frame to use it summary statistics
trans_genus_data_t1 <- melt(genus_data_t1)
trans_genus_data_t1
trans_genus_data_t2 <- melt(genus_data_t2)
trans_genus_data_t2
trans_genus_data_t3 <- melt(genus_data_t3)
trans_genus_data_t3
trans_genus_data_t4 <- melt(genus_data_t4)
trans_genus_data_t4
#Changing the column name from 'variable' to 'Genus'
colnames(trans_genus_data_t1)[4] <- "Genus"
trans_genus_data_t1
colnames(trans_genus_data_t2)[4] <- "Genus"
trans_genus_data_t2
colnames(trans_genus_data_t3)[4] <- "Genus"
trans_genus_data_t3
colnames(trans_genus_data_t4)[4] <- "Genus"
trans_genus_data_t4

##Obtain summary for genus by group
summary_genus_t3 <- Summarize(trans_genus_data_t1$value ~ trans_genus_data_t1$Genus*trans_genus_data_t1$Group)
summary_genus_t2 <- Summarize(trans_genus_data_t2$value ~ trans_genus_data_t2$Genus*trans_genus_data_t2$Group)
summary_genus_t3 <- Summarize(trans_genus_data_t3$value ~ trans_genus_data_t3$Genus*trans_genus_data_t3$Group)
summary_genus_t4 <- Summarize(trans_genus_data_t4$value ~ trans_genus_data_t4$Genus*trans_genus_data_t4$Group)

##Exporting these files as CSV
write.csv(summary_genus_t1, "summary_genus_t1.csv")
write.csv(summary_genus_t2, "summary_genus_t2.csv")
write.csv(summary_genus_t3, "summary_genus_t3.csv")
write.csv(summary_genus_t4, "summary_genus_t4.csv")


###Extracting genus level data now by MOM groups
genus_mom_data <- read.csv("otu_table_genus_rm_MOM_clr_new.csv", sep = ",", header = T)
genus_mom_data

#Subset by time point
genus_data_t1 <- genus_mom_data[ which(genus_mom_data$TimePoint=='T1'), ]
genus_data_t2 <- genus_mom_data[ which(genus_mom_data$TimePoint=='T2'), ]
genus_data_t3 <- genus_mom_data[ which(genus_mom_data$TimePoint=='T3'), ]
genus_data_t4 <- genus_mom_data[ which(genus_mom_data$TimePoint=='T4'), ]

#transpose data frame and add a column variable header 'Genus'
#Melting the data frame to use it summary statistics
trans_genus_data_t1 <- melt(genus_data_t1)
trans_genus_data_t1
trans_genus_data_t2 <- melt(genus_data_t2)
trans_genus_data_t2
trans_genus_data_t3 <- melt(genus_data_t3)
trans_genus_data_t3
trans_genus_data_t4 <- melt(genus_data_t4)
trans_genus_data_t4
#Changing the column name from 'variable' to 'Genus'
colnames(trans_genus_data_t1)[4] <- "Genus"
trans_genus_data_t1
colnames(trans_genus_data_t2)[4] <- "Genus"
trans_genus_data_t2
colnames(trans_genus_data_t3)[4] <- "Genus"
trans_genus_data_t3
colnames(trans_genus_data_t4)[4] <- "Genus"
trans_genus_data_t4

##Obtain summary for genus by group
summary_genus_t1 <- Summarize(trans_genus_data_t1$value ~ trans_genus_data_t1$Genus*trans_genus_data_t1$MOM_Group)
summary_genus_t2 <- Summarize(trans_genus_data_t2$value ~ trans_genus_data_t2$Genus*trans_genus_data_t2$MOM_Group)
summary_genus_t3 <- Summarize(trans_genus_data_t3$value ~ trans_genus_data_t3$Genus*trans_genus_data_t3$MOM_Group)
summary_genus_t4 <- Summarize(trans_genus_data_t4$value ~ trans_genus_data_t4$Genus*trans_genus_data_t4$MOM_Group)

##Exporting these files as CSV
write.csv(summary_genus_t1, "summary_genus_mom_t1.csv")
write.csv(summary_genus_t2, "summary_genus_mom_t2.csv")
write.csv(summary_genus_t3, "summary_genus_mom_t3.csv")
write.csv(summary_genus_t4, "summary_genus_mom_t4.csv")


#####The chunk above was used to extract the CLR transformed abundance values and to apply summary statistics: Median, IQR etc.
#####The Following code is used to extract the relative abundance values and to apply summary statistics. Please note the differential abundance testing and longitudinal analysis (trend analysis) was performed on CLR transformed abundance.

###Extracting relative abundance data by RCT groups

psnoncontam_tru_final_rare_rm

##Agglomerate at phyla level and extract the OTU table
genus_data_total_rm <- tax_glom(psnoncontam_tru_final_rare_rm,taxrank = "Genus")
genus_data_total_rm

#Relative abundances
genus_data_total_rm.rel <- transform(genus_data_total_rm, "compositional")
genus_data_total_rm.rel

#replacing ASVs with short string
taxa_names(genus_data_total_rm.rel) <- paste0("SV", seq(ntaxa(genus_data_total_rm.rel)))
#View(tax_table(genus_data_rm_zclr))

#Extracting short strings and their respective taxonomy
genus_list<- tax_table(genus_data_total_rm.rel)@.Data
write.csv(genus_list, "genus_total_list_rel_rct.csv")

#extracting as csv and replacing string names with phyla names, transposed and added metadata
write_phyloseq(genus_data_total_rm.rel, type= 'all')


#transposed, added additional metadata and importing the file
genus_rel_rm_rct <- read.table("otu_table_rm_rel_rct_new.csv", header = TRUE, sep= ",", check.names = F)
genus_rel_rm_rct

#Subset by time point
genus_data_t1 <- genus_rel_rm_rct[ which(genus_rel_rm_rct$TimePoint=='T1'), ]
genus_data_t2 <- genus_rel_rm_rct[ which(genus_rel_rm_rct$TimePoint=='T2'), ]
genus_data_t3 <- genus_rel_rm_rct[ which(genus_rel_rm_rct$TimePoint=='T3'), ]
genus_data_t4 <- genus_rel_rm_rct[ which(genus_rel_rm_rct$TimePoint=='T4'), ]

#transpose data frame and add a column variable header 'Genus'
#Melting the data frame to use it summary statistics
trans_genus_data_t1 <- melt(genus_data_t1)
trans_genus_data_t1
trans_genus_data_t2 <- melt(genus_data_t2)
trans_genus_data_t2
trans_genus_data_t3 <- melt(genus_data_t3)
trans_genus_data_t3
trans_genus_data_t4 <- melt(genus_data_t4)
trans_genus_data_t4
#Changing the column name from 'variable' to 'Genus'
colnames(trans_genus_data_t1)[4] <- "Genus"
trans_genus_data_t1
colnames(trans_genus_data_t2)[4] <- "Genus"
trans_genus_data_t2
colnames(trans_genus_data_t3)[4] <- "Genus"
trans_genus_data_t3
colnames(trans_genus_data_t4)[4] <- "Genus"
trans_genus_data_t4

##Obtain summary for genus by group
summary_genus_t1 <- Summarize(trans_genus_data_t1$value ~ trans_genus_data_t1$Genus*trans_genus_data_t1$Group)
summary_genus_t2 <- Summarize(trans_genus_data_t2$value ~ trans_genus_data_t2$Genus*trans_genus_data_t2$Group)
summary_genus_t3 <- Summarize(trans_genus_data_t3$value ~ trans_genus_data_t3$Genus*trans_genus_data_t3$Group)
summary_genus_t4 <- Summarize(trans_genus_data_t4$value ~ trans_genus_data_t4$Genus*trans_genus_data_t4$Group)

##Exporting these files as CSV
write.csv(summary_genus_t1, "summary_genus_rel_rct_t1.csv")
write.csv(summary_genus_t2, "summary_genus_rel_rct_t2.csv")
write.csv(summary_genus_t3, "summary_genus_rel_rct_t3.csv")
write.csv(summary_genus_t4, "summary_genus_rel_rct_t4.csv")


###Extracting genus level data now by MOM groups


genus_mom_rel_data <- read.csv("otu_table_rm_rel_mom_new.csv", sep = ",", header = T)
genus_mom_rel_data

#Subset by time point
genus_data_t1 <- genus_mom_rel_data[ which(genus_mom_rel_data$TimePoint=='T1'), ]
genus_data_t2 <- genus_mom_rel_data[ which(genus_mom_rel_data$TimePoint=='T2'), ]
genus_data_t3 <- genus_mom_rel_data[ which(genus_mom_rel_data$TimePoint=='T3'), ]
genus_data_t4 <- genus_mom_rel_data[ which(genus_mom_rel_data$TimePoint=='T4'), ]

#transpose data frame and add a column variable header 'Genus'
#Melting the data frame to use it summary statistics
trans_genus_data_t1 <- melt(genus_data_t1)
trans_genus_data_t1
trans_genus_data_t2 <- melt(genus_data_t2)
trans_genus_data_t2
trans_genus_data_t3 <- melt(genus_data_t3)
trans_genus_data_t3
trans_genus_data_t4 <- melt(genus_data_t4)
trans_genus_data_t4
#Changing the column name from 'variable' to 'Genus'
colnames(trans_genus_data_t1)[4] <- "Genus"
trans_genus_data_t1
colnames(trans_genus_data_t2)[4] <- "Genus"
trans_genus_data_t2
colnames(trans_genus_data_t3)[4] <- "Genus"
trans_genus_data_t3
colnames(trans_genus_data_t4)[4] <- "Genus"
trans_genus_data_t4

##Obtain summary for genus by group
summary_genus_t1 <- Summarize(trans_genus_data_t1$value ~ trans_genus_data_t1$Genus*trans_genus_data_t1$MOM_Group)
summary_genus_t2 <- Summarize(trans_genus_data_t2$value ~ trans_genus_data_t2$Genus*trans_genus_data_t2$MOM_Group)
summary_genus_t3 <- Summarize(trans_genus_data_t3$value ~ trans_genus_data_t3$Genus*trans_genus_data_t3$MOM_Group)
summary_genus_t4 <- Summarize(trans_genus_data_t4$value ~ trans_genus_data_t4$Genus*trans_genus_data_t4$MOM_Group)

##Exporting these files as CSV
write.csv(summary_genus_t1, "summary_genus_mom_rel_t1.csv")
write.csv(summary_genus_t2, "summary_genus_mom_rel_t2.csv")
write.csv(summary_genus_t3, "summary_genus_mom_rel_t3.csv")
write.csv(summary_genus_t4, "summary_genus_mom_rel_t4.csv")


#####Multivaribale regression using MaAsLin2 package (https://huttenhower.sph.harvard.edu/maaslin/). Performed regression cross-sectional and longitudinal (mixed effects model) using Maasalin2 package

###Installing Maasalin2

#library(devtools)
#devtools::install_github("biobakery/maaslin2")


#library(Maaslin2)

##Extracting phyloseq object by time points for regression analysis
pstotal_rm <- readRDS("phyloseq_shreyas.rds")


##Apply relative abundance
#ps.glom.prop <- transform_sample_counts(ps.glom, function(otu) otu/sum(otu))

###Agglomerate at genus level for all time points
pstotal_glom <- tax_glom(physeq = pstotal_rm,taxrank = 'Genus')
pstotal_glom

#extracting otu table and metadata from the phyloseq object
df <- as.data.frame.matrix(otu_table(pstotal_glom))
tax <- as.matrix(tax_table(pstotal_glom))
colnames(df) <- as.character(tax[,'Genus'])


#Apply CLR transformation and replace the original counts with clr abundance
abund_table <- t(cmultRepl((df), method="CZM", output="prop",))
abund_clr <- t(apply(abund_table, 2, function(x){log(x) - mean (log(x))}))
samdf <- as.matrix(sample_data(pstotal_glom))
write.table(x = abund_clr,file = 'abundance_genus_total_clr.tsv',quote = F,sep = '\t',col.names = NA)
write.table(x = samdf,file = 'sample_data_genus_total_clr.tsv',sep = '\t',quote = F,col.names = NA)

abund_clr <- read.table("abundance_genus_total_clr.tsv", sep = '\t', header = T, row.names = 1)
samdf <- read.table("sample_data_genus_total_clr.tsv", sep = '\t', header = T, row.names = 1)

reg_t1 <- abund_clr[as.character(samdf$Description)== "T1",]
reg_t2 <- abund_clr[as.character(samdf$Description)== "T2",]
reg_t3 <- abund_clr[as.character(samdf$Description)== "T3",]
reg_t4 <- abund_clr[as.character(samdf$Description)== "T4",]

samdf_t1 <- samdf[as.character(samdf$Description)== "T1",]
samdf_t2 <- samdf[as.character(samdf$Description)== "T2",]
samdf_t3 <- samdf[as.character(samdf$Description)== "T3",]
samdf_t4 <- samdf[as.character(samdf$Description)== "T4",]

write.table(x = reg_t1,file = 'abundance_genus_t1_final.tsv',quote = F,sep = '\t',col.names = NA)
write.table(x = samdf_t1,file = 'sample_data_t1_final.tsv',sep = '\t',quote = F,col.names = NA)

write.table(x = reg_t2,file = 'abundance_genus_t2_final.tsv',quote = F,sep = '\t',col.names = NA)
write.table(x = samdf_t2,file = 'sample_data_t2_final.tsv',sep = '\t',quote = F,col.names = NA)

write.table(x = reg_t3,file = 'abundance_genus_t3_final.tsv',quote = F,sep = '\t',col.names = NA)
write.table(x = samdf_t3,file = 'sample_data_t3_final.tsv',sep = '\t',quote = F,col.names = NA)

write.table(x = reg_t4,file = 'abundance_genus_t4_final.tsv',quote = F,sep = '\t',col.names = NA)
write.table(x = samdf_t4,file = 'sample_data_t4_final.tsv',sep = '\t',quote = F,col.names = NA)


library(Maaslin2)
#Regression at each time point

#For all time points
genus_data_maasalin_t1 = Maaslin2(
  input_data = 'abundance_genus_t1_final.tsv', 
  input_metadata = 'sample_data_t1_final.tsv', 
  output = "genus_reg_final_t1", 
  fixed_effects = c("Group", "MOM_Group"), normalization = 'NONE',transform = 'NONE',correction = 'BH',analysis_method = 'LM', min_abundance = -50, min_prevalence = 0.0, min_variance = 0, max_significance = 0.05)

genus_data_maasalin_t2 = Maaslin2(
  input_data = 'abundance_genus_t2_final.tsv', 
  input_metadata = 'sample_data_t2_final.tsv', 
  output = "genus_reg_final_t2", 
  fixed_effects = c("Group", "MOM_Group"), normalization = 'NONE',transform = 'NONE',correction = 'BH',analysis_method = 'LM', min_abundance = -50, min_prevalence = 0.0, min_variance = 0, max_significance = 0.05)


genus_data_maasalin_t3 = Maaslin2(
  input_data = 'abundance_genus_t3_final.tsv', 
  input_metadata = 'sample_data_t3_final.tsv', 
  output = "genus_reg_final_t3", 
  fixed_effects = c("Group", "MOM_Group"), normalization = 'NONE',transform = 'NONE',correction = 'BH',analysis_method = 'LM', min_abundance = -50, min_prevalence = 0.0, min_variance = 0, max_significance = 0.05)

genus_data_maasalin_t4 = Maaslin2(
  input_data = 'abundance_genus_t4_final.tsv', 
  input_metadata = 'sample_data_t4_final.tsv', 
  output = "genus_reg_final_t4", 
  fixed_effects = c("Group", "MOM_Group"), normalization = 'NONE',transform = 'NONE',correction = 'BH',analysis_method = 'LM', min_abundance = -50, min_prevalence = 0.0, min_variance = 0, max_significance = 0.05)


###At phylum level


#library(Maaslin2)
###Agglomerate at Phylum level
ps.glom <- tax_glom(physeq = psnoncontam_tru_final_rare_rm,taxrank = 'Phylum')
ps.glom
#extracting otu table and metadata from the phyloseq object
df <- as.data.frame.matrix(otu_table(ps.glom))
tax <- as.matrix(tax_table(ps.glom))
colnames(df) <- as.character(tax[,'Phylum'])

#Apply CLR transformation
abund_table_z <- t(cmultRepl((df), method="CZM", output="prop"))
abund_clr <- t(apply(abund_table_z, 2, function(x){log(x) - mean (log(x))}))


samdf <- as.matrix(sample_data(ps.glom))


write.table(x = abund_clr,file = 'abundance_phylum.tsv',quote = F,sep = '\t',col.names = NA)
write.table(x = samdf,file = 'sample_data_phylum.tsv',sep = '\t',quote = F,col.names = NA)

#library(Maaslin2)
phylum_data_maasalin = Maaslin2(
  input_data = 'abundance_phylum.tsv', 
  input_metadata = 'sample_data_phylum.tsv', 
  output = "phylum_level_output", 
  fixed_effects = c("Group", "MOM_Group"),
  random_effects = c("SubjectID"),normalization = 'NONE',transform = 'NONE',correction = 'BH',analysis_method = 'LM')


##At Family level


#library(Maaslin2)
###Agglomerate at Family level
ps.glom <- tax_glom(physeq = psnoncontam_tru_final_rare_rm,taxrank = 'Family')
ps.glom
#extracting otu table and metadata from the phyloseq object
df <- as.data.frame.matrix(otu_table(ps.glom))
tax <- as.matrix(tax_table(ps.glom))
colnames(df) <- as.character(tax[,'Family'])

#Apply CLR transformation
abund_table_z <- t(cmultRepl((df), method="CZM", output="prop"))
abund_clr <- t(apply(abund_table_z, 2, function(x){log(x) - mean (log(x))}))


samdf <- as.matrix(sample_data(ps.glom))


write.table(x = abund_clr,file = 'abundance_family.tsv',quote = F,sep = '\t',col.names = NA)
write.table(x = samdf,file = 'sample_data_family.tsv',sep = '\t',quote = F,col.names = NA)

#library(Maaslin2)
family_data_maasalin = Maaslin2(
  input_data = 'abundance_family.tsv', 
  input_metadata = 'sample_data_family.tsv', 
  output = "family_level_output", 
  fixed_effects = c("Group", "MOM_Group"),
  random_effects = c("SubjectID"),normalization = 'NONE',transform = 'NONE',correction = 'BH',analysis_method = 'LM')


######################Complex HeatMap##########################
rm(list = ls())
library(ComplexHeatmap)
library(circlize)

####Importing abundance table and results from regression output for each time point
#Defining the color scale for heatmap
col_fun = colorRamp2(c(min(dat.t), 0, max(dat.t)), c("#275da1", "white", "#FF0000"))


####Preparing data frames for using to plot heatmaps at each time point

dat = read.table(file = "abundance_genus_total_clr.tsv", header = T, sep = "\t", check.names = F, row.names = 1)
met = read.table(file = "sample_data_genus_total_clr.tsv", header = T, sep = "\t", check.names = F, row.names = 1)
dat.t = t(dat)
met$MOM_Group = factor(x = met$MOM_Group, levels = c("LMOM", "HMOM"))

View(dat.t)
get_reg_obj = function(filename) {
  reg = read.table(file = filename, header = T, sep = "\t", check.names = F)
  reg.momgroup = reg[as.character(reg$metadata) == "MOM_Group",]
  reg.group = reg[as.character(reg$metadata) == "Group",]
  rownames(reg.momgroup) = as.character(reg.momgroup$feature)
  rownames(reg.group) = as.character(reg.group$feature)
  reg.momgroup = reg.momgroup[,c("coef", "stderr", "qval")]
  reg.group = reg.group[,c("coef", "stderr", "qval")]
  return(list("momgroup" = reg.momgroup, "group" = reg.group))
}

for (i in list.files(path = ".", pattern = "all_results*")) {
  sampleName = gsub(pattern = "all_results_t(\\d+)_final.tsv",  replacement = "T\\1", x = i)
  temp.obj = get_reg_obj(filename = i)
  momgroup.reg = temp.obj$momgroup
  group.reg = temp.obj$group
  group.reg = group.reg[rownames(dat.t),]
  momgroup.reg = momgroup.reg[rownames(dat.t),]
  
  if(sampleName == "T1"){
    group.ht = Heatmap(dat.t[,as.character(met$Description) == sampleName], name = "group", col= col_fun, show_column_names = F, cluster_columns = F, show_parent_dend_line = F, row_names_side = "left", column_order = order(met$Group[as.character(met$Description) == sampleName]), cluster_rows = F, show_row_dend = F, show_row_names = T, border = T, top_annotation = HeatmapAnnotation("Group"= met$Group[as.character(met$Description) == sampleName], show_annotation_name = FALSE, col = list("Group" = c("Control_group" = "#00AFBB", "Treatment_group" = "#E7B800"))), right_annotation = rowAnnotation(width = unit(1.5, "cm"), gap = unit(0.2, "cm"), "coef" = anno_barplot(group.reg$coef), "qval" = anno_text(format(round(group.reg$qval, 2), nsmall = 2), location = 0.5, just = "center", gp = gpar(col = "black", fontsize = 10))), heatmap_legend_param = list(title = "Abundance (CLR)", legend_height = unit(4, "cm"), title_position = "topcenter", direction = "horizontal"))
    
    momgroup.ht = Heatmap(dat.t[,as.character(met$Description) == sampleName], name = "group", col= col_fun, show_column_names = F, cluster_columns = F, show_parent_dend_line = F, row_names_side = "left", column_order = order(met$MOM_Group[as.character(met$Description) == sampleName]), cluster_rows = F, show_row_dend = F, show_row_names = T, border = T, top_annotation = HeatmapAnnotation("MOM_Group"= met$MOM_Group[as.character(met$Description) == sampleName], show_annotation_name = FALSE, col = list("MOM_Group" = c("LMOM" = "#A9BFCC", "HMOM" = "#DC354B"))), right_annotation = rowAnnotation(width = unit(1.5, "cm"), gap = unit(0.2, "cm"), "coef" = anno_barplot(momgroup.reg$coef), "qval" = anno_text(format(round(momgroup.reg$qval, 2), nsmall = 2), location = 0.5, just = "center", gp = gpar(col = "black", fontsize = 10))), heatmap_legend_param = list(title = "Abundance (CLR)", legend_height = unit(4, "cm"), title_position = "topcenter", direction = "horizontal"))
    
  } else {
    group.ht = group.ht + Heatmap(dat.t[,as.character(met$Description) == sampleName], name = "group", col= col_fun, show_column_names = F, cluster_columns = F, show_parent_dend_line = F, row_names_side = "left", column_order = order(met$Group[as.character(met$Description) == sampleName]), cluster_rows = F, show_row_dend = F, show_row_names = F, border = T, top_annotation = HeatmapAnnotation("Group"= met$Group[as.character(met$Description) == sampleName], show_annotation_name = FALSE, col = list("Group" = c("Control_group" = "#00AFBB", "Treatment_group" = "#E7B800"))), right_annotation = rowAnnotation(width = unit(1.5, "cm"), gap = unit(0.2, "cm"), "coef" = anno_barplot(group.reg$coef), "qval" = anno_text(format(round(group.reg$qval, 2), nsmall = 2), location = 0.5, just = "center", gp = gpar(col = "black", fontsize = 10))), heatmap_legend_param = list(title = "Abundance (CLR)", legend_height = unit(4, "cm"), title_position = "topcenter", direction = "horizontal"))
    
    momgroup.ht = momgroup.ht + Heatmap(dat.t[,as.character(met$Description) == sampleName], name = "group", col= col_fun, show_column_names = F, cluster_columns = F, show_parent_dend_line = F, row_names_side = "left", column_order = order(met$MOM_Group[as.character(met$Description) == sampleName]), cluster_rows = F, show_row_dend = F, show_row_names = F, border = T, top_annotation = HeatmapAnnotation("MOM_Group"= met$MOM_Group[as.character(met$Description) == sampleName], show_annotation_name = FALSE, col = list("MOM_Group" = c("LMOM" = "#A9BFCC", "HMOM" = "#DC354B"))), right_annotation = rowAnnotation(width = unit(1.5, "cm"), gap = unit(0.2, "cm"), "coef" = anno_barplot(momgroup.reg$coef), "qval" = anno_text(format(round(momgroup.reg$qval, 2), nsmall = 2), location = 0.5, just = "center", gp = gpar(col = "black", fontsize = 10))), heatmap_legend_param = list(title = "Abundance (CLR)", legend_height = unit(4, "cm"), title_position = "topcenter", direction = "horizontal"))
  }
}

draw(object = group.ht, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
draw(object = momgroup.ht, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")

#######End#############