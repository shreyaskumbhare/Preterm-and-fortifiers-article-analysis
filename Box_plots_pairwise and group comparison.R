##Plotting box plots for phyla at four different time points by RCT

##Using the rarified RM phyloseq file for plotting phlya trend


psnoncontam_tru_final_rare_rm

##Agglomerate at phyla level and extract the OTU table
phylum_data_total_rm <- tax_glom(psnoncontam_tru_final_rare_rm,taxrank = "Phylum")
phylum_data_total_rm

#Relative abundances
#phylum_data_total_rm.rel <- transform(phylum_data_total_rm, "compositional")
#phylum_data_total_rm.rel

#using CLR transformed abundances to plot

#Using CLR transformed abundances here
#Write the phyloseq object as a CSV table with ASVs and counts
phylum_data_total_rm
write_phyloseq(phylum_data_total_rm)

#Read the CSV file
abund_table <- read.csv("phylum_rm_rare_clr_new.csv", row.names=1, check.names= FALSE) #Added header "Features" to the first clumn of ASVs
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
phylum_data_rm_zclr <- phylum_data_total_rm
otu_table(phylum_data_rm_zclr) <- otu_table(abund_clr, taxa_are_rows = F)
phylum_data_rm_zclr

#replacing ASVs with short string
taxa_names(phylum_data_rm_zclr) <- paste0("SV", seq(ntaxa(phylum_data_rm_zclr)))
#View(tax_table(phylum_data_rm_zclr))

#Extracting short strings and their respective taxonomy
phylum_list<- tax_table(phylum_data_rm_zclr)@.Data
write.csv(phylum_list, "phylum_total_list.csv")

#extracting as csv and replacing string names with phyla names, transposed and added metadata
write_phyloseq(phylum_data_rm_zclr, type= 'all')


#transposed, added additional metadata and importing the file
phylum_zibr_rm <- read.table("otu_table_phyla_rm_clr_new.csv", header = TRUE, sep= ",", check.names = F)
phylum_zibr_rm

#Melting the data frame to use it for boxplot and statistics
trans_phylum_zibr_rm <- melt(phylum_zibr_rm)
trans_phylum_zibr_rm

#Changing the column name from 'variable' to 'Phylum'
colnames(trans_phylum_zibr_rm)[4] <- "Phylum"
trans_phylum_zibr_rm

#Rstatix
stat.test.phylum <- trans_phylum_zibr_rm %>%
  group_by(Phylum, TimePoint) %>%
  wilcox_test(value ~ Group) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()
stat.test.phylum

#Adding XY position
stat.test.phylum <- stat.test.phylum %>% add_xy_position(x = "TimePoint")
stat.test.phylum

#Box plot
bxp.phylum.zibr <- ggboxplot(data= trans_phylum_zibr_rm, x = "TimePoint", y= "value", facet.by = "Phylum", xlab = "Time Points", ylab="Relative abundance", fill = "Group", width = 0.8, outlier.shape= NA)

#modifying the plot to use individual y-scale for each facet
bxp.phylum.zibr <- facet(p = bxp.phylum.zibr,facet.by = "Phylum",scales = "free_y") + scale_fill_manual(values = c("#80D7DD", "#F3DB80")) +
  geom_point(aes(color = Group), position = position_dodge(0.7)) + stat_summary(
    fun = median,
    geom = 'line',
    aes(group = Group, colour = Group),
    position = position_dodge(width = 0.7) #this has to be added
  ) +
  scale_color_manual(values = c("#354360", "#7D6608")) + xlab("Time points") +
  ylab("Abundance (CLR)")

bxp.phylum.zibr + theme_bw() + theme(panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank()) + stat_pvalue_manual(stat.test.phylum,label= "p.adj.signif", tip.length=0.02, hide.ns = T)



####Box plots for genera at four time points by RCT groups


psnoncontam_tru_final_rare_rm

##Agglomerate at phyla level and extract the OTU table
genus_data_total_rm <- tax_glom(psnoncontam_tru_final_rare_rm,taxrank = "Genus")
genus_data_total_rm

#Relative abundances
#genus_data_total_rm.rel <- transform(genus_data_total_rm, "compositional")
#genus_data_total_rm.rel

#Using CLR abundances for plot
#Using CLR transformed abundances here
#Write the phyloseq object as a CSV table with ASVs and counts
genus_data_total_rm
write_phyloseq(genus_data_total_rm)

#Read the CSV file
abund_table <- read.csv("genus_rm_rare_clr_new.csv", row.names=1, check.names= FALSE) #Added header "Features" to the first clumn of ASVs
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
genus_data_rm_zclr <- genus_data_total_rm
otu_table(genus_data_rm_zclr) <- otu_table(abund_clr, taxa_are_rows = F)
genus_data_rm_zclr

#replacing ASVs with short string
taxa_names(genus_data_rm_zclr) <- paste0("SV", seq(ntaxa(genus_data_rm_zclr)))
#View(tax_table(genus_data_rm_zclr))

#Extracting short strings and their respective taxonomy
genus_list<- tax_table(genus_data_rm_zclr)@.Data
write.csv(genus_list, "genus_total_list.csv")

#extracting as csv and replacing string names with phyla names, transposed and added metadata
write_phyloseq(genus_data_rm_zclr, type= 'all')


#transposed, added additional metadata and importing the file
genus_zibr_rm <- read.table("otu_table_genus_rm_clr_new.csv", header = TRUE, sep= ",", check.names = F)
genus_zibr_rm

#Melting the data frame to use it for boxplot and statistics
trans_genus_zibr_rm <- melt(genus_zibr_rm)
trans_genus_zibr_rm

#Changing the column name from 'variable' to 'Genus'
colnames(trans_genus_zibr_rm)[4] <- "Genus"
trans_genus_zibr_rm

#Rstatix
stat.test.genus <- trans_genus_zibr_rm %>%
  group_by(Genus, TimePoint) %>%
  wilcox_test(value ~ Group) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()
stat.test.genus

#Exporting the post-hoc analysis as csv file
#write.csv(stat.test.genus, "genus_summary_posthoc_statistics_rct.csv")


#Adding XY position
stat.test.genus <- stat.test.genus %>% add_xy_position(x = "TimePoint")
stat.test.genus

#Box plot
bxp.genus.zibr <- ggboxplot(data= trans_genus_zibr_rm, x = "TimePoint", y= "value", facet.by = "Genus", xlab = "Time Points", ylab="Abundance (CLR)", fill = "Group", width = 0.8, outlier.shape= NA)

#modifying the plot to use individual y-scale for each facet
bxp.genus.zibr <- facet(p = bxp.genus.zibr,facet.by = "Genus",scales = "free_y") + scale_fill_manual(values = c("#80D7DD", "#F3DB80")) +
  geom_point(aes(color = Group), position = position_dodge(0.7)) + stat_summary(
    fun = median,
    geom = 'line',
    aes(group = Group, colour = Group),
    position = position_dodge(width = 0.7) #this has to be added
  ) +
  scale_color_manual(values = c("#354360", "#7D6608")) + xlab("Time points") +
  ylab("Abundance (CLR)")

bxp.genus.zibr + theme_bw() + theme(panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank()) + stat_pvalue_manual(stat.test.genus,label= "p.adj.signif", tip.length=0.02, hide.ns = T)


####Box plots for taxa at four time points by MOM groups

##At phylum level


#psnoncontam_tru_final_rare_rm

##Agglomerate at phyla level and extract the OTU table
#phylum_data_total_rm <- tax_glom(psnoncontam_tru_final_rare_rm,taxrank = "Phylum")
#phylum_data_total_rm

#Relative abundances
#phylum_data_total_rm.rel <- transform(phylum_data_total_rm, "compositional")
#phylum_data_total_rm.rel

#using CLR transformed abundances to plot

#Using CLR transformed abundances here
#Write the phyloseq object as a CSV table with ASVs and counts
#phylum_data_total_rm
#write_phyloseq(phylum_data_total_rm)

#Read the CSV file
#abund_table <- read.csv("phylum_rm_rare_clr.csv", row.names=1, check.names= FALSE) #Added header "Features" to the first clumn of ASVs
#abund_table_t<-t(abund_table)
#View(abund_table_t)

#install.packages("zCompositions")
#library (zCompositions)
#abund_table_r <- t(cmultRepl((abund_table_t), method="CZM", output="prop")) #prop calculates the imputed proportion is similar to "counts" and then taking proportions in the previous version; CZM: count zero multiplicative

#Check the before and after files
#View(abund_table)
#View(abund_table_r)

#Apply CLR transformation
#abund_clr <- t(apply(abund_table_r, 2, function(x){log(x) - mean (log(x))}))
#View(abund_clr)

#Replacing original counts with CLR transformed values
#phylum_data_rm_zclr <- phylum_data_total_rm
#otu_table(phylum_data_rm_zclr) <- otu_table(abund_clr, taxa_are_rows = F)
#phylum_data_rm_zclr

#replacing ASVs with short string
#taxa_names(phylum_data_rm_zclr) <- paste0("SV", seq(ntaxa(phylum_data_rm_zclr)))
#View(tax_table(phylum_data_rm_zclr))

#Extracting short strings and their respective taxonomy
#phylum_list<- tax_table(phylum_data_rm_zclr)@.Data
#write.csv(phylum_list, "phylum_total_list.csv")

#extracting as csv and replacing string names with phyla names, transposed and added metadata
write_phyloseq(phylum_data_rm_zclr, type= 'all')


#transposed, added additional metadata and importing the file
phylum_zibr_rm_MOM <- read.table("otu_table_phyla_rm_MOM_clr_new.csv", header = TRUE, sep= ",", check.names = F)
phylum_zibr_rm_MOM

#Melting the data frame to use it for boxplot and statistics
trans_phylum_zibr_rm_MOM <- melt(phylum_zibr_rm_MOM)
trans_phylum_zibr_rm_MOM

#Changing the column name from 'variable' to 'Phylum'
colnames(trans_phylum_zibr_rm_MOM)[4] <- "Phylum"
trans_phylum_zibr_rm_MOM

#Rstatix
stat.test.phylum.MOM <- trans_phylum_zibr_rm_MOM %>%
  group_by(Phylum, TimePoint) %>%
  wilcox_test(value ~ MOM_Group) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()
stat.test.phylum.MOM

#Adding XY position
stat.test.phylum.MOM <- stat.test.phylum.MOM %>% add_xy_position(x = "TimePoint")
stat.test.phylum.MOM

#Box plot
bxp.phylum.zibr.MOM <- ggboxplot(data= trans_phylum_zibr_rm_MOM, x = "TimePoint", y= "value", facet.by = "Phylum", xlab = "Time Points", ylab="Relative abundance", fill = "MOM_Group", width = 0.8, outlier.shape= NA)

#modifying the plot to use individual y-scale for each facet
bxp.phylum.zibr.MOM <- facet(p = bxp.phylum.zibr.MOM,facet.by = "Phylum",scales = "free_y") + scale_fill_manual(values = c("#E56878", "#D4DFE6")) +
  geom_point(aes(color = MOM_Group), position = position_dodge(0.7)) + stat_summary(
    fun = median,
    geom = 'line',
    aes(group = MOM_Group, colour = MOM_Group),
    position = position_dodge(width = 0.7) #this has to be added
  ) +
  scale_color_manual(values=c("#DC354B", "#A9BFCC")) + xlab("Time points") +
  ylab("Abundance (CLR)")

bxp.phylum.zibr.MOM + theme_bw() + theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank()) + stat_pvalue_manual(stat.test.phylum.MOM,label= "p.adj.signif", tip.length=0.02, hide.ns = T)

###For genus level

psnoncontam_tru_final_rare_rm

##Agglomerate at phyla level and extract the OTU table
genus_data_total_rm <- tax_glom(psnoncontam_tru_final_rare_rm,taxrank = "Genus")
genus_data_total_rm
write_phyloseq(genus_data_total_rm)

#Read the CSV file
abund_table <- read.csv("otu_table_genus_rm_clr_april_2021.csv", row.names=1, check.names= FALSE) #Added header "Features" to the first clumn of ASVs
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
genus_data_rm_zclr <- genus_data_total_rm
otu_table(genus_data_rm_zclr) <- otu_table(abund_clr, taxa_are_rows = F)
genus_data_rm_zclr

#replacing ASVs with short string
#taxa_names(genus_data_rm_zclr) <- paste0("SV", seq(ntaxa(genus_data_rm_zclr)))
#View(tax_table(genus_data_rm_zclr))

#Extracting short strings and their respective taxonomy
#genus_list<- tax_table(genus_data_rm_zclr)@.Data
#write.csv(genus_list, "genus_total_list.csv")

#extracting as csv and replacing string names with phyla names, transposed and added metadata
write_phyloseq(genus_data_rm_zclr, type= 'all')


#transposed, added additional metadata and importing the file
genus_zibr_rm.MOM <- read.table("otu_table_genus_rm_MOM_clr_new.csv", header = TRUE, sep= ",", check.names = F)
genus_zibr_rm.MOM

#Melting the data frame to use it for boxplot and statistics
trans_genus_zibr_rm.MOM <- melt(genus_zibr_rm.MOM)
trans_genus_zibr_rm.MOM

#Changing the column name from 'variable' to 'Genus'
colnames(trans_genus_zibr_rm.MOM)[4] <- "Genus"
trans_genus_zibr_rm.MOM

#Rstatix
stat.test.genus.MOM <- trans_genus_zibr_rm.MOM %>%
  group_by(Genus, TimePoint) %>%
  wilcox_test(value ~ MOM_Group) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()
stat.test.genus.MOM

#Exporting this post-hoc analysis data as a csv file
#write.csv(stat.test.genus.MOM, "genus_summary_posthoc_statistics_MOM.csv")

#Adding XY position
stat.test.genus.MOM <- stat.test.genus.MOM %>% add_xy_position(x = "TimePoint")
stat.test.genus.MOM

#Box plot
bxp.genus.zibr <- ggboxplot(data= trans_genus_zibr_rm.MOM, x = "TimePoint", y= "value", facet.by = "Genus", xlab = "Time Points", ylab="Abundance (CLR)", fill = "MOM_Group", width = 0.8, outlier.shape= NA)

#modifying the plot to use individual y-scale for each facet
bxp.genus.zibr <- facet(p = bxp.genus.zibr,facet.by = "Genus",scales = "free_y") + scale_fill_manual(values = c("#E56878", "#D4DFE6")) +
  geom_point(aes(color = MOM_Group), position = position_dodge(0.7)) + stat_summary(
    fun = median,
    geom = 'line',
    aes(group = MOM_Group, colour = MOM_Group),
    position = position_dodge(width = 0.7) #this has to be added
  ) +
  scale_color_manual(values=c("#DC354B", "#A9BFCC")) + xlab("Time points") +
  ylab("Abundance (CLR)")

bxp.genus.zibr + theme_bw() + theme(panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank()) + stat_pvalue_manual(stat.test.genus.MOM,label= "p.adj.signif", tip.length=0.02, hide.ns = T)
