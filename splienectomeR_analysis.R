###############################Longitudinal analysis using splinectomeR##########################################
### For repeated measures analysis on taxa, using the rarified phyloseq object and then removing the subjects for those paired samples at all time points are unavialble
psnoncontam_tru_final_rare

#Removing subjects that don't have samples at all time points
psnoncontam_tru_final_rare_rm <- subset_samples(psnoncontam_tru_final_rare, sample_data(psnoncontam_tru_final_rare)$SubjectID!= "W28" & sample_data(psnoncontam_tru_final_rare)$SubjectID!= "W30" & sample_data(psnoncontam_tru_final_rare)$SubjectID!= "W06")

##No need to rerarify the dataset after removing these subjects, as the sample with least number of reads is none of these subject
psnoncontam_tru_final_rare_rm

#filtering out taxa not present in the true samples
psnoncontam_tru_final_rare_rm <- prune_taxa(taxa_sums(psnoncontam_tru_final_rare_rm)>0, psnoncontam_tru_final_rare_rm)
psnoncontam_tru_final_rare_rm

psnoncontam_tru_final_rare_rm

##Agglomerate at phyla level and extract the OTU table
phylum_data_total_rm <- tax_glom(psnoncontam_tru_final_rare_rm,taxrank = "Phylum")
phylum_data_total_rm

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


##Extracting the abundance and taxa info as separate data frames
otus <- as.data.frame(otu_table(phylum_data_rm_zclr))
colnames(otus) <- as.character(tax_table(phylum_data_rm_zclr)[,'Phylum'])
tax_list <- colnames(otus)
otus$Sample_name <- rownames(otus) 

#Extracting the metadata as a separate data frame
metadata <- as.matrix(sample_data(phylum_data_rm_zclr))
results <- matrix(ncol = 2)

###Performing splinectomeR based longitudinal analysis
##Actino for RCT group
temp_df_rct_phylum <- data.frame(taxa=otus[,'Actinobacteria'],Sample_name=otus[,'Sample_name'])
colnames(temp_df_rct_phylum) <- c('Actinobacteria','Sample_name')
temp_df_rct_phylum <- merge(x = temp_df_rct_phylum,y = metadata,by.x = 'Sample_name',by.y = 'SampleName')

temp_df_rct_phylum$Description <- as.numeric(gsub(pattern = 'T',replacement = '',x = temp_df_rct_phylum$Description))

actino.rct.result <- permuspliner(data = temp_df_rct_phylum, xvar = 'Description',yvar = 'Actinobacteria', perms = 999,category = 'Group', cases = 'Sample_name',quiet = T)

p.actino.rct <- permuspliner.plot.permsplines(data = actino.rct.result,
                                              xvar = 'timepoint',
                                              yvar = 'Abundance (CLR)')
p.actino.rct <- p.actino.rct + scale_colour_manual(values = c("#80D7DD", "#F3DB80")) + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

p.actino.rct

actino_rct_dist <- permuspliner.plot.permdistance(actino.rct.result, xlabel = 'timepoint') + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

actino_rct_dist

actino.rct.result$pval


####Actino for MOM group



actino.mom.result <- permuspliner(data = temp_df_rct_phylum, xvar = 'Description',yvar = 'Actinobacteria', perms = 999,category = 'MOM_Group', cases = 'Sample_name',quiet = T)


p.actino.mom <- permuspliner.plot.permsplines(data = actino.mom.result,
                                              xvar = 'timepoint',
                                              yvar = 'Abundance (CLR)')
p.actino.mom<- p.actino.mom + scale_color_manual(values=c("#DC354B", "#A9BFCC")) + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))
p.actino.mom

actino_mom_dist <- permuspliner.plot.permdistance(actino.mom.result, xlabel = 'timepoint') + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

actino_mom_dist


actino.mom.result$pval



##Firmicutes for RCT group


temp_df_rct_phylum <- data.frame(taxa=otus[,'Firmicutes'],Sample_name=otus[,'Sample_name'])
colnames(temp_df_rct_phylum) <- c('Firmicutes','Sample_name')
temp_df_rct_phylum <- merge(x = temp_df_rct_phylum,y = metadata,by.x = 'Sample_name',by.y = 'SampleName')

temp_df_rct_phylum$Description <- as.numeric(gsub(pattern = 'T',replacement = '',x = temp_df_rct_phylum$Description))

firmi.rct.result <- permuspliner(data = temp_df_rct_phylum, xvar = 'Description',yvar = 'Firmicutes', perms = 999,category = 'Group', cases = 'Sample_name',quiet = T)

p.firmi.rct <- permuspliner.plot.permsplines(data = firmi.rct.result,
                                             xvar = 'timepoint',
                                             yvar = 'Abundance (CLR)')
p.firmi.rct <- p.firmi.rct + scale_colour_manual(values = c("#80D7DD", "#F3DB80")) + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

p.firmi.rct

firmi_rct_dist <- permuspliner.plot.permdistance(firmi.rct.result, xlabel = 'timepoint') + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

firmi_rct_dist

firmi.rct.result$pval



##Firmicutes for MOM group


firmi.mom.result <- permuspliner(data = temp_df_rct_phylum, xvar = 'Description',yvar = 'Firmicutes', perms = 999,category = 'MOM_Group', cases = 'Sample_name',quiet = T)


p.firmi.mom <- permuspliner.plot.permsplines(data = firmi.mom.result,
                                             xvar = 'timepoint',
                                             yvar = 'Abundance (CLR)')
p.firmi.mom<- p.firmi.mom + scale_color_manual(values=c("#DC354B", "#A9BFCC")) + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))
p.firmi.mom

firmi_mom_dist <- permuspliner.plot.permdistance(firmi.mom.result, xlabel = 'timepoint') + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

firmi_mom_dist


firmi.mom.result$pval



#Proteobacteria for RCT


temp_df_rct_phylum <- data.frame(taxa=otus[,'Proteobacteria'],Sample_name=otus[,'Sample_name'])
colnames(temp_df_rct_phylum) <- c('Proteobacteria','Sample_name')
temp_df_rct_phylum <- merge(x = temp_df_rct_phylum,y = metadata,by.x = 'Sample_name',by.y = 'SampleName')

temp_df_rct_phylum$Description <- as.numeric(gsub(pattern = 'T',replacement = '',x = temp_df_rct_phylum$Description))

proteo.rct.result <- permuspliner(data = temp_df_rct_phylum, xvar = 'Description',yvar = 'Proteobacteria', perms = 999,category = 'Group', cases = 'Sample_name',quiet = T)

p.proteo.rct <- permuspliner.plot.permsplines(data = proteo.rct.result,
                                              xvar = 'timepoint',
                                              yvar = 'Abundance (CLR)')
p.proteo.rct <- p.proteo.rct + scale_colour_manual(values = c("#80D7DD", "#F3DB80")) + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

p.proteo.rct

proteo_rct_dist <- permuspliner.plot.permdistance(proteo.rct.result, xlabel = 'timepoint') + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

proteo_rct_dist

proteo.rct.result$pval

##Proteobacteria for MOM group


proteo.mom.result <- permuspliner(data = temp_df_rct_phylum, xvar = 'Description',yvar = 'Proteobacteria', perms = 999,category = 'MOM_Group', cases = 'Sample_name',quiet = T)


p.proteo.mom <- permuspliner.plot.permsplines(data = proteo.mom.result,
                                              xvar = 'timepoint',
                                              yvar = 'Abundance (CLR)')
p.proteo.mom<- p.proteo.mom + scale_color_manual(values=c("#DC354B", "#A9BFCC")) + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))
p.proteo.mom

proteo_mom_dist <- permuspliner.plot.permdistance(proteo.mom.result, xlabel = 'timepoint') + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

proteo_mom_dist


proteo.mom.result$pval

##Combining the trend and distance plots at phylum level for RCT groups

ggarrange(plotlist = list(p.firmi.rct, p.proteo.rct, p.actino.rct, firmi_rct_dist, proteo_rct_dist, actino_rct_dist), common.legend = TRUE, ncol = 3, nrow = 2)

###Combining phyla plots for MOM groups

ggarrange(plotlist = list(p.firmi.mom, p.proteo.mom, p.actino.mom, firmi_mom_dist, proteo_mom_dist, actino_mom_dist), common.legend = TRUE, ncol = 3, nrow = 2)


######Simlar analysis was carried out at genus level and individual splinectome plots were generated and combined as presented in the figures in the manuscript.


