##Univariable analysis: covariate testing (Ref: https://jkzorz.github.io/2020/04/04/NMDS-extras.html)

###Performed univariate analysis at 4 time points by using the bray curtis distances on clr transformed abundances
#For first time point
##Testing only on important covriates for which known association with gut microbiome is known

#reading the metadata T1 file
df_t1 <- read.csv("metadata_T1_final_copy_revised.csv", sep=",", header = TRUE)
bf_env <- subset_samples(ps_bf_total_zclr_pseudo, SubjectID!= "W34" & SubjectID!= "W36" & SubjectID!= "W22") 

#ordination
bf_bray_ord <- ordinate(bf_env, "NMDS", "bray")

#EnvFit function
enT3 = envfit(bf_bray_ord, df_t3, permutations = 9999)
enT3

#View(enT3$vectors$pvals)

##Function for adjusting p-values from envfit

p.adjust.envfit <- function (x, method = 'fdr', n)
{
  x.new <- x
  if (!is.null (x$vectors)) pval.vectors <- x$vectors$pvals else pval.vectors <- NULL
  if (!is.null (x$factors)) pval.factors <- x$factors$pvals else pval.factors <- NULL
  if (missing (n)) n <- length (pval.vectors) + length (pval.factors)
  if (!is.null (x$vectors)) x.new$vectors$pvals <- p.adjust (x$vectors$pvals, method = method, n = n)
  if (!is.null (x$factors)) x.new$factors$pvals <- p.adjust (x$factors$pvals, method = method, n = n)
  cat ('Adjustment of significance by', method, 'method')
  return (x.new)
}

#adjusting p-values using the above function
adjusted.enT1 <- p.adjust.envfit(enT1, method = 'fdr')
#Extracting the data
ent1.factors <- adjusted.enT3$factors$r
write.csv(ent1.factors, "ent1_factors.csv")
ent1.vectors <- adjusted.enT1$vectors$r
write.csv(ent1.vectors, "ent1_vectors.csv")

adjusted.enT1

#Combined the data of vectors and factors, importing the file
#ent1.covariates <- read.csv("ent1_covariates.csv", header = TRUE, sep = ",")
#ent1.covariates
#Bar plot
#p<-ggplot(data=ent1.covariates, aes(x=Covariate, y=R2)) +
# geom_bar(stat="identity") + scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, #0.5, by = 0.3))
#p + coord_flip() +  theme_minimal() + ylab("Amount of variance (r2)")

##Time point T2

#reading the metadata T3 file
df_t2 <- read.csv("metadata_T2_final_copy_revised.csv", sep=",", header = TRUE)

df_env <- subset_samples(ps_df_total_zclr_pseudo, SubjectID!= "W20" & SubjectID!= "W22") 

#ordination
df_bray_ord <- ordinate(df_env, "NMDS", "bray")
#EnvFit function
enT2 = envfit(df_bray_ord, df_t2, permutations = 9999, na.rm = TRUE)
enT2

#adjusting p-values using the above function
adjusted.enT2 <- p.adjust.envfit(enT2, method = 'fdr')
adjusted.enT2

#Extracting the data
ent2.factors <- adjusted.enT2$factors$r
write.csv(ent2.factors, "ent2_factors.csv")
ent2.vectors <- adjusted.enT2$vectors$r
write.csv(ent2.vectors, "ent2_vectors.csv")


##########################Use only to plot separartely#################################
#Combined the data of vectors and factors, importing the file
ent2.covariates <- read.csv("ent2_covariates.csv", header = TRUE, sep = ",")
ent2.covariates
#Bar plot
p<-ggplot(data=ent2.covariates, aes(x=Covariate, y=R2)) +
  geom_bar(stat="identity") + scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, by = 0.3))
p + coord_flip() +  theme_minimal() + ylab("Amount of variance (r2)")
#######################################################################################


##Time point T3
#reading the metadata T3 file
df_t3 <- read.csv("metadata_T3_final_copy_revised.csv", sep=",", header = TRUE)
#ordination
ef_bray_ord <- ordinate(ps_ef_total_zclr_pseudo, "NMDS", "bray")
#EnvFit function
enT3 = envfit(ef_bray_ord, df_t3, permutations = 9999, na.rm = TRUE)
enT3

#adjusting p-values using the above function
adjusted.enT3 <- p.adjust.envfit(enT3, method = 'fdr')
adjusted.enT3

#Extracting the data
ent3.factors <- adjusted.enT3$factors$r
write.csv(ent3.factors, "ent3_factors.csv")
ent3.vectors <- adjusted.enT3$vectors$r
write.csv(ent3.vectors, "ent3_vectors.csv")

#######################Use only for plotting separately################################
#Combined the data of vectors and factors, importing the file
ent3.covariates <- read.csv("ent3_covariates.csv", header = TRUE, sep = ",")
ent3.covariates
#Bar plot
p<-ggplot(data=ent3.covariates, aes(x=Covariate, y=R2)) +
  geom_bar(stat="identity") + scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, by = 0.3))
p + coord_flip() +  theme_minimal() + ylab("Amount of variance (r2)")
#######################################################################################


##Time point T4

#reading the metadata T3 file
df_t4 <- read.csv("metadata_T4_final_copy_revised.csv", sep=",", header = TRUE)

ff_env <- subset_samples(ps_ff_total_zclr_pseudo, SubjectID!= "W04" & SubjectID!= "W39") 

#ordination
ff_bray_ord <- ordinate(ff_env, "NMDS", "bray")
#EnvFit function
enT4 = envfit(ff_bray_ord, df_t4, permutations = 9999, na.rm = TRUE)
enT4

#adjusting p-values using the above function
adjusted.enT4 <- p.adjust.envfit(enT4, method = 'fdr')
adjusted.enT4

#Extracting the data
ent4.factors <- adjusted.enT4$factors$r
write.csv(ent4.factors, "ent4_factors.csv")
ent4.vectors <- adjusted.enT4$vectors$r
write.csv(ent4.vectors, "ent4_vectors.csv")


########################Use only to plot separately####################################
#Combined the data of vectors and factors, importing the file
ent4.covariates <- read.csv("ent4_covariates.csv", header = TRUE, sep = ",")
ent4.covariates
#Bar plot
p<-ggplot(data=ent4.covariates, aes(x=Covariate, y=R2)) +
  geom_bar(stat="identity") + scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, by = 0.3))
p + coord_flip() +  theme_minimal() + ylab("Amount of variance (r2)")
#######################################################################################


###Exported results and combined the file
###Plotting envfit results

envresults <- read.csv("envfit_results.csv", header = TRUE, sep= ",")
envresults

#Ordering the dataframe
envresults_reorder <- envresults
envresults_reorder$Covariate <- factor(envresults_reorder$Covariate, levels= unique(envresults_reorder$Covariate[order(envresults_reorder$Factor)]))

#Bar plot
p<-ggplot(data=envresults_reorder, aes(x=Covariate, y=R2, fill=Factor)) +
  geom_bar(stat="identity", colour="black") + scale_y_continuous(limits = c(0, 0.4), breaks = seq(0, 0.4, by = 0.3)) 
p3 <- p + coord_flip() +  theme_minimal() + ylab("Amount of variance (r2)")
p2 <-facet(p3, facet.by = "TimePoint", nrow = 3) + scale_fill_brewer(palette = "Spectral") + theme(panel.background = element_rect(fill = "white", colour = "grey50"), panel.grid.minor = element_blank()) + theme(legend.position = "bottom")
p2
