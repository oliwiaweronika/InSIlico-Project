install.packages('BiocManager')
BiocManager::install("phyloseq")
library("phyloseq")
library("ggplot2")
library("scales")
install.packages("ggpubr")
library("ggpubr")
install.packages("vegan")
library("vegan")
library("dplyr")
BiocManager::install("DESeq2")
library(DESeq2)
install.packages("remotes")
remotes::install_github("vmikk/metagMisc")
library(metagMisc)


#define colors for each tretament to be used in most charts
genotypec<-c("#DC9C03","#965499","#00A087FF")
treatmentc<-c("#4DBBD5FF", "#E64B35FF")

otu1 <- otu_table(ASV, taxa_are_rows=TRUE)

sam1 <- sample_data(metadata) 

tax1 <- tax_table(taxonomy)
## tax1 <- data.table::as.data.table(taxonomy)

#make sure rows of taxa match the rows of OTU table
colnames(tax1) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
rownames(tax1)<- rownames(otu1)

#build a class data phyloseq object (named FLNA) that contains OTU table, taxonomy info and metadata info for all samples
FLNA <- phyloseq(otu1,tax1,sam1)

#check data to see if all data makes sense/readable
most_abundant_taxa <- sort(taxa_sums(FLNA), TRUE)[1:25]
ex2 <- prune_taxa(names(most_abundant_taxa), FLNA)
topFamilies <- tax_table(ex2)[, "Family"]
as(topFamilies, "vector")

#######################################################################################
#As a first analysis, we will look at the distribution of read counts from our samples
sample_sum_df <- data.frame(sum = sample_sums(FLNA))
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

# output mean, max and min of sample read counts
smin <- min(sample_sums(FLNA))
smin
smean <- mean(sample_sums(FLNA))
smean
ssd<- sd(sample_sums(FLNA))
ssd
smax <- max(sample_sums(FLNA))
smax

#calculate coverage for each sample - note: you can do this before or after rarefying the data
phyloseq_coverage(FLNA,correct_singletons = FALSE, add_attr = T)


#rarefy data (if desired; if not jump this step and normalize the easy way(below))

FLNA_rare<-rarefy_even_depth(FLNA, sample.size = min(sample_sums(FLNA)),
                                 rngseed = 771, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

# Convert phyloseq to a data frame
df <- psmelt(FLNA_rare)

# Calculate min, max, mean, and sd
min_rare <- min(df$Abundance, na.rm = TRUE)
max_rare <- max(df$Abundance, na.rm = TRUE)
mean_rare <- mean(df$Abundance, na.rm = TRUE)
sd_rare <- sd(df$Abundance, na.rm = TRUE)

# Print the results
print(paste("Min: ", min_rare))
print(paste("Max: ", max_rare))
print(paste("Mean: ", mean_rare))
print(paste("SD: ", sd_rare))

#quick look into the new OTU table
otu_FLNA_rare<-otu_table(FLNA_rare)
otu_FLNA_rare

write.csv(otu_FLNA_rare, file = "otu_FLNA_rare.csv")

#fix the order at which samples appear on charts: first control. etc..
sample_data(FLNA_rare)$genotype<- factor(sample_data(FLNA_rare)$genotype,levels = c("wt", "Q","R"))
sample_data(FLNA_rare)$treatment<- factor(sample_data(FLNA_rare)$treatment, levels=c("antibiotics", "none"))

#subset samples by genotype
FLNA_rare_wt<-subset_samples(FLNA_rare, genotype=="wt")
FLNA_rare_Q<-subset_samples(FLNA_rare, genotype=="Q")
FLNA_rare_R<-subset_samples(FLNA_rare, genotype=="R")

#alpha_diversity metrics by genotype
alpha_all<-plot_richness(FLNA_rare, x="genotype",color="genotype",measures=c("InvSimpson", "Shannon","Observed"))+scale_color_manual(values=genotypec)
alpha_all+geom_boxplot(lwd=1,width=0.5)+geom_point()+stat_compare_means(method = 'wilcox',label = 'p.signif',bracket.size = 0.1,label.y.npc = "top",tip.length=0.01,hide.ns = TRUE,label.x.npc = "center")+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust=1,size=9))

#alpha_diversity metrics by treatment
alpha_all<-plot_richness(FLNA_rare, x="treatment",color="treatment",measures=c("InvSimpson", "Shannon","Observed"))+scale_color_manual(values=treatmentc)
alpha_all+geom_boxplot(lwd=1,width=0.5)+geom_point()+stat_compare_means(method = 'wilcox',label = 'p.signif',bracket.size = 0.1,label.y.npc = "top",tip.length=0.01,hide.ns = TRUE,label.x.npc = "center",ref.group = "antibiotics")+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust=0.7,size=9))

#alpha_diversity metrics by treatment but only Q samples
alpha_all<-plot_richness(FLNA_rare_Q, x="treatment",color="treatment",measures=c("InvSimpson", "Shannon","Observed"))+scale_color_manual(values=treatmentc)
alpha_all+geom_boxplot(lwd=1,width=0.5)+geom_point()+stat_compare_means(method = 'wilcox',label = 'p.signif',bracket.size = 0.1,label.y.npc = "top",tip.length=0.01,hide.ns = TRUE,label.x.npc = "center",ref.group = "antibiotics")+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust=0.7,size=9))

#alpha_diversity metrics by treatment but only R samples
alpha_all<-plot_richness(FLNA_rare_R, x="treatment",color="treatment",measures=c("InvSimpson", "Shannon","Observed"))+scale_color_manual(values=treatmentc)
alpha_all+geom_boxplot(lwd=1,width=0.5)+geom_point()+stat_compare_means(method = 'wilcox',label = 'p.signif',bracket.size = 0.1,label.y.npc = "top",tip.length=0.01,hide.ns = TRUE,label.x.npc = "center",ref.group = "antibiotics")+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust=0.7,size=9))


#extract important data for a new alpha div simple table - treatment
richness.treatment.rare<- cbind(estimate_richness(FLNA_rare, 
                                        measures = c('Observed','Shannon',"Simpson","InvSimpson",'treatment')),
                      sample_data(FLNA)$treatment)

write.csv(richness.treatment.rare, file = "alphaDivFLNA_treat.csv")

colnames(richness.treatment.rare) <- c('Observed','Shannon',"Simpson","InvSimpson",'treatment')
richness.treatment.rare



#extract important data for a new alpha div simple table - genotype
richness.genotype.rare<- cbind(estimate_richness(FLNA_rare, 
                                                  measures = c('Observed','Shannon',"Simpson","InvSimpson",'genotype')),
                                sample_data(FLNA)$genotype)

write.csv(richness.genotype.rare, file = "alphaDivFLNA_geno.csv")


colnames(richness.genotype.rare) <- c('Observed','Shannon',"Simpson","InvSimpson",'genotype')
richness.genotype.rare

### summary statistics for alpha diversity measures - treatment 

summary_stats_treatment <- richness.treatment.rare %>%
  group_by(treatment) %>%
  summarise(
    min_observed = min(Observed, na.rm = TRUE),
    max_observed = max(Observed, na.rm = TRUE),
    median_observed = median(Observed, na.rm = TRUE),
    lower_quartile_observed = quantile(Observed, 0.25, na.rm = TRUE),
    upper_quartile_observed = quantile(Observed, 0.75, na.rm = TRUE),
    
    min_shannon = min(Shannon, na.rm = TRUE),
    max_shannon = max(Shannon, na.rm = TRUE),
    median_shannon = median(Shannon, na.rm = TRUE),
    lower_quartile_shannon = quantile(Shannon, 0.25, na.rm = TRUE),
    upper_quartile_shannon = quantile(Shannon, 0.75, na.rm = TRUE),
    
    min_InvSimpson = min(InvSimpson, na.rm = TRUE),
    max_InvSimpson = max(InvSimpson, na.rm = TRUE),
    median_InvSimpson = median(InvSimpson, na.rm = TRUE),
    lower_quartile_InvSimpson = quantile(InvSimpson, 0.25, na.rm = TRUE),
    upper_quartile_InvSimpson = quantile(InvSimpson, 0.75, na.rm = TRUE),
 
  )

print(summary_stats_treatment)


### summary statistics for alpha diversity measures - genotype

summary_stats_genotype <- richness.genotype.rare %>%
  group_by(genotype) %>%
  summarise(
    min_observed = min(Observed, na.rm = TRUE),
    max_observed = max(Observed, na.rm = TRUE),
    median_observed = median(Observed, na.rm = TRUE),
    lower_quartile_observed = quantile(Observed, 0.25, na.rm = TRUE),
    upper_quartile_observed = quantile(Observed, 0.75, na.rm = TRUE),
    
    min_shannon = min(Shannon, na.rm = TRUE),
    max_shannon = max(Shannon, na.rm = TRUE),
    median_shannon = median(Shannon, na.rm = TRUE),
    lower_quartile_shannon = quantile(Shannon, 0.25, na.rm = TRUE),
    upper_quartile_shannon = quantile(Shannon, 0.75, na.rm = TRUE),
    
    min_InvSimpson = min(InvSimpson, na.rm = TRUE),
    max_InvSimpson = max(InvSimpson, na.rm = TRUE),
    median_InvSimpson = median(InvSimpson, na.rm = TRUE),
    lower_quartile_InvSimpson = quantile(InvSimpson, 0.25, na.rm = TRUE),
    upper_quartile_InvSimpson = quantile(InvSimpson, 0.75, na.rm = TRUE),
    
  )

print(summary_stats_genotype)

#additional stats: note if data is normally distributed, you can use a t-test instead of wilcox - treat 
compare_means(Shannon ~ treatment,  data = richness.treatment.rare,
              ref.group = "antibiotics", method = "wilcox")
compare_means(InvSimpson ~ treatment,  data = richness.treatment.rare,
              ref.group = "antibiotics", method = "wilcox")
compare_means(Observed ~ treatment,  data = richness.treatment.rare,
              ref.group = "antibiotics", method = "wilcox")

#additional stats: note if data is normally distributed, you can use a t-test instead of wilcox - geno (Q)
compare_means(Shannon ~ genotype,  data = richness.genotype.rare,
              ref.group = "Q", method = "wilcox")
compare_means(InvSimpson ~ genotype,  data = richness.genotype.rare,
              ref.group = "Q", method = "wilcox")
compare_means(Observed ~ genotype,  data = richness.genotype.rare,
              ref.group = "Q", method = "wilcox")

#additional stats: note if data is normally distributed, you can use a t-test instead of wilcox - geno (R)
compare_means(Shannon ~ genotype,  data = richness.genotype.rare,
              ref.group = "R", method = "wilcox")
compare_means(InvSimpson ~ genotype,  data = richness.genotype.rare,
              ref.group = "R", method = "wilcox")
compare_means(Observed ~ genotype,  data = richness.genotype.rare,
              ref.group = "R", method = "wilcox")

###################################Calculate rel abundance and filter out some data###########################

#calculate relative abundance 

FLNA_rare_Rel  = transform_sample_counts(FLNA_rare, function(x) x / sum(x) )

ASV_Rel = as(otu_table(FLNA_rare_Rel), "matrix")
sample_data(FLNA_rare_Rel)$ID<- factor(sample_data(FLNA_rare_Rel)$Mouse_ID,levels = c("4003","4004","4038","4039","4037","4040","4044",
                               "4046","4022","4020","4021","4023","7659","7703","7643","7647","7653","7657","7656","7654","7658",
                               "7704","7644","7646","7701","7700","7702","7757","7623","7622","7650","7651","7707","7706","7758",
                               "7636","7633","7769","7767","7766","7649","7648"))

######################################STACKED BARPLOTS##########################################################

# Barplots/stacked plots of phylum
genus_colors <- c(
  "#CBD588", "#5F7FC7", "orange", "#508578", "#DA5724","#D14285","#CD9BCD",
  "#AD6F3B","#5F7FC7", "#a258b0", "#652926", "#C84248", 
  "#5E738F","#D1A33D", "#8569D5","#599861","#5fb3c7","Grey40"
)

family_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)

# Find out what our top 10 genus and family are 

genus_count <- table(taxonomy$Genus)
sorted_genus <- sort(genus_count, decreasing = TRUE)
top_10_genus <- head(sorted_genus, 10)
print(top_10_genus)

family_count <- table(taxonomy$Family)
sorted_family <- sort(family_count, decreasing = TRUE)
top_10_family <- head(sorted_family, 10)
print(top_10_family)

# Assign colours to the 10 most common genus and family
# i thought thats what it was but now not sure and tbh not even sure if i need it?

colours_htot<-c("Bacillaceae"="#66C1A5", "Beijerinckiaceae"="coral3", "Sericytochromatia_unclassified"="#D4D943","Brevibacteriaceae"="#CF9797","Burkholderiaceae"="#8C9FCB", "Comamonadaceae"="#BD96C6", "Corynebacteriaceae"="cadetblue3","Enterobacteriaceae"="coral2","Lactobacillaceae"="#E68AC3", "Mitochondria"="#CDB48F" , "Moraxellaceae"="#A5D853","Other"="darkcyan", "Paenibacillaceae"="darkgoldenrod", "Pseudomonadaceae"="chocolate","Sphingomonadaceae"="coral1","Staphylococcaceae"="#CDBBA3","Streptococcaceae"="#B2B2B2","Unclassified_unclassified"="#F2CD6B","Vibrionaceae"= "#F2CD6B","Xanthomonadaceae"= "cadetblue2","unclassified_Bacteroidales"="chocolate","Coriobacteriaceae"="darkcyan","Clostridiales_Incertae_Sedis_XIII"="darkgoldenrod","Coriobacteriaceae"="azure2","Micrococcaceae"="#F2CD6B","unclassified_Firmicutes"="darkcyan", "Micrococcales_unclassified"="azure4")
colours_htot_phyla<-c("Actinobacteriota"="#66C1A5", "Firmicutes"="coral2", "Proteobacteria"="darkcyan","Bacteroidota"="coral3","Other"="#8C9FCB", "Patescibacteria"="#F2CD6B","Cyanobacteria"="darkgoldenrod","Spirochaetota"="#E68AC3","Unclassified_unclassified"="#CDB48F")

######################################Family level########################################
FLNA_Rel_family <- tax_glom(FLNA_rare_Rel, taxrank = 'Family',NArm=FALSE) # agglomerate taxa
FamilyRel <- psmelt(FLNA_Rel_family ) # create dataframe from phyloseq object
FamilyRel$Family <- as.character(FamilyRel$Family) #convert to character
FamilyRel$Family[FamilyRel$Abundance < 0.05] <- "Other" #rename taxa with < 5% abundance

ggplot(FamilyRel, aes(x =ID, y = Abundance, fill = Family))+geom_bar(stat="identity") +ylab("Relative abundance")+xlab("")+theme_bw()+theme(legend.text = element_text(size=8),legend.title = element_text(size=10),axis.title.y=element_text(size=12),axis.title.x = element_text(size=12),axis.text.y = element_text(size=8),axis.text.x = element_text(angle = 45, hjust=1,size=7),strip.text = element_text(size = 10))+facet_wrap(genotype ~ treatment,scales= "free_x",nrow=1)#+scale_fill_manual(values=colours_htot)

#stats to analyse
analyze_family <- function(Family) {
  avg_abundance <- aggregate(Abundance ~ Family + genotype + treatment, FamilyRel, mean)
  avg_abundance[avg_abundance$Family == Family, ]
}

families_present <- unique(FamilyRel$Family)
results_family <- lapply(families_present, analyze_family)
results_family

##################################Phylum Level###################################

FLNA_Rel_phyla <- tax_glom(FLNA_rare_Rel, taxrank = 'Phylum') # agglomerate taxa
PhylaRel<- psmelt(FLNA_Rel_phyla) # create dataframe from phyloseq object
PhylaRel$Phylum<- as.character(PhylaRel$Phylum) #convert to character
PhylaRel$Phylum[PhylaRel$Abundance < 0.01] <- "Other" #rename taxa with < 1% abundance

ggplot(PhylaRel, aes(x =ID, y = Abundance, fill = Phylum))+geom_bar(stat="identity") +ylab("Relative abundance")+xlab("")+theme_bw()+theme(legend.text = element_text(size=8),legend.title = element_text(size=10),axis.title.y=element_text(size=12),axis.title.x = element_text(size=12),axis.text.y = element_text(size=8),axis.text.x = element_text(angle = 45, hjust=1,size=7),strip.text = element_text(size = 10))+facet_wrap(genotype ~ treatment,scales= "free_x",nrow=1)#+scale_fill_manual(values=colours_htot_phyla)

#stats to analyse
analyze_phylum <- function(Phylum) {
  avg_abundance <- aggregate(Abundance ~ Phylum + genotype + treatment, PhylaRel, mean)
  avg_abundance[avg_abundance$Phylum == Phylum, ]
}

phylum_present <- unique(PhylaRel$Phylum)
results_phylum <- lapply(phylum_present, analyze_phylum)
results_phylum

######################################Genus level########################################
FLNA_Rel_genus <- tax_glom(FLNA_rare_Rel, taxrank = 'Genus',NArm=FALSE) # agglomerate taxa
GenusRel <- psmelt(FLNA_Rel_genus ) # create dataframe from phyloseq object
GenusRel$Genus <- as.character(GenusRel$Genus) #convert to character
GenusRel$Genus[GenusRel$Abundance < 0.05] <- "Other" #rename taxa with < 5% abundance

ggplot(GenusRel, aes(x =ID, y = Abundance, fill = Genus))+geom_bar(stat="identity") +ylab("Relative abundance")+xlab("")+theme_bw()+theme(legend.text = element_text(size=8),legend.title = element_text(size=10),axis.title.y=element_text(size=12),axis.title.x = element_text(size=12),axis.text.y = element_text(size=8),axis.text.x = element_text(angle = 45, hjust=1,size=7),strip.text = element_text(size = 10))+facet_wrap(genotype ~ treatment,scales= "free_x",nrow=1)#+scale_fill_manual(values=colours_htot)

#stats to analyse
analyze_genus <- function(Genus) {
  avg_abundance <- aggregate(Abundance ~ Genus + genotype + treatment, GenusRel, mean)
  avg_abundance[avg_abundance$Genus == Genus, ]
}

genus_present <- unique(GenusRel$Genus)
results_genus <- lapply(genus_present, analyze_genus)
results_genus

####################################BETA DIVERSITY ANALYSIS#########################################################
#calculate Bray Curtis dissimilarity matrix
bc_dist = phyloseq::distance(FLNA_rare, method="bray", weighted=F)
#PERMANOVA test to see if there is a cluster by group (in this case the treatment)
adonis2(bc_dist ~ sample_data(FLNA_rare)$treatment)
bc_dist


#ordination plot (NMDS) for mucosa vs digesta with centroids and elipses
## im not sure if this is neant to be sites for treatment or genotype and species for treatment or genotype 

ord.meta <- data.frame(sample_data(FLNA_rare))
ord.asvs<-otu_table(FLNA_rare)
ord.asvst<-t(ord.asvs)

bc.nmds <- metaMDS(ord.asvst,distance = "bray",
                   autotransform = T)


### based on treatment ###
site.scrs <- as.data.frame(scores(bc.nmds, display = "sites"))
site.scrs <- cbind(site.scrs, sites = ord.meta$treatment)
head(site.scrs)

NMDS.mean=aggregate(site.scrs[,c("NMDS1", "NMDS2")], 
                    list(group = site.scrs$sites), mean)

NMDS.mean
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell.dune.treatment <- data.frame() #sets up a data frame before running the function.
for(g in levels(site.scrs$sites)){
  df_ell.dune.treatment <- rbind(df_ell.dune.treatment, cbind(as.data.frame(with(site.scrs [site.scrs$site==g,],
                                                                                 veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))) ,site=g))
}

site.scrs$sites <- as.factor(site.scrs$sites)
sites <- c("#4DBBD5FF", "#E64B35FF")

#draw the NMDS ordination plot with elipses for each dataset
NMDS_plot<-ggplot(data = site.scrs, aes(NMDS1, NMDS2)) + geom_point(aes(color = sites),show.legend=NA,size=1.5)+scale_color_manual(values=sites)
NMDS_plot

NMDS_plot_elipse <- ggplot(data = site.scrs, aes(NMDS1, NMDS2)) +
  geom_point(aes(color = sites), show.legend=NA, size=1.5) +
  scale_color_manual(values=sites) +
  stat_ellipse(data = site.scrs, aes(NMDS1, NMDS2, color = sites), type = "t", level = 0.95)
NMDS_plot_elipse

### NMDS_plot_elipse<-NMDS_plot+geom_path(data= df_ell.dune.treatment, aes(x=NMDS1,y=NMDS2,colour=sites), size=1, linetype=2)+theme_bw()+annotate("text",x=NMDS.mean$MDS1,y=NMDS.mean$MDS2,label=NMDS.mean$sites)

### based on genotype ###
geno.scrs <- as.data.frame(scores(bc.nmds, display = "sites"))
geno.scrs <- cbind(geno.scrs, sites = ord.meta$genotype)
head(geno.scrs)

NMDS.mean=aggregate(geno.scrs[,c("NMDS1", "NMDS2")], 
                    list(group = geno.scrs$sites), mean)

NMDS.mean
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell.dune.treatment <- data.frame() #sets up a data frame before running the function.
for(g in levels(geno.scrs$sites)){
  df_ell.dune.treatment <- rbind(df_ell.dune.treatment, cbind(as.data.frame(with(site.scrs [site.scrs$site==g,],
                                                                                 veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))) ,site=g))
}

geno.scrs$sites <- as.factor(geno.scrs$sites)
sites <- c("#DC9C03","#965499","#00A087FF")

#draw the NMDS ordination plot with elipses for each dataset
NMDS_plot<-ggplot(data = geno.scrs, aes(NMDS1, NMDS2)) + geom_point(aes(color = sites),show.legend=NA,size=1.5)+scale_color_manual(values=sites)
NMDS_plot

NMDS_plot_elipse <- ggplot(data = geno.scrs, aes(NMDS1, NMDS2)) +
  geom_point(aes(color = sites), show.legend=NA, size=1.5) +
  scale_color_manual(values=sites) +
  stat_ellipse(data = geno.scrs, aes(NMDS1, NMDS2, color = sites), type = "t", level = 0.95)
NMDS_plot_elipse

###########DESEQ2 anaysis#################


###################Comparison treatments: Q vs R in None###################

FLNA_prune = prune_taxa(taxa_sums(FLNA) > 10, FLNA)

FLNA_prune_treat<- subset_samples(FLNA_prune, treatment %in% c("none")) 

sample_data(FLNA_prune_treat)$genotype<- factor(sample_data(FLNA_prune_treat)$genotype,levels = c("Q", "wt","R"))

treatdss = phyloseq_to_deseq2(FLNA_prune_treat, ~genotype)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(treatdss), 1, gm_mean)
treatdss = estimateSizeFactors(treatdss, geoMeans = geoMeans)
treatdss = DESeq(treatdss, fitType="local")

res = results(treatdss)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab
sigtab = cbind(as(sigtab, "data.frame"), 
               as(tax_table(FLNA_prune)[rownames(sigtab), ], "matrix"))
sigtabgen = subset(sigtab, !is.na(Genus))

write.csv(sigtab, file = "DeSeq2_Q_R_none.csv")

x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))

x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))

x = tapply(sigtabgen$log2FoldChange, sigtabgen$Order, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Order = factor(as.character(sigtabgen$Order), levels=names(x))


color_vector <- rainbow(11)

qr1<-ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Family)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+theme_bw()+scale_color_manual(values=color_vector)
qr1

qr2<-ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+theme_bw()+scale_color_manual(values=color_vector)
qr2

###################Comparison treatments: wt vs R in None###################

FLNA_prune = prune_taxa(taxa_sums(FLNA) > 10, FLNA)

FLNA_prune_treat<- subset_samples(FLNA_prune, treatment %in% c("none")) 

sample_data(FLNA_prune_treat)$genotype<- factor(sample_data(FLNA_prune_treat)$genotype,levels = c("wt", "Q","R"))

treatdss = phyloseq_to_deseq2(FLNA_prune_treat, ~genotype)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(treatdss), 1, gm_mean)
treatdss = estimateSizeFactors(treatdss, geoMeans = geoMeans)
treatdss = DESeq(treatdss, fitType="local")

res = results(treatdss)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab
sigtab = cbind(as(sigtab, "data.frame"), 
               as(tax_table(FLNA_prune)[rownames(sigtab), ], "matrix"))
sigtabgen = subset(sigtab, !is.na(Genus))

write.csv(sigtab, file = "DeSeq2_WT_R_None.csv")

x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))

x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))

x = tapply(sigtabgen$log2FoldChange, sigtabgen$Order, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Order = factor(as.character(sigtabgen$Order), levels=names(x))


color_vector <- rainbow(17)

wtr1<-ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Family)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+theme_bw()+scale_color_manual(values=color_vector)
wtr1

color_vector <- rainbow(5)

wtr2<-ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+theme_bw()+scale_color_manual(values=color_vector)
wtr2

###################Comparison treatments: wt vs Q in None###################

FLNA_prune = prune_taxa(taxa_sums(FLNA) > 10, FLNA)

FLNA_prune_treat<- subset_samples(FLNA_prune, treatment %in% c("none")) 

sample_data(FLNA_prune_treat)$genotype<- factor(sample_data(FLNA_prune_treat)$genotype,levels = c("wt", "R","Q"))

treatdss = phyloseq_to_deseq2(FLNA_prune_treat, ~genotype)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(treatdss), 1, gm_mean)
treatdss = estimateSizeFactors(treatdss, geoMeans = geoMeans)
treatdss = DESeq(treatdss, fitType="local")

res = results(treatdss)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab
sigtab = cbind(as(sigtab, "data.frame"), 
               as(tax_table(FLNA_prune)[rownames(sigtab), ], "matrix"))
sigtabgen = subset(sigtab, !is.na(Genus))

write.csv(sigtab, file = "DeSeq2_WT_Q_None.csv")

x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))

x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))

x = tapply(sigtabgen$log2FoldChange, sigtabgen$Order, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Order = factor(as.character(sigtabgen$Order), levels=names(x))


color_vector <- rainbow(8)

wtq1<-ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Family)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+theme_bw()+scale_color_manual(values=color_vector)
wtq1

color_vector <- rainbow(3)

wtq2<-ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+theme_bw()+scale_color_manual(values=color_vector)
wtq2

###################Comparison treatments: Q vs R in Antibiotics###################
################### don't include figures in diss but include data in main text ###########

FLNA_prune = prune_taxa(taxa_sums(FLNA) > 10, FLNA)

FLNA_prune_treat<- subset_samples(FLNA_prune, treatment %in% c("antibiotics")) 

sample_data(FLNA_prune_treat)$genotype<- factor(sample_data(FLNA_prune_treat)$genotype,levels = c("Q", "wt","R"))

treatdss = phyloseq_to_deseq2(FLNA_prune_treat, ~genotype)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(treatdss), 1, gm_mean)
treatdss = estimateSizeFactors(treatdss, geoMeans = geoMeans)
treatdss = DESeq(treatdss, fitType="local")

res = results(treatdss)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab
sigtab = cbind(as(sigtab, "data.frame"), 
               as(tax_table(FLNA_prune)[rownames(sigtab), ], "matrix"))
sigtabgen = subset(sigtab, !is.na(Genus))

write.csv(sigtab, file = "DeSeq2_Q_R_Treat.csv")

x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))

x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))

x = tapply(sigtabgen$log2FoldChange, sigtabgen$Order, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Order = factor(as.character(sigtabgen$Order), levels=names(x))


color_vector <- rainbow(2)

qra1<-ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Family)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+theme_bw()+scale_color_manual(values=color_vector)
qra1

qra2<-ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+theme_bw()+scale_color_manual(values=color_vector)
qra2



