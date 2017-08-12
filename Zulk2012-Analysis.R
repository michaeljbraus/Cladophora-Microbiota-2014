# Zulkifly 16S Amplicon RE-Anysis from 2012
# Last Run 2017-03-11 by MJB

# Make code output to pdf figures: 
#pdf("filename.pdf")
#{code, plots, etc.}
#dev.off()

# Turn console output into a text file. 
#sink("filename.txt")
#{code, calculations, etc.}
#sink()

# Write a table to a text file. 
#write.table(dataframe,"Figs/dataframe.txt",sep="\t",row.names=FALSE)

# Install Phyloseq: 
#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')

library(phyloseq)
library(ggplot2)
library(Rmisc)
library(vegan)
library(gridExtra)
library(dplyr)

setwd("~/Dropbox/Cladophora/Cladophora-Microbiota-2014/")

# Load biom file. 
biom <- import_biom("otu_table-Zulk2012.biom", parseFunction=parse_taxonomy_greengenes); biom 
sample_names(biom)
rank_names(biom)

# Subset only bacteria, remove singletons, remove plastids, and normalize by relative abundance. 
bact.withsingletons <- subset_taxa(biom, Kingdom=="Bacteria"); bact.withsingletons
bact <- subset_taxa(bact.withsingletons, taxa_sums(bact.withsingletons) > 1); bact
bact <- subset_taxa(bact, Class!="D_3__Chloroplast"); bact
bact.relabund <- transform_sample_counts(bact, function(x) x / sum(x)); bact.relabund

# Plot richness metrics. 
pdf("Zulk2012-Figs/plot_richness_raw.pdf")
plot_richness(bact) 
dev.off()

# Plot the classes.  
pdf("Zulk2012-Figs/class_relabund.pdf")
p = plot_bar(bact.relabund, fill="Class")
p + geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack") 
dev.off()

# Find classes_zulk of interest and subset bact.relabund. 
bact.relabund.melt <- psmelt(bact.relabund)
bact.relabund.melt.classzulk <- bact.relabund.melt%>%
  group_by(Sample, Class)%>%
  mutate(ClassAbundance = sum(Abundance))%>%
  distinct(Sample, ClassAbundance, Class)
#View(bact.relabund.melt.classzulk)
write.csv(bact.relabund.melt.classzulk, file = "Zulk2012.bact.relabund.class.csv", row.names = FALSE)
