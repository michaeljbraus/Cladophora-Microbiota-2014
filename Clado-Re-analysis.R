# Cladophora 16S Amplicon Anysis, 2014
# Last Run 2017-03-11 by MJB

# Make code output to pdf figures: 
#pdf("filename.pdf")
#{code, plots, etc.}
#dev.off()

# Turn console output into a text file. 
#sink("filename.txt")
#{code, calculations, etc.}
#sink()

# Install Phyloseq: 
#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')

library(phyloseq)
library(ggplot2)
library(Rmisc)
library(vegan)

setwd("~/Dropbox/Cladophora/Cladophora-Microbiota-2014/")

# Load biom file. 
biom <- import_biom("clado_otu_table.biom", parseFunction=parse_taxonomy_default); biom 
sample_names(biom)
rank_names(biom)

# Load and merge sample metadata with read data (biom).
sam.data <- read.csv(file="sample.data.csv", row.names=1, header=TRUE)
sample_data(biom) <- sam.data; sample_data(biom) 

# Subset only bacteria, remove singletons, remove plastids, and normalize by relative abundance. 
bact.withsingletons <- subset_taxa(biom, Rank1=="k__Bacteria"); bact.withsingletons
bact <- subset_taxa(bact.withsingletons, taxa_sums(bact.withsingletons) > 1); bact
bact <- subset_taxa(bact, Rank3!="c__Chloroplast"); bact
bact.relabund <- transform_sample_counts(bact, function(x) x / sum(x)); bact.relabund

# Plot richness metrics. 
plot_richness(bact, color = "Date", shape = "Site") 
bact.rich.est <- estimate_richness(bact, measures = NULL)
bact.rich.est$SampleID.1 <- row.names(bact.rich.est)
bact.rich.est <- merge(bact.rich.est, sam.data, by = "SampleID.1")
bact.rich.est$Site <- as.factor(bact.rich.est$Site)
bact.rich.est$Date <- as.numeric(bact.rich.est$Date)
head(bact.rich.est)  

# Plot Shannon index richness. 
pdf("Clado-Re-analysis-Figs/richness_shannon.pdf", width=6, height=3)
bact.rich.est.sha <- summarySE(bact.rich.est, measurevar="Shannon", groupvars=c("Date","Site")); bact.rich.est.sha
p.sha <- ggplot(bact.rich.est.sha, aes(x=Date, y=Shannon)) + facet_grid(~Site) +
  geom_point(size = 2) +  geom_errorbar(aes(ymin=Shannon-se, ymax=Shannon+se)) 
p.sha + theme_bw() +  
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust=1),axis.text.y = element_text(size = 10))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #+ 
# geom_line(lm(data = bact.rich.est.sha), color = ) + 
# geom_line(lm(data = bact.rich.est.sha)) + 
# geom_line(lm(data = bact.rich.est.sha))
dev.off()

# PERMANOVA. 
df = as(sample_data(bact.relabund), "data.frame")
d = phyloseq::distance(bact.relabund, "bray") 
clado.adonis = adonis(d ~ Site*Date, df)
sink("Clado-Re-analysis-Figs/permanova.txt")
clado.adonis
sink()

# Ordination plot. Maybe tweak k value. 
pdf("Clado-Re-analysis-Figs/nmds_plot.pdf")
ordNMDS <- ordinate(bact.relabund, method="NMDS", distance="bray", k=2)  
ord <- plot_ordination(bact.relabund, ordNMDS, color = "Date", shape = "Site") + geom_point(size=3)
ord +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())   
dev.off()

# Plot the phyla of the top 5000 OTUs. 
pdf("Clado-Re-analysis-Figs/phyla_top5000_otus.pdf")
TopNOTUs <- names(sort(taxa_sums(bact.relabund), TRUE)[1:5000])
bact5000   <- prune_taxa(TopNOTUs, bact.relabund)
p = plot_bar(bact5000, fill="Rank2")
p + geom_bar(aes(color=Rank2, fill=Rank2), stat="identity", position="stack") 
dev.off()

# Plot the class of the top 500 OTUs. 
pdf("Clado-Re-analysis-Figs/class_top500_otus.pdf")
TopNOTUs <- names(sort(taxa_sums(bact.relabund), TRUE)[1:500])
bact500   <- prune_taxa(TopNOTUs, bact.relabund)
p = plot_bar(bact500, fill="Rank3")
p + geom_bar(aes(color=Rank3, fill=Rank3), stat="identity", position="stack") 
dev.off()

# Plot the phyla of the top 1000 OTUs. 
pdf("Clado-Re-analysis-Figs/phyla_top1000_otus.pdf")
TopNOTUs <- names(sort(taxa_sums(bact.relabund), TRUE)[1:1000])
bact1000   <- prune_taxa(TopNOTUs, bact.relabund)
p = plot_bar(bact1000, fill="Rank2", facet_grid=~Site)
p + geom_bar(aes(color=Rank2, fill=Rank2), stat="identity", position="stack") 
dev.off()

# Plot the phyla of the top 100 OTUs. 
pdf("Clado-Re-analysis-Figs/phyla_top100_otus.pdf")
TopNOTUs <- names(sort(taxa_sums(bact.relabund), TRUE)[1:100])
bact100   <- prune_species(TopNOTUs, bact.relabund)
p = plot_bar(bact100, fill="Rank2", facet_grid=~Site)
p + geom_bar(aes(color=Rank2, fill=Rank2), stat="identity", position="stack") 
dev.off()

# Plot the genera of the top 1000 OTUs. 
pdf("Clado-Re-analysis-Figs/genera_top1000_otus.pdf", width = 12, height = 8)
TopNOTUs <- names(sort(taxa_sums(bact.relabund), TRUE)[1:1000])
bact1000 <- prune_taxa(TopNOTUs, bact.relabund)
p = plot_bar(bact1000, "Date", fill="Rank6", facet_grid=~Site)
p + geom_bar(aes(color=Rank6, fill=Rank6), stat="identity", position="stack") 
dev.off()
