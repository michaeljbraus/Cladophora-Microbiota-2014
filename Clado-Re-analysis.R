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

setwd("~/Dropbox/Github/Cladophora-Microbiota-2014/")

# Load biom file. 
biom <- import_biom("clado_otu_table_SILVA.biom", parseFunction=parse_taxonomy_greengenes); biom 
sample_names(biom)
rank_names(biom)

# Load and merge sample metadata with read data (biom).
sam.data <- read.csv(file="sample.data.csv", row.names=1, header=TRUE)
sam.data$Date <- as.factor(sam.data$Date)
sample_data(biom) <- sam.data; str(sample_data(biom))

# Subset only bacteria, remove singletons, remove plastids, and normalize by relative abundance. 
bact.withsingletons <- subset_taxa(biom, Rank1=="D_0__Bacteria"); bact.withsingletons
bact <- subset_taxa(bact.withsingletons, taxa_sums(bact.withsingletons) > 1); bact
bact <- subset_taxa(bact, Rank3!="D_3__Chloroplast"); bact
bact.relabund <- transform_sample_counts(bact, function(x) x / sum(x)); bact.relabund

# Plot richness metrics. 
pdf("Figs/plot_richness_raw.pdf")
plot_richness(bact, color = "Date", shape = "Site") 
dev.off()

# Dataframe bacterial richness metrics
bact.rich.est <- estimate_richness(bact, measures = NULL)
bact.rich.est$SampleID.1 <- row.names(bact.rich.est)
bact.rich.est <- merge(bact.rich.est, sam.data, by = "SampleID.1")
bact.rich.est$Site <- as.factor(bact.rich.est$Site)
bact.rich.est$Date <- as.numeric(as.character(bact.rich.est$Date)); bact.rich.est$Date
write.table(bact.rich.est,"Figs/bact.rich.est-CHECK.txt",sep="\t",row.names=FALSE)
# Linear regressions of diversity over time at each site. 
bact.rich.est.sha <- summarySE(bact.rich.est, measurevar="Shannon", groupvars=c("Date","Site")); bact.rich.est.sha
north <- subset(bact.rich.est.sha, Site=="North")
point <- subset(bact.rich.est.sha, Site=="Point")
south <- subset(bact.rich.est.sha, Site=="South")
sink("Figs/diversity-linear-regressions-CHECK.txt")
lm.north <- lm(north$Shannon ~ north$Date); lm.north; summary(lm.north)
lm.point <- lm(point$Shannon ~ point$Date); lm.point; summary(lm.point)
lm.south <- lm(south$Shannon ~ south$Date); lm.south; summary(lm.south)
sink()

# Plot Shannon index richness. 
pdf("Figs/richness_shannon-CHECK.pdf", width=6, height=3)
p.sha <- ggplot(bact.rich.est.sha, aes(x=Date, y=Shannon)) + facet_grid(~Site) +
  geom_point(size = 2) +  geom_errorbar(aes(ymin=Shannon-se, ymax=Shannon+se)) 
p.sha + theme_bw() + geom_smooth(method="lm", se = FALSE) +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust=1),axis.text.y = element_text(size = 10))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

# PERMANOVA. 
df = as(sample_data(bact.relabund), "data.frame")
d = phyloseq::distance(bact.relabund, "bray") 
clado.adonis = adonis(d ~ Site*Date, df)
sink("Figs/permanova.txt")
clado.adonis
sink()

# Ordination. Maybe tweak k value. 
sink("Figs/ordination_stress_k2.txt")
ordNMDS.k2 <- ordinate(bact.relabund, method="NMDS", distance="bray", k=2)  
sink()
sink("Figs/ordination_stress_k3.txt")
ordNMDS.k3 <- ordinate(bact.relabund, method="NMDS", distance="bray", k=3)  
sink()
# Plot ordination of stress k = 3
pdf("Figs/nmds_plot_messy_k2.pdf", width = 5, height = 4)
ord.k2 <- plot_ordination(bact.relabund, ordNMDS.k2, shape="Site", color = "Date") + geom_point(size=5)
ord.k2 + theme_bw() + scale_colour_hue(h=c(300, 500))+
  geom_point(colour="white", size = 3)+
  geom_point(colour="black", size = 1)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())   
dev.off()
# Plot ordination of stress k = 3.
pdf("Figs/nmds_plot_messy_k3.pdf", width = 5, height = 4)
ord.k3 <- plot_ordination(bact.relabund, ordNMDS.k3, shape="Site", color = "Date") + geom_point(size=5)
ord.k3 + theme_bw() + scale_colour_hue(h=c(300, 500))+
  geom_point(colour="white", size = 3)+
  geom_point(colour="black", size = 1)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())   
dev.off()

# Plot the phyla of the top 5000 OTUs. 
pdf("Figs/phyla_top5000_otus.pdf")
TopNOTUs <- names(sort(taxa_sums(bact.relabund), TRUE)[1:5000])
bact5000   <- prune_taxa(TopNOTUs, bact.relabund)
p = plot_bar(bact5000, fill="Rank2")
p + geom_bar(aes(color=Rank2, fill=Rank2), stat="identity", position="stack") 
dev.off()

# Plot the class of the top 500 OTUs. 
pdf("Figs/class_top500_otus.pdf")
TopNOTUs <- names(sort(taxa_sums(bact.relabund), TRUE)[1:500])
bact500   <- prune_taxa(TopNOTUs, bact.relabund)
p = plot_bar(bact500, fill="Rank3")
p + geom_bar(aes(color=Rank3, fill=Rank3), stat="identity", position="stack") 
dev.off()

# Plot the phyla of the top 1000 OTUs. 
pdf("Figs/phyla_top1000_otus.pdf")
TopNOTUs <- names(sort(taxa_sums(bact.relabund), TRUE)[1:1000])
bact1000   <- prune_taxa(TopNOTUs, bact.relabund)
p = plot_bar(bact1000, fill="Rank2", facet_grid=~Site)
p + geom_bar(aes(color=Rank2, fill=Rank2), stat="identity", position="stack") 
dev.off()

# Plot the phyla of the top 100 OTUs. 
pdf("Figs/phyla_top100_otus.pdf")
TopNOTUs <- names(sort(taxa_sums(bact.relabund), TRUE)[1:100])
bact100   <- prune_species(TopNOTUs, bact.relabund)
p = plot_bar(bact100, fill="Rank2", facet_grid=~Site)
p + geom_bar(aes(color=Rank2, fill=Rank2), stat="identity", position="stack") 
dev.off()

# Plot the genera of the top 1000 OTUs. 
pdf("Figs/genera_top1000_otus.pdf", width = 12, height = 8)
TopNOTUs <- names(sort(taxa_sums(bact.relabund), TRUE)[1:1000])
bact1000 <- prune_taxa(TopNOTUs, bact.relabund)
p = plot_bar(bact1000, "Date", fill="Rank6", facet_grid=~Site)
p + geom_bar(aes(color=Rank6, fill=Rank6), stat="identity", position="stack") 
dev.off()

# Find top 30 genera and subset bact.relabund. 
sort.genera <- sort(tapply(taxa_sums(bact.relabund), tax_table(bact.relabund)[, "Rank6"], sum), TRUE)
top.genera <- sort.genera[1:30]
top.genera.list <- names(top.genera)
bact.relabund.subset = subset_taxa(bact.relabund, Rank6 %in% top.genera.list)
bact.relabund.subset.taxa <- subset_taxa(bact.relabund.subset, Rank6 %in% as.factor(top.genera.list))
bact.relabund.subset.taxa
relabund.top.genera <- psmelt(bact.relabund.subset.taxa)
relabund.top.genera.genus <- relabund.top.genera%>%
  group_by(Sample, Rank6)%>%
  mutate(GenusAbundance = sum(Abundance))%>%
  distinct(Sample, GenusAbundance, TreatmentGroup, Site, Date, Rank2, Rank6, Rank6)
head(relabund.top.genera.genus) 
# Summary of genus abundance of top 30 genera. 
relabund.top.genera.genus.est <- summarySE(relabund.top.genera.genus, measurevar="GenusAbundance", groupvars=c("Site","Date", "Rank6"))
head(relabund.top.genera.genus.est)
relabund.top.genera.genus.est$Date <- as.character(relabund.top.genera.genus.est$Date)
relabund.top.genera.genus.est$Date <- as.numeric(relabund.top.genera.genus.est$Date) 
# Plot summary of genus abundance of top 30 genera. 
pdf("Figs/top30genera.pdf", height = 20, width = 11)
p <- ggplot(relabund.top.genera.genus.est, aes(x=Date, y=GenusAbundance, color = Site, shape = Site)) + geom_point(size = 2) +  geom_errorbar(aes(ymin=GenusAbundance-se, ymax=GenusAbundance+se)) + facet_wrap(~Rank6, ncol = 3, scales="free_y") + scale_colour_hue(h=c(400, 120))
p + theme_bw() + theme(axis.text.x = element_text(size = 10, angle = 45, hjust=1),axis.text.y = element_text(size = 10))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()

# Make a big taxa table. 
write.table(bact.relabund@tax_table, "Figs/taxa_table.txt",sep="\t",row.names=FALSE)

# Find methanotrophic bacteria by genus. 
methanolist <- read.table(file ="taxa-of-interest/methanos.txt")
methanolist <- as.vector(methanolist$V1)
bact.relabund.methanos <- subset_taxa(bact.relabund, Rank6 %in% as.factor(methanolist))
bact.relabund.methanos
relabund.methanos <- psmelt(bact.relabund.methanos)
relabund.methanos.genus <- relabund.methanos%>%
  group_by(Sample, Rank6)%>%
  mutate(GenusAbundance = sum(Abundance))%>%
  distinct(Sample, GenusAbundance, TreatmentGroup, Site, Date, Rank2, Rank6, Rank6)
relabund.methanos.genus.est <- summarySE(relabund.methanos.genus, measurevar="GenusAbundance", groupvars=c("Site","Date", "Rank6"))
head(relabund.methanos.genus.est)
relabund.methanos.genus.est$Date <- as.character(relabund.methanos.genus.est$Date)
relabund.methanos.genus.est$Date <- as.numeric(relabund.methanos.genus.est$Date) 

# Plot summary of genus abundance of methanotrophic genera. 
pdf("Figs/methanos.pdf", width = 7, height = 6)
p <- ggplot(relabund.methanos.genus.est, aes(x=Date, y=GenusAbundance, color = Site, shape = Site)) + geom_point(size = 2) +  geom_errorbar(aes(ymin=GenusAbundance-se, ymax=GenusAbundance+se)) + facet_wrap(~Rank6, ncol = 2, scales="free_y") + scale_colour_hue(h=c(400, 120))
p + theme_bw() + theme(axis.text.x = element_text(size = 10, angle = 45, hjust=1),axis.text.y = element_text(size = 10))+ labs(x="Date",y="Mean Relative Abundance")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(strip.text = element_text(face = "italic")) 
dev.off()

# We want to calculate the total relative abundance of each phylum. 
# "Melt" the phyloseq data into a dataframe and then take the top X most abundant phyla. 
bact.melt <- psmelt(bact.relabund)
bact.melt.sorted <- bact.melt %>%
  group_by(Sample,Rank2) %>%
  summarize(PhyAbund = sum(Abundance))%>%
  group_by(Rank2)%>%
  summarize(MeanPhyAbund = mean(PhyAbund))%>%
  arrange(-MeanPhyAbund) 
# List of nPhyla top most abundant phyla. 
nPhyla =16
PhylumList <- bact.melt.sorted[1:nPhyla,1]
PhylumList <- PhylumList[is.na(PhylumList)==FALSE,]
PhylumList <- levels(droplevels(as.factor(PhylumList$Rank2)))
PhylumList
# Subset bact.melt for phyla. 
bact.subset <- subset_taxa(bact.relabund, Rank2 %in% PhylumList)
bact.subset.melt <- psmelt(bact.subset)
bact.subset.melt.sorted = bact.subset.melt %>%
  group_by(Sample,Site,Date,Rank2) %>%
  summarize(PhyAbund = sum(Abundance)) 
# Summarize phylum abundances. 
bact.subset.melt.sorted.est <- summarySE(bact.subset.melt.sorted, measurevar="PhyAbund", groupvars=c("Site","Date", "Rank2"))
head(bact.subset.melt.sorted.est)
bact.subset.melt.sorted.est$Date <- as.character(bact.subset.melt.sorted.est$Date)
bact.subset.melt.sorted.est$Date <- as.numeric(bact.subset.melt.sorted.est$Date) 
# Plot top phyla abundances. 
pdf("Figs/phyla_top.pdf", height = 6, width = 12)
p <- ggplot(bact.subset.melt.sorted.est, aes(x=Date, y=PhyAbund, color = Site, shape = Site)) + geom_point(size = 2) +  geom_errorbar(aes(ymin=PhyAbund-se, ymax=PhyAbund+se)) + facet_wrap(~Rank2, ncol = 4, scales="free_y") + scale_colour_hue(h=c(400,120))
p + theme_bw() + theme(axis.text.x = element_text(size = 10, angle = 45, hjust=1),axis.text.y = element_text(size = 10)) + labs(x="Date",y="Mean Relative Abundance") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(strip.text = element_text(face = "italic")) 
dev.off()




# Find classes_zulk of interest and subset bact.relabund. 
classzulk <- read.csv(file = "taxa-of-interest/classes-zulkifly.csv", header = T)
classzulk.list <- as.vector(classzulk$Rank3)
bact.relabund.subset.classzulk = subset_taxa(bact.relabund, Rank3 %in% classzulk.list)
bact.relabund.subset.classzulk.taxa <- subset_taxa(bact.relabund.subset.classzulk, Rank3 %in% as.factor(classzulk.list))
bact.relabund.subset.classzulk.taxa
relabund.classzulk <- psmelt(bact.relabund.subset.classzulk.taxa)
relabund.classzulk.class <- relabund.classzulk%>%
  group_by(Sample, Rank3)%>%
  mutate(Rank3Abundance = sum(Abundance))%>%
  distinct(Sample, Rank3Abundance, TreatmentGroup, Site, Date, Rank3, Rank6, Rank3)
length(relabund.classzulk.class$Rank3); head(relabund.classzulk.class)
relabund.classzulk.class <- subset(relabund.classzulk.class, Rank3!="Chloroplast")
length(relabund.classzulk.class$Rank3); head(relabund.classzulk.class) 

# Summary of class abundance of classes_zulk of interest. 
relabund.classzulk.class.est <- summarySE(relabund.classzulk.class, measurevar ="Rank3Abundance", groupvars=c("Site","Date","Rank3"))
length(relabund.classzulk.class.est$Rank3)
relabund.classzulk.class.est$Date <- as.character(relabund.classzulk.class.est$Date)
relabund.classzulk.class.est$Date <- as.numeric(relabund.classzulk.class.est$Date) 

# Plot MERGED summary of class abundance of classes_zulk of interest and Zulkifly et al. class abundances. 
relabund.classzulk.class.est.merge <- merge(relabund.classzulk.class.est, classzulk, by = "Rank3")
head(relabund.classzulk.class.est.merge)
write.table(relabund.classzulk.class.est.merge,"Figs/relabund.classzulk.class.est.merge.txt",sep="\t",row.names=FALSE)
pdf("Figs/relabund.classzulk.class.est.merge.pdf", height = 7, width = 8)
p <- ggplot(relabund.classzulk.class.est.merge, aes(x=Date, y=Rank3Abundance, color = Site, shape = Site)) + geom_point(size = 2) +  geom_errorbar(aes(ymin=Rank3Abundance-se, ymax=Rank3Abundance+se)) + facet_wrap(~Rank3, ncol = 3, scales="free_y") + scale_colour_hue(h=c(400, 120)) + geom_hline(aes(yintercept = ClassAbundance_zulkifly), linetype="dashed", color="red")+ scale_colour_hue(h=c(400, 120))
p + theme_bw() + theme(axis.text.x = element_text(size = 10, angle = 45, hjust=1),axis.text.y = element_text(size = 10))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(strip.text = element_text(face = "italic")) + ylab("Mean Relative Abundance")
dev.off()

####### SANDBOX ##########

# Summary of class abundance of classes_zulk of interest. 
relabund.classzulk.class.est.nosplit <- summarySE(relabund.classzulk.class.est, measurevar ="Rank3Abundance", groupvars=c("Rank3","Rank2"))
head(relabund.classzulk.class.est.nosplit)
relabund.classzulk.class.est.merge.nosplit <- merge(relabund.classzulk.class.est.nosplit, classzulk, by ="Rank3")
relabund.classzulk.class.est.merge.nosplit 

df <- relabund.classzulk.class.est.merge.nosplit
limits <- aes(ymin = Rank3Abundance-se, ymax= Rank3Abundance+se)
p <- ggplot(df, aes(Rank3))
p <- p + 
  ylab("Mean Relative Abundance")+
  facet_grid(~Rank3, scales = "free", space="free") +
  guides(color=guide_legend(title="Growth Season")) + 
  #theme_bw() + 
  theme(legend.position="top")+
  theme(strip.text.x = element_text(size = 8, angle = 90, face="italic")) +
  theme(axis.text.x = element_text(size = 10, angle = 55, hjust=1),axis.text.y = element_text(size = 10))+
  geom_point(aes(y = Rank3Abundance, color = "2014 (this study)  "), size = 3, alpha = 0.5) + 
  geom_errorbar(limits, color = "black", width=0.5)+
  geom_point(aes(y = ClassAbundance_zulkifly, color = "2011 (Zulkifly et al., 2012)"), size = 3, alpha = 0.5)+
  scale_colour_manual(values=cbPalette)
p 

df <- read.csv(file = "braus-zulkifly.csv")
limits <- aes(ymin = ClassAbundance-se, ymax= ClassAbundance+se)
p <- ggplot(df, aes(Class, color = Year))
p <- p + 
  geom_point(aes(y = ClassAbundance), size = 2) + 
  geom_errorbar(limits, color = "black", width=0.5) +
  ylab("Mean Relative Abundance") +
  facet_grid(~Rank3, scales = "free", space="free") +
  guides(color=guide_legend(title="Growth Season")) + 
  theme_bw() + 
  theme(legend.position="top") +
  theme(strip.text.x = element_text(size = 8, angle = 90, face="italic")) +
  theme(axis.text.x = element_text(size = 10, angle = 55, hjust=1),axis.text.y = element_text(size = 10)) +
  scale_colour_manual(values=c("#b5b5b5", "#212121"))
p+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


####### NOT REFACTORED ##########

# Plot summary of genus abundance of methanotrophic genera. 
p <- ggplot(relabund.methanos.genus.est, aes(x=Date, y=GenusAbundance, color = Site, shape = Site)) + geom_point(size = 2) +  geom_errorbar(aes(ymin=GenusAbundance-se, ymax=GenusAbundance+se)) + facet_wrap(~Rank6, ncol = 5, scales="free_y") + scale_colour_manual(values=cbPalette)
p + theme_bw() + theme(axis.text.x = element_text(size = 10, angle = 45, hjust=1),axis.text.y = element_text(size = 10))+ labs(x="Date",y="Mean Relative\nAbundance")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(strip.text = element_text(face = "italic")) 

# Find all genera
all.genera <- sort(get_taxa_unique(bact.relabund, "Rank6"), decreasing=FALSE)
bact.relabund.all.genera <- subset_taxa(bact.relabund, Rank6 %in% as.factor(all.genera))
bact.relabund.all.genera <- psmelt(bact.relabund.all.genera)
bact.relabund.all.genera.genus <- bact.relabund.all.genera%>%
  group_by(Sample, Rank6)%>%
  mutate(GenusAbundance = sum(Abundance))%>%
  distinct(Sample, GenusAbundance, TreatmentGroup, Site, Date, Rank6, Rank6) 

bact.relabund.all.genera.genus.est <- summarySE(bact.relabund.all.genera.genus, measurevar="GenusAbundance", groupvars=c("Site","Date", "Rank6"))
head(bact.relabund.all.genera.genus.est)
bact.relabund.all.genera.genus.est$Date <- as.character(bact.relabund.all.genera.genus.est$Date)
bact.relabund.all.genera.genus.est$Date <- as.numeric(bact.relabund.all.genera.genus.est$Date)
p <- ggplot(bact.relabund.all.genera.genus.est, aes(x=Date, y=GenusAbundance, color = Site, shape = Site)) + geom_point(size = 2) +  geom_errorbar(aes(ymin=GenusAbundance-se, ymax=GenusAbundance+se)) + facet_wrap(~Rank6, ncol = 3, scales="free_y") + scale_colour_hue(h=c(400, 120))
p + theme_bw() + theme(axis.text.x = element_text(size = 10, angle = 45, hjust=1),axis.text.y = element_text(size = 10))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

# Plot top phyla abundances. 
p <- ggplot(bact.subset.melt.sorted.est, aes(x=Date, y=PhyAbund, color = Site, shape = Site)) + geom_point(size = 2) +  geom_errorbar(aes(ymin=PhyAbund-se, ymax=PhyAbund+se)) + facet_wrap(~Rank2, ncol = 3, scales="free_y") + scale_colour_manual(values=cbPalette)
p + theme_bw() + theme(axis.text.x = element_text(size = 10, angle = 45, hjust=1),axis.text.y = element_text(size = 10)) + labs(x="Date",y="Mean Relative Abundance") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(strip.text = element_text(face = "italic")) 

# Find top 30 classes and subset bact.relabund. 
sort.classes <- sort(tapply(taxa_sums(bact.relabund), tax_table(bact.relabund)[, "Class"], sum), TRUE)
top.classes <- sort.classes[1:30]
top.classes.list <- names(top.classes)
bact.relabund.subset = subset_taxa(bact.relabund, Class %in% top.classes.list)
bact.relabund.subset.taxa <- subset_taxa(bact.relabund.subset, Class %in% as.factor(top.classes.list))
bact.relabund.subset.taxa
relabund.top.classes <- psmelt(bact.relabund.subset.taxa)
relabund.top.classes.class <- relabund.top.classes%>%
  group_by(Sample, Class)%>%
  mutate(ClassAbundance = sum(Abundance))%>%
  distinct(Sample, ClassAbundance, TreatmentGroup, Site, Date, Rank2, Rank6, Class)
relabund.top.classes.class.print <- summarySE(relabund.top.classes.class, measurevar="ClassAbundance", groupvars=c("Class"))
relabund.top.classes.class.print[,c("Class","N","ClassAbundance")] 

# Summary of class abundance of top 30 classes. 
relabund.top.classes.class.est <- summarySE(relabund.top.classes.class, measurevar="ClassAbundance", groupvars=c("Site","Date", "Class"))
head(relabund.top.classes.class.est)
relabund.top.classes.class.est$Date <- as.character(relabund.top.classes.class.est$Date)
relabund.top.classes.class.est$Date <- as.numeric(relabund.top.classes.class.est$Date) 

# Plot summary of class abundance of top 30 classes. 
p <- ggplot(relabund.top.classes.class.est, aes(x=Date, y=ClassAbundance, color = Site, shape = Site)) + geom_point(size = 2) +  geom_errorbar(aes(ymin=ClassAbundance-se, ymax=ClassAbundance+se)) + facet_wrap(~Class, ncol = 3, scales="free_y") + scale_colour_hue(h=c(400, 120))
p + theme_bw() + theme(axis.text.x = element_text(size = 10, angle = 45, hjust=1),axis.text.y = element_text(size = 10))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

# Find all classes
all.classes <- sort(get_taxa_unique(bact.relabund, "Class"), decreasing=FALSE)
bact.relabund.all.classes <- subset_taxa(bact.relabund, Class %in% as.factor(all.classes))
bact.relabund.all.classes <- psmelt(bact.relabund.all.classes)
bact.relabund.all.classes.class <- bact.relabund.all.classes%>%
  group_by(Sample, Class)%>%
  mutate(ClassAbundance = sum(Abundance))%>%
  distinct(Sample, ClassAbundance, TreatmentGroup, Site, Date, Rank6, Class) 

bact.relabund.all.classes.class.est <- summarySE(bact.relabund.all.classes.class, measurevar="ClassAbundance", groupvars=c("Site","Date", "Class"))
head(bact.relabund.all.classes.class.est)
bact.relabund.all.classes.class.est$Date <- as.character(bact.relabund.all.classes.class.est$Date)
bact.relabund.all.classes.class.est$Date <- as.numeric(bact.relabund.all.classes.class.est$Date) 

p <- ggplot(bact.relabund.all.classes.class.est, aes(x=Date, y=ClassAbundance, color = Site, shape = Site)) + geom_point(size = 2) +  geom_errorbar(aes(ymin=ClassAbundance-se, ymax=ClassAbundance+se)) + facet_wrap(~Class, ncol = 3, scales="free_y") + scale_colour_hue(h=c(400, 120))
p + theme_bw() + theme(axis.text.x = element_text(size = 10, angle = 45, hjust=1),axis.text.y = element_text(size = 10))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

# Data downloaded from 
data <- read.csv("north_temperate_lakes_lter__daily_water_temperature_-_lake_mendota.csv", header=T)
head(data)
nolegend <- theme(legend.position="none")
Jdays <- c(172,178,185,199,206,214) 

# Subset all wtemp at or below 1.5 meters (sample depth up to shore). 
clado.depth <- subset(data, depth<=1.5, select = year4:wtemp) 

# All years (2006-2016), full year. 
p <- qplot(clado.depth$daynum, y = clado.depth$wtemp)
p <- p + geom_point(aes(colour = clado.depth$year4), size=1, alpha = 0.8) + 
  scale_y_continuous(limits = c(4,30))+
  xlab("Date") + ylab("Temperature (°C)")+ 
  geom_vline(xintercept=Jdays, colour="darkgrey", linetype = "longdash") +
  guides(color=guide_legend(title="Year")) + 
  scale_colour_gradient(low = "darkred")
p+ theme_bw() + theme(axis.text.x = element_text(size = 10, angle = 45, hjust=1),axis.text.y = element_text(size = 10))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

# All years (2006-2016), sample dates ± 25 days
clado.depth.zoom <- subset(clado.depth, daynum>=150 & daynum<=250, select = year4:wtemp)
p <- qplot(clado.depth.zoom$daynum, y = clado.depth.zoom$wtemp) + 
  scale_y_continuous(limits = c(10,30))+
  geom_point(aes(colour = clado.depth.zoom$year4), size=1, alpha = 0.8) + 
  labs(title = "Water Temperature (°C) <1.5m Depth, Lake Mendota", x="Date", y="Temperature °C") + 
  geom_vline(xintercept=Jdays, colour="darkgreen", linetype = "longdash") + 
  guides(color=guide_legend(title="Year")) + scale_colour_gradient(low = "darkred")
p+ theme_bw() + theme(axis.text.x = element_text(size = 10, angle = 45, hjust=1),axis.text.y = element_text(size = 10))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

# Smarter way to make grid of plots. 
p <- ggplot(clado.depth.zoom, aes(daynum, wtemp))
p <- p + geom_point(size = 0.5) + 
  ylim(15,30) + 
  facet_wrap(~year4, ncol = 3) + 
  geom_vline(xintercept=Jdays, colour="grey30", linetype = "longdash") +
  guides(color=guide_legend(title="Year")) +
  labs(title = "Water Temperature (°C) <1.5m Depth, Lake Mendota", x="Date", y="Temperature °C")
p + theme_bw() + theme(axis.text.x = element_text(size = 10, angle = 45, hjust=1),axis.text.y = element_text(size = 10))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

clado.depth.not14 <- subset(clado.depth, year4!=2014)
clado.depth.est.not14 <- summarySE(clado.depth.not14, measurevar= "wtemp", groupvars=c("daynum"))
head(clado.depth.est.not14)
p <- ggplot(clado.depth.est.not14, aes(x=daynum, y=wtemp)) + 
  geom_point(size = 0.5) +  
  geom_errorbar(aes(ymin=wtemp-se, ymax=wtemp+se))+
  geom_vline(xintercept=Jdays, colour="darkgrey", linetype = "longdash") + 
  labs(x="Date", y="Temperature (°C) depths ≤1.5m")
p + theme_bw() + 
  scale_y_continuous(limits = c(10,26))+
  scale_x_continuous(limits = c(120,320))+
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust=1),axis.text.y = element_text(size = 10))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

clado.depth.deep <- subset(clado.depth, depth>0)
clado.depth.not14 <- subset(clado.depth.deep, year4!=2014)
clado.depth.est.not14 <- summarySE(clado.depth.not14, measurevar= "wtemp", groupvars=c("daynum", "year4"))
clado.depth.est.not14.est <- summarySE(clado.depth.est.not14, measurevar= "wtemp", groupvars=c("daynum"))
head(clado.depth.est.not14.est)
p <- ggplot(clado.depth.est.not14.est, aes(x=daynum, y=wtemp)) + 
  geom_point(size = 0.5) +  
  geom_errorbar(aes(ymin=wtemp-se, ymax=wtemp+se))+
  geom_vline(xintercept=Jdays, colour="darkgreen", linetype = "longdash") + 
  labs(x="Date", y="Temperature (°C), 0.5-1.5m")
p + theme_bw() + 
  scale_y_continuous(limits = c(10,26))+
  scale_x_continuous(limits = c(120,320))+
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust=1),axis.text.y = element_text(size = 10))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

clado.depth.deep <- subset(clado.depth, depth>0)
clado.depth.not14 <- subset(clado.depth.deep, year4!=2014)
clado.depth.est.not14 <- summarySE(clado.depth.not14, measurevar= "wtemp", groupvars=c("daynum", "year4"))
clado.depth.est.not14.est <- summarySE(clado.depth.est.not14, measurevar= "wtemp", groupvars=c("daynum"))
head(clado.depth.est.not14.est)
p <- ggplot(clado.depth.est.not14.est, aes(x=daynum, y=wtemp)) + 
  geom_point(size = 0.5) +  
  geom_errorbar(aes(ymin=wtemp-se, ymax=wtemp+se))+
  geom_vline(xintercept=Jdays, colour="grey10", linetype = "longdash") + 
  labs(x="Date", y="Temperature (°C), 0.5-1.5m")
p + theme_bw() + 
  scale_y_continuous(limits = c(10,26))+
  scale_x_continuous(limits = c(120,320))+
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust=1),axis.text.y = element_text(size = 10))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

# Find genera of interest and subset bact.relabund. 
genofint <- read.table(file = "taxa-of-interest/genera-of-interest.txt")
genofint.list <- as.vector(genofint$V1)
bact.relabund.subset.genofint = subset_taxa(bact.relabund, Rank6 %in% genofint.list)
bact.relabund.subset.genofint.taxa <- subset_taxa(bact.relabund.subset.genofint, Rank6 %in% as.factor(genofint.list))
bact.relabund.subset.genofint.taxa 
relabund.genofint <- psmelt(bact.relabund.subset.genofint.taxa)
relabund.genofint.genus <- relabund.genofint%>%
  group_by(Sample, Rank6)%>% 
  mutate(GenusAbundance = sum(Abundance))%>%
  distinct(Sample, GenusAbundance, TreatmentGroup, Site, Date, Rank6, Rank6)
head(relabund.genofint.genus) 

# Summary of genus abundance of genera of interest. 
relabund.genofint.genus.est <- summarySE(relabund.genofint.genus, measurevar="GenusAbundance", groupvars=c("Site","Date","Rank6")) 
head(relabund.genofint.genus.est)
relabund.genofint.genus.est$Date <- as.character(relabund.genofint.genus.est$Date)
relabund.genofint.genus.est$Date <- as.numeric(relabund.genofint.genus.est$Date) 

# Plot summary of genus abundance of genera of interest. 
p <- ggplot(relabund.genofint.genus.est, aes(x=Date, y=GenusAbundance, color = Site, shape = Site)) + geom_point(size = 2) +  geom_errorbar(aes(ymin=GenusAbundance-se, ymax=GenusAbundance+se)) + facet_wrap(~Rank6, ncol = 3, scales="free_y") + scale_colour_hue(h=c(400, 120))
p + theme_bw() + theme(axis.text.x = element_text(size = 10, angle = 45, hjust=1),axis.text.y = element_text(size = 10)) + ylab("Mean Relative Abundance")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(strip.text = element_text(face = "italic"))  

# Plot summary of genus abundance of genera of interest. 
p <- ggplot(relabund.genofint.genus.est, aes(x=Date, y=GenusAbundance, color = Site, shape = Site)) + geom_point(size = 2) +  geom_errorbar(aes(ymin=GenusAbundance-se, ymax=GenusAbundance+se)) + facet_wrap(~Rank6, ncol = 3, scales="free_y") + scale_colour_manual(values=cbPalette) 
p + theme_bw() + theme(axis.text.x = element_text(size = 10, angle = 45, hjust=1),axis.text.y = element_text(size = 10)) + ylab("Mean Relative Abundance")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(strip.text = element_text(face = "italic")) 
