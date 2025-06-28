# 8. Visualize abundance plots #

# Load environmental variables.
load('figure3-out.RData')

### Plot abundance without controls, with unknowns ###

# Remove negative and positive controls.
ps_prok_noc <- subset_samples(ps_noncontam, cultureID != "Bio_neg_control")
ps_prok_noc <- subset_samples(ps_prok_noc, cultureID != "PCR_neg_control")
ps_prok_noc <- subset_samples(ps_prok_noc, cultureID != "Zymo_standard")
ps_prok_noc
df_prok_noc <- psmelt(ps_prok_noc)

# Relative abundance

# Plot top 50 abundance
top50 <- df_prok_noc %>%
  dplyr::select(OTU, sampleID, cultureID, Genus, Abundance) %>%
  group_by(sampleID) %>%
  mutate(totalSum = sum(Abundance)) %>%
  ungroup() %>%
  group_by(sampleID, cultureID, Genus, totalSum) %>%
  reframe(
    Abundance = sum(Abundance),
    Genus = ifelse(Abundance < 220, ".Others", Genus)) %>% 
  group_by(sampleID, cultureID, Genus, totalSum) %>%
  reframe(
    Abundance = sum(Abundance),
    RelAb = Abundance/totalSum) %>%
  unique()
# Check number of unique genera. The number returned should = 51 because top 20 families + 1 "Others"
length(unique(top50$Genus))

length(unique(df_prok_noc$Genus))

# Reorder x-axis labels
level_order <- c("S.mic_KB8",
                 "S.mic_CCMP2464",
                 "B.aen_MAC-04-234",
                 "D.tre_CCMP3408",
                 "C203",
                 "E.vor_RT-383")
level_label <- c("S.mic_KB8",
                 "S.mic_CCMP2464",
                 "B.aen_MAC-04-234",
                 "D.tre_CCMP3408",
                 "C203",
                 "E.vor_RT-383")

top50_plot <- ggplot(data=top50, aes(cultureID, x = factor(cultureID,level=level_order))) + 
  geom_bar(aes(color=Genus, fill=Genus,y=Abundance), stat="identity",position="fill") + 
  theme(axis.text.x=element_text(angle=90)) + xlab("Culture IDs") + ylab("Relative Abundance") +
  scale_fill_manual(values=mypalette(51)) + 
  scale_color_manual(values=mypalette(51)) +
  guides(fill=guide_legend(ncol=2,bycol=T)) + 
  theme(text=element_text(size=20))
top50_plot

pdf(file = "/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/results/figures/rel_abund_top50.pdf",   # The directory you want to save the file in
    width = 16, # The width of the plot in inches
    height = 10) # The height of the plot in inches

top50_plot

dev.off()

# Absolute abundance
abs_plot = ggplot(top50, aes(cultureID, x = factor(cultureID,level=level_order))) + 
  geom_bar(aes(color=Genus, fill=Genus,y=Abundance), stat="identity",position="stack") + 
  theme(axis.text.x=element_text(angle=90)) + xlab("Culture IDs") + ylab("Absolute Abundance") +
  scale_fill_manual(values=mypalette(51)) + 
  scale_color_manual(values=mypalette(51)) +
  guides(fill=guide_legend(ncol=2,bycol=T)) + 
  theme(text=element_text(size=20))
abs_plot

pdf(file = "/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/results/figures/abs_abund_top50.pdf",   # The directory you want to save the file in
    width = 16, # The width of the plot in inches
    height = 10) # The height of the plot in inches

abs_plot

dev.off()

# Identify core microbiome members.

### UpSet plot ###
library("MicrobiotaProcess")

metadata_noncontam=sample_data(ps_noncontam)
unique_gen=get_upset(ps_prok_noc,factorNames="Symbiodiniaceae_genus")

library(UpSetR)
library(ComplexHeatmap)
#library(ComplexUpset)

upset=upset(unique_gen, nsets=30, order.by="degree", decreasing=T, text.scale=1.5, set_size.show=T)
upset

pdf(file = "/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/results/figures/upset_indiv.pdf",   # The directory you want to save the file in
    width = 14, # The width of the plot in inches
    height = 10) # The height of the plot in inches

upset

dev.off()

# Determine identities of genera shared between all cultures
library(phylosmith)
species_shared_all=common_taxa(ps_noncontam, treatment="cultureID")
species_shared_all
sink(file="/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/results/figures/species_shared_all.txt")
species_shared_all
sink(file = NULL)

species_shared_noc=common_taxa(ps_prok_noc, treatment="cultureID")
species_shared_noc
sink(file="/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/results/figures/species_shared_noc.txt")
species_shared_noc
sink(file = NULL)

species_shared_allcultures=common_taxa(ps_noncontam, treatment="sampleID")
species_shared_allcultures
sink(file="/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/results/figures/species_shared_allcultures.txt")
species_shared_allcultures
sink(file = NULL)


# Abundance of core microbiome members.
df_core <- df_prok_noc[df_prok_noc$Genus == c("Muricauda", "Labrenzia", "Marinobacter"),]
core_plot = ggplot(df_core, aes(sampleID)) + 
  geom_bar(aes(color=Genus, fill=Genus,y=Abundance), stat="identity",position="fill") + 
  theme(axis.text.x=element_text(angle=90)) + xlab("Culture IDs") + ylab("Relative Abundance") +
  scale_fill_manual(values=mypalette(6)) + 
  scale_color_manual(values=mypalette(6)) +
  guides(fill=guide_legend(ncol=1,bycol=T)) + 
  theme(text=element_text(size=20))
core_plot

pdf(file = "/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/results/figures/rel_abund_core.pdf",   # The directory you want to save the file in
    width = 20, # The width of the plot in inches
    height = 14) # The height of the plot in inches

core_plot

dev.off()

core_abs = ggplot(df_core, aes(sampleID)) + 
  geom_bar(aes(color=Genus, fill=Genus,y=Abundance), stat="identity",position="stack") + 
  theme(axis.text.x=element_text(angle=90)) + xlab("Culture IDs") + ylab("Relative Abundance") +
  scale_fill_manual(values=mypalette(6)) + 
  scale_color_manual(values=mypalette(6)) +
  guides(fill=guide_legend(ncol=1,bycol=T)) + 
  theme(text=element_text(size=20))
core_abs

pdf(file = "/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/results/figures/abs_abund_core.pdf",   # The directory you want to save the file in
    width = 20, # The width of the plot in inches
    height = 14) # The height of the plot in inches

core_abs

dev.off()

# Plot NMDS for CLR transformed abundances.
ps_noc_clr <- subset_samples(ps_clr, cultureID != "Bio_neg_control")
ps_noc_clr <- subset_samples(ps_noc_clr, cultureID != "PCR_neg_control")
ps_noc_clr <- subset_samples(ps_noc_clr, cultureID != "Zymo_standard")
ps_noc_clr
df_noc_clr <- psmelt(ps_noc_clr)

ait_clr <- phyloseq::distance(ps_noc_clr, "euclidean")
ait.nmds_clr = ordinate(ps_noc_clr, method="NMDS", distance=ait_clr)
nmds_clr <- plot_ordination(ps_clr, ait.nmds_clr, type = "samples", color="cultureID") + geom_point(size = 3) + scale_color_manual(values=setpalette(6), name="cultureID") 
nmds_clr
nmds_clr = nmds_clr + stat_ellipse(aes(group=cultureID,color=cultureID)) + 
  theme(text=element_text(size=20))
nmds_clr

ggsave("/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/results/figures/nmds_clr_cultureID.pdf",nmds_clr)

# Stat tests all.
library(veganEx)
library(vegan)

set.seed(100)
meta_clr=sample_data(ps_noc_clr)
(anosim = vegan::anosim(ait_clr, meta_clr$cultureID)) # If you just want to know whether within-group similarities are greater than between-group similarities (without pinpointing how or where): Use ANOSIM
(permanova = vegan::adonis2(ait_clr ~ meta_clr$cultureID)) # If you're interested in whether group centroids differ in multivariate space (i.e., actual multivariate "location" differences): Use PERMANOVA
# Always check for dispersion differences (via PERMDISP) before trusting PERMANOVA results.
# If possible, report both PERMANOVA and ANOSIM (and their respective statistics) to support robust interpretation.

# Pairwise
anosim_cultures=anosim.pairwise(ait_clr, meta_clr$cultureID)
anosim_cultures
library(RVAideMemoire)
pairwise_permanova_cultures=pairwise.perm.manova(ait_clr,meta_clr$cultureID,nperm=999)
pairwise_permanova_cultures

# Average number of taxa per sample.
## Make separate lists of Genera present in each culture.
df_clr_noc=psmelt(ps_prok_noc)

df_3 = subset(df_clr_noc, cultureID=="S.mic_KB8")
df_3 = subset(df_3,Abundance>0)
n_distinct(df_3$Genus)
#gen_3 = unique(df_3$Genus)

df_4 = subset(df_clr_noc, cultureID=="S.mic_CCMP2464")
df_4 = subset(df_4,Abundance>0)
n_distinct(df_4$Genus)
#gen_4 = unique(df_4$Genus)

df_6 = subset(df_clr_noc, cultureID=="B.aen_MAC-04-234")
df_6 = subset(df_6,Abundance>0)
n_distinct(df_6$Genus)
#gen_6 = unique(df_6$Genus)

df_1 = subset(df_clr_noc, cultureID=="D.tre_CCMP3408")
df_1 = subset(df_1,Abundance>0)
n_distinct(df_1$Genus)
#gen_1 = unique(df_1$Genus)

df_12 = subset(df_clr_noc, cultureID=="C203")
df_12 = subset(df_12,Abundance>0)
n_distinct(df_12$Genus)
#gen_12 = unique(df_12$Genus)

df_2 = subset(df_clr_noc, cultureID=="E.vor_RT-383")
df_2 = subset(df_2,Abundance>0)
n_distinct(df_2$Genus)
#gen_2 = unique(df_2$Genus)