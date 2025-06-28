# 1. Import data and create phyloseq object #

# Import ASV abundance table.
library(tidyverse)
abund = as.data.frame(read_tsv("/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/results/pilot_silva_q20_95percentsimilarity_gen/abundance_table_genus.tsv"))
row.names(abund)=abund$tax
abund=as.matrix(abund[,-1])

# Import tax table.
taxa <- as.data.frame(read_tsv("/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/results/pilot_silva_q20_95percentsimilarity_gen/tax_table_genus.tsv"))
row.names(taxa)=taxa$tax
taxa=as.matrix(taxa[,-1])

# Import the metadata.
metadata <- as.data.frame(read_csv("/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/results/pilot_silva_q20_95percentsimilarity_gen/metadata.csv"))
row.names(metadata)=metadata$barcode
metadata=data.frame(metadata[,-1])
metadata$bio_rep <- as.character(metadata$bio_rep)
metadata$seq_run <- as.character(metadata$seq_run)
#metadata$quant_reading=metadata$quant_reading + 0.04

# Create phyloseq object.
library(phyloseq)
ps <- phyloseq(otu_table(abund, taxa_are_rows=T), 
               phyloseq::tax_table(taxa),
               sample_data(metadata))
# Check phyloseq object.
ps
# Create a data frame from the phyloseq object.
df= psmelt(ps)

# It looks like EPI2ME already filters out the mitochondria, chloroplast, and zero reads. 
# However, we still need to remove eukaryotic OTUs.
ps_prok <- subset_taxa(ps, 
                       Superkingdom  != "Eukaryota")
df_prok=psmelt(ps_prok)

ps_prok_nozero = prune_taxa(taxa_sums(ps_prok) > 0, ps_prok)
df_prok_nozero = psmelt(ps_prok_nozero)

# Remove singletons. Only if needed.
ps_prok_nosing <- prune_taxa(taxa_sums(ps_prok) > 10, ps_prok)
df_prok_nosing=psmelt(ps_prok_nosing)


# Remove <1% abundance taxa per sample. 
FSr  = transform_sample_counts(ps_prok_nosing, function(x) x / sum(x))
physeqrF = filter_taxa(FSr, function(x) sum(x) < .01, TRUE)
rmtaxa = taxa_names(physeqrF)
alltaxa = taxa_names(ps_prok_nosing)

myTaxa = alltaxa[!alltaxa %in% rmtaxa]

ps_filtered <- prune_taxa(myTaxa,ps_prok_nosing)

df_filtered=psmelt(ps_filtered)


# Define color palettes.
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(12,"Paired"))
setpalette <- colorRampPalette(brewer.pal(8,"Dark2"))

# 2. Decontaminate dataset #

library(decontam)

# Let’s take a quick first look at the library sizes (i.e. the number of reads) in each sample, as a function of whether that sample was a true positive sample or a negative control:

df <- as.data.frame(sample_data(ps_filtered)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps_filtered)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_control)) + geom_point()

# Combined method.
contamdf <- isContaminant(ps_filtered, method="combined", conc="DNA_quant_reading", neg="is.neg", threshold=0.5)

table(contamdf$contaminant)

# Now that we have identified likely contaminants, let’s remove them from the phyloseq object:

ps_filtered
ps_noncontam <- prune_taxa(!contamdf$contaminant, ps_filtered)
ps_noncontam
df_noncontam <- psmelt(ps_noncontam)

# Remove singletons.
ps_noncontam = prune_taxa(taxa_sums(ps_noncontam) > 1, ps_noncontam)

# 3. Centered log-ratio transformation #

# Remove positive and negative controls so as not to skew data.
#ps_prok_noc <- subset_samples(ps_noncontam, cultureID != "Zymo_standard")
#ps_prok_noc <- subset_samples(ps_prok_noc, cultureID != "Bio_neg_control")
#ps_prok_noc <- subset_samples(ps_prok_noc, cultureID != "PCR_neg_control")
#ps_prok_noc
#df_prok_noc <- psmelt(ps_prok_noc)

# Add pseudocount to remove 0 values. 
otu_prok = data.frame(ps_noncontam@otu_table) # Extract otu table from ps object.
otu_pseudocount = otu_prok + 1 # Add pseudocount 1 to all counts.

library(compositions)

otu_clr = as.data.frame(clr(otu_pseudocount)) # Perform centered log-ratio transformation.

# Make new ps object with transformed otu table.
ps_clr <- phyloseq(otu_table(otu_clr, taxa_are_rows=T), 
                   phyloseq::tax_table(taxa),
                   sample_data(metadata))

# 4. Statistical tests, PCoA, and NMDS #
library('ggplot2')

# Plot PCoA for CLR transformed abundances.
ait_clr <- phyloseq::distance(ps_clr, "euclidean")
ait.pcoa_clr = ordinate(ps_clr, method="PCoA", distance=ait_clr)
plot_scree(ait.pcoa_clr, "Screen plot")
theme_set(theme_gray())
pcoa_clr<-0
pcoa_clr <- plot_ordination(ps_clr, ait.pcoa_clr, type = "samples", color="seq_run", shape="cultureID") + theme_bw() + geom_point(size = 4, stroke=1.5) + scale_color_manual(values=mypalette(6), name="Sequencing run") + scale_shape_manual(values=1:9)
pcoa_clr = pcoa_clr + stat_ellipse(aes(group=seq_run,color=seq_run), linewidth=1.5) + 
  theme(text=element_text(size=20))
pcoa_clr

ggsave("/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/results/figures/pcoa_clr_samples.pdf",pcoa_clr)

# Plot NMDS for CLR transformed abundances.
ait_clr <- phyloseq::distance(ps_clr, "euclidean")
ait.nmds_clr = ordinate(ps_clr, method="NMDS", distance=ait_clr)
nmds_clr <- plot_ordination(ps_clr, ait.nmds_clr, type = "samples", color="seq_run", shape="cultureID") + theme_bw() + geom_point(size = 4, stroke=1.5) + scale_color_manual(values=mypalette(6), name="Sequencing run") + scale_shape_manual(values=1:9)
nmds_clr
nmds_clr = nmds_clr + stat_ellipse(aes(group=seq_run,color=seq_run), linewidth=1.5) + 
  theme(text=element_text(size=20))
nmds_clr

ggsave("/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/results/figures/nmds_clr_samples.pdf",nmds_clr)

# Stat tests all.
library(veganEx)
library(vegan)

set.seed(100)
meta_clr=sample_data(ps_clr)
(anosim = vegan::anosim(ait_clr, meta_clr$seq_run)) # If you just want to know whether within-group similarities are greater than between-group similarities (without pinpointing how or where): Use ANOSIM
(permanova = vegan::adonis2(ait_clr ~ meta_clr$seq_run)) # If you're interested in whether group centroids differ in multivariate space (i.e., actual multivariate "location" differences): Use PERMANOVA
# Always check for dispersion differences (via PERMDISP) before trusting PERMANOVA results.
# If possible, report both PERMANOVA and ANOSIM (and their respective statistics) to support robust interpretation.

# Save environmental variables.
save.image(file='figure3-out.RData')
