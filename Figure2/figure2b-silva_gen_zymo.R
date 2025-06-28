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
ps_prok_nosing <- prune_taxa(taxa_sums(ps_prok) > 1, ps_prok)
df_prok_nosing=psmelt(ps_prok_nosing)


# Remove taxa per sample with <10 reads. Only if needed.
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
setpalette <- colorRampPalette(brewer.pal(8,"Set2"))

# 2. Decontaminate dataset #

library(decontam)

# Let’s take a quick first look at the library sizes (i.e. the number of reads) in each sample, as a function of whether that sample was a true positive sample or a negative control:

df <- as.data.frame(sample_data(ps_filtered)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps_filtered)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_control)) + geom_point()

# Either frequency or prevalence method.
contamdf <- isContaminant(ps_filtered, method="combined", conc="DNA_quant_reading", neg="is.neg", threshold=0.5)

# Check decontam scores distribution
hist(contamdf$p)
pdf(file="/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/results/figures/suppfigure1a-decontam_hist_zymo_gen.pdf", width = 6, height = 6)
hist=hist(contamdf$p)
dev.off()

table(contamdf$contaminant)

# Now that we have identified likely contaminants, let’s remove them from the phyloseq object:

ps_filtered
ps_noncontam <- prune_taxa(!contamdf$contaminant, ps_filtered)
ps_noncontam
df_noncontam <- psmelt(ps_noncontam)

# Remove singletons.
ps_noncontam = prune_taxa(taxa_sums(ps_noncontam) > 1, ps_noncontam)

# 3 - Check positive controls #
ps_zymo <- subset_samples(ps_noncontam, cultureID == "Zymo_standard")
ps_zymo = prune_taxa(taxa_sums(ps_zymo) > 10, ps_zymo)
df_zymo=psmelt(ps_zymo)

# Plot relative abundance with raw abundances.
zymo_gen=ggplot(data=df_zymo, aes(cultureID), width=0.5) + theme_bw() +
  geom_bar(aes(color=Genus, fill=Genus,y=Abundance), stat="identity",position="fill", width=0.5) +
  #facet_grid(. ~ factor(seq_run), scales="free_x", space="free_x") +
  #theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=c("Bacillus"="#a5cee2", 
                             "Enterococcus"="#3f8faa",
                             "Escherichia-Shigella"="#79c360",
                             "Limosilactobacillus"="#f88a89",
                             "Listeria"="#f06c45",
                             "Pseudomonas"="#fc870f",
                             "Salmonella"="#b295c7",
                             "Staphylococcus"="#c7b699", 
                             "Unknown" = "#7f7f7f")) + 
  scale_color_manual(values = c("Bacillus"="#a5cee2", 
                                "Enterococcus"="#3f8faa",
                                "Escherichia-Shigella"="#79c360",
                                "Limosilactobacillus"="#f88a89",
                                "Listeria"="#f06c45",
                                "Pseudomonas"="#fc870f",
                                "Salmonella"="#b295c7",
                                "Staphylococcus"="#c7b699", 
                                "Unknown" = "#7f7f7f")) + 
  theme(text=element_text(size=20))
zymo_gen

ggsave("/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/results/figures/zymo_silva_gen_all.pdf",zymo_gen, width = 6, height = 10)

# Identify percentage of Unknown reads
ps_zymo_prop  = transform_sample_counts(ps_zymo, function(x) x / sum(x))
df_zymo_prop = psmelt(ps_zymo_prop)

# 4. Check negative controls #
# Keep only the negative controls.
ps_neg <- subset_samples(ps_noncontam, Sample_or_control == "Control")
df_neg=psmelt(ps_neg)

# Remove singletons.
ps_neg_pruned <- prune_taxa(taxa_sums(ps_neg) > 1, ps_neg)
df_neg_pruned=psmelt(ps_neg_pruned)

neg_sp=ggplot(data=df_neg_pruned, aes(sampleID)) +
  geom_bar(aes(color=Genus, fill=Genus,y=Abundance), stat="identity",position="fill") + 
  #facet_grid(. ~ factor(seq_run), scales="free_x", space="free_x") +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=mypalette(41)) + 
  scale_color_manual(values=mypalette(41)) +
  scale_x_discrete(drop=F)
neg_sp