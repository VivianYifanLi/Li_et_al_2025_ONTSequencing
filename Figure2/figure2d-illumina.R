# 1. Import data and create phyloseq object #

# Import ASV abundance table.
library(tidyverse)
abund = as.data.frame(read_csv("/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/minion_manuscript_illumina/output/seqtab-nochim.csv"))
abund = t(abund)
colnames(abund)=abund[1,]
abund=as.matrix(abund[-1,])
class(abund) <- "numeric"

# Import tax table.
taxa <- as.data.frame(read_csv("/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/minion_manuscript_illumina/output/tax-table.csv"))
row.names(taxa)=taxa[,1]
taxa=as.matrix(taxa[,-1])

# Import the metadata.
metadata <- as.data.frame(read_csv("/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/minion_manuscript_illumina/metadata_pooled.csv"))
row.names(metadata)=metadata$sampleID
#metadata=data.frame(metadata[,-1])
#metadata$bio_rep <- as.character(metadata$bio_rep)
#metadata$seq_run <- as.character(metadata$seq_run)
#metadata$quant_reading=metadata$quant_reading + 0.04

# Create phyloseq object.
library(phyloseq)
ps <- phyloseq(phyloseq::otu_table(abund, taxa_are_rows=T), 
               phyloseq::tax_table(taxa),
               sample_data(metadata))
# Check phyloseq object.
ps
# Create a data frame from the phyloseq object.
df= psmelt(ps)

# Filter out the mitochondria, chloroplast, and zero reads. 
# Remove eukaryotic OTUs.
ps <- subset_taxa(ps, 
                  Family  != "Mitochondria")
ps_prok <- subset_taxa(ps, 
                       Kingdom  != "Eukaryota")
df_prok=psmelt(ps_prok)

ps_prok_nozero = prune_taxa(taxa_sums(ps_prok) > 0, ps_prok)
df_prok_nozero = psmelt(ps_prok_nozero)

# Remove singletons. Only if needed.
ps_prok_nosing <- prune_taxa(taxa_sums(ps_prok_nozero) > 1, ps_prok_nozero)
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

# 2. Decontamination #
library(decontam)

contamdf <- isContaminant(ps_filtered, method="combined", conc="quant_reading", neg="is.neg", threshold=0.5)

hist(contamdf$p)
pdf(file="/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/results/figures/suppfigure1c-decontam_hist_illumina_gen.pdf", width = 6, height = 6)
hist=hist(contamdf$p)
dev.off()

table(contamdf$contaminant)

# Examine contaminant sequences
ps_contaminants <- prune_taxa(contamdf$contaminant=="TRUE", ps_filtered)
df_contaminants <- psmelt(ps_contaminants)
distinct = distinct(df_contaminants, Genus, .keep_all=TRUE)

ps_filtered
ps_noncontam <- prune_taxa(!contamdf$contaminant, ps_filtered)
ps_noncontam
df_noncontam <- psmelt(ps_noncontam)

# Remove singletons.
ps_noncontam = prune_taxa(taxa_sums(ps_noncontam) > 1, ps_noncontam)

# 3. Check positive controls #
# Keep only the positive controls.
ps_zymo <- subset_samples(ps_noncontam, original_cultureID == "Zymo_standard")
ps_zymo = prune_taxa(taxa_sums(ps_zymo) > 10, ps_zymo)
df_zymo=psmelt(ps_zymo)

# Plot relative abundance with raw abundances.
zymo_gen=ggplot(data=df_zymo, aes(original_cultureID), width=0.5) + theme_bw() +
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
  scale_color_manual(values=c("Bacillus"="#a5cee2", 
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

ggsave("/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/results/figures/zymo_gen_all_illumina.pdf",zymo_gen, width = 6, height = 10)

# Identify percentage of Unknown reads in positive controls
ps_zymo_prop  = transform_sample_counts(ps_zymo, function(x) x / sum(x))
df_zymo_prop = psmelt(ps_zymo_prop)


# 4. Check negative controls #
# Keep only the negative controls.
ps_neg <- subset_samples(ps_noncontam, sample_or_control == "Control")
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

ggsave("/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/allseqs_thesisChp1/results/figures/neg_species.pdf", plot=neg_sp)

# Check number of genera without neg controls.
ps_prok_noc <- subset_samples(ps_noncontam, original_cultureID != "Bio_neg_control")
ps_prok_noc <- subset_samples(ps_prok_noc, original_cultureID != "PCR_neg_control")
ps_prok_noc
df_prok_noc <- psmelt(ps_prok_noc)

pos_genera = distinct(df_noncontam, Genus, .keep_all=TRUE)

# Identify percentage of Unknown reads in all samples excluding negatives.
ps_prop = transform_sample_counts(ps_prok_noc, function(x) x / sum(x))
df_prop = psmelt(ps_prop)

df_prok_noc[is.na(df_prok_noc)] = "NA"
df_prok_NA = subset(df_prok_noc,Genus=="NA")
sum(df_prok_NA$Abundance)
sum(df_prok_noc$Abundance)
Unknown_prop = sum(df_prok_NA$Abundance)/sum(df_prok_noc$Abundance)

# Save environmental variables.
save.image(file='/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/code copy/figure5-out.RData')
