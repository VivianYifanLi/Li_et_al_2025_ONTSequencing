########## Plot summary of reads distribution ##########

library(tidyr)
library(tidyverse)
library(ggplot2)
library(decontam)

##### ONT #####

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
df_noncontam <- psmelt(ps_noncontam)

reads=as.data.frame(sample_sums(ps_noncontam))
metadata_reads_ONT = merge(metadata, reads, by = 'row.names', all = TRUE)
metadata_reads_ONT = na.omit(metadata_reads_ONT)

##### Illumina #####

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


# Remove <1% abundance taxa per sample. Only if needed.
# Convert counts to proportions
ps_prop <- transform_sample_counts(ps_prok_nosing, function(x) x*100 / sum(x))

# Keep only taxa with at least 0.25% abundance in at least one sample
keep <- filter_taxa(ps_prop, function(x){max(x) > 0.01})

# Filter out all other taxa
ps_filtered <- prune_taxa(taxa = keep, x = ps_prok_nosing)

df_filtered=psmelt(ps_filtered)

## Decontam ##

library(decontam)

# Let’s take a quick first look at the library sizes (i.e. the number of reads) in each sample, as a function of whether that sample was a true positive sample or a negative control:

df <- as.data.frame(sample_data(ps_filtered)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps_filtered)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=sample_or_control)) + geom_point()

# Either freq or prev.
contamdf <- isContaminant(ps_filtered, method="either", conc="quant_reading", neg="is.neg", threshold=c(0.1,0.5))

ps_noncontam <- prune_taxa(!contamdf$contaminant, ps_filtered)
ps_noncontam

# Remove singletons.
ps_noncontam = prune_taxa(taxa_sums(ps_noncontam) > 1, ps_noncontam)
df_noncontam <- psmelt(ps_noncontam)

reads=as.data.frame(sample_sums(ps_noncontam))
metadata_reads_Illumina = merge(metadata, reads, by = 'row.names', all = TRUE)
metadata_reads_Illumina = na.omit(metadata_reads_Illumina)

#reads_tidy <- gather(reads_all_runs,
#                      key = "attributes",
#                      value = "values",
#                      -c(Platform, Run_ID),
#                      factor_key=T)

metadata_reads_ONT=metadata_reads_ONT[c("Row.names","seq_run","Platform","sample_sums(ps_noncontam)")]
metadata_reads_Illumina=metadata_reads_Illumina[c("Row.names","seq_run","Platform","sample_sums(ps_noncontam)")]

combined_reads = rbind(metadata_reads_ONT, metadata_reads_Illumina)

options(scipen = 999)
set.seed=6
combined_plot = ggplot(combined_reads, aes(x=seq_run, y=`sample_sums(ps_noncontam)`, color=Platform)) + 
            geom_boxplot() +
            scale_color_manual(values=c("#F8766D","#00BFC4")) +
            geom_jitter(width=0.1, height=0) +
            #facet_wrap(~attributes, scales="free_y") +
            scale_y_continuous(limits = c(0,NA)) +
            theme_bw() +
            stat_summary(fun=mean, geom="point", shape=23, size=3, color="black", fill="black") +
            theme(text = element_text(size = 16), 
                  axis.title.x = element_text(margin = margin(t = 20)),
                  axis.title.y = element_text(margin = margin(r = 20)), 
                  axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
            labs(y = "Number of reads per sample", x = "Sequencing runs and platforms") +
            scale_x_discrete(labels = c("ONT Run 1", "ONT Run 2", "ONT Run 3", "Illumina"))
combined_plot

ggsave("~/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/results/figures/combined_reads.pdf", height = 4, width = 7, units = "in")

########## Plot summary of total reads and genera detected ##########

total_reads = read_csv("~/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/reads_distribution.csv")
total_reads=as.data.frame(total_reads)

ggplot(total_reads, aes(x=Run_ID, y=Total_number_of_reads, color=Platform, fill=Platform), levels=c("ONT","Illuma")) + geom_bar(stat="identity") + theme_bw() +
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)), 
        axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
  labs(y = "Total number of reads per run", x = "Sequencing runs and platforms") +
  scale_x_discrete(labels = c("ONT Run 1", "ONT Run 2", "ONT Run 3", "Illumina"))

ggsave("~/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/results/figures/total_reads.pdf", height = 4, width = 7, units = "in")

ggplot(total_reads, aes(x=Run_ID, y=`Number_of_genera_detected_(SILVA)`, color=Platform, fill=Platform), levels=c("ONT","Illuma")) + geom_bar(stat="identity") + theme_bw() +
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)), 
        axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
  labs(y = "Total number of bacteria genera detected per run (SILVA)", x = "Sequencing runs and platforms") +
  scale_x_discrete(labels = c("ONT Run 1", "ONT Run 2", "ONT Run 3", "Illumina"))

ggsave("~/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/results/figures/total_genera.pdf", height = 4, width = 7, units = "in")
