# 2. Summary statistics #

# Load environmental variables.
load('myEnvironment.RData')

# All sequences before and after data processing.

ps_prok_noc <- subset_samples(ps_noncontam, cultureID != "Bio_neg_control")
ps_prok_noc <- subset_samples(ps_prok_noc, cultureID != "PCR_neg_control")
ps_prok_noc
df_prok_noc <- psmelt(ps_prok_noc)

sum(sample_sums(ps_prok_noc))
mean(sample_sums(ps_prok_noc))
median(sample_sums(ps_prok_noc))
min(sample_sums(ps_prok_noc))
max(sample_sums(ps_prok_noc))

ps_prok_noc1 = subset_samples(ps_prok_noc, seq_run==1)
df_prok_noc1=psmelt(ps_prok_noc1)
df_prok_noc1=subset(df_prok_noc1, Abundance>1)

ps_prok_noc2 = subset_samples(ps_prok_noc, seq_run==2)
df_prok_noc2=psmelt(ps_prok_noc2)
df_prok_noc2=subset(df_prok_noc2, Abundance>1)

ps_prok_noc3 = subset_samples(ps_prok_noc, seq_run==3)
df_prok_noc3=psmelt(ps_prok_noc3)
df_prok_noc3=subset(df_prok_noc3, Abundance>1)

sum(sample_sums(ps_prok_noc1))
mean(sample_sums(ps_prok_noc1))
median(sample_sums(ps_prok_noc1))
min(sample_sums(ps_prok_noc1))
max(sample_sums(ps_prok_noc1))

sum(sample_sums(ps_prok_noc2))
mean(sample_sums(ps_prok_noc2))
median(sample_sums(ps_prok_noc2))
min(sample_sums(ps_prok_noc2))
max(sample_sums(ps_prok_noc2))

sum(sample_sums(ps_prok_noc3))
mean(sample_sums(ps_prok_noc3))
median(sample_sums(ps_prok_noc3))
min(sample_sums(ps_prok_noc3))
max(sample_sums(ps_prok_noc3))

# Remove taxa with low abundance according to neg controls.
#ps_prok_pruned <- prune_taxa(taxa_sums(ps_prok_nosing) > 9, ps_prok_nosing)
#df_prok_pruned=psmelt(ps_prok_pruned)
#sum(sample_sums(ps_prok_pruned))

# Plot number of reads per sample
reads = as.data.frame(sample_sums(ps_prok_noc))
reads$barcodes = row.names(reads)
reads$barcode_seqrun = row.names(reads)
colnames(reads) = c("reads","barcodes","barcodes_seqrun")
reads <- reads %>% separate(barcodes_seqrun, c("barcode", "seq_run"), sep= "_") %>% mutate(barcode=str_remove(barcode, "barcode"))

total_reads=ggplot(reads, aes(x=barcode, y=reads, fill=barcode)) + 
  theme_grey(base_size=20) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=mypalette(24)) +
  theme(axis.text.x=element_text(angle=90, vjust=1, hjust=1, size=10)) +
  facet_grid(. ~ factor(seq_run)) 
total_reads

# Save plot as PDF.
ggsave("/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/results/figures/total_reads.pdf",total_reads)

# Plot number of reads per culture after pooling (everything; NOT nolow).
ps_per_culture=merge_samples(ps_prok_nosing, "cultureID", fun=sum)
df_per_culture=psmelt(ps_per_culture)

reads = as.data.frame(sample_sums(ps_per_culture))
reads$cultureID = row.names(reads)
colnames(reads) = c("reads","cultureID")

library(scales)
reads_per_culture=ggplot(reads, aes(x=cultureID, y=reads, fill=cultureID)) + 
  theme_grey(base_size=20) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=mypalette(25)) +
  theme(axis.text.x=element_text(angle=90, vjust=1, hjust=1, size=10)) +
  scale_y_continuous(labels = label_comma())
reads_per_culture

# Save plot as PDF.
ggsave("/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/allseqs_thesisChp1/results/figures/total_reads_perculture.pdf",reads_per_culture)

# Save environmental variables.
save.image(file='myEnvironment.RData')