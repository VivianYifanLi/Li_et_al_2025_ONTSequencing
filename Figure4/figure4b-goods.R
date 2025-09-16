# Calculate Good's coverage indices #

# Load environmental variables.
load('figure3-out.RData')
library(ggplot2)
library(dplyr)
library(phyloseq)

# Remove negative controls.
ps_prok_noc <- subset_samples(ps_noncontam, cultureID != "Bio_neg_control")
ps_prok_noc <- subset_samples(ps_prok_noc, cultureID != "PCR_neg_control")
ps_prok_noc
df_prok_noc <- psmelt(ps_prok_noc)

# For run 1
ps_prok_noc1=subset_samples(ps_prok_noc, seq_run==1)
df_prok_noc1=psmelt(ps_prok_noc1)

goods_run1=df_prok_noc1 %>% 
  group_by(sampleID) %>% 
  summarize(n_seqs=sum(Abundance),
            n_sing=sum(Abundance==1),
            goods = 100 * (1 - n_sing/n_seqs)) %>%
  ggplot(aes(x=n_seqs, y=goods, label=sampleID)) + theme_bw() + geom_point(aes(color = sampleID, shape = sampleID, size=3)) + scale_color_manual(values=setpalette(18)) + scale_shape_manual(values=1:18) + geom_text(vjust=0, hjust=-0.3, size=5) + scale_y_continuous(name="Good's coverage index (%)", limits=c(99.9, 100), n.breaks=8) + theme(text=element_text(size=15)) + xlab("Number of reads") + theme(text=element_text(size=20)) + guides(color=FALSE, shape=FALSE, size=FALSE)
goods_run1
ggsave("/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/results/figures/goods_run1.pdf",goods_run1, width=16, height=6)
