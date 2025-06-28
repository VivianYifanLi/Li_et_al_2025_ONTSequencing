# Idenitfy shared vs different taxa between sequencing runs. #

# Load environmental variables.
load('figure6-out.RData')

## Venn diagrams after CLR transformation.
library(ggVennDiagram)

ps_clr1 = subset_samples(ps_noncontam, seq_run==1)
df_clr1=psmelt(ps_clr1)
df_clr1=subset(df_clr1, Abundance>1)

ps_clr2 = subset_samples(ps_noncontam, seq_run==2)
df_clr2=psmelt(ps_clr2)
df_clr2=subset(df_clr2, Abundance>1)

ps_clr3 = subset_samples(ps_noncontam, seq_run==3)
df_clr3=psmelt(ps_clr3)
df_clr3=subset(df_clr3, Abundance>1)

vennsp = list("Run 1"=df_clr1$Genus, "Run 2"=df_clr2$Genus, "Run 3"=df_clr3$Genus)
vennsp_seqruns = ggVennDiagram((vennsp), edge_size=1, label_size = 10, set_color = c("#a5cde2", "#97d174", "#f06667")) +
                  scale_color_brewer(palette = "Paired") +
                  scale_fill_distiller(palette = "Blues", direction = 1) +
                  theme(text=element_text(size=20))
vennsp_seqruns
ggsave("/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/results/figures/vennsp_seqruns1.pdf",vennsp_seqruns)

ps_clr1 = subset_samples(ps_clr_noc, bio_rep==1)
df_clr1=psmelt(ps_clr1)
ps_clr2 = subset_samples(ps_clr_noc, bio_rep==2)
df_clr2=psmelt(ps_clr2)
ps_clr3 = subset_samples(ps_clr_noc, bio_rep==3)
df_clr3=psmelt(ps_clr3)

vennsp = list("Rep 1"=df_clr1$Species, "Rep 2"=df_clr2$Species, "Rep 3"=df_clr3$Species)
vennsp_bioreps = ggVennDiagram(vennsp, set_color=c("blue","green","red"), edge_size=0.3)
ggsave("/Users/vivianli/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/results/w60-allruns/figures/vennsp_bioreps.png",vennsp_bioreps)