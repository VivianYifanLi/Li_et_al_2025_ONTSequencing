library(tidyverse)
abund = as.data.frame(read_csv("/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/zymo_theoretical.csv"))

zymo_gen=ggplot(data=abund, aes(Sample), width=0.5) + theme_bw() +
  geom_bar(aes(color=Species, fill=Species,y=Abundance), stat="identity",position="fill", width=0.5) +
  #facet_grid(. ~ factor(seq_run), scales="free_x", space="free_x") +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=c("Bacillus subtilis"="#a5cee2", 
                             "Enterococcus faecalis"="#3f8faa",
                             "Escherichia coli"="#79c360",
                             "Lactobacillus fermentum"="#f88a89",
                             "Listeria monocytogenes"="#f06c45",
                             "Pseudomonas aeruginosa"="#fc870f",
                             "Salmonella enterica"="#b295c7",
                             "Staphylococcus aureus"="#c7b699")) + 
  scale_color_manual(values = c("Bacillus subtilis"="#a5cee2", 
                                "Enterococcus faecalis"="#3f8faa",
                                "Escherichia coli"="#79c360",
                                "Lactobacillus fermentum"="#f88a89",
                                "Listeria monocytogenes"="#f06c45",
                                "Pseudomonas aeruginosa"="#fc870f",
                                "Salmonella enterica"="#b295c7",
                                "Staphylococcus aureus"="#c7b699")) + 
  theme(text=element_text(size=20))
zymo_gen
ggsave("/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/results/figures/zymo_theoretical.pdf",zymo_gen, width = 6, height = 10)
