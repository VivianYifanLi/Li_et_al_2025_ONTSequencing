# Plot summary of reads distribution

library(tidyr)

reads_all_runs = read.table("/Users/yifanli/Documents/PhD/Projects/Algal_microbiome_characterization/Nanopore_pilot/reads_distribution.csv", sep=",", header=T)

reads_tidy <- gather(reads_all_runs,
                      key = "attributes",
                      value = "values",
                      -c(Platform, Run_ID),
                      factor_key=T)

options(scipen = 999)
ggplot(reads_tidy, aes(x=Platform, y=values, fill=Platform)) +
  geom_boxplot() +
  facet_wrap(~attributes, scales="free_y") +
  scale_y_continuous(limits = c(0,NA)) +
  theme_bw()
