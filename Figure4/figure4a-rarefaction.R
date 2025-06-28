# Plot rarefaction curves #

library(iNEXT)
library(ggplot2)
load('figure3-out.RData')

# Remove negative controls.
ps_prok_noc <- subset_samples(ps_noncontam, cultureID != "Bio_neg_control")
ps_prok_noc <- subset_samples(ps_prok_noc, cultureID != "PCR_neg_control")
ps_prok_noc
df_prok_noc <- psmelt(ps_prok_noc)

# For run 1
ps_prok_noc1=subset_samples(ps_prok_noc, seq_run==1)
df_prok_noc1=psmelt(ps_prok_noc1)

otu_rare_noc1=data.frame(otu_table(ps_prok_noc1))
#otu_rare_nolow1=t(otu_rare_nolow1)
#rarecurve_run1=rarecurve(otu_rare1, label=F)

otu_rare1_inext=iNEXT(otu_rare_noc1, endpoint=50000) # set x-axis endpoint depending on data
g=ggiNEXT(otu_rare1_inext, type=1, se=F) + theme_bw() + ylim(0,65) + xlab("Number of reads") + scale_color_manual(values=setpalette(18)) + scale_shape_manual(values=1:18) + guides(color=FALSE)
gb1=ggplot_build(g + theme(text=element_text(size=20), legend.text = element_text(size = 20))) 
gb1$data[[2]]$size <- 1
gb1$data[[2]]$linewidth <- 1
gt1 <- ggplot_gtable(gb1)
library(grid)
gt1_plot=grid.draw(gt1)

# For run 2
ps_prok_noc2=subset_samples(ps_prok_noc, seq_run==2)
df_prok_noc2=psmelt(ps_prok_noc2)

otu_rare_noc2=data.frame(otu_table(ps_prok_noc2))

otu_rare2_inext=iNEXT(otu_rare_noc2, endpoint = 22000) # set x-axis endpoint depending on data
g2=ggiNEXT(otu_rare2_inext, type=1, se=F) + theme_bw() + ylim(0, 65) + xlab("Number of reads")
gb2=ggplot_build(g2 + theme(legend.text = element_text(size = 20)))
gb2$data[[2]]$size <- 1
gb2$data[[2]]$linewidth <- 1
gt2 <- ggplot_gtable(gb2)
gt2_plot=grid.draw(gt2)

# For run 3
ps_prok_noc3=subset_samples(ps_prok_noc, seq_run==3)
df_prok_noc3=psmelt(ps_prok_noc3)

otu_rare_noc3=data.frame(otu_table(ps_prok_noc3))

otu_rare3_inext=iNEXT(otu_rare_noc3, endpoint = 22000) # set x-axis endpoint depending on data
g3=ggiNEXT(otu_rare3_inext, type=1, se=F) + theme_bw() + ylim(0, 65) + xlab("Number of reads")
gb3=ggplot_build(g3 + theme(legend.text = element_text(size = 20)))
gb3$data[[2]]$size <- 1
gb3$data[[2]]$linewidth <- 1
gt3 <- ggplot_gtable(gb3)
gt3_plot=grid.draw(gt3)
