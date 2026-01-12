#### Jinxin Meng, 20250502, 20250922 ####
setwd('/data/mengjx/tmp_proj/05.20250905_apex_EP_from_BGC/')
pacman::p_load(tidyverse, ggpubr)

#### AA composition ####
data <- rbind(
  read.delim('final_AMPs.AA.pct.tsv') %>% add_column(class = 'EPs'),
  read.delim('amps.pep.ref.removeXBZ.AA.pct.tsv') %>% add_column(class = 'Database')
  ) %>% 
  mutate(AA = factor(AA))

ggbarplot(data, 'AA', 'pct', fill = 'class', position = position_dodge2(), color = NA,
          palette = c('#B291B5','#3BA997'), xlab = 'Amino Acid', ylab = 'Normalized frequency',
          legend = 'right') +
  scale_y_continuous(expand = c(0, 0), position = 'right') +
  theme(axis.ticks.length = unit(2, 'mm'))
ggsave('Figure/AA_pct.pdf', width = 8, height = 4)

#### peptide properties ####

data <- rbind(
  data.table::fread('amps.pep.ref.removeXBZ.globalDescriptor.tsv') %>% 
    add_column(class = 'Database'),
  data.table::fread('final_AMPs.globalDescriptor.tsv') %>% 
    add_column(class = 'EPs')
)

# Charge
p1 <- ggboxplot(data, 'class', 'Charge', fill = 'class', outliers = F, legend = 'none',
          xlab = '', ylab = 'Net Charge', palette = c('#B291B5','#3BA997')) +
  ggsignif::geom_signif(comparisons = list(c('Database', 'EPs')), 
                        y_position = 9, tip_length = .005) +
  theme(aspect.ratio = 3/2, 
        axis.ticks.length = unit(2, 'mm'))

# InstabilityInd
p2 <- ggboxplot(data, 'class', 'InstabilityInd', fill = 'class', outliers = F, legend = 'none',
                xlab = '', ylab = 'Instability Index', palette = c('#B291B5','#3BA997')) +
  ggsignif::geom_signif(comparisons = list(c('Database', 'EPs')), 
                        y_position = 125, tip_length = .01) +
  theme(aspect.ratio = 3/2, 
        axis.ticks.length = unit(2, 'mm'))

# AliphaticInd
p3 <- ggboxplot(data, 'class', 'AliphaticInd', fill = 'class', outliers = F, legend = 'none',
          xlab = '', ylab = 'Aliphatic Index', palette = c('#B291B5','#3BA997')) +
  ggsignif::geom_signif(comparisons = list(c('Database', 'EPs')), 
                        y_position = 220, tip_length = .015) +
  theme(aspect.ratio = 3/2, 
        axis.ticks.length = unit(2, 'mm'))

# BomanInd
p4 <- ggboxplot(data, 'class', 'BomanInd', fill = 'class', outliers = F, legend = 'none',
                xlab = '', ylab = 'Boman Index', palette = c('#B291B5','#3BA997')) +
  ggsignif::geom_signif(comparisons = list(c('Database', 'EPs')), 
                        y_position = 6, tip_length = .015) +
  theme(aspect.ratio = 3/2, 
        axis.ticks.length = unit(2, 'mm'))

data <- rbind(
  data.table::fread('amps.pep.ref.removeXBZ.H_uH.tsv') %>% 
    add_column(class = 'Database'),
  data.table::fread('final_AMPs.H_uH.tsv') %>% 
    add_column(class = 'EPs')
)

# H_Eisenberg
p5 <- ggboxplot(data, 'class', 'H_Eisenberg', fill = 'class', outliers = F, legend = 'none',
          xlab = '', ylab = 'Eisenberg Hydrophobicity', palette = c('#B291B5','#3BA997')) +
  ggsignif::geom_signif(comparisons = list(c('Database', 'EPs')), 
                        y_position = .8, tip_length = .02) +
  theme(aspect.ratio = 3/2, 
        axis.ticks.length = unit(2, 'mm'))

cowplot::plot_grid(p1, p2, p3, p4, p5, nrow = 1, align = 'v')
ggsave('Figure/Pep.properties.pdf', width = 14, height = 4)

