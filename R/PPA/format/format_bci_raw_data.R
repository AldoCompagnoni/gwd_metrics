library(tidyverse)
library(gridExtra)

# read stem and taxonomic information
# At https://datadryad.org/stash/dataset/doi:10.15146/5xcp-0d46
load("C:/CODE/gwd_metrics/data/bci_raw_data_2/bci.spptable.rdata")
load("C:/CODE/gwd_metrics/data/bci_raw_data_2/bci.tree/bci.tree8.rdata")

# Basal areas of each species 
bci_ba_df <- bci.tree8 %>% 
                  # remove dead trees
                  subset( !(ba == 0) ) %>% 
                  group_by( sp ) %>% 
                  summarise( ba = sum(ba) ) %>% 
                  ungroup %>% 
                  left_join( select(bci.spptable,sp,Latin) )

# Number of individual of each species
bci_n_df  <- bci.tree8 %>% 
                  # remove dead trees
                  subset( !(ba == 0) ) %>% 
                  count( sp ) %>% 
                  left_join( select(bci.spptable,sp,Latin) )

# relative abundances in BCI
bci_ra_df <- left_join( bci_ba_df, bci_n_df ) %>% 
              mutate( prop_ba = ba / sum(ba),
                      prop_n  = n  / sum(n) ) %>% 
              mutate( adj_ba  = 1/prop_ba,
                      adj_n   = 1/prop_n )

# store results
write.csv( bci_ra_df, 
           'R/PPA/input/bci_adj_fec.csv',
           row.names = F )

# plots ----------------------


p1 <- ggplot(bci_ra_df) +
        geom_point( aes(prop_ba, prop_n),
                    size = 4, alpha = 0.5) +
        theme_minimal() +
        labs( x = 'Proportions based on basal area',
              y = 'Proportions based on abundance')
  
p2 <- ggplot(bci_ra_df) +
        geom_point( aes(log(prop_ba), log(prop_n)),
                    size = 4, alpha = 0.5 ) +
        theme_minimal() +
    
        labs( x = 'log(Proportions based on basal area)',
              y = 'log(Proportions based on abundance)')

tiff( 'results/PPA/plots/summary/community_props.tiff',
      unit = 'in', res = 600, 
      width = 6.3, height = 6.3, compression = 'lzw') 

grid.arrange(p1,p2, ncol=1)

dev.off()

