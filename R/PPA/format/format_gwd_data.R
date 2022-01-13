library(tidyverse)

# format
bci_dat  <- read.csv( 'data/gwd_db/bci_demog_means.csv' )
bci_tax  <- read.csv( 'data/gwd_db/bci_taxa.csv' ) %>% 
              select( LCVP_Accepted_Taxon, Submitted_Name ) %>% 
              rename( Latin = Submitted_Name )
spVitals <- read.table( 'C:/CODE/gwd_metrics/R/PPA/input/PPA_FG5_filtered.csv', 
                        sep = ",", header = TRUE)
bci_adj  <- read.csv( 'R/PPA/input/bci_adj_fec.csv' )

# select 
n_sel <- c( 'LCVP_Accepted_Taxon',
            paste0('growth_layer',1:4,'_mean'),
            paste0('survival_layer',1:4,'_mean'),
            'recruits_per_year_mean' )

# write down properly formatted file
bci_dat[,n_sel] %>% 
  mutate( wd = 0.55 ) %>%
  setNames( names(spVitals) ) %>% 
  mutate( LCVP_Accepted_Taxon   = SpeciesID,
          max_dbh               = bci_dat$max_dbh_6_largest_ind_mean ) %>% 
  arrange( SpeciesID ) %>% 
  mutate( SpeciesID = as.factor(SpeciesID) %>% as.numeric ) %>% 
  # REMOVE "double". Temporary. This need be addressed by Ollie
  subset( !(LCVP_Accepted_Taxon == 'Piper aequale Vahl') ) %>% 
  # Add taxonomic data
  left_join( bci_tax ) %>% 
  # Add "adjusted fecundity"
  left_join( bci_adj ) %>% 
  # REMOVE adj == NA. Temporary. Should be resolved by Ollie
  subset( !(is.na(adj_ba) | is.na(adj_n)) ) %>% 
  # calculate "death rates"
  mutate( mu1 = 1 - mu1,
          mu2 = 1 - mu2 ,
          mu3 = 1 - mu3,
          mu4 = 1 - mu4 ) %>% 
  # calculate growth in mm
  mutate( G1 = G1/10,
          G2 = G2/10,
          G3 = G3/10,
          G4 = G4/10 ) %>% 
  # adjust fecundity! 
  mutate( F = ifelse( rep(ADJUST_BY_BA, nrow(.)),
                      `F` * adj_ba,
                      `F` * adj_n) ) %>% 
  write.csv( 'R/PPA/input/ppa_bci.csv',
             row.names = F )
