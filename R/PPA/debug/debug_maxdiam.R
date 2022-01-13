library(tidyverse)

# Raw BCI data
# At https://datadryad.org/stash/dataset/doi:10.15146/5xcp-0d46
load("C:/CODE/gwd_metrics/data/doi_10.15146_5xcp-0d46__v2/bci.spptable.rdata")
load("C:/CODE/gwd_metrics/data/doi_10.15146_5xcp-0d46__v2/bci.tree/bci.tree8.rdata")

# Kambach data formatted
ppa_bci <- read.csv( "C:/CODE/gwd_metrics/R/PPA/input/ppa_bci.csv" )

# Kambach data
kamback <- read.table('C:/CODE/gwd/data/Demographies_with_dbh_means_Kambach.txt',
                      header = T ) 

# Maximum DBH  abundances in BCI
bci_max_dbh <- bci.tree8 %>% 
                group_by( sp ) %>% 
                summarise( dbh_max_bci = max(dbh, na.rm=T) ) %>% 
                ungroup %>% 
                # remove species with no DBH information
                subset( !(dbh_max_bci == -Inf) ) 

# Kambach's estimates
kam_max_dbh <- kamback %>% 
                subset( site == 'bci' ) %>% 
                select( sp, max_dbh, 
                        max_dbh_6_largest_ind_mean ) %>% 
                rename( dbh_max_absolute_max_kam       = max_dbh, 
                        dbh_max_6_largest_ind_mean_kam = max_dbh_6_largest_ind_mean )
  

max_dbh_df  <- left_join( bci_max_dbh, kam_max_dbh ) %>% 
                 gather( dbh_estimate, dbh_max_kam, 
                         dbh_max_absolute_max_kam:dbh_max_6_largest_ind_mean_kam ) %>% 
                 mutate( dbh_estimate = gsub('dbh_max_','',dbh_estimate) ) %>% 
                 mutate( dbh_estimate = gsub('_kam','',dbh_estimate) ) 
  
# Check two estimates of max DBH
ggplot( max_dbh_df ) +
  geom_point( aes(dbh_max_bci, 
                  dbh_max_kam) ) + 
  facet_wrap( ~ dbh_estimate ) +
  theme_minimal() +
  labs( x = 'Max. DBH - BCI raw data (mm)',
        y = 'Max. DBH - Kambach data (mm)' ) +
  theme( axis.text.x = element_text( angle = 90,
                                     hjust = 1,
                                     vjust = 0.5),
         axis.text.y = element_text( angle = 0,
                                     hjust = 0,
                                     vjust = 0.5),
         strip.text  = element_text( size  = 15 ) )
  ggsave( 'C:/CODE/gwd_metrics/results/PPA/plots/debug/dbh_max_check.tiff',
          width = 6.3, height = 4, compression = 'lzw')

# max_dbh_df %>% 
#   subset( dbh_estimate == 'absolute_max' ) %>% 
#   ggplot(  ) +
#     geom_point( aes(dbh_max_kam,
#                     dbh_max_bci) ) + 
#     facet_wrap( ~ dbh_estimate ) +
#     theme_minimal() +
#     labs( x = 'Max. DBH - BCI raw data (mm)',
#           y = 'Max. DBH - Kambach data (mm)' ) +
#     theme( axis.text.x = element_text( angle = 90,
#                                        hjust = 1,
#                                        vjust = 0.5),
#            axis.text.y = element_text( angle = 0,
#                                        hjust = 0,
#                                        vjust = 0.5),
#            strip.text  = element_text( size  = 15 ) )  
  

cia <- 
left_join( select( ppa_bci, sp, max_dbh ), 
           select( kamback, sp, max_dbh_6_largest_ind_mean ) %>% 
             rename( max_dbh_kam = max_dbh_6_largest_ind_mean ) ) 

ggplot(cia) +
  geom_point( aes(max_dbh, max_dbh_kam) )

plot(cia$max_dbh, cia$max_dbh_kam)
