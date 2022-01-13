library(tidyverse)
library(ggthemes)

# read and set up the data
sim_df    <- read.table( "data/out_skratch_FG5.txt", sep="\t", header=TRUE)
orig_df   <- read.table( "data/out_secondary_FG5.txt", sep="\t", header=TRUE)
fgdata    <- read.table( 'data/PPA_FG5.txt',sep="\t",header=TRUE)
cohort    <- readRDS( 'data/out_skratch_FG5.rds' )
# fgdata_d  <- read.table( 'data/PPA_FG5_maxD.txt',sep="\t",header=TRUE)
fgdata_d  <- read.table( 'data/out_skratch_FG5.txt',sep="\t",header=TRUE)
group_n   <- dplyr::select(fgdata, fg, fgname ) %>% 
              rename( `Functional group` = fgname ) 

# Age distribution -------------------------------------------------------------

# main matrix with columns: 
# (1) cohort diameter 
# (2) number of individuals 
# (3) crown class 
# (4) functional group
# (5) age

# debugging ------
final_df <- cohort[[200]][[1]] %>% 
  as.data.frame() %>% 
  setNames( c('diam','n','layer','fg','age') )

fgdata_d %>%
  dplyr::select( time, max_d_1:max_d_5 ) %>% 
  subset( time == 1000 )

final_df %>% 
  subset( fg == 1 ) %>% 
  # subset( layer == 1 ) %>% 
  subset( diam > (109.5/2)  ) %>%
  .$n %>% 
  sum


final_df <- cohort[[200]][[1]] %>% 
              as.data.frame() %>% 
              setNames( c('diam','n','layer','fg','age') ) %>% 
              mutate( n = n*1000 ) 
       
final_maxD <- subset( fgdata_d, time == 1000 ) %>% 
                dplyr::select( max_d_1:max_d_5 ) %>% 
                gather( fg, max_d, max_d_1:max_d_5) %>% 
                mutate( fg = gsub('max_d_','',fg) ) %>% 
                mutate( fg = as.numeric(fg) ) %>% 
                left_join( group_n ) %>% 
                mutate( `Functional group` = gsub('Intermediate',
                                                  'Intermed.',
                                                  `Functional group`) 
                )

final_df %>% 
  subset( layer == 1 ) %>%
  group_by( fg ) %>% 
  summarise( min_age = min(age) ) %>% 
  ungroup

# repicate data frame rows
rep_df_rows <- function( ii ){
  
  out <- final_df[ii,]
  out <- out[rep(1, round(final_df[ii,]$n)),]
  # replicate( , 
  #            final_df[ii,], simplify=F ) %>% 
  #   bind_rows
  out
}

# replicated 
final_rep_df <- lapply( 1:nrow(final_df), 
                        rep_df_rows) %>% 
                  bind_rows 

final_rep_df %>% 
  left_join( group_n ) %>% 
  mutate( `Functional group` = gsub('Intermediate',
                                    'Intermed.',
                                    `Functional group`) 
          ) %>% 
  mutate( layer = paste0('layer_',layer) ) %>% 
  ggplot( ) +
  geom_density( aes(age) ) +
  facet_grid( layer ~ `Functional group`,
              scale = 'free_y') +
  theme_minimal() +
  labs( x = 'Age',
        y = 'Kernel density estimate' ) +
  theme( strip.text  = element_text( size  = 15 ),
         axis.title  = element_text( size  = 20 ),
         axis.text.x = element_text( size  = 10,
                                     angle = 90,
                                     hjust = 1,
                                     vjust = 0.5) ) +
  ggsave( 'results/rugers/hist_age.tiff',
          width = 6.3, height = 6.3, 
          compression = 'lzw' )
  
final_rep_df %>% 
    left_join( group_n ) %>% 
    mutate( `Functional group` = gsub('Intermediate',
                                      'Intermed.',
                                      `Functional group`) 
    ) %>% 
    mutate( layer = paste0('layer_',layer) ) %>% 
    ggplot( ) +
    geom_density( aes(diam) ) +
    geom_vline( data = final_maxD,
                aes(xintercept = (max_d/2) ),
                lty = 2, color = 'red' ) +
    facet_grid( layer ~ `Functional group`,
                scale = 'free_y') +
    theme_minimal() +
    labs( x = 'Size',
          y = 'Kernel density estimate' ) +
    theme( strip.text  = element_text( size  = 15 ),
           axis.title  = element_text( size  = 20 ),
           axis.text.x = element_text( size  = 10,
                                       angle = 90,
                                       hjust = 1,
                                       vjust = 0.5) ) +
    ggsave( 'results/rugers/hist_size.tiff',
            width = 6.3, height = 6.3, 
            compression = 'lzw' )
  

final_rep_df %>% 
  left_join( group_n ) %>% 
  subset( `Functional group` == 'Intermediate' & 
           layer == 1 ) %>% 
  ggplot() +
  geom_histogram( aes(age) ) 



# Abundance dynamics -----------------------------------------------------------
fgdata_d %>%
# orig_df %>% 
  dplyr::select( time, n_Slow:agb_Intermediate ) %>% 
  gather( type, value, n_Slow:agb_Intermediate ) %>% 
  separate( type, c('metric','fg'),sep='_') %>% 
  mutate( metric = gsub('agb','Biomass',metric) ) %>% 
  mutate( metric = gsub('ba', 'Basal area',metric) ) %>% 
  mutate( metric = gsub('n',  'Number',metric) ) %>% 
  ggplot(  ) +
  geom_point( aes(time, value, 
                  group = fg,
                  color = fg),
              alpha = 0.3) +
  facet_wrap(~metric, scale = 'free') +
  theme_minimal() +
  scale_color_colorblind() +
  labs( x = 'Time',
        y = 'Abundance',
        color = 'Functional Group') +
  theme( axis.text.x = element_text(angle = 70) ) +
  ggsave( 'results/rugers/abundance_orig.tiff',
          width = 6.3, height = 3, 
          compression = 'lzw' )


# Age dynamics -----------------------------------------------------------------
age_by_time <- function( tt ){
  cohort[[tt]] %>% 
    as.data.frame %>% 
    group_by( cohorts.4 ) %>% 
    summarise(  = weighted.mean(cohorts.5, cohorts.2) ) %>% 
    ungroup %>% 
    rename( fg = cohorts.4 ) %>% 
    left_join( group_n ) %>% 
    mutate( time = tt*5,
            fg   = as.factor(fg) ) 
}

lapply( 1:200, age_by_time) %>% 
  bind_rows %>% 
  ggplot() +
  geom_point( aes(time, avg_age,
                  group = `Functional group`,
                  color = `Functional group`) ) + 
  theme_minimal() +
  scale_color_colorblind() +
  labs( y = 'Average age of functional group',
        x = 'Time' ) +
  theme( axis.title   = element_text( size = 15),
         legend.title = element_text( size = 15),
         legend.text = element_text( size = 15) ) +
  ggsave( 'results/rugers/avg_age.tiff',
          width = 5, height = 4, 
          compression = 'lzw')
  

# generation time based on canopy trees ----------------------------------------
fgdata_d %>% 
  dplyr::select( time, gt_1:gt_5 ) %>% 
  gather( fg, gt, gt_1:gt_5) %>% 
  mutate( fg = gsub('gt_','',fg) ) %>% 
  mutate( fg = as.numeric(fg) ) %>% 
  left_join( group_n ) %>% 
  ggplot() +
  geom_point( aes( time, gt, 
                   color = `Functional group`) ) +
  scale_color_colorblind() +
  theme_minimal( ) +
  ylab( 'Average age of parents' ) +
  ggtitle( 'Parents as canopy trees' ) +
  theme( plot.title = element_text( hjust = 0.5),
         axis.title   = element_text( size = 15),
         legend.title = element_text( size = 15),
         legend.text = element_text( size = 15) ) +
  ggsave( 'results/rugers/generation_time_canopy.tiff',
          width = 5, height = 4, 
          compression = 'lzw' )

# # generation time based on canopy trees
# # But age calculated via diam / growth
# sim_df %>% 
#   dplyr::select( time, gt_g_1:gt_g_5 ) %>% 
#   gather( fg, gt, gt_g_1:gt_g_5 ) %>% 
#   mutate( fg = gsub('gt_g_','',fg) ) %>% 
#   mutate( fg = as.numeric(fg) ) %>% 
#   left_join( group_n ) %>% 
#   ggplot() +
#   geom_point( aes( time, gt, 
#                    color = `Functional group`) ) +
#   scale_color_colorblind() +
#   theme_minimal( ) +
#   ylab( 'Average age of parents' ) +
#   ggtitle( 'Parents as canopy trees' ) +
#   theme( plot.title = element_text( hjust = 0.5) ) 

  
# generation time based on half of maximum diameter
fgdata_d %>% 
  dplyr::select( time, gt_diam_1:gt_diam_5 ) %>%
  gather( fg, gt, gt_diam_1:gt_diam_5 ) %>% 
  mutate( fg = gsub('gt_diam_','',fg) ) %>% 
  mutate( fg = as.numeric(fg) ) %>% 
  left_join( group_n ) %>% 
  ggplot() +
  geom_point( aes( time, gt, 
                   color = `Functional group`) ) +
  scale_color_colorblind() +
  theme_minimal( ) +
  ylab( 'Average age of parents' ) + 
  ggtitle( 'Parents at half the maximum diameter' ) +
  theme( plot.title   = element_text( hjust = 0.5),
         axis.title   = element_text( size = 15),
         legend.title = element_text( size = 15),
         legend.text  = element_text( size = 15) ) +
  ggsave( 'results/rugers/generation_time_maxDiam.tiff',
          width = 5, height = 4, 
          compression = 'lzw' )


# Maximum diameter
fgdata_d %>% 
  dplyr::select( time, max_d_1:max_d_5 ) %>% 
  gather( fg, max_d, max_d_1:max_d_5 ) %>% 
  left_join( mutate(group_n, fg = paste0('max_d_',fg)) ) %>% 
  ggplot() +
  geom_point( aes( time, max_d, 
                   color = `Functional group` ) ) +
  scale_color_colorblind() +
  theme_minimal() +
  labs( x = 'Time',
        y = 'Maximum diameter' ) +
  theme( plot.title   = element_text( hjust = 0.5),
         axis.title   = element_text( size = 15),
         axis.text.x  = element_text( size  = 10,
                                      angle = 80),
         legend.title = element_text( size = 15),
         legend.text  = element_text( size = 15) ) +
  ggsave( 'results/rugers/maxDiam_vs_time.tiff',
          width = 5, height = 4, compression = 'lzw' )


# Store maximum diameters
sim_df %>% 
  subset( time == 1000 ) %>% 
  dplyr::select( max_d_1:max_d_5 ) %>% 
  gather( fg, max_d, max_d_1:max_d_5 ) %>% 
  mutate( fg = gsub('max_d_','',fg) ) %>% 
  mutate( fg = as.numeric(fg) ) %>% 
  right_join( fgdata ) %>% 
  dplyr::select( fg, fgname:wd, max_d ) %>% 
  write.table( paste(mainfolder,"PPA_FG5_maxD.txt",sep=""), sep="\t" )
