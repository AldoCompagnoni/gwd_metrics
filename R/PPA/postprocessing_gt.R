# TODO Attach species PCA traits to the processed data.frames

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(ggthemes)
library(GGally)

# Directories and data sets -----------------------------------------------

base_directory    <- dirname(rstudioapi::getActiveDocumentContext()$path)
output_directory  <- paste0( getwd(), "/results/PPA/output/")
input_directory   <- paste0(base_directory, "/input/")
# species_file      <- paste0(base_directory, "/input/PPA_FG5_filtered.csv")
species_file      <- paste0(base_directory, "/input/ppa_bci.csv")

cohorts           <- readRDS(paste0(output_directory, "PPA_output_raw_cohort_ba_1.rds"))
mortality         <- readRDS(paste0(output_directory, "PPA_output_raw_cohort_mortality_ba_1.rds"))
spVitals          <- read.table(species_file, sep = ",", header = TRUE)
rugers            <- readxl::read_xlsx( 'data/rugers_et_al_2020/aaz4797_Ruger_Data_S1.xlsx', sheet = 2 )
load("C:/CODE/gwd_metrics/data/bci_raw_data_2/bci.tree/bci.tree8.rdata")

# Generate generation time estimates --------------------------------------

# Plot and store generation times-by-year estimates
gt_estimate_plot <- function( spp_i, adj_type ){
    
    print( spp_i )
    
    # read up on the cohort of interest
    cohorts       <- readRDS( paste0(output_directory, 
                              "PPA_output_raw_cohort_",
                              adj_type,
                              spp_i,".rds") )
    
    # Species name
    spp           <- spVitals$LCVP_Accepted_Taxon[spp_i]
    
    # Obeseved MaxDiam
    maxdiam_obs   <- spVitals$max_dbh[spp_i] / 10
    
    # Maximum diameter X year
    maxdiam_df    <- cohorts %>%
                        group_by( Year ) %>% 
                        summarise( Max_Diam = quantile(Diameter, 
                                                       probs = 0.95) )
    
    # generation time based on 1/2 of SIMULATED max diameter
    gt_maxdiam_df <- cohorts %>%
                        group_by( Year ) %>% 
                        mutate( Max_Diam = quantile(Diameter, 
                                                    probs = 0.95) ) %>% 
                        ungroup %>% 
                        mutate( Half_Diam = Max_Diam / 2 ) %>% 
                        subset( Diameter > Half_Diam ) %>% 
                        mutate( w = N * (phi*Diameter^theta) ) %>% 
                        group_by( Year ) %>% 
                        summarise( MaxDiam_sim = weighted.mean(Age, 
                                                           w) ) %>% 
                        ungroup
    
    # generation time based on 1/2 max observed diameter
    gt_maxdiam_obs <- cohorts %>%
                        subset( Diameter > (maxdiam_obs/2) ) %>% 
                        mutate( w = N * (phi*Diameter^theta) ) %>% 
                        group_by( Year ) %>% 
                        summarise( MaxDiam_obs = weighted.mean(Age, 
                                                           w) ) %>% 
                        ungroup
    
    # generation time based on canopy individuals
    gt_layer_df <- cohorts %>%
                        subset( Layer == 1 ) %>% 
                        group_by( Year ) %>% 
                        mutate( w = N * (phi*Diameter^theta) ) %>% 
                        group_by( Year ) %>% 
                        summarise( Layer = weighted.mean(Age, 
                                                            w) ) %>% 
                        ungroup
    
    gt_df <- left_join(gt_maxdiam_df, gt_maxdiam_obs) %>% 
                left_join( gt_layer_df) %>% 
                gather( Method, gt, MaxDiam_sim:Layer )
    
    # store plot
    p <- ggplot( gt_df ) +
            geom_point( aes(Year, gt,
                            color = Method) ) +
            labs( y        = 'Generation Time',
                  title    = paste0('Species n. ',spp_i),
                  subtitle = spp ) +
            theme_minimal() +
            scale_color_colorblind() +
            theme( plot.title    = element_text( hjust = 0.5,
                                                 size  = 20 ),
                   plot.subtitle = element_text( hjust = 0.5,
                                                 size  = 17 ),
                   axis.title    = element_text( size = 15),
                   axis.text     = element_text( size = 15),
                   legend.text   = element_text( size = 15),
                   legend.title  = element_text( size = 15) ) 
    ggsave( paste0('results/PPA/plots/species_',adj_type,spp_i,'.tiff'),
            plot = p, width = 6.3, height = 5, compression = 'lzw' )

    p <- ggplot(maxdiam_df) +
            geom_point( aes(Year,Max_Diam) ) +
            theme_minimal() +
            theme( plot.title    = element_text( hjust = 0.5,
                                                size  = 20 ),
                   plot.subtitle = element_text( hjust = 0.5,
                                                 size  = 17 ),
                   axis.title    = element_text( size = 15),
                   axis.text     = element_text( size = 15),
                   legend.text   = element_text( size = 15),
                   legend.title  = element_text( size = 15) ) +
            labs( title    = paste0('Species n. ',spp_i),
                  subtitle = spp )
    ggsave( paste0('results/PPA/plots/maxdiam_species_',adj_type,spp_i,'.tiff'),
            plot = p, width = 6.3, height = 5, compression = 'lzw' )
    
    # store results
    write.csv( mutate( gt_df, spp = spp ), 
               paste0('results/PPA/output/gen_time_spp_',adj_type,spp_i,'.csv'),
               row.names = F )
    
}

# plots
sapply(1:247, gt_estimate_plot, 'ba_')
sapply(1:247, gt_estimate_plot, 'n_')

# store maximum diameter from simulations
calc_max_diam <- function( spp_i ){
    
    # read up on the cohort of interest
    cohorts       <- readRDS( paste0(output_directory, 
                                     "PPA_output_raw_cohort_ba_",
                                     spp_i,".rds") )
    
    # output
    cohorts %>%
        subset( Year == 2000 ) %>% 
        group_by( Year, SpeciesID ) %>% 
        summarise( Max_Diam = quantile(Diameter, 
                                       probs = 0.95) ) %>% 
        ungroup 
    
}

# store maximum diameter
maxdiam_sim_df <- lapply(1:247, calc_max_diam) %>% bind_rows

# store generation times estimates
gt_extract <- function( spp_i, adj_type ){
 
    read.csv( paste0('results/PPA/output/gen_time_spp_',
                     adj_type,spp_i,'.csv') ) %>% 
        subset( Year == 2000 )
    
}

# store generation time
gt_all_df   <- lapply(1:247, gt_extract, 'ba_') %>% bind_rows
gt_all_n_df <- lapply(1:247, gt_extract, 'n_') %>% bind_rows


# Plots -----------------------------------------------------------

# generation time stats
gt_stats <- gt_all_df %>% 
    group_by( Method ) %>% 
    summarise( gt_mean   = mean(gt,na.rm=T),
               gt_median = median(gt,na.rm=T),
               gt_min    = min(gt,na.rm=T),
               gt_max    = max(gt,na.rm=T) ) %>% 
    mutate( abs_range    = gt_max - gt_min ) %>% 
    ungroup

# Density plot for generation time
gt_all_df %>% 
    ggplot() +
    geom_density( aes(x = gt, fill = Method),
                  alpha = 0.5) +
    scale_fill_colorblind() +
    labs( x = 'Generation time estimate',
          y = 'Kernel density estimate' ) +
    theme_minimal() + 
    theme( axis.text.x  = element_text( angle = 90,
                                        vjust = 0,
                                        hjust = 1
                                        ),
           axis.title   = element_text( size = 20),
           legend.text  = element_text( size = 15),
           legend.title = element_text( size = 15),
           axis.text    = element_text( size = 10)
           ) +
    ggsave( 'C:/CODE/gwd_metrics/results/PPA/plots/summary/gt_kernel_densities.tiff',
            width = 6.3, height = 4, compression = 'lzw' )

# Density plot for generation time; two adjustments
bind_rows( mutate( gt_all_df,   adj = 'Basal area'),
           mutate( gt_all_n_df, adj = 'Abundance') ) %>% 
    ggplot() +
    geom_density( aes(x = gt, fill = Method),
                  alpha = 0.5) +
    scale_fill_colorblind() +
    labs( x = 'Generation time estimate',
          y = 'Kernel density estimate' ) +
    facet_wrap( ~adj ) +
    theme_minimal() + 
    theme( axis.text.x  = element_text( angle = 90,
                                        vjust = 0,
                                        hjust = 1),
    axis.title   = element_text( size = 20),
    legend.text  = element_text( size = 15),
    legend.title = element_text( size = 15),
    strip.text   = element_text( size = 15),
    axis.text    = element_text( size = 10)
    ) +
    ggsave( 'C:/CODE/gwd_metrics/results/PPA/plots/summary/gt_baVSn_kernel_densities.tiff',
            width = 6.3, height = 4, compression = 'lzw' )

# Pairs plot among generation time estimates
p <- gt_all_df %>% 
    pivot_wider( names_from  = Method,
                 values_from = gt ) %>% 
    select( -Year, -spp ) %>% 
    setNames( c('MaxDiam (sim.)',
                'MaxDiam (obs.)',
                'Layers') ) %>% 
    ggpairs(  ) +
    theme_minimal() + 
    theme( axis.text.x = element_text( angle = 90,
                                       vjust = 0.5,
                                       hjust = 1,
                                       size  = 10),
           axis.text.y = element_text( size  = 10),
           strip.text  = element_text( size  = 15) )

ggsave( 'C:/CODE/gwd_metrics/results/PPA/plots/summary/gt_est_corrs.tiff',
        plot = p,
        width = 6.3, height = 6.3, compression = 'lzw' )


# plot observed max. diameter vs. simulated max. diameter
maxdiam_sim_df %>% 
    select(-Year) %>% 
    left_join( select(spVitals, SpeciesID, max_dbh) ) %>% 
    rename( max_dbh_obs = max_dbh,
            max_dbh_sim = Max_Diam ) %>% 
    mutate( max_dbh_obs = max_dbh_obs / 10 ) %>% 
    ggplot() +
    geom_point( aes(max_dbh_sim, max_dbh_obs) ) +
    theme_minimal() +
    labs( y = 'Maximum dbh, observed (cm)',
          x = 'Maximum dbh, simulated (cm)' ) +
    theme( axis.title = element_text( size = 20 ),
           axis.text  = element_text( size = 15 ) ) +
    ggsave( 'C:/CODE/gwd_metrics/results/PPA/plots/summary/maxdbh_obsVssim.tiff',
            width = 6.3, height = 6.3, compression = 'lzw' )

# compare observed max diameter vs. max_dbh in spVitals
bci.tree8 %>% 
    subset( !(ba == 0) ) %>% 
    group_by( sp ) %>% 
    summarise( max_dbh_8 = max(dbh) ) %>% 
    ungroup %>% 
    left_join( select(spVitals, sp, max_dbh) ) %>% 
    ggplot( ) +
    geom_point( aes(max_dbh, max_dbh_8) ) +
    theme_minimal()
    
# test: do we share all 247 species?
intersect( bci.tree8$sp %>% unique,
           spVitals$sp %>% unique ) %>% 
    length %>% 
    testthat::expect_equal( 247 )

# format to compare GT adjusted by basal area and abundance
gt_n_ba_df <- bind_rows( mutate( gt_all_df,   
                                 adj = 'Adjusted by basal area'),
                         mutate( gt_all_n_df, 
                                 adj = 'Adjusted by abundance') ) %>% 
                pivot_wider( values_from = gt, 
                             names_from  = adj )

# plot comparison of GT adjusted by basal area and abundance
gt_n_ba_df %>% 
    ggplot() +
    geom_point( aes(x=`Adjusted by basal area`, 
                    y=`Adjusted by abundance`) ) +
    facet_wrap( ~Method ) +
    theme_minimal() +
    theme( axis.text.x = element_text( angle = 90,
                                       vjust = 0.5,
                                       hjust = 1,
                                       size  = 10 ) ) +
    ggsave( 'C:/CODE/gwd_metrics/results/PPA/plots/summary/gt_baVSn_corr.tiff',
            width = 6.3, height = 6.3, compression = 'lzw' )
    
# calculate correlations between GT estimates
gt_n_ba_df %>% 
    select( c("Method", "Adjusted by basal area", "Adjusted by abundance") ) %>% 
    drop_na %>% 
    group_by( Method ) %>% 
    summarise( corr = cor(`Adjusted by basal area`,
                          `Adjusted by abundance`) ) %>% 
    ungroup

# correlations with demographic rates --------------------------------

# demographic correlation data frame
gt_demog_cor_df <- gt_all_df %>% 
    subset( Method == 'MaxDiam_obs' ) %>% 
    select( -Year ) %>% 
    rename( LCVP_Accepted_Taxon = spp ) %>% 
    left_join( spVitals ) %>% 
    select( Method:F ) %>% 
    pivot_longer( G1:F, 
                  names_to  = 'demog_par',
                  values_to = 'demog_par_val') %>% 
    mutate( demog_par_val_log = log(demog_par_val),
            demog_par_log     = paste0('log(',demog_par,')'),
            gt_log            = log(gt) )

# actual correlation values
gt_demog_cor_vals <- gt_demog_cor_df %>% 
    drop_na %>% 
    group_by( demog_par_log ) %>% 
    summarise( gt_demog_corr     = cor( gt,     demog_par_val_log ),
               gt_log_demog_corr = cor( gt_log, demog_par_val_log ) ) %>% 
    ungroup

# GT versus logged demo. rages
p <- gt_demog_cor_df %>% 
        left_join( gt_demog_cor_vals ) %>% 
        mutate( demog_par_log = paste0( demog_par_log,
                                        '; r=', 
                                        round(gt_demog_corr, 2)) ) %>%
        ggplot() +
        geom_point( aes(gt, demog_par_val_log),
                    alpha = 0.5, size = 2 ) +
        facet_wrap( ~ demog_par_log,
                    scale = 'free' ) +
        theme_minimal() +
        labs( x = 'Generation time',
              y = 'Logged demographic rate') +
        theme( axis.title = element_text( size = 15 ),
               strip.text = element_text( size = 12.5 ) )

ggsave( 'C:/CODE/gwd_metrics/results/PPA/plots/summary/gt_vs_logDemogrates.tiff',
        plot = p,
        width = 6.3, height = 6.3, compression = 'lzw' )
    
    
# LOGGED GT versus logged demo. rages
p <- gt_demog_cor_df %>% 
        left_join( gt_demog_cor_vals ) %>% 
        mutate( demog_par_log = paste0( demog_par_log,
                                        '; r=', 
                                        round(gt_log_demog_corr, 2)) ) %>%
        ggplot() +
        geom_point( aes(gt_log, demog_par_val_log),
                    alpha = 0.5, size = 2 ) +
        facet_wrap( ~ demog_par_log,
                    scale = 'free' ) +
        theme_minimal() +
        labs( x = 'Generation time',
              y = 'Logged demographic rate') +
        theme( axis.title = element_text( size = 15 ),
               strip.text = element_text( size = 12.5 ) )

ggsave( 'C:/CODE/gwd_metrics/results/PPA/plots/summary/logGt_vs_logDemogrates.tiff',
        plot = p,
        width = 6.3, height = 6.3, compression = 'lzw' )


# format Rügers' PC values
rugers <- rugers %>% 
        mutate( sp    = tolower(sp),
                Latin = paste0(Genus,' ',Species) ) %>% 
        select( sp, Latin, PC1score, PC2score )

# gt versus PC scores
gt_pc_df <- gt_all_df %>% 
    subset( Method == 'MaxDiam_obs' ) %>% 
    select( -Year ) %>% 
    rename( LCVP_Accepted_Taxon = spp ) %>% 
    left_join( spVitals ) %>% 
    left_join( rugers ) %>% 
    select( Method:SpeciesID, PC1score:PC2score ) %>% 
    pivot_longer( PC1score:PC2score, 
                  names_to  = 'demog_par',
                  values_to = 'demog_par_val') %>% 
    mutate( gt_log    = log(gt),
            demog_par = gsub('score','',demog_par) ) 

# correlations between GT and PC values
gt_pc_stats <- gt_pc_df %>% 
    select( demog_par, gt_log, demog_par_val ) %>% 
    drop_na %>% 
    group_by( demog_par ) %>% 
    summarise( corr = cor(gt_log, demog_par_val) ) %>% 
    ungroup


# plot log_GT versus PC axes 
p <- gt_pc_df %>% 
        left_join( gt_pc_stats ) %>% 
        mutate( demog_par = paste0( demog_par,
                                    '; r=', 
                                    round(corr, 2)) ) %>%

        ggplot() +
        geom_point( aes(gt_log, demog_par_val),
                    alpha = 0.5, size = 2 ) +
        facet_wrap( ~ demog_par,
                    scale = 'free' ) +
        theme_minimal() +
        labs( x = 'Logged generation time',
              y = 'Logged demographic rate') +
        theme( axis.title = element_text( size = 15 ),
               strip.text = element_text( size = 12.5 ) )

ggsave( 'C:/CODE/gwd_metrics/results/PPA/plots/summary/logGt_vs_PCaxes.tiff',
        plot = p,
        width = 6.3, height = 4, compression = 'lzw' )


# Regression model ------------------------------------------

# compare most important demographic rates to GT
gt_demog_mod_df <- gt_demog_cor_df %>% 
        select( SpeciesID, gt_log, demog_par_log, 
                demog_par_val_log) %>% 
        pivot_wider( values_from = demog_par_val_log,
                     names_from  = demog_par_log ) 

# Predict GT with mu1 and G4
mod <- lm( gt_log ~ `log(mu1)` + `log(G4)`, 
           data = gt_demog_mod_df )

mod %>% summary

# Predict GT with PC axes
gt_pc_df %>% 
    select( SpeciesID, gt_log, demog_par, 
            demog_par_val) %>% 
    pivot_wider( values_from = demog_par_val,
                 names_from  = demog_par ) %>% 
    lm( gt_log ~ PC1 + PC2, 
        data = . ) %>% 
    summary
    