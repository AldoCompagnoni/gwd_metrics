# wrapper file 


# Michael Crawford's code ----------------------------

# Use raw BCI data to calculate relative proportions of single species
# to produce an adjustment factor for fecundity
source('R/PPA/format/format_bci_raw_data.R')

# Run year-specific simulations of size and age structure, one species at a time
source('R/PPA/runscript_age.R')

# calculate year-specific generation times, maximum diameter, and store the results
source('R/PPA/postprocessing_gt.R')


# # Farrior's code -------------------------------------
# 
# # Farrior's stochastic simulation that tracks age
# source('R/farrior/ppa_workscript_age.r')
# 
# 
# # Rugers ---------------------------------------------
# 
# source('R/rugers/PPA_5PFTs_4layers_age.R')
# source('R/rugers/plot_data.R')
