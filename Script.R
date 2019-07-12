
#####################First, load the packages.


library(ape)
library(phytools)
library(BioGeoBEARS)
library(parallel)
library(minqa)
library(optimx)

#NOTE: The BioGeoBEARS from GitHub must be installed (above), version 1.1.1


#####################
#######################################################
# Stochastic mapping on DEC with matrix of distances
#######################################################
load("CyperaceaeMadagascar_DEC_M0_constrained_x.Rdata")
res = res
#######################################################

clado_events_tables = NULL
ana_events_tables = NULL
lnum = 0

#######################################################

BSM_inputs_fn = "BSM_inputs_file.Rdata"
runInputsSlow = TRUE
if (runInputsSlow)
{
  stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res=res)
  save(stochastic_mapping_inputs_list, file=BSM_inputs_fn)
} else {
  # Loads to "stochastic_mapping_inputs_list"
  load(BSM_inputs_fn)
} # END if (runInputsSlow)

# Check inputs (doesn't work the same on unconstr)
names(stochastic_mapping_inputs_list)
stochastic_mapping_inputs_list$phy2
stochastic_mapping_inputs_list$COO_weights_columnar
stochastic_mapping_inputs_list$unconstr
set.seed(seed=as.numeric(Sys.time()))

runBSMslow = TRUE
if (runBSMslow == TRUE)
{
  # Saves to: RES_clado_events_tables.Rdata
  # Saves to: RES_ana_events_tables.Rdata
  BSM_output = runBSM(res, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, maxnum_maps_to_try=100, nummaps_goal=50, maxtries_per_branch=40000, save_after_every_try=TRUE, savedir=getwd(), seedval=12345, wait_before_save=0.01)
  
  RES_clado_events_tables = BSM_output$RES_clado_events_tables
  RES_ana_events_tables = BSM_output$RES_ana_events_tables
} else {
  # Load previously saved...
  
  # Loads to: RES_clado_events_tables
  load(file="RES_clado_events_tables.Rdata")
  # Loads to: RES_ana_events_tables
  load(file="RES_ana_events_tables.Rdata")
  BSM_output = NULL
  BSM_output$RES_clado_events_tables = RES_clado_events_tables
  BSM_output$RES_ana_events_tables = RES_ana_events_tables
} # END if (runBSMslow == TRUE)

save(BSM_output, file = "BSM_output.RData")