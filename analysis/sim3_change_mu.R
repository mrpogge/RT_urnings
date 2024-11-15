#-------------------------------Description----------------------------------
#Simulation how well using the multiple urn Urnings algorithm  based system can track
#linear changes on the logit scale
#--------------------------------------------
################################################################################
# loading C script, wrapper, settings and utils
################################################################################
dyn.load("CRT.so")
library(tidyverse)
source("crt_wrapper.R")
source("analysis/fixed_factors_change.R")
source("util.R")

MSE_baseline = readRDS("output/MSE_baseline_change.Rds")
################################################################################
# Design
################################################################################
Design = matrix(0, nrow = length(student_urn_size), ncol = 4)
Design[,1] = 999 # if we would like to change something else
Design[,2:4] = student_urn_size_mu

################################################################################
# Simulation
################################################################################
#container = array(0, dim = c(n_students, n_games, n_reps))
mse_container = matrix(0, nrow = n_reps, ncol = n_games)
result_list = vector(mode = "list", length = nrow(Design))
mse_change_us_container = matrix(0, nrow = length(change_amounts) * length(student_urn_size), ncol = n_reps + 2)
mse_change_us_container[,1] = rep(change_amounts, times = length(student_urn_size))
mse_change_us_container[,2] = rep(student_urn_size, each = length(change_amounts))

for(i in 1:nrow(Design)){
  print(paste0("We are at case no. ", i))
  for(r in 1:n_reps){
    #---------------------------------------------------------------------------
    #setting the simulation seed 
    set.seed(seeds_for_rep[r])
    
    #---------------------------------------------------------------------------
    #initialise starting values 
    
    #student: 0.5 urnings value (0 logit)
    student_starting_MU = matrix(c(create_starting(Theta[,1], n_students, Design[i,2], items = FALSE, new = 0),
                                create_starting(Theta[,1], n_students, Design[i,3], items = FALSE, new = 0),
                                create_starting(Theta[,1], n_students, Design[i,4], items = FALSE, new = 0)), 
                              nrow = 3, ncol = n_students, byrow = TRUE)
      
    
    #---------------------------------------------------------------------------
    # running the system 
    temp = urnings_separate_RT(student_starting_MU,
                               item_starting_MU,
                               Theta,
                               Delta,
                               n_students,
                               n_items,
                               n_games,
                               Design[i,2:4], # urn size
                               item_urn_size_MU,
                               adaptive = 1, # adaptivity on
                               m_adapt = mu_p_baseline,
                               sd_adapt = sigma_p,
                               paired = 1, # paired update on
                               return = "simple", # just the necessary outcomes
                               OS = "MAC") # for other OS you need to compile the C 
    # code accordingly and set this to 
    # "WINDOWS" or "LINUX"      
    #---------------------------------------------------------------------------
    # storing the results
    combined_urnings = (temp[[1]] + temp[[2]] + temp[[3]]) / sum(Design[i,2:4])
    mse_s = (combined_urnings - transform_logit(Theta))^2
    mse_container[r,] = colMeans(mse_s)
    #---------------------------------------------------------------------------
    mse_change = cbind(change_assignments, rowMeans(mse_s[,302:501]))
    colnames(mse_change)[2] = "MSE"
    mse_change_mean = mse_change %>%
      as.data.frame() %>%
      group_by(change_assignments) %>%
      summarise(mMSE = mean(MSE)) %>%
      as.matrix()
    mse_change_us_container[mse_change_us_container[,2] == sum(Design[i,2:4]), r + 2] = mse_change_mean[,2]
    
  }
  #---------------------------------------------------------------------------
  # assembling return object
  result_list[[i]]$MSE = colMeans(mse_container) #mse each repeated system
  result_list[[i]]$coverage = coverage(result_list[[i]]$MSE, MSE_baseline[MSE_baseline[,2] == sum(Design[i,2:4]), 4], 
                                       MSE_baseline[MSE_baseline[,2] == sum(Design[i,2:4]), 5])
  
}
mse_change_res = cbind(mse_change_us_container[,1:2], rowMeans(mse_change_us_container[,3:ncol(mse_change_us_container)]))

#-------------------------------------------------------------------------------
# saving the results
saveRDS(mse_change_res, "output/sim3_change_mu_mmse.rds")
saveRDS(result_list, "output/sim3_change_mu.rds")
