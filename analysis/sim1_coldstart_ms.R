#-------------------------------Description----------------------------------
#Simulation how fast using the multiple update Urnings algorithm  based system can converge
#from a total cold start scenario to an equilibrium
#-------------------------------------------------------------------------------

################################################################################
# loading C script, wrapper, settings and utils
################################################################################
dyn.load("CRT.so")
source("crt_wrapper.R")
source("analysis/fixed_factors.R")
source("util.R")

MSE_baseline = readRDS("output/MSE_baseline.rds")

################################################################################
# Design
################################################################################
Design = matrix(0, nrow = 3, ncol = 2)
Design[,1] = 999 # total cold start indicator in case we want to have other cold start levels
Design[,2] = student_urn_size


################################################################################
# Simulation
################################################################################
#container = array(0, dim = c(n_students, n_games, n_reps))
mse_container = matrix(0, nrow = n_reps, ncol = n_games)
result_list = vector(mode = "list", length = nrow(Design))

for(i in 1:nrow(Design)){
  print(paste0("We are at case no. ", i))
  for(r in 1:n_reps){
    #---------------------------------------------------------------------------
    #setting the simulation seed 
    set.seed(seeds_for_rep[r])
    
    #---------------------------------------------------------------------------
    #initialise starting values 
    
    #student: 0.5 urnings value (0 logit)
    student_starting = rep(Design[i,2]/2, times = n_students)
    
    # item: 0.5 urnings value (0 logit)
    item_starting = rep(item_urn_size/2, times = n_items)
    
    #---------------------------------------------------------------------------
    # running the system 
    temp = urnings_combined_RT(student_starting,
                              item_starting,
                              Theta,
                              Delta,
                              n_students,
                              n_items,
                              n_games,
                              Design[i,2], # urn size
                              item_urn_size,
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
    #container[,,r] = temp / Design[i,2] #saving rating instead of the urnings value
    mse_container[r,] = colMeans((temp / Design[i,2] - transform_logit(Theta))^2)
    #---------------------------------------------------------------------------
  }
  #---------------------------------------------------------------------------
  # assembling return object
  #result_list[[i]]$mean = apply(container, c(1,2), mean)
  HT_i = apply(mse_container, MARGIN = 1, calculate_ht_MSE, MSE_baseline[i, "upper"]) 
  result_list[[i]]$MSE = colMeans(mse_container) #mse each repeated system 
  result_list[[i]]$HT = mean(HT_i) #HT for the whole system
  result_list[[i]]$HT_sd = quantile(HT_i, c(0.025,0.975))
}

#-------------------------------------------------------------------------------
# saving the results

saveRDS(result_list, "output/sim1_coldstart_ms.rds")

