#-------------------------------Description----------------------------------
#Simulation how fast using the original Urnings algorithm  based system can converge
#from a total cold start scenario to an equilibrium
#-------------------------------------------------------------------------------

################################################################################
# loading C script, wrapper, settings and utils
################################################################################
dyn.load("CRT.so")
source("crt_wrapper.R")
source("extra/fixed_factors.R")
source("util.R")


urnings_baseline = readRDS("output/urnings_baseline_extra.rds")

################################################################################
# Design
################################################################################
Design = matrix(0, nrow = 3, ncol = 2)
Design[,1] = 8 # number of new student in case we want change this
Design[,2] = student_urn_size


################################################################################
# Simulation
################################################################################
container = array(0, dim = c(n_new, n_games, n_reps))
result_list = vector(mode = "list", length = nrow(Design))

for(i in 1:nrow(Design)){
  print(paste0("We are at case no. ", i))
  for(r in 1:n_reps){
    #---------------------------------------------------------------------------
    #setting the simulation seed 
    set.seed(seeds_for_rep[r])
    
    #---------------------------------------------------------------------------
    #initialise starting values 
    
    #student: start from their invariant, except the new ones
    student_starting = create_starting(Theta[,1], n_students, Design[i,2], items = FALSE, new = n_new, cold_value = cold_value)
    
    #---------------------------------------------------------------------------
    # running the system 
    temp = urnings_simple(student_starting,
                          item_starting,
                          Theta_new,
                          Delta,
                          n_students + n_new,
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
    container[,,r] = temp[1001:1008, ] / Design[i,2] #saving rating instead of the urnings value
    
    #---------------------------------------------------------------------------
  }
  #---------------------------------------------------------------------------
  # assembling return objec
  urnings_temp = apply(container, c(1,2), mean)
  result_list[[i]]$mean = urnings_temp
  ht_v = numeric(n_new)
  for(nw in 1:n_new){
    ht_v[nw] = calculate_ht(urnings_temp[nw,], 
                            urnings_baseline[urnings_baseline[,1] == Design[i,2], "lower"][nw],
                            urnings_baseline[urnings_baseline[,1] == Design[i,2], "upper"][nw])
  }
  result_list[[i]]$ht = ht_v
}

#-------------------------------------------------------------------------------
# saving the results

saveRDS(result_list, "output/sim2_newstudents_acc_extra.rds")


plot(result_list[[1]]$mean[1,], type = "l", col = "red", ylim = c(0,1))
lines(rep(transform_logit(new_student_thetas[1]), times = n_games), lty = "dotted")
for(i in 2:n_new){
  lines(result_list[[1]]$mean[i,], col = "red")
  lines(rep(transform_logit(new_student_thetas[i]), times = n_games), lty = "dotted")
}