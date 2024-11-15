#------------------------Description--------------------------------------------
# This script is used to generate the baseline MSE error of a converged system
# and the corresponding 95% central interval of the sampling distribution.
#
# The second part is the the 95% central interval of the sampling distribution of
# given true ability levels and urn sizes
#-------------------------------------------------------------------------------

################################################################################
# Dependencies
################################################################################
source("extra/fixed_factors.R")
source("util.R")

################################################################################
# MSE baseline
################################################################################
n_samples = 10000
pi_theta = transform_logit(Theta[,1])
#-------------------------------------------------------------------------------
# Simulate 500 responses for each student from the corresponding limiting dist
# repeat this process 10000 times

#-------------------------------------------------------------------------------
# Creating containers for the interval and the mean MSE (the latter is just a check
# whether we get the same from theoretical calcs)

MSE_baseline = matrix(0, nrow = length(student_urn_size), ncol = 5)
colnames(MSE_baseline) = c("urn_size", "mean", "lower", "upper", "mean_theoretical")
MSE_baseline[,1] = student_urn_size

for(us in 1:length(student_urn_size)){
  print(paste0(us, " out of ", length(student_urn_size)))
  repeated_system_container = numeric(length=n_samples)
  for(r_sample in 1:n_samples){
    mock_system_container = numeric(length=n_students)
    for(i in 1:n_students){

      #-------------------------------------------------------------------------
      # Simulate 500 binomials with the given urn sizes and true ability levels
      
      #student_ratings = rbinom(n_games, student_urn_size[us], pi_theta[i]) / student_urn_size[us]
      student_ratings = rbinom(1, student_urn_size[us], pi_theta[i]) / student_urn_size[us]
      mock_system_container[i] = (student_ratings - pi_theta[i])^2
      
      #-------------------------------------------------------------------------
    }
    #---------------------------------------------------------------------------
    # Getting system level MSE for each sample
    repeated_system_container[r_sample] = mean(mock_system_container)
    
    #---------------------------------------------------------------------------
  }
  #-----------------------------------------------------------------------------
  # Getting the mean and the 95% central interval of the sampling distribution
  MSE_baseline[us, 2] = mean(repeated_system_container)
  MSE_baseline[us, 3:4] = quantile(repeated_system_container, c(0.025, 0.975))
  
  #-----------------------------------------------------------------------------
  # Calculating theoretical mean MSE
  # calculate the variance of each person divide by the urn size squared and 
  #average across students
  MSE_baseline[us, 5] = mean(pi_theta * (1-pi_theta) / student_urn_size[us])
}

#-------------------------------------------------------------------------------
# save the resulting baseline values to the output folder

saveRDS(MSE_baseline, "output/MSE_baseline.Rds")

################################################################################
# baseline for urn sizes and true ability levels
################################################################################
#-------------------------------------------------------------------------------
# rescaling the new student's theta values
pi_theta_new = transform_logit(new_student_thetas)

#-------------------------------------------------------------------------------
# Creating containers for the urn sizes, new true values and lower and upper bounds
# of the 95% central interval of the sampling distribution

urnings_baseline = matrix(0, nrow = n_new * length(student_urn_size), ncol = 4)
colnames(urnings_baseline) = c("urn_size", "theta", "lower", "upper")
urnings_baseline[,1] = rep(student_urn_size, times = n_new)
urnings_baseline[,2] = rep(pi_theta_new, each = length(student_urn_size))

#-------------------------------------------------------------------------------
for(i in 1:nrow(urnings_baseline)){
  print(paste0(i, " out of ", nrow(urnings_baseline)))
  repeated_system_container = numeric(length=n_samples)
  for(r_sample in 1:n_samples){
    #---------------------------------------------------------------------------
    # simulate n_reps amount of ratings to mimic the simulation study and average it
    repeated_system_container[r_sample] = mean(rbinom(n_reps, urnings_baseline[i,1], urnings_baseline[i,2]) / urnings_baseline[i,1])
    
    #---------------------------------------------------------------------------
  }
  #-----------------------------------------------------------------------------
  # saving the 2.5 and 97.5 quantiles of the sampling distribution
  urnings_baseline[i, 3:4] = quantile(repeated_system_container, c(0.025, 0.975))
}

#-------------------------------------------------------------------------------
# saving the baseline for urnings

saveRDS(urnings_baseline, "output/urnings_baseline_extra.Rds")
