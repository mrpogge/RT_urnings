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
source("analysis/fixed_factors_change.R")
source("util.R")

################################################################################
# MSE baseline
################################################################################
n_samples = 1000
pi_theta = transform_logit(Theta)

#-------------------------------------------------------------------------------
# Simulate 500 responses for each student from the corresponding limiting dist
# repeat this process 10000 times

#-------------------------------------------------------------------------------
# Creating containers for the interval and the mean MSE (the latter is just a check
# whether we get the same from theoretical calcs)

MSE_baseline = matrix(0, nrow = length(student_urn_size) * n_games, ncol = 5)
colnames(MSE_baseline) = c("iter","urn_size", "mean", "lower", "upper")
MSE_baseline[,1] = rep(1:n_games, times = length(student_urn_size))
MSE_baseline[,2] = rep(student_urn_size, each = n_games)

for(us in 1:length(student_urn_size)){
  print(paste0(us, " out of ", length(student_urn_size)))
  repeated_system_container = matrix(0, nrow=n_samples, ncol = 501)
  for(r_sample in 1:n_samples){
    mock_system_container = matrix(0, nrow = n_students, ncol = 501)
    for(i in 1:n_students){
      
      #-------------------------------------------------------------------------
      # Simulate 500 binomials with the given urn sizes and true ability levels
      
      #student_ratings = rbinom(n_games, student_urn_size[us], pi_theta[i]) / student_urn_size[us]
      student_ratings = rbinom(501, student_urn_size[us], pi_theta[i,]) / student_urn_size[us]
      mock_system_container[i,] = (student_ratings - pi_theta[i,])^2
      
      #-------------------------------------------------------------------------
    }
    #---------------------------------------------------------------------------
    # Getting system level MSE for each sample
    repeated_system_container[r_sample, ] = colMeans(mock_system_container)
    
    #---------------------------------------------------------------------------
  }
  #-----------------------------------------------------------------------------
  # Getting the mean and the 95% central interval of the sampling distribution
  MSE_baseline[MSE_baseline[,2] == student_urn_size[us], 3] = colMeans(repeated_system_container)
  #calculating quantiles for each iteration
  MSE_baseline[MSE_baseline[,2] == student_urn_size[us], 4:5] = t(apply(repeated_system_container, 2, quantile, probs = c(0.025, 0.975)))
}

saveRDS(MSE_baseline, "output/MSE_baseline_change.Rds")