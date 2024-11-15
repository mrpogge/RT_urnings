################################################################################
# loading C script, wrapper, and settings
################################################################################
dyn.load("CRT.so")
source("crt_wrapper.R")
source("simulation/fixed_factors.R")

################################################################################
# fix values
################################################################################
n_students = n_students-n_new
item_urn_size = item_urn_size_mu

# true item difficulties
delta=qnorm(seq(1/(n_items+1),n_items/(n_items+1),length=n_items),mu_delta,sigma_delta)
Delta=matrix(0, ncol=n_games,nrow=n_items)
for(j in 1:n_games){
  Delta[,j]=delta
}


item_starting = matrix(c(create_starting(delta, n_items, item_urn_size[1], items = TRUE),
                         create_starting(delta, n_items, item_urn_size[2], items = TRUE),
                         create_starting(delta, n_items, item_urn_size[3], items = TRUE)),
                       3,n_items, byrow = TRUE)
rm(delta)
################################################################################
# creating design for the simulation 
################################################################################
Design = matrix(0, 12, 4)
Design[,1] = rep(1:3, each = 4)
Design[,2] = rep(student_urn_size_mu[,1], times = 4)
Design[,3] = rep(mu_p, times = 3)
Design[,4] = rep(mu_theta, times = 3) 

################################################################################
# simulations
################################################################################
container_1 = array(0, dim = c(n_students + 8, n_games, n_reps))
container_2 = array(0, dim = c(n_students + 8, n_games, n_reps))
container_3 = array(0, dim = c(n_students + 8, n_games, n_reps))
result_list = vector(mode = "list", length = nrow(Design))

for(i in 1:nrow(Design)){
  print(Design[i,])
  for(r in 1:n_reps){
    Theta = create_theta(Design[i,4], sigma_theta, n_games, n_students)
    student_uss = c(Design[i,2], Design[i,2]/2, Design[i,2]/4)
    student_starting = matrix(c(create_starting(Theta[,1], n_students, student_uss[1], items = FALSE, new = 8),
                                create_starting(Theta[,1], n_students, student_uss[2], items = FALSE, new = 8),
                                create_starting(Theta[,1], n_students, student_uss[3], items = FALSE, new = 8)),
                              3, n_students+8, byrow = "TRUE")
    temp = urnings_separate_RT(student_starting,
                               item_starting,
                               Theta,
                               Delta,
                               n_students + 8,
                               n_items,
                               n_games,
                               student_uss, #different urn sizes
                               item_urn_size,
                               adaptive = 1,
                               m_adapt = Design[i,3], #different adaptivity levels
                               sd_adapt = sigma_p,
                               paired = 1,
                               return = "simple",
                               OS = "MAC")
    container_1[,,r] = temp[[1]] 
    container_2[,,r] = temp[[2]] 
    container_3[,,r] = temp[[3]] 
    
  }
  result_list[[i]]$mean_1 = apply(container_1, c(1,2), mean)
  result_list[[i]]$var_1 = apply(container_1, c(1,2), var)
  result_list[[i]]$mean_2 = apply(container_2, c(1,2), mean)
  result_list[[i]]$var_2 = apply(container_2, c(1,2), var)
  result_list[[i]]$mean_3 = apply(container_3, c(1,2), mean)
  result_list[[i]]$var_3 = apply(container_3, c(1,2), var)
  result_list[[i]]$true = Theta[,1]
}

saveRDS(result_list, "output/MUResultUrnXAdapt.Rds")