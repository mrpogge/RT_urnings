################################################################################
# loading C script, wrapper, and settings
################################################################################
dyn.load("CRT.so")
source("crt_wrapper.R")
source("simulation/fixed_factors.R")

################################################################################
# fix values
################################################################################
#global settings
n_students = 1000
item_urn_size = item_urn_size_mu
# true item difficulties
delta=qnorm(seq(1/(n_items+1),n_items/(n_items+1),length=n_items),0,1)
Delta=matrix(0, ncol=n_games,nrow=n_items)
for(j in 1:n_games){
  Delta[,j]=delta
}

Theta = create_theta(mu_theta_baseline, sigma_theta, n_games, n_students, new = 0)

################################################################################
# creating design for the simulation 
################################################################################
Design = matrix(0, 3, 2)
Design[,1] = 999
Design[,2] = student_urn_size_mu[,1]
################################################################################
# simulations
################################################################################
container_1 = array(0, dim = c(n_students, n_games, n_reps))
container_2 = array(0, dim = c(n_students, n_games, n_reps))
container_3 = array(0, dim = c(n_students, n_games, n_reps))
result_list = vector(mode = "list", length = nrow(Design))
for(i in 1:nrow(Design)){
  print(Design[i,])
  for(r in 1:n_reps){
    if(Design[i,1] == 100){
      item_starting = matrix(c(create_starting(delta, n_items, item_urn_size[1], items = TRUE, new = 0),
                               create_starting(delta, n_items, item_urn_size[2], items = TRUE, new = 0),
                               create_starting(delta, n_items, item_urn_size[3], items = TRUE, new = 0)),
                             3, n_items, byrow = "TRUE")
    } else {
      item_starting = matrix(c(rep(item_urn_size[1]/2, times = n_items), 
                               rep(item_urn_size[2]/2, times = n_items),
                               rep(item_urn_size[3]/2, times = n_items)),
                             3, n_items, byrow = "TRUE")
    }
    student_uss = c(Design[i,2], Design[i,2]/2, Design[i,2]/4)
    student_starting = matrix(c(rep(student_uss[1]/2, times = n_students), 
                                rep(student_uss[2]/2, times = n_students),
                                rep(student_uss[3]/2, times = n_students)),
                              3, n_students, byrow = "TRUE")
    #check and create true and starting values
    temp = urnings_separate_RT(student_starting,
                               item_starting,
                               Theta,
                               Delta,
                               n_students,
                               n_items,
                               n_games,
                               student_uss,
                               item_urn_size,
                               adaptive = 1,
                               m_adapt = mu_p_baseline,
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

saveRDS(result_list, "output/MUResultCsXUrn_2.Rds")
