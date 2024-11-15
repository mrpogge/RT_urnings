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
Design[,1] = rep(999, times = 3)
Design[,2] = student_urn_size
################################################################################
# simulations
################################################################################
container = array(0, dim = c(n_students, n_games, n_reps))
result_list = vector(mode = "list", length = nrow(Design))
for(i in 1:nrow(Design)){
  print(Design[i,])
  for(r in 1:n_reps){
    if(Design[i,1] == 100){
      item_starting = create_starting(delta, n_items, item_urn_size[1], items = TRUE)
    } else {
      item_starting = rep(item_urn_size[1]/2, times = n_items)
    }
    student_starting = rep(Design[i,2]/2, times = n_students)
    #check and create true and starting values
    temp = urnings_simple(student_starting,
                          item_starting,
                          Theta,
                          Delta,
                          n_students,
                          n_items,
                          n_games,
                          Design[i,2],
                          item_urn_size[1],
                          adaptive = 1,
                          m_adapt = mu_p_baseline,
                          sd_adapt = sigma_p,
                          paired = 1,
                          return = "simple",
                          OS = "MAC")
    container[,,r] = temp / Design[i,2]
  }
  result_list[[i]]$mean = apply(container, c(1,2), mean)
  result_list[[i]]$var = apply(container, c(1,2), var)
  result_list[[i]]$true = Theta[,1]
}

saveRDS(result_list, "output/AccResultCsXUrn_2.Rds")
