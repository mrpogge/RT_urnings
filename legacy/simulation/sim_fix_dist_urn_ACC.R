################################################################################
# loading C script, wrapper, and settings
################################################################################
dyn.load("CRT.so")
source("crt_wrapper.R")
source("simulation/fixed_factors.R")
################################################################################
# fix values
################################################################################
item_urn_size = item_urn_size[1]

# true item difficulties
delta=qnorm(seq(1/(n_items+1),n_items/(n_items+1),length=n_items),0,1)
Delta=matrix(0, ncol=n_games,nrow=n_items)

for(j in 1:n_games){
  Delta[,j]=delta
}
item_starting = create_starting(delta, n_items, item_urn_size, items = TRUE)
rm(delta)
################################################################################
# creating design for the simulation 
################################################################################
Design = matrix(0, 9, 2)
Design[,1] = rep(c(log(0.6/0.4),log(0.7/0.3),log(0.8/0.2)), times = 3)
Design[,2] = rep(student_urn_size, each = 3)


################################################################################
# simulations
################################################################################
container = array(0, dim = c(n_students + 8, n_games, n_reps))
result_list = vector(mode = "list", length = nrow(Design))
for(i in 1:nrow(Design)){
  print(Design[i,])
  for(r in 1:n_reps){
    #check and create true and starting values
    mu_theta = Design[i,1]
    Theta = create_theta(mu_theta, sigma_theta, n_games, n_students)
    student_starting = create_starting(Theta[,1], n_students, Design[i,2], items = FALSE, new = 8)
    temp = urnings_simple(student_starting,
                               item_starting,
                               Theta,
                               Delta,
                               n_students + 8,
                               n_items,
                               n_games,
                               Design[i,2],
                               item_urn_size,
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

saveRDS(result_list, "output/ACCResultDistXUrn.Rds")

