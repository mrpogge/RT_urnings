################################################################################
# loading C script, wrapper, and settings
################################################################################
dyn.load("CRT.so")
source("crt_wrapper.R")
source("simulation/fixed_factors.R")
################################################################################
# fix values
################################################################################
n_students = 1000 - 40
item_urn_size = item_urn_size[1]

# true item difficulties
delta=qnorm(seq(1/(n_items+1),n_items/(n_items+1),length=n_items),0,1)

Delta=matrix(0, ncol=n_games,nrow=n_items)
for(j in 1:n_games){
  Delta[,j]=delta
}
item_starting = create_starting(delta, n_items, item_urn_size, items = TRUE)

amount_of_change = rep(c(-0.5, 0, 0.5, 1, 2)/n_games, times = 200)
theta  = qnorm(seq(1/(n_students+1),n_students/(n_students+1),by=1/(n_students+1)),mu_theta_baseline,sigma_theta)
theta = sample(theta, n_students)
theta = c(theta, rep(seq(-1,1, length.out = 8), each = 5))
n_students_add = n_students + 40
Theta=matrix(0, ncol=n_games,nrow=n_students_add)
Theta[,1] = theta + amount_of_change
for(j in 2:n_games){
  Theta[,j]=Theta[,j-1] + amount_of_change
}


################################################################################
# creating design for the simulation 
################################################################################
Design = matrix(0, 3, 1)
Design[,1] = student_urn_size

################################################################################
# simulations
################################################################################
container = array(0, dim = c(n_students_add, n_games, n_reps))
result_list = vector(mode = "list", length = nrow(Design))
for(i in 1:nrow(Design)){
  print(Design[i,])
  for(r in 1:n_reps){
    true_probs = 1/(1+exp(-theta))
    student_starting=rbinom(n_students_add,Design[i,1],true_probs)
    temp = urnings_combined_RT(student_starting,
                          item_starting,
                          Theta,
                          Delta,
                          n_students_add,
                          n_items,
                          n_games,
                          Design[i,1],
                          item_urn_size,
                          adaptive = 1,
                          m_adapt = mu_p_baseline,
                          sd_adapt = sigma_p,
                          paired = 1,
                          return = "simple",
                          OS = "MAC")
    container[,,r] = temp / Design[i,1]
  }
  result_list[[i]]$mean = apply(container, c(1,2), mean)
  result_list[[i]]$var = apply(container, c(1,2), var)
  result_list[[i]]$true = Theta
}

saveRDS(result_list, "output/MSChange.Rds")
