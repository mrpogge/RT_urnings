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
item_urn_size = item_urn_size[1]

# true item difficulties
delta=qnorm(seq(1/(n_items+1),n_items/(n_items+1),length=n_items),mu_delta,sigma_delta)
Delta=matrix(0, ncol=n_games,nrow=n_items)
for(j in 1:n_games){
  Delta[,j]=delta
}

item_starting = create_starting(delta, n_items, item_urn_size, items = TRUE)
rm(delta)
################################################################################
# creating design for the simulation 
################################################################################
Design = matrix(0, 12, 4)
Design[,1] = rep(1:3, each = 4)
Design[,2] = rep(student_urn_size, times = 4)
Design[,3] = rep(mu_p, times = 3)
Design[,4] = rep(mu_theta, times = 3) 
################################################################################
# simulations
################################################################################
container = array(0, dim = c(n_students + 8, n_games, n_reps))
result_list = vector(mode = "list", length = nrow(Design))

for(i in 1:nrow(Design)){
  print(Design[i,])
  for(r in 1:n_reps){
    Theta = create_theta(Design[i,4], sigma_theta, n_games, n_students)
    student_urn_size = Design[i,2]
    student_starting = create_starting(Theta[,1], n_students, student_urn_size, items = FALSE, new = 8)
    temp = urnings_combined_RT(student_starting,
                          item_starting,
                          Theta,
                          Delta,
                          n_students + 8,
                          n_items,
                          n_games,
                          student_urn_size, #different urn sizes
                          item_urn_size,
                          adaptive = 1,
                          m_adapt = Design[i,3], #different adaptivity levels
                          sd_adapt = sigma_p,
                          paired = 1,
                          return = "simple",
                          OS = "MAC")
    container[,,r] = temp / student_urn_size
  }
  result_list[[i]]$mean = apply(container, c(1,2), mean)
  result_list[[i]]$var = apply(container, c(1,2), var)
  result_list[[i]]$true = Theta[,1]
}

saveRDS(result_list, "output/MUpResultUrnXAdapt.Rds")