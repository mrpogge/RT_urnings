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
item_urn_size_curr = item_urn_size[1]

# true item difficulties
delta=qnorm(seq(1/(n_items+1),n_items/(n_items+1),length=n_items),mu_delta,sigma_delta)
Delta=matrix(0, ncol=n_games,nrow=n_items)
for(j in 1:n_games){
  Delta[,j]=delta
}

Theta = create_theta(mu_theta_baseline, sigma_theta, n_games, n_students)

item_starting = create_starting(delta, n_items, item_urn_size_curr, items = TRUE)
rm(delta)

################################################################################
# simulations
################################################################################
container = array(0, dim = c(n_students + 8, n_games, n_reps))
result_list = vector(mode = "list", length = length(student_urn_size))

for(i in 1:length(student_urn_size)){
  print(student_urn_size[i])
  for(r in 1:n_reps){
    student_starting = create_starting(Theta[,1], n_students, student_urn_size[i], items = FALSE, new = 8)
    temp = urnings_simple(student_starting,
                          item_starting,
                          Theta,
                          Delta,
                          n_students + 8,
                          n_items,
                          n_games,
                          student_urn_size[i], #different urn sizes
                          item_urn_size_curr,
                          adaptive = 1,
                          m_adapt = mu_p_baseline,
                          sd_adapt = sigma_p,
                          paired = is_pu,
                          return = "simple",
                          OS = "MAC")
    container[,,r] = temp / student_urn_size[i]
  }
  result_list[[i]]$mean = apply(container, c(1,2), mean)
  result_list[[i]]$var = apply(container, c(1,2), var)
  result_list[[i]]$true = Theta[,1]
}

saveRDS(result_list, "output/AccResultUrn.Rds")