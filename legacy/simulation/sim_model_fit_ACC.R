set.seed(123457)
################################################################################
# loading C script, wrapper, and settings
################################################################################
dyn.load("CRT.so")
source("crt_wrapper.R")

################################################################################
# fix values
################################################################################
n_students = 1000
n_items = 200
n_games = 2000
n_reps = 1
mu_theta = 0
sd_theta = 1
sd_P = 0.5
weight = 4
item_urn_size = 200
student_urn_size = 80

# true item difficulties
delta=qnorm(seq(1/(n_items+1),n_items/(n_items+1),length=n_items),0,1)
theta = rnorm(n_students, 0,1)

Delta=matrix(0, ncol=n_games,nrow=n_items)
Theta=matrix(0, ncol=n_games,nrow=n_students)
for(j in 1:n_games){
  Delta[,j]=delta
  Theta[,j]=theta
}
item_starting = create_starting(delta, n_items, item_urn_size, items = TRUE)
true_probs = 1/(1+exp(-theta))
student_starting=rbinom(n_students,student_urn_size,true_probs)

################################################################################
# creating design for the simulation 
################################################################################
Design = matrix(0, 15, 2)
Design[,1] = rep(1:3, each = 5)
Design[,2] = rep(1:5, times = 3)

Design = Design[1:5, ]

################################################################################
# simulations
################################################################################
result_list = vector(mode = "list", length = nrow(Design))

for(i in 1:nrow(Design)){
    print(Design[i,])
    if(Design[i,2] > 2){
      rho = switch(as.character(Design[i,2]),
                   "3" = 0.5,
                   "4" = -0.5,
                   "5" = 0)
      params = MASS::mvrnorm(n_students, c(0,0), Sigma = matrix(c(1,rho,rho,1), 2,2))
      theta = params[,1]
      theta_speed = params[,2]
      Theta=matrix(0, ncol=n_games,nrow=n_students)
      for(j in 1:n_games){
        Theta[,j]=theta
      }
      true_probs = 1/(1+exp(-theta))
      student_starting=rbinom(n_students,student_urn_size,true_probs)
    }
    ordering = matrix(sample(0:(n_items-1), n_students*n_games, replace = TRUE), n_students, n_games)
    #generate response array
    misfit = switch(as.character(Design[i,2]),
                    "1" = sim_observed_SRT(theta, delta, n_games, ordering)[,,1],
                    "2" = sim_observed_CISRT(theta, delta, n_games, ordering)[,,1],
                    "3" = sim_observed_LSRT(theta, theta_speed, delta, n_games, ordering)[,,1],
                    "4" = sim_observed_LSRT(theta, theta_speed, delta, n_games, ordering)[,,1],
                    "5" = sim_observed_LSRT(theta, theta_speed, delta, n_games, ordering)[,,1])
    
    temp = urnings_simple_misfit(misfit,
                                 ordering,
                          student_starting,
                          item_starting,
                          Theta,
                          Delta,
                          n_students,
                          n_items,
                          n_games,
                          student_urn_size, #different urn sizes
                          item_urn_size,
                          adaptive = 1,
                          m_adapt = 0.5, 
                          sd_adapt = sd_P,
                          paired = 1,
                          return = "simple",
                          OS = "MAC")
    
  result_list[[i]]$est = temp / student_urn_size
  if(Design[i,2] < 3){
    result_list[[i]]$true = Theta[,1]
  } else {
    result_list[[i]]$true = matrix(c(Theta[,1], theta_speed), 2, n_students, byrow = TRUE)
  }
}
saveRDS(result_list, "output/AccModelFit.Rds")
dyn.unload("CRT.so")

