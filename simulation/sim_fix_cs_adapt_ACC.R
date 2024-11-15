################################################################################
# loading C script, wrapper, and settings
################################################################################
dyn.load("CRT.so")
source("crt_wrapper.R")

################################################################################
# fix values
################################################################################
#global settings
n_students = 1000
n_items = 200
n_games = 500
n_reps = 200
mu_theta = 0
sd_theta = 1
sd_P = 0.5
weight = 4
item_urn_size = 200
student_urn_size = 80

# true item difficulties
delta=qnorm(seq(1/(n_items+1),n_items/(n_items+1),length=n_items),0,1)
Delta=matrix(0, ncol=n_games,nrow=n_items)
for(j in 1:n_games){
  Delta[,j]=delta
}
item_starting = rep(item_urn_size/2, times = n_items)

Theta = create_theta(mu_theta, sd_theta, n_games, n_students)
student_starting = rep(student_urn_size/2, times = n_students)
################################################################################
# creating design for the simulation 
################################################################################
Design = matrix(0, 12, 3)
Design[,1] = rep(1:3, each = 4)
Design[,2] = rep(c(100,100,999,999), times = 3)
Design[,3] = rep(c(0.5,0.7), times = 6)

Design = Design[1:4, ]
################################################################################
# simulations
################################################################################
container = array(0, dim = c(n_students, n_games, n_reps))
result_list = vector(mode = "list", length = nrow(Design))
for(i in 1:nrow(Design)){
  print(Design[i,])
  for(r in 1:n_reps){
    if(Design[i,2] == 100){
      item_starting = create_starting(delta, n_items, item_urn_size, items = TRUE)
    }
    #check and create true and starting values
    temp = urnings_simple(student_starting,
                          item_starting,
                          Theta,
                          Delta,
                          n_students,
                          n_items,
                          n_games,
                          student_urn_size,
                          item_urn_size,
                          adaptive = 1,
                          m_adapt = Design[i,3],
                          sd_adapt = sd_P,
                          paired = 1,
                          return = "simple",
                          OS = "MAC")
    container[,,r] = temp / student_urn_size
  }
  result_list[[i]]$mean = apply(container, c(1,2), mean)
  result_list[[i]]$var = apply(container, c(1,2), var)
  result_list[[i]]$true = Theta[,1]
}

saveRDS(result_list, "output/AccResultCsXAdapt.Rds")
