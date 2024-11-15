################################################################################
# loading C script, wrapper, and settings
################################################################################
dyn.load("CRT.so")
source("crt_wrapper.R")

################################################################################
# fix values
################################################################################
n_students = 1000 - 8
n_items = 200
n_games = 500
n_reps = 200
mu_theta = 0
sd_theta = 1
sd_P = 0.5
weight = 4
item_urn_size = 200

# true item difficulties
delta=qnorm(seq(1/(n_items+1),n_items/(n_items+1),length=n_items),0,1)
Delta=matrix(0, ncol=n_games,nrow=n_items)
for(j in 1:n_games){
  Delta[,j]=delta
}

Theta = create_theta(mu_theta, sd_theta, n_games, n_students)

item_starting = create_starting(delta, n_items, item_urn_size, items = TRUE)
rm(delta)
################################################################################
# creating design for the simulation 
################################################################################
Design = matrix(0, 18, 3)
Design[,1] = rep(1:3, each = 6)
Design[,2] = rep(c(40,40,80,80,160,160), times = 3)
Design[,3] = rep(c(0.5,0.7), times = 9)

Design = Design[1:6,]

################################################################################
# simulations
################################################################################
container = array(0, dim = c(n_students + 8, n_games, n_reps))
result_list = vector(mode = "list", length = nrow(Design))

for(i in 1:nrow(Design)){
  print(Design[i,])
  for(r in 1:n_reps){
    student_urn_size = Design[i,2]
    student_starting = create_starting(Theta[,1], n_students, student_urn_size, items = FALSE, new = 8)
    temp = urnings_simple(student_starting,
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

saveRDS(result_list, "output/AccResultUrnXAdapt.Rds")