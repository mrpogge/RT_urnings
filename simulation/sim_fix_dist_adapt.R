################################################################################
# loading C script, wrapper, and settings
################################################################################
dyn.load("CRT.so")
source("crt_wrapper.R")
source("simulation/base_settings.R")

################################################################################
# fix values
################################################################################
n_students = 1000
n_items = 200
n_games = 500
sd_theta = 1
sd_P = 0.5

################################################################################
# creating design for the simulation 
################################################################################
Design = matrix(0, 12, 3)
Design[,1] = rep(1:3, each = 4)
Design[,2] = rep(c(0,0,-1,-1), times = 3)
Design[,3] = rep(c(0.5,0.7), times = 6)

