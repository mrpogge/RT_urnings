################################################################################
# dependencies to create true values and starting values
################################################################################
set.seed(18131219)
n_reps = 1000
n_reps_p = 10000
seeds_for_rep = sample(18131221:19131237, n_reps, replace = FALSE)
seeds_for_rep_p = sample(18131221:19131237, n_reps_p, replace = FALSE)
source("crt_wrapper.R")

################################################################################
# agents and systems settings
################################################################################
n_students = 1000
n_new = 8
n_items = 200
n_games = 500

################################################################################
# Urnings fix params
################################################################################
item_urn_size_MU = c(200,100,50) # urn sizes for the multiple urns algorithm
item_urn_size = sum(item_urn_size_MU) # urn sizes for the rest

#parameters of the normal kernel method
mu_p_baseline = 0.7
sigma_p = 0.5

is_pu = 1 #paired update always on

weight = 4

################################################################################
# latent variable params
################################################################################
#student ability hyperparams
mu_theta = log(0.7/0.3)
sigma_theta = 1/4 #the original param is 1 but Urnings needs rescaling with the weight see above

new_student_mean = 0
new_student_thetas = (new_student_mean + seq(-2,2,length=n_new)) / 4

#true student abilities
Theta = create_theta(mu_theta, sigma_theta, n_games, n_students, new = 0, seed = 18131219)
Theta_new = create_theta(mu_theta, sigma_theta, n_games, n_students, new = 8, seed = 18131219, new_mean = new_student_mean)
#items
mu_delta = 0
sigma_delta = 2/4 #same as before

# true item difficulties
delta=qnorm(seq(1/(n_items+1),n_items/(n_items+1),length=n_items),mu_delta,sigma_delta)
Delta=matrix(0, ncol=n_games,nrow=n_items)
for(j in 1:n_games){
  Delta[,j]=delta
}

# starting values for items
item_starting = create_starting(delta, n_items, item_urn_size, items = TRUE)
item_starting_MU = matrix(c(create_starting(delta, n_items, item_urn_size_MU[1], items = TRUE),
                            create_starting(delta, n_items, item_urn_size_MU[2], items = TRUE),
                            create_starting(delta, n_items, item_urn_size_MU[3], items = TRUE)),
                          3,n_items, byrow = TRUE)
rm(delta)

# starting point for students
cold_value = 0.5
################################################################################
# simulation factors
################################################################################
student_urn_size_mu = matrix(c(40,20,10,80,40,20,160,80,40), nrow = 3, ncol = 3, byrow = TRUE) #70 140 280 accuracy and multiple stakes
student_urn_size = c(70,140,280)

