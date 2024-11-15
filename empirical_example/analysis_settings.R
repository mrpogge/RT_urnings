################################################################################
# Settings for the empirical analysis
################################################################################
#nreps
nreps = 100

# urn sizes
accuracy_urn_sizes = seq(20, 160, by = 4)
pseudoRT_urn_sizes = accuracy_urn_sizes / 2

combined_urn_sizes = accuracy_urn_sizes + pseudoRT_urn_sizes

item_urn_sizes_MU = c(40,20)
item_urn_size = 60

# Design
Design_mu = matrix(0, length(accuracy_urn_sizes) * 3, 3)
Design_mu[,1] = rep(1:3, times = length(accuracy_urn_sizes))
Design_mu[,2] = rep(accuracy_urn_sizes, each = 3)
Design_mu[,3] = rep(pseudoRT_urn_sizes, each = 3)

Design = matrix(0, length(accuracy_urn_sizes) * 3, 2)
Design[,1] = rep(1:3, times = length(accuracy_urn_sizes))
Design[,2] = rep(combined_urn_sizes, each = 3)