################################################################################
# Loading libraries
################################################################################
library(tidyverse)
library(patchwork)
source("simulation/fixed_factors.R")
source("util.R")
set.seed(19186553)
################################################################################
# loading the data
################################################################################
#design
Design = matrix(0, 3, 2)
Design[,1] = rep(999, times = 3)
Design[,2] = student_urn_size


accuracy = readRDS("output/AccResultCsXUrn_2.Rds")

acc_dat = map_dfr(1:length(accuracy), 
                  ~(accuracy[[.x]]$mean - transform_logit(accuracy[[.x]]$true))^2 %>% 
                    as.data.frame() %>%
                    mutate(person = 1:1000) %>%
                    pivot_longer(-person, names_to = "game", values_to = "value") %>%
                    mutate(game = as.numeric(str_remove(game, "V")),
                           urn_size = Design[.x,2],
                           cs = Design[.x,1],
                           algo = "Accuracy",
                           type = "est"))

multiple_stakes = readRDS("output/MSResultCsXUrn_2.Rds")

ms_dat = map_dfr(1:length(multiple_stakes), 
                 ~(multiple_stakes[[.x]]$mean - transform_logit(multiple_stakes[[.x]]$true))^2 %>% 
                   as.data.frame() %>%
                   mutate(person = 1:1000) %>%
                   pivot_longer(-person, names_to = "game", values_to = "value") %>%
                   mutate(game = as.numeric(str_remove(game, "V")),
                          urn_size = Design[.x,2],
                          cs = Design[.x,1],
                          algo = "Multiple Updates",
                          type = "est"))

multiple_urns = readRDS("output/MUResultCsXUrn_2.Rds")
mu_dat = map_dfr(1:length(multiple_urns), 
                 ~((multiple_urns[[.x]]$mean_1+multiple_urns[[.x]]$mean_2+multiple_urns[[.x]]$mean_3)/sum(student_urn_size_mu[.x, ]) - transform_logit(multiple_urns[[.x]]$true))^2  %>% 
                   as.data.frame() %>%
                   mutate(person = 1:1000) %>%
                   pivot_longer(-person, names_to = "game", values_to = "value") %>%
                   mutate(game = as.numeric(str_remove(game, "V")),
                          urn_size = Design[.x,2],
                          cs = Design[.x,1],
                          algo = "Multiple Urns",
                          type = "est"))

dat = rbind(acc_dat, ms_dat, mu_dat)


dat = dat %>%
  group_by(urn_size, cs, algo, game) %>%
  summarise(value = mean(value))

################################################################################
# hitting times
################################################################################
#hitting time
true_vals =  qnorm(seq(1/(n_students+1),n_students/(n_students+1),by=1/(n_students+1)),mu_theta_baseline,sigma_theta)
true_vals = sample(true_vals, n_students)

HT_baseline = matrix(0,nrow = 3*length(true_vals), ncol = 5)
HT_baseline[,1] = rep(student_urn_size, each = length(true_vals))
HT_baseline[,2] = rep(transform_logit(true_vals), times = 3)
true_dist_holder = numeric(length=500)
for(i in 1:nrow(HT_baseline)){
  for(j in 1:500){
    true_dist_holder[j] = rbinom(1, HT_baseline[i,1], HT_baseline[i,2]) / HT_baseline[i,1]
  }
  HT_baseline[i, 3:4] = t.test(true_dist_holder)$conf.int[1:2] #we might not need this
  HT_baseline[i, 5] = (mean(true_dist_holder) - HT_baseline[i,2])^2
}

colnames(HT_baseline) = c("urn_size", "true_value", "lCI", "uCI", "MSE")
baseline_MSE = HT_baseline %>%
  as.data.frame() %>%
  group_by(urn_size) %>%
  summarise(lower = t.test(MSE)$conf.int[1],
            upper = t.test(MSE)$conf.int[2],
            mMSE = mean(MSE))


# HT_container = matrix(0,nrow=3*length(true_vals),ncol = 5)
# HT_container[,1:2] = HT_baseline[,1:2]
# ht_dat contains both algorithm type, and all the mean matricies
# ht_dat = map_dfr(1:3, ~accuracy[[.x]]$mean %>%
#                   as.data.frame())
# 
# ht_dat_ms = map_dfr(1:3, ~multiple_stakes[[.x]]$mean %>%
#                      as.data.frame())
# ht_dat_mu = map_dfr(1:3, ~as.data.frame(multiple_urns[[.x]]$mean_1+multiple_urns[[.x]]$mean_2+multiple_urns[[.x]]$mean_3)/sum(student_urn_size_mu[.x, ]))
#
# for(i in 1:nrow(HT_container)){
#   HT_container[i,3] = calculate_ht(as.numeric(ht_dat[i,]), HT_baseline[i,3], HT_baseline[i,4])
#   HT_container[i,4] = calculate_ht(as.numeric(ht_dat_ms[i,]), HT_baseline[i,3], HT_baseline[i,4])
#   HT_container[i,5] = calculate_ht(as.numeric(ht_dat_mu[i,]), HT_baseline[i,3], HT_baseline[i,4])
# }
# colnames(HT_container) = c("urn_size", "true_value", "ht_acc", "ht_ms", "ht_mu")

HT_container = matrix(0,nrow=9,ncol = 5)
HT_container[,1] = rep(student_urn_size, times = 3)
HT_container[,3:4] = rep(as.matrix(baseline_MSE[,2:3]), times = 3)
HT_container = as.data.frame(HT_container)
HT_container[,2] = rep(c("Accuracy", "Multiple Updates", "Multiple Stakes"), each = 3)
colnames(HT_container) = c("urn_size", "algo", "lCI", "uCI", "HT")

mse_chain = dat %>%
  left_join(HT_container, by = c("urn_size", "algo")) %>%
  as.data.frame() %>%
  group_by(urn_size, algo) %>%
  mutate(HT_bool = uCI > value) %>% #get the game value when HT_bool is first true in each group
  filter(HT_bool == TRUE)

################################################################################
# mse convergence plot
################################################################################
dat = dat %>%
  left_join(baseline_MSE, by = "urn_size")

plt_mse = dat %>% 
  ggplot(aes(x=game,
             y=value,
             color = algo)) +
  geom_line() +
  geom_hline(aes(yintercept = MSE), linetype = "dashed") +
  facet_wrap(vars(urn_size), ncol = 3) + 
  labs(x = "Items Answered", y = "MSE") +
  jtools::theme_apa(legend.font.size = 10) + #change the names of the legends
  scale_color_manual(values = c("black","red","green"))

################################################################################
# given a hitting time how more precise can we be
################################################################################





