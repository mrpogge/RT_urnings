################################################################################
# loading libraries
################################################################################
library(tidyverse)
library(patchwork)
source("simulation/fixed_factors.R")
source("util.R")

################################################################################
# loading data
################################################################################
multiple_stakes = readRDS("output/MSResultUrnXAdapt.Rds")
accuracy = readRDS("output/AccResultUrnXAdapt.Rds")

#recreating design matrix for simplicty
Design = matrix(0, 12, 4)
Design[,1] = rep(1:3, each = 4)
Design[,2] = rep(student_urn_size, times = 4)
Design[,3] = rep(mu_p, times = 3)
Design[,4] = rep(mu_theta, times = 3) 

################################################################################
# preparing data for trace
################################################################################
# accuracy
acc_dat = map_dfr(1:length(multiple_stakes), 
                  ~accuracy[[.x]]$mean[993:1000,] %>% 
                    as.data.frame() %>%
                    mutate(person = 1:8) %>% 
                    pivot_longer(-person, names_to = "game", values_to = "value") %>%
                    mutate(game = as.numeric(str_remove(game, "V")),
                           urn_size = Design[.x,2],
                           adapt = Design[.x,3],
                           algo = "accuracy",
                           type = "est"))
# multiple stakes
ms_dat = map_dfr(1:length(multiple_stakes), 
                 ~multiple_stakes[[.x]]$mean[993:1000,] %>% 
                   as.data.frame() %>%
                   mutate(person = 1:8) %>% 
                   pivot_longer(-person, names_to = "game", values_to = "value") %>%
                   mutate(game = as.numeric(str_remove(game, "V")),
                          urn_size = Design[.x,2],
                          adapt = Design[.x,3],
                          algo = "multiple_stakes",
                          type = "est"))
dat = rbind(acc_dat, ms_dat)

# true values
trues = seq(-1,1, length.out = 8)
true_vals = rep(exp(trues)/(1+exp(trues)),each = 500)
dat_true = dat
dat_true$value = rep(true_vals, times = 24)
dat_true$type = "true"
dat_true$algo = "accuracy"
dat = rbind(dat,dat_true)
rm(acc_dat, ms_dat,dat_true)

################################################################################
# trace plot (looks terrible)
################################################################################
plt_trace = dat %>% 
  ggplot(aes(x=game,
             y=value,
             group = interaction(factor(person),type, algo, adapt),
             color = interaction(adapt),
             linetype = type)) +
  geom_line() +
  facet_wrap(vars(urn_size,algo), ncol = 2) +
  labs(x = "Items Answered", y = "Rating") +
  jtools::theme_apa(legend.font.size = 10) 

################################################################################
# hitting time
################################################################################
HT_baseline = matrix(0,nrow = 12*8, ncol = 5)
HT_baseline[,1] = rep(Design[,2], each = 8)
HT_baseline[,2] = rep(exp(trues)/(1+exp(trues)), times = 12)
HT_baseline[,3] = rep(Design[,3], each = 8)
true_dist_holder = numeric(length=1000)

for(i in 1:nrow(HT_baseline)){
  for(j in 1:1000){
    true_dist_holder[j] = rbinom(1, HT_baseline[i,1], HT_baseline[i,2]) / HT_baseline[i,1]
  }
  HT_baseline[i, 4:5] = t.test(true_dist_holder)$conf.int[1:2]
}

HT_container = matrix(0,nrow=12*8,ncol = 5)
HT_container[,1:3] = HT_baseline[,1:3]

#ht_dat contains both algorithm type, and all the mean matricies 
ht_dat = map_dfr(1:12, ~accuracy[[.x]]$mean[993:1000,] %>% 
                   as.data.frame())

ht_dat_ms = map_dfr(1:12, ~multiple_stakes[[.x]]$mean[993:1000,] %>% 
                      as.data.frame())

for(i in 1:(nrow(HT_container))){
  HT_container[i,4] = calculate_ht(as.numeric(ht_dat[i,]), HT_baseline[i,4], HT_baseline[i,5])
  HT_container[i,5] = calculate_ht(as.numeric(ht_dat_ms[i,]), HT_baseline[i,4], HT_baseline[i,5])
}

colnames(HT_container) = c("urn_size","true_val", "adapt","acc","ms")

HT_container_person = HT_container %>% 
  as.data.frame() %>%
  pivot_longer(cols = c("acc","ms"), names_to = "algo", values_to = "ht") %>%
  mutate(algo = ifelse(algo == "acc", "accuracy", "multiple_stakes")) 

HT_container_plt = HT_container_person %>%
  group_by(urn_size, algo, adapt) %>%
  summarise(mean = mean(ht),
            sd = sd(ht))

HT_container_adapt = HT_container_person %>%
  group_by(algo, adapt) %>%
  summarise(mean = mean(ht),
            sd = sd(ht))
#create point plot where x is urn size, y is hitting time, and color is algorithm
#connect points with line and add sd as an errorbar
plt_urnXavgHT = 
  HT_container_plt %>% ggplot(aes(x = factor(urn_size), y = mean, color = algo)) +
  geom_point() +
  geom_line(aes(group = algo)) +
  facet_wrap(~adapt) +
  labs(x = "Urn Size", y = "Hitting Time") +
  jtools::theme_apa(legend.font.size = 10) 

plt_urnXavgHT = 
  HT_container_adapt %>% ggplot(aes(x = factor(adapt), y = mean, color = algo)) +
  geom_point() +
  geom_line(aes(group = algo)) +
  labs(x = "Probability Correct", y = "Hitting Time") +
  jtools::theme_apa(legend.font.size = 10) 

#create a similar point plot, where x are the true vals, y is hitting time, and color is algorithm
#make facets for the different urn sizes
plt_trueXurnHT = HT_container_person %>% mutate_at(vars(true_val), round, 2) %>%
  ggplot(aes(x = factor(true_val), y = ht, color = algo)) +
  geom_point() +
  geom_line(aes(group = algo)) +
  facet_wrap(adapt~urn_size, ncol = 3) +
  labs(x = "True Value", y = "Hitting Time") +
  jtools::theme_apa(legend.font.size = 10) + 
  scale_color_manual(values = c("black","red"))

################################################################################
#MSE plots
################################################################################
#MSE of the first 250 iterations
MSE_container = HT_container
ht_dat = map_dfr(1:12, ~(accuracy[[.x]]$mean[993:1000,1:250] - matrix(rep(unique(true_vals), times = 250), nrow = 8, ncol = 250))^2%>% 
                   as.data.frame())
ht_dat = rowMeans(ht_dat)
ht_dat_ms = map_dfr(1:12, ~(multiple_stakes[[.x]]$mean[993:1000,1:250] - matrix(rep(unique(true_vals), times = 250), nrow = 8, ncol = 250))^2%>% 
                      as.data.frame())
ht_dat_ms = rowMeans(ht_dat_ms)
MSE_container[,4] = ht_dat
MSE_container[,5] = ht_dat_ms

MSE_container_plt_person = MSE_container %>% 
  as.data.frame() %>%
  pivot_longer(cols = c("acc","ms"), names_to = "algo", values_to = "mse") %>%
  mutate(algo = ifelse(algo == "acc", "accuracy", "multiple_stakes"))

MSE_container_plt = MSE_container_plt_person %>%
  group_by(adapt, urn_size, algo) %>%
  summarise(mean = mean(mse),
            sd = sd(mse))

#plot urn size x mse, color is algorithm
plt_urnXavgMSE = 
  MSE_container_plt %>% ggplot(aes(x = factor(urn_size), y = mean, color = algo)) +
  geom_point() +
  geom_line(aes(group = algo)) +
  facet_wrap(~adapt) +
  labs(x = "Urn Size", y = "MSE") +
  jtools::theme_apa(legend.font.size = 10) + 
  scale_color_manual(values = c("black","red"))

#plot true val x mse, color is algorithm facet by urn size
plt_trueXurnMSE = MSE_container_plt_person %>% mutate_at(vars(true_val), round, 2) %>%
  ggplot(aes(x = factor(true_val), y = mse, color = algo)) +
  geom_point() +
  geom_line(aes(group = algo)) +
  facet_wrap(adapt~urn_size, ncol = 3) +
  labs(x = "True Value", y = "MSE") +
  jtools::theme_apa(legend.font.size = 10) + 
  scale_color_manual(values = c("black","red"))




