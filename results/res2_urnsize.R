library(tidyverse)
library(patchwork)
source("simulation/fixed_factors.R")
source("util.R")
multiple_stakes = readRDS("output/MSResultUrn.Rds")
accuracy = readRDS("output/AccResultUrn.Rds")
multiple_urns = readRDS("output/MUResultUrn.Rds")

#data wrangling
urn_sizes = c(70,140,280)
dat_acc = map_dfr(1:3, ~accuracy[[.x]]$mean[993:1000,] %>% 
                    as.data.frame() %>%
                    mutate(person = 1:8) %>% 
                    pivot_longer(-person, names_to = "game", values_to = "value") %>%
                    mutate(game = as.numeric(str_remove(game, "V")),
                    urn_size = urn_sizes[.x],
                    algo = "accuracy",
                    type = "est"))

dat_ms = map_dfr(1:3, ~multiple_stakes[[.x]]$mean[993:1000,] %>% 
                    as.data.frame() %>%
                    mutate(person = 1:8) %>% 
                    pivot_longer(-person, names_to = "game", values_to = "value") %>%
                    mutate(game = as.numeric(str_remove(game, "V")),
                    urn_size = urn_sizes[.x],
                    algo = "multiple_stakes",
                    type = "est"))

dat_mu = map_dfr(1:3, ~((multiple_urns[[.x]]$mean_1[993:1000,] + 
                        multiple_urns[[.x]]$mean_2[993:1000,] +
                        multiple_urns[[.x]]$mean_3[993:1000,])/sum(student_urn_size_mu[.x, ])) %>% 
                    as.data.frame() %>%
                    mutate(person = 1:8) %>% 
                    pivot_longer(-person, names_to = "game", values_to = "value") %>%
                    mutate(game = as.numeric(str_remove(game, "V")),
                    urn_size = urn_sizes[.x],
                    algo = "multiple_urns",
                    type = "est"))

true_vals = rep(exp(seq(-1,1, length.out = 8))/(1+exp(seq(-1,1, length.out = 8))),each = 500)
#construct a true value tibble

dat = rbind(dat_acc, dat_ms, dat_mu)

dat_true = dat
dat_true$value = rep(true_vals, times = 9)
dat_true$type = "true"
dat_true$algo = "accuracy"
dat = rbind(dat,dat_true)
rm(dat_acc, dat_ms,dat_true)

#plotting 
#3 facets based on urn size, black lines accuracy, red lines multiple stakes, dashed lines true values
#each person is plotted seperately with 1 line.

plt_trace1 = dat %>% 
  ggplot(aes(x=game,
             y=value,
             group = interaction(factor(person),algo,type),
             color = algo,
             linetype = type)) +
  geom_line() +
  facet_wrap(~urn_size) +
  labs(x = "Items Answered", y = "Rating") +
  jtools::theme_apa(legend.font.size = 10) + #change lengend setting
  scale_color_manual(values = c("black","red", "green"))
  
  
#hitting time
#calculate the true distribution of the true values based on the urn size
HT_baseline = matrix(0,nrow = 24, ncol = 4)
HT_baseline[,1] = rep(urn_sizes, each = 8)
HT_baseline[,2] = rep(unique(true_vals), times = 3)
true_dist_holder = numeric(length=500)
for(i in 1:nrow(HT_baseline)){
  for(j in 1:500){
    true_dist_holder[j] = rbinom(1, HT_baseline[i,1], HT_baseline[i,2]) / HT_baseline[i,1]
  }
  HT_baseline[i, 3:4] = t.test(true_dist_holder)$conf.int[1:2]
}

HT_container = matrix(0,nrow=24,ncol = 5)
HT_container[,1:2] = HT_baseline[,1:2]

#ht_dat contains both algorithm type, and all the mean matricies 
ht_dat = map_dfr(1:3, ~accuracy[[.x]]$mean[993:1000,] %>% 
                   as.data.frame())
                 
ht_dat_ms = map_dfr(1:3, ~multiple_stakes[[.x]]$mean[993:1000,] %>% 
                                      as.data.frame())
ht_dat_mu = map_dfr(1:3, ~as.data.frame(multiple_urns[[.x]]$mean_1[993:1000,]+multiple_urns[[.x]]$mean_2[993:1000,]+multiple_urns[[.x]]$mean_3[993:1000,])/sum(student_urn_size_mu[.x, ]))

for(i in 1:nrow(HT_container)){
  HT_container[i,3] = calculate_ht(as.numeric(ht_dat[i,]), HT_baseline[i,3], HT_baseline[i,4])
  HT_container[i,4] = calculate_ht(as.numeric(ht_dat_ms[i,]), HT_baseline[i,3], HT_baseline[i,4])
  HT_container[i,5] = calculate_ht(as.numeric(ht_dat_mu[i,]), HT_baseline[i,3], HT_baseline[i,4])
}

colnames(HT_container) = c("urn_size","true_val","acc","ms","mu")

HT_container_person = HT_container %>% 
  as.data.frame() %>%
  pivot_longer(cols = c("acc","ms","mu"), names_to = "algo", values_to = "ht") %>%
  mutate(algo = ifelse(algo == "acc", "accuracy", ifelse(algo == "ms", "multiple_stakes", "multiple_urns"))) 

HT_container_plt = HT_container_person %>%
  group_by(urn_size, algo) %>%
  summarise(mean = mean(ht),
            sd = sd(ht))

#create point plot where x is urn size, y is hitting time, and color is algorithm
#connect points with line and add sd as an errorbar
plt_urnXavgHT = 
  HT_container_plt %>% ggplot(aes(x = factor(urn_size), y = mean, color = algo)) +
  geom_point() +
  geom_line(aes(group = algo)) +
  labs(x = "Urn Size", y = "Hitting Time") +
  jtools::theme_apa(legend.font.size = 10) + 
  scale_color_manual(values = c("black","red","green"))

#create a similar point plot, where x are the true vals, y is hitting time, and color is algorithm
#make facets for the different urn sizes
plt_trueXurnHT = HT_container_person %>% mutate_at(vars(true_val), round, 2) %>%
  ggplot(aes(x = factor(true_val), y = ht, color = algo)) +
  geom_point() +
  geom_line(aes(group = algo)) +
  facet_wrap(~urn_size) +
  labs(x = "True Value", y = "Hitting Time") +
  jtools::theme_apa(legend.font.size = 10) + 
  scale_color_manual(values = c("black","red","green"))

plt_patch_urnsize = plt_urnXavgHT + plt_trueXurnHT + plot_layout(guides = "collect", ncol = 1, nrow = 2)

#MSE of the first 250 iterations
MSE_container = HT_container
ht_dat = map_dfr(1:3, ~(accuracy[[.x]]$mean[993:1000,1:250] - matrix(rep(unique(true_vals), times = 250), nrow = 8, ncol = 250))^2%>% 
                   as.data.frame())
ht_dat = rowMeans(ht_dat)
ht_dat_ms = map_dfr(1:3, ~(multiple_stakes[[.x]]$mean[993:1000,1:250] - matrix(rep(unique(true_vals), times = 250), nrow = 8, ncol = 250))^2%>% 
                   as.data.frame())
ht_dat_ms = rowMeans(ht_dat_ms)
MSE_container[,3] = ht_dat
MSE_container[,4] = ht_dat_ms

MSE_container_plt_person = MSE_container %>% 
  as.data.frame() %>%
  pivot_longer(cols = c("acc","ms"), names_to = "algo", values_to = "mse") %>%
  mutate(algo = ifelse(algo == "acc", "accuracy", "multiple_stakes"))

MSE_container_plt = MSE_container_plt_person %>%
  group_by(urn_size, algo) %>%
  summarise(mean = mean(mse),
            sd = sd(mse))

#plot urn size x mse, color is algorithm
plt_urnXavgMSE = 
  MSE_container_plt %>% ggplot(aes(x = factor(urn_size), y = mean, color = algo)) +
  geom_point() +
  geom_line(aes(group = algo)) +
  labs(x = "Urn Size", y = "MSE") +
  jtools::theme_apa(legend.font.size = 10) + 
  scale_color_manual(values = c("black","red"))

#plot true val x mse, color is algorithm facet by urn size
plt_trueXurnMSE = MSE_container_plt_person %>% mutate_at(vars(true_val), round, 2) %>%
  ggplot(aes(x = factor(true_val), y = mse, color = algo)) +
  geom_point() +
  geom_line(aes(group = algo)) +
  facet_wrap(~urn_size) +
  labs(x = "True Value", y = "MSE") +
  jtools::theme_apa(legend.font.size = 10) + 
  scale_color_manual(values = c("black","red"))

#final plots
plt_trace1
plt_patch_urnXmeasure = plt_trueXurnHT + plt_trueXurnMSE + plot_layout(guides = "collect", ncol = 1, nrow = 2)
plt_patch_urnXavg = plt_urnXavgHT + plt_urnXavgMSE + plot_layout(guides = "collect", ncol = 1, nrow = 2)

ggsave("figures/fig2_trace.png", plt_trace1)

#checking item params
plot(colMeans(multiple_stakes[[3]]$item_mean), type = "l")

