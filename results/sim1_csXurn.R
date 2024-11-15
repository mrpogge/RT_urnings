################################################################################
# Loading libraries
################################################################################
library(tidyverse)
library(patchwork)
source("simulation/fixed_factors.R")
source("util.R")

################################################################################
# loading the data
################################################################################
#design
Design = matrix(0, 6, 2)
Design[,1] = rep(c(100,999), times = 3)
Design[,2] = rep(student_urn_size, times = 2)

accuracy = readRDS("output/AccResultCsXUrn_2.Rds")

acc_dat = map_dfr(1:length(accuracy), 
                  ~(accuracy[[.x]]$mean - transform_logit(accuracy[[.x]]$true))^2 %>% 
                    as.data.frame() %>%
                    mutate(person = 1:1000) %>%
                    pivot_longer(-person, names_to = "game", values_to = "value") %>%
                    mutate(game = as.numeric(str_remove(game, "V")),
                           urn_size = Design[.x,2],
                           cs = Design[.x,1],
                           algo = "accuracy",
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
                          algo = "multiple_stakes",
                          type = "est"))

multiple_urns = readRDS("output/MUResultCsXUrn_2.Rds")

mu_dat = map_dfr(1:length(multiple_urns), 
                 ~((multiple_urns[[.x]]$mean_1+multiple_urns[[.x]]$mean_2+multiple_urns[[.x]]$mean_3)/sum(student_urn_size_mu[.x, ])^2 - transform_logit(multiple_urns[[.x]]$true))^2  %>% 
                   as.data.frame() %>%
                   mutate(person = 1:1000) %>%
                   pivot_longer(-person, names_to = "game", values_to = "value") %>%
                   mutate(game = as.numeric(str_remove(game, "V")),
                          urn_size = Design[.x,2],
                          cs = Design[.x,1],
                          algo = "multiple_urns",
                          type = "est"))

dat = rbind(acc_dat, ms_dat, mu_dat)

dat = dat %>%
  group_by(urn_size, cs, algo, game) %>%
  summarise(value = mean(value))

################################################################################
# mse convergence plot
################################################################################
plt_mse = dat %>% 
  ggplot(aes(x=game,
             y=value,
             color = algo)) +
  geom_line() +
  facet_wrap(urn_size~cs, ncol = 2) +
  labs(x = "Items Answered", y = "MSE") +
  jtools::theme_apa(legend.font.size = 10) + #change lengend setting
  scale_color_manual(values = c("black","red"))


################################################################################
# loading the data
################################################################################
#design
Design = matrix(0, 6, 2)
Design[,1] = rep(c(100,999), times = 3)
Design[,2] = rep(student_urn_size, times = 2)

accuracy = readRDS("output/AccResultCsXUrn_6.Rds")

acc_dat = map_dfr(1:length(accuracy), 
                  ~(accuracy[[.x]]$mean - transform_logit(accuracy[[.x]]$true))^2 %>% 
                    as.data.frame() %>%
                    mutate(person = 1:1000) %>%
                    pivot_longer(-person, names_to = "game", values_to = "value") %>%
                    mutate(game = as.numeric(str_remove(game, "V")),
                           urn_size = Design[.x,2],
                           cs = Design[.x,1],
                           algo = "accuracy",
                           type = "est"))

multiple_stakes = readRDS("output/MSResultCsXUrn_6.Rds")

ms_dat = map_dfr(1:length(multiple_stakes), 
                 ~(multiple_stakes[[.x]]$mean - transform_logit(multiple_stakes[[.x]]$true))^2 %>% 
                   as.data.frame() %>%
                   mutate(person = 1:1000) %>%
                   pivot_longer(-person, names_to = "game", values_to = "value") %>%
                   mutate(game = as.numeric(str_remove(game, "V")),
                          urn_size = Design[.x,2],
                          cs = Design[.x,1],
                          algo = "multiple_stakes",
                          type = "est"))
dat = rbind(acc_dat, ms_dat)

dat = dat %>%
  group_by(urn_size, cs, algo, game) %>%
  summarise(value = mean(value))

################################################################################
# mse convergence plot
################################################################################
plt_mse_6 = dat %>% 
  ggplot(aes(x=game,
             y=value,
             color = algo)) +
  geom_line() +
  facet_wrap(urn_size~cs, ncol = 2) +
  labs(x = "Items Answered", y = "MSE") +
  jtools::theme_apa(legend.font.size = 10) + #change lengend setting
  scale_color_manual(values = c("black","red"))

################################################################################
# loading the data
################################################################################
#design
Design = matrix(0, 6, 2)
Design[,1] = rep(c(100,999), times = 3)
Design[,2] = rep(student_urn_size, times = 2)

accuracy = readRDS("output/AccResultCsXUrn_8.Rds")

acc_dat = map_dfr(1:length(accuracy), 
                  ~(accuracy[[.x]]$mean - transform_logit(accuracy[[.x]]$true))^2 %>% 
                    as.data.frame() %>%
                    mutate(person = 1:1000) %>%
                    pivot_longer(-person, names_to = "game", values_to = "value") %>%
                    mutate(game = as.numeric(str_remove(game, "V")),
                           urn_size = Design[.x,2],
                           cs = Design[.x,1],
                           algo = "accuracy",
                           type = "est"))

multiple_stakes = readRDS("output/MSResultCsXUrn_8.Rds")

ms_dat = map_dfr(1:length(multiple_stakes), 
                 ~(multiple_stakes[[.x]]$mean - transform_logit(multiple_stakes[[.x]]$true))^2 %>% 
                   as.data.frame() %>%
                   mutate(person = 1:1000) %>%
                   pivot_longer(-person, names_to = "game", values_to = "value") %>%
                   mutate(game = as.numeric(str_remove(game, "V")),
                          urn_size = Design[.x,2],
                          cs = Design[.x,1],
                          algo = "multiple_stakes",
                          type = "est"))
dat = rbind(acc_dat, ms_dat)

dat = dat %>%
  group_by(urn_size, cs, algo, game) %>%
  summarise(value = mean(value))

################################################################################
# mse convergence plot
################################################################################
plt_mse_8 = dat %>% 
  ggplot(aes(x=game,
             y=value,
             color = algo)) +
  geom_line() +
  facet_wrap(urn_size~cs, ncol = 2) +
  labs(x = "Items Answered", y = "MSE") +
  jtools::theme_apa(legend.font.size = 10) + #change lengend setting
  scale_color_manual(values = c("black","red"))




plt_mse / plt_mse_6 / plt_mse_8
