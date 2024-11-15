################################################################################
# libraries and deps
################################################################################

library(tidyverse)
library(patchwork)
source("util.R")
source("simulation/fixed_factors.R")

################################################################################
# loading data
################################################################################
multiple_stakes = readRDS("output/MSChange.Rds")
accuracy = readRDS("output/AccChange.Rds")

Design = matrix(student_urn_size, nrow = 3, ncol = 1)

acc_dat = map_dfr(1:length(accuracy), 
                   ~accuracy[[.x]]$mean[961:1000,] %>% 
                     as.data.frame() %>%
                     mutate(person = 1:40) %>% 
                     mutate(person_starting = rep(1:8, each = 5)) %>%
                     mutate(person_changing = rep(1:5, times = 8)) %>%
                     pivot_longer(-starts_with("person"), names_to = "game", values_to = "value") %>%
                     mutate(game = as.numeric(str_remove(game, "V")),
                            urn_size = Design[.x,1],
                            algo = "accuracy",
                            type = "est"))

acc_dat_true = map_dfr(1:length(accuracy), 
                        ~transform_logit(accuracy[[.x]]$true[961:1000,]) %>% 
                          as.data.frame() %>%
                          mutate(person = 1:40) %>% 
                          mutate(person_starting = rep(1:8, each = 5)) %>%
                          mutate(person_changing = rep(1:5, times = 8)) %>%
                          pivot_longer(-starts_with("person"), names_to = "game", values_to = "value") %>%
                          mutate(game = as.numeric(str_remove(game, "V")),
                                 urn_size = Design[.x,1],
                                 algo = "accuracy",
                                 type = "true"))

ms_dat = map_dfr(1:length(multiple_stakes), 
                 ~multiple_stakes[[.x]]$mean[961:1000,] %>% 
                   as.data.frame() %>%
                   mutate(person = 1:40) %>% 
                   mutate(person_starting = rep(1:8, each = 5)) %>%
                   mutate(person_changing = rep(1:5, times = 8)) %>%
                   pivot_longer(-starts_with("person"), names_to = "game", values_to = "value") %>%
                   mutate(game = as.numeric(str_remove(game, "V")),
                          urn_size = Design[.x,1],
                          algo = "multiple_stakes",
                          type = "est"))

ms_dat_true = map_dfr(1:length(accuracy), 
                      ~transform_logit(multiple_stakes[[.x]]$true[961:1000,]) %>% 
                        as.data.frame() %>%
                        mutate(person = 1:40) %>% 
                        mutate(person_starting = rep(1:8, each = 5)) %>%
                        mutate(person_changing = rep(1:5, times = 8)) %>%
                        pivot_longer(-starts_with("person"), names_to = "game", values_to = "value") %>%
                        mutate(game = as.numeric(str_remove(game, "V")),
                               urn_size = Design[.x,1],
                               algo = "accuracy",
                               type = "true"))

dat = rbind(acc_dat, ms_dat, acc_dat_true, ms_dat_true)

################################################################################
# trace
################################################################################
plt_trace_1 = dat %>% 
  filter(person_starting == 1) %>%
  ggplot(aes(x=game,
             y=value,
             group = interaction(factor(person),algo,type),
             color = algo,
             linetype = type)) +
  geom_line() +
  facet_wrap(~urn_size) +
  labs(x = "Items Answered", y = "Rating") +
  jtools::theme_apa(legend.font.size = 10) + #change lengend setting
  scale_color_manual(values = c("black","red"))

plt_trace_4 = dat %>% 
  filter(person_starting == 4) %>%
  ggplot(aes(x=game,
             y=value,
             group = interaction(factor(person),algo,type),
             color = algo,
             linetype = type)) +
  geom_line() +
  facet_wrap(~urn_size) +
  labs(x = "Items Answered", y = "Rating") +
  jtools::theme_apa(legend.font.size = 10) + #change lengend setting
  scale_color_manual(values = c("black","red"))

plt_trace_8 = dat %>% 
  filter(person_starting == 8) %>%
  ggplot(aes(x=game,
             y=value,
             group = interaction(factor(person),algo,type),
             color = algo,
             linetype = type)) +
  geom_line() +
  facet_wrap(~urn_size) +
  labs(x = "Items Answered", y = "Rating") +
  jtools::theme_apa(legend.font.size = 10) + #change lengend setting
  scale_color_manual(values = c("black","red"))

################################################################################
# Mean difference
################################################################################
dat_est = rbind(acc_dat, ms_dat)
ms_dat_true$algo = "multiple_stakes"
dat_true = rbind(acc_dat_true, ms_dat_true)

dat_md = dat_est %>% 
  left_join(dat_true, by = c("person", "person_starting", "person_changing", "game", "urn_size", "algo")) %>%
  mutate(md = (value.x - value.y)^2) %>%
  select(person, person_starting, person_changing, game, urn_size, algo, md)

dat_md_aggregate = dat_md %>% 
  group_by(person_starting, person_changing, urn_size, algo) %>%
  summarise(md = mean(md))

#starting value x md
plot_starting = dat_md_aggregate %>% 
                  group_by(person_starting, urn_size, algo) %>% 
                  summarise(md = mean(md)) %>%
                  ggplot(aes(x = person_starting, y = md, color = algo)) +
                  geom_point(aes(group = algo)) +
                  geom_line() +
                  facet_wrap(~urn_size) +
                  labs(x = "Starting Value", y = "Mean Difference") +
                  jtools::theme_apa(legend.font.size = 10) #change lengend setting

plot_changing = dat_md_aggregate %>% 
  group_by(person_changing, urn_size, algo) %>% 
  summarise(md = mean(md)) %>%
  ggplot(aes(x = person_changing, y = md, color = algo)) +
  geom_point(aes(group = algo)) +
  geom_line() +
  facet_wrap(~urn_size) +
  labs(x = "Amount of Change", y = "MSE") +
  jtools::theme_apa(legend.font.size = 10) #change lengend setting

plot_urn = dat_md_aggregate %>% 
  group_by(urn_size, algo) %>% 
  summarise(md = mean(md)) %>%
  ggplot(aes(x = as.factor(urn_size), y = md, color = algo)) +
  geom_point() +
  geom_line(aes(group = algo)) +
  labs(x = "Urn Size", y = "MSE") +
  jtools::theme_apa(legend.font.size = 10) #change lengend setting


