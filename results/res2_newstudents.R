library(tidyverse)
library(patchwork)
library(kableExtra)
source("util.R")
source("analysis/fixed_factors.R")

multiple_stakes = readRDS("output/sim2_newstudents_ms.rds")
accuracy = readRDS("output/sim2_newstudents_acc.rds")
multiple_urns = readRDS("output/sim2_newstudents_mu.rds")


dat_acc = map_dfr(1:3, ~(transform_expit(accuracy[[.x]]$mean)*4) %>% 
                    as.data.frame() %>%
                    mutate(person = 1:8) %>% 
                    pivot_longer(-person, names_to = "game", values_to = "value") %>%
                    mutate(game = as.numeric(str_remove(game, "V")),
                           urn_size = student_urn_size[.x],
                           algo = "Accuracy",
                           type = "Estimate"))

dat_ms = map_dfr(1:3, ~(transform_expit(multiple_stakes[[.x]]$mean)*4) %>% 
                   as.data.frame() %>%
                   mutate(person = 1:8) %>% 
                   pivot_longer(-person, names_to = "game", values_to = "value") %>%
                   mutate(game = as.numeric(str_remove(game, "V")),
                          urn_size = student_urn_size[.x],
                          algo = "Multiple Updates",
                          type = "Estimate"))

dat_mu = map_dfr(1:3, ~(transform_expit(multiple_urns[[.x]]$mean)*4) %>% 
                   as.data.frame() %>%
                   mutate(person = 1:8) %>% 
                   pivot_longer(-person, names_to = "game", values_to = "value") %>%
                   mutate(game = as.numeric(str_remove(game, "V")),
                          urn_size = student_urn_size[.x],
                          algo = "Multiple Urn",
                          type = "Estimate"))
dat = rbind(dat_acc, dat_ms, dat_mu) 

true_dat = dat
true_dat$value = rep(rep(rep((new_student_thetas*4), each = 250), times = 3), times = 3)
true_dat$type = "True value"
true_dat$algo = "Accuracy"
dat = rbind(dat, true_dat)

plt_trace1 = dat %>% 
  ggplot(aes(x=game,
             y=value,
             group = interaction(factor(person),algo,type),
             color = algo,
             linetype = type)) +
  geom_line() +
  facet_wrap(~urn_size) +
  labs(x = "Items Answered", y = "Estimate") +
  scale_color_manual(values = c("black","red", "green"))+
  jtools::theme_apa() + #change lengend setting
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
ggsave("figures/fig2_trace.png", plt_trace1, dpi = 300)

#--------------------------------------------------------------------------------
# adding hitting time/ urnings plot

acc_ht = map_dfr(1:3, ~accuracy[[.x]]$ht %>% 
                   as.data.frame() %>%
                   mutate(person = new_student_thetas*4, 
                          urn_size = student_urn_size[.x],
                          algo = "accuracy"))
ms_ht = map_dfr(1:3, ~multiple_stakes[[.x]]$ht %>% 
                   as.data.frame() %>%
                   mutate(person = new_student_thetas*4, 
                          urn_size = student_urn_size[.x],
                          algo = "multiple_stakes"))
mu_ht = map_dfr(1:3, ~multiple_urns[[.x]]$ht %>% 
                   as.data.frame() %>%
                   mutate(person = new_student_thetas*4, 
                          urn_size = student_urn_size[.x],
                          algo = "multiple_urns"))
ht_dat = rbind(acc_ht, ms_ht, mu_ht)
colnames(ht_dat)[1] = "HT"
#--------------------------------------------------------------------------------
#plotting mean hitting time over urn size

plt_ht = ht_dat %>% 
  group_by(urn_size, algo) %>%
  summarise(mHT = mean(HT)) %>%
  ggplot(aes(x=factor(urn_size),
             y=mHT,
             color = algo)) +
  geom_point() +
  geom_line(aes(group = algo)) +
  jtools::theme_apa()+ #change lengend setting
  scale_color_manual(values = c("black","red", "green"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

#--------------------------------------------------------------------------------
#plotting mean hitting time over distance of true vals from 0 

plt_ht2 = ht_dat %>%
  group_by(algo, person) %>%
  mutate(person = round(person, 2)) %>%
  summarise(mHT = mean(HT)) %>%
  ggplot(aes(x=factor(person),
             y=mHT,
             color = algo)) +
  geom_point() +
  geom_line(aes(group = algo)) +
  jtools::theme_apa()+ #change lengend setting
  scale_color_manual(values = c("black","red", "green"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

ht_dat %>% 
  mutate(person = round(person, 2)) %>%
  pivot_wider(names_from = person, values_from = HT) %>%
  kable("latex",
        booktabs = T,
        longtable = F,
        digits=3,
        col.names = c("Urn size", "Algorithm",  -2, -1.43, -0.86, -0.29, 0.29, 0.86, 1.43, 2),
        caption = "Hitting times of the new students using different urn sizes and different true values") %>%
  kable_classic(full_width = F, html_font = "helvetica")


#--------------------------------------------------------------------------------
# descriptive stats

ht_dat_desc = ht_dat %>% 
  mutate(person = round(person, 2)) 

ht_dat_desc %>%
  group_by(algo, urn_size) %>%
  summarise(mHT = mean(HT))

ht_dat_desc2 = ht_dat %>% 
  mutate(person = round(person, 2)) %>%
  pivot_wider(names_from = person, values_from = HT)
  
colMeans(ht_dat_desc2[,3:10])
