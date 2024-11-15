#-------------------------------------------------------------------------------
# Description: This script is used to plot the results of the first simulation
#-------------------------------------------------------------------------------

################################################################################
# Dependencies
################################################################################
library(tidyverse)
library(patchwork)
library(kableExtra)
source("util.R")
source("analysis/fixed_factors.R")

################################################################################
# Loading data
################################################################################
accuracy = readRDS("output/sim1_coldstart_acc.rds")
multiple_stakes = readRDS("output/sim1_coldstart_ms.rds")
multiple_urns = readRDS("output/sim1_coldstart_mu.rds")

MSE_baselines = readRDS("output/MSE_baseline.rds")

################################################################################
# Data wrangling and plotting
################################################################################

#------------------------------------------------------------------------------
# Creating a suitable data format for plotting

acc_dat = map_dfr(1:length(accuracy), 
                  ~accuracy[[.x]]$MSE %>% 
                    as.data.frame() %>%
                    mutate(iter = 1:n_games) %>% 
                    mutate(urn_size = student_urn_size[.x],
                    algo = "Accuracy"))

ms_dat = map_dfr(1:length(multiple_stakes), 
                  ~multiple_stakes[[.x]]$MSE %>% 
                    as.data.frame() %>%
                    mutate(iter = 1:n_games) %>% 
                    mutate(urn_size = student_urn_size[.x],
                    algo = "Multiple Updates"))

mu_dat = map_dfr(1:length(multiple_urns), 
                  ~multiple_urns[[.x]]$MSE %>% 
                    as.data.frame() %>%
                    mutate(iter = 1:n_games) %>% 
                    mutate(urn_size = student_urn_size[.x],
                    algo = "Multiple Urns"))

dat = rbind(acc_dat, ms_dat, mu_dat)
colnames(dat)[1] = "value"
rm(acc_dat, ms_dat, mu_dat)

#------------------------------------------------------------------------------
# Plot MSE for each urn size facet wrapped by algorithm

p_mse = dat %>%
  left_join(as.data.frame(MSE_baselines), by = "urn_size") %>%
  ggplot(aes(x = iter, y = value, color = algo, linetype = algo)) +
  geom_line() +
  geom_hline(aes(yintercept = mean_theoretical), linetype = "dashed") +
  facet_wrap(~urn_size, scales = "free", ncol = 1) +
  labs(x = "Game",
       y = "MSE") +
  scale_color_manual(values = c("black","red","green")) + 
  jtools::theme_apa(legend.font.size = 10) +
  theme(legend.position = "none") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 


#------------------------------------------------------------------------------
# plot Hitting times with central intervals

acc_ht = map_dfr(1:length(accuracy), 
                  ~accuracy[[.x]]$HT %>% 
                    as.data.frame() %>%
                    mutate(urn_size = student_urn_size[.x],
                           algo = "Accuracy",
                           lower = accuracy[[.x]]$HT_sd[1],
                           upper = accuracy[[.x]]$HT_sd[2]))

ms_ht = map_dfr(1:length(multiple_stakes),
                 ~multiple_stakes[[.x]]$HT %>% 
                   as.data.frame() %>%
                   mutate(urn_size = student_urn_size[.x],
                          algo = "Multiple Updates",
                          lower = multiple_stakes[[.x]]$HT_sd[1],
                          upper = multiple_stakes[[.x]]$HT_sd[2]))

mu_ht = map_dfr(1:length(multiple_urns),
                 ~multiple_urns[[.x]]$HT %>% 
                   as.data.frame() %>%
                   mutate(urn_size = student_urn_size[.x],
                          algo = "Multiple Urns",
                          lower = multiple_urns[[.x]]$HT_sd[1],
                          upper = multiple_urns[[.x]]$HT_sd[2]))


dat_ht = rbind(acc_ht, ms_ht, mu_ht)
rm(acc_ht, ms_ht, mu_ht)
colnames(dat_ht)[1] = "HT"

#-------------------------------------------------------------------------------
# Plot hitting times as points with central intervals as error bars

p_ht = dat_ht %>%
  ggplot(aes(x = as.factor(urn_size), y = HT, color = algo)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) + # add a line plot to connect the points of the same color
  geom_line(aes(group = algo)) +
  labs(x = "Urn Size",
       y = "Hitting Time") +
  scale_color_manual(values = c("black","red","green")) + 
  jtools::theme_apa() + #change lengend setting
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

# making a table in apa format which can be exported to a latex file
dat_ht = dat_ht %>%
          select(algo, urn_size, HT, lower, upper)
dat_ht_raw = dat_ht
colnames(dat_ht) = c("Algorithm", "Urn Size", "Hitting Time", "2.5%", "97.5%")

dat_ht %>%
  kbl(caption = "Hitting times and [2.5%, 97.5%] quantiles of the respective sampling distribution 
      for each algorithm and each urn sizes 95% central intervals",
      format = "latex") %>%
  kable_classic(full_width = TRUE, html_font = "helvetica")


#-------------------------------------------------------------------------------
# Combine the two plots

patch_mse_ht = p_mse + p_ht + plot_layout(guides = "auto", ncol = 2)

#ggsave("figures/fig1_coldstart.png", patch_mse_ht, dpi = 300)

tikz(file = "figures/fig1_coldstart.tex", width = 5, height = 5)
patch_mse_ht = p_mse + p_ht + plot_layout(guides = "auto", ncol = 2) +
  theme(axis.title.y = element_text(size = rel(0.8)),
        axis.title.x = element_text(size = rel(0.8)),
        legend.text = element_text(size = rel(0.8)),
        plot.title = element_text(size = rel(0.8)))
print(patch_mse_ht)
dev.off()

#-------------------------------------------------------------------------------
# Descriptive values for the text
dat_desc = dat %>%
  left_join(as.data.frame(MSE_baselines), by = "urn_size") 
dat_desc %>%
  filter(iter<100) %>%
  group_by(algo) %>%
  summarise(mMSE = mean(value),
            sdMSE = sd(value),
            mean_theoretical = mean(mean_theoretical),
            lower = mean(lower),
            upper = mean(upper))

dat_desc %>%
  filter(iter == 50)


#comparison of HT-s
dat_ht_raw %>%
  mutate(range = upper - lower) 

dat_ht_raw %>%
  group_by(algo) %>%
  summarise(mHT = mean(HT))
