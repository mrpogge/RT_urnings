library(tidyverse)
library(patchwork)
library(kableExtra)
source("util.R")
source("analysis/fixed_factors_change.R")

#reading the results lists
multiple_stakes = readRDS("output/sim3_change_ms.rds")
accuracy = readRDS("output/sim3_change_acc.rds")
multiple_urns = readRDS("output/sim3_change_mu.rds")
MSE_baseline = readRDS("output/MSE_baseline_change.Rds")

acc_mse = map_dfr(1:length(accuracy), ~accuracy[[.x]]$MSE %>% 
                    as.data.frame() %>%
                    mutate(game = 1:501,
                           urn_size = student_urn_size[.x],
                           algo = "Accuracy",
                           is_est = "1est"))


ms_mse = map_dfr(1:length(multiple_stakes), ~multiple_stakes[[.x]]$MSE %>% 
                    as.data.frame() %>%
                    mutate(game = 1:501,
                           urn_size = student_urn_size[.x],
                           algo = "Multiple Updates",
                           is_est = "1est"))

mu_mse = map_dfr(1:length(multiple_urns), ~multiple_urns[[.x]]$MSE %>% 
                    as.data.frame() %>%
                    mutate(game = 1:501,
                           urn_size = student_urn_size[.x],
                           algo = "Multiple Urns",
                           is_est = "1est"))

acc_bl =  MSE_baseline %>% 
  as.data.frame() %>%
  select(mean, iter, urn_size) %>%
  mutate(algo = "Accuracy",
         is_est = "2baseline")

ms_bl = MSE_baseline %>% 
  as.data.frame() %>%
  select(mean, iter, urn_size) %>%
  mutate(algo = "Multiple Updates",
         is_est = "2baseline")
  
mu_bl = MSE_baseline %>% 
  as.data.frame() %>%
  select(mean, iter, urn_size) %>%
  mutate(algo = "Multiple Urns",
         is_est = "2baseline")

colnames(acc_bl)[1:2] = colnames(ms_bl)[1:2] = colnames(mu_bl)[1:2] = c("MSE", "game")
colnames(acc_mse)[1] = colnames(ms_mse)[1] = colnames(mu_mse)[1] = "MSE"
mse_dat = rbind(acc_mse, ms_mse, mu_mse, acc_bl, ms_bl, mu_bl)
mse_dat$is_est = factor(mse_dat$is_est, levels = c("1est", "2baseline"), labels = c("Estimate", "Baseline"))

mse_change = mse_dat %>%
  ggplot(aes(x = game, y = MSE, color = factor(urn_size), linetype = is_est)) +
  geom_line() +
  facet_wrap(vars(algo)) + 
  jtools::theme_apa(legend.font.size = 10) +
  labs(x = "Game",
       y = "MSE") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

ggsave("figures/fig3_changeMSE.png", mse_change)

acc_mean = map_dfr(c(1,3,7), ~accuracy[[.x]]$MSE %>% 
                     as.data.frame() %>%
                     mutate(game = 1:501,
                            urn_size = student_urn_size[.x],
                            algo = "accuracy"))

#--------------------------------------------------------------------------------
# coverage

acc_coverage = map_dfr(1:length(accuracy), ~accuracy[[.x]]$coverage %>% 
                         as.data.frame() %>%
                         mutate(urn_size = student_urn_size[.x],
                                algo = "accuracy"))

ms_coverage = map_dfr(1:length(multiple_stakes), ~multiple_stakes[[.x]]$coverage %>%
                         as.data.frame() %>%
                         mutate(urn_size = student_urn_size[.x],
                                algo = "multiple_stakes"))

mu_coverage = map_dfr(1:length(multiple_urns), ~multiple_urns[[.x]]$coverage %>%
                         as.data.frame() %>%
                         mutate(urn_size = student_urn_size[.x],
                                algo = "multiple_urns"))

coverage_dat = rbind(acc_coverage, ms_coverage, mu_coverage)
colnames(coverage_dat)[1] = "coverage"

coverage_change = 0.84 %>%
  ggplot(aes(y = coverage, x = factor(urn_size), group = factor(algo))) +
  geom_line(aes(color = factor(algo))) +
  jtools::theme_apa(legend.font.size = 10) +
  labs(x = "Student Urn Size",
       y = "Coverage") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

#coverage table
coverage_table = coverage_dat %>%
  pivot_wider(names_from = algo, values_from = coverage) %>%
  kable("latex",
        booktabs = T,
        longtable = F,
        digits=3,
        col.names = c("Urn Size", "Accuracy","Multiple Stakes", "Multiple Urns"),
        caption = "Coverage of the 95% confidence interval for each algorithm type and urn size.") %>%
  kable_classic(full_width = F, html_font = "helvetica")
#-------------------------------------------------------------------------------
# best urn size per change and algorithm

acc_mmse = readRDS("output/sim3_change_acc_mmse.rds")
mu_mmse = readRDS("output/sim3_change_mu_mmse.rds")
ms_mmse = readRDS("output/sim3_change_ms_mmse.rds")

algo_mmse = as.data.frame(rbind(acc_mmse, mu_mmse, ms_mmse))
algo_mmse = cbind(algo_mmse, rep(c("Accuracy", "Multiple Urns", "Multiple Updates"), 
                                 each = nrow(acc_mmse)))

colnames(algo_mmse) = c("change_amounts","urn_size","mse_200_val", "algo")
best_urn = algo_mmse %>%
            group_by(change_amounts, algo, urn_size) %>%
            summarise(mse_m = mean(mse_200_val)) %>%
            summarise(best_urn = urn_size[which.min(mse_m)]) %>%
            pivot_wider(names_from = algo, values_from = best_urn) %>%
            mutate(change_amounts = change_amounts * 501) %>%
            kable("latex",
                  booktabs = T,
                  longtable = F,
                  digits=3,
                  col.names = c("Amount of change", "Accuracy","Multiple Stakes", "Multiple Urns"),
                  caption = "Urn sizes that provided the lowest average MSE over the last 200 iterations for each algorithm type and amount of change.") %>%
            kable_classic(full_width = F, html_font = "helvetica")
