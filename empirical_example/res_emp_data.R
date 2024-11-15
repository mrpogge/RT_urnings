#-------------------------------------------------------------------------------
# dependencies

source("CRT_wrapper.R")
source("util.R")
source("empirical_example/analysis_settings.R")
library(tidyverse)
library(patchwork)

#-------------------------------------------------------------------------------
# reading the results
accuracy = readRDS("empirical_data/emp_analysis_acc_pred.rds")
multiple_stakes = readRDS("empirical_data/emp_analysis_ms_pred.rds")
multiple_urns = readRDS("empirical_data/emp_analysis_mu_pred.rds")

pred_accuracies = as.data.frame(rbind(accuracy, multiple_stakes, multiple_urns))
pred_accuracies[,4] = c(rep("Accuracy", nrow(accuracy)), 
                        rep("Multiple Updates", nrow(multiple_stakes)), 
                        rep("Multiple Urns", nrow(multiple_urns)))

colnames(pred_accuracies) = c("Fast_Indicator", "Urn_Size", "Prediction_accuracy", "Algorithm")
#-------------------------------------------------------------------------------
# fit plots
pred_plot = pred_accuracies %>%
  mutate(Fast_Indicator = recode(Fast_Indicator, "1" = "grand median", "2" = "item median", "3" = "90th quantile")) %>%
  ggplot(aes(x = Urn_Size, y = Prediction_accuracy, color = Algorithm)) +
  geom_line() + 
  facet_wrap(~Fast_Indicator) +
  labs(x = "Urn size",
       y = "Mean squared prediction error") +
  jtools::theme_apa(legend.font.size = 10)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

ggsave("figures/fig4_pred.png", pred_plot)

pred_accuracies %>%
  group_by(Algorithm, Fast_Indicator) %>% #get the minimum prediction error and the urn size where this prediction error happens
  summarise(min_pred = min(Prediction_accuracy), 
            min_urn = Urn_Size[which.min(Prediction_accuracy)],
            max_pred = max(Prediction_accuracy)) %>%
  mutate(min_urn = ifelse(Algorithm != "Accuracy", min_urn + min_urn/2, min_urn)) #add the pRT urn size to the urn size for the multiple urns algorithm

min(pred_accuracies[pred_accuracies[,"Algorithm"] == "Accuracy",2:3])
min(pred_accuracies[pred_accuracies[,"Algorithm"] == "Multiple Urns",2:3])

data = readRDS("empirical_data/emp_dat_analysis.rds")
res = urnings_combined_RT_data(data,
                               student_urn_size = 60,
                               item_urn_size = 60,
                               fast_indicator = "half_90_quantile")
p_multiple_stakes_ms = evaluate_fit_ms(res)
ggsave("figures/fig5_fit.png", p_multiple_stakes_ms)


