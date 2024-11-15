#-------------------------------------------------------------------------------
# dependencies

source("CRT_wrapper.R")
source("util.R")
source("empirical_example/analysis_settings.R")

#-------------------------------------------------------------------------------
# read data

data = readRDS("empirical_data/emp_dat_analysis.rds")

#-------------------------------------------------------------------------------
pred_accuracy_matrix = matrix(0, nrow = nrow(Design), ncol = 3)
pred_accuracy_matrix[,1:2] = Design_mu[,1:2]

for(i in 1:nrow(Design)){
  print(Design[i, ])
  pred_container = numeric(nreps)
  fast_indicator = switch(as.character(Design[i,1]),
                          "1" = "grand_median",
                          "2" = "item_median",
                          "3" = "half_90_quantile")
  for(r in 1:nreps){
    res = urnings_combined_RT_data(data,
                                   student_urn_size = Design[i,2],
                                   item_urn_size = item_urn_size,
                                   fast_indicator = fast_indicator)
    pred_container[r] = prediction_accuracy(res, is.multiple_urn = FALSE, level = "STUDENT")
  }
  pred_accuracy_matrix[i, 3] = mean(pred_container)
}

saveRDS(pred_accuracy_matrix, "empirical_data/emp_analysis_ms_pred.rds")