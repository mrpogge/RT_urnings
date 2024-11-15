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
pred_accuracy_matrix_2 = matrix(0, nrow = nrow(Design), ncol = 3)
pred_accuracy_matrix_2[,1:2] = Design_mu[,1:2]
for(i in 1:nrow(Design_mu)){
  pred_container = numeric(nreps)
  pred_container_2 = numeric(nreps)
  print(Design_mu[i, ])
  fast_indicator = switch(as.character(Design_mu[i,1]),
                          "1" = "grand_median",
                          "2" = "item_median",
                          "3" = "half_90_quantile")
  for(r in 1:nreps){
    res = urnings_separate_RT_data(data,
                                   student_urn_size = Design_mu[i,2],
                                   item_urn_size = item_urn_sizes_MU[1],
                                   pRT_student_urn_size = Design_mu[i,3],
                                   pRT_item_urn_size =item_urn_sizes_MU[2],
                                   fast_indicator = fast_indicator)
    
    pred_container[r] = prediction_accuracy(res, is.multiple_urn = FALSE, level = "STUDENT")
    pred_container_2[r] = prediction_accuracy(res, is.multiple_urn = TRUE, level = "STUDENT")
  }
  pred_accuracy_matrix[i, 3] = mean(pred_container)
  pred_accuracy_matrix_2[i, 3] = mean(pred_container_2)
}

saveRDS(pred_accuracy_matrix, "empirical_data/emp_analysis_acc_pred.rds")
saveRDS(pred_accuracy_matrix_2, "empirical_data/emp_analysis_mu_pred.rds")
