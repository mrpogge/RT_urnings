#-------------------------------------------------------------------------------
# dependencies

source("CRT_wrapper.R")
source("util.R")
source("empirical_example/analysis_settings.R")

#-------------------------------------------------------------------------------
# read data

data = readRDS("empirical_data/emp_dat_analysis.rds")

#-------------------------------------------------------------------------------
res_list = list()
for(i in 1:nrow(Design_mu)){
  print(Design[i, ])
  fast_indicator = switch(as.character(Design_mu[i,1]),
                          "1" = "grand_median",
                          "2" = "item_median",
                          "3" = "half_90_quantile")
  res = urnings_separate_RT_data(data,
                                 student_urn_size = Design_mu[i,2],
                                 item_urn_size = Design_mu[i,2],
                                 pRT_student_urn_size = Design_mu[i,3],
                                 pRT_item_urn_size = Design_mu[i,3],
                                 fast_indicator = fast_indicator)
  res_list[[i]] = res
}

saveRDS(res_list, "empirical_data/emp_analysis_mu.rds")