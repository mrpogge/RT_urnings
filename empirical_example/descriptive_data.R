library(tidyverse)

dat = readRDS("empirical_data/total_emp_dat.rds")

#-------------------------------------------------------------------------------
# prop correct

dat_prop = dat %>%
  summarise(prop_correct = mean(response),
            sample_size = n())

#n_students
length(unique(dat$student_id))
length(unique(dat$item_id))

#average response per student
resps = dat %>%
  group_by(student_id) %>%
  summarise(n = n()) %>%
  summarise(mean_response = mean(n),
            sd_response = sd(n),
            median_response = median(n),
            iqr_response = IQR(n),
            min_response = min(n),
            max_response = max(n))

resps_it = dat %>%
  group_by(item_id) %>%
  summarise(n = n()) %>%
  summarise(mean_response = mean(n),
            sd_response = sd(n),
            median_response = median(n),
            iqr_response = IQR(n),
            min_response = min(n),
            max_response = max(n))

#response times
RT_it = dat %>%
  summarise(median_rt = median(response_time),
            iqr_rt = IQR(response_time),
            min_rt = min(response_time),
            max_rt = max(response_time)) %>%
  mutate_all(function(x){x/1000})

#pseudo-RT-s
median(dat$fast_1)


