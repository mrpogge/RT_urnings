library(tidyverse)

# load the data
raw_dat = read.csv("empirical_data/answers.csv")


#--------------------------------------------------------------------------------
# filtering and wrangling rules

# steps: 
# 1. remove rows with missing user or missing response
# 2. select t2d item type and responses from the chech republic
# 3. create a response variable that is 1 if the item_asked is the same as the item_answered
# 4. convert the time variable to a date-time object
# 5. sort the data by time
# 6. group the data by user and create a user_answers variable that counts the number of answers per user
# 7. filter out users with less than 10 answers
# 8. group the data by item_asked and create an item_answers variable that counts the number of answers per item
# 9. filter out items with less than 100 answers
# 10. create a new variable which is one if the response time is less than the median of all response times and 0 otherwise
# 11. create a median response time for each item and create a fast2 variable which is one if the response time is less than the median of the given item and 0 otherwise 

#--------------------------------------------------------------------------------
options(digits.secs = 6)
dat = raw_dat %>%
  filter(!is.na(user)) %>%
  filter(!is.na(item_answered)) %>% 
  filter(type == "t2d" & ip_country == "CZ") %>%      
  mutate(response = 1*(item_asked == item_answered),
         time = ymd_hms(time)) %>%
  arrange(time) %>%
  group_by(item_asked) %>%
  mutate(item_answers = n()) %>%
  ungroup() %>%
  filter(item_answers >= 100) %>%
  group_by(user) %>%
  mutate(user_answers = n()) %>%
  ungroup() %>%
  filter(user_answers >= 10) %>%
  mutate(fast_1 = 1*(response_time < median(response_time)))

#-------------------------------------------------------------------------------
# create a median response time for each item
median_time = dat %>%
  group_by(item_asked) %>%
  summarise(median_time = median(response_time), 
            quantile_twice = quantile(response_time, probs = 0.9)/2) %>%
  mutate(item_id = 1:nrow(.))

student_ids = dat %>%
  group_by(user) %>%
  summarise(freq = n()) %>%
  mutate(student_id = 1:nrow(.)) %>%
  select(-freq)

#-------------------------------------------------------------------------------
# add the descriptive stat variables to the dat and then create the other two pseudoresponse time vars
# rename the items from 1:n_items
dat = dat %>%
  left_join(median_time, by = "item_asked") %>%
  left_join(student_ids, by = "user") %>%
  mutate(fast_2 = 1*(response_time < median_time),
         fast_3 = 1*(response_time < quantile_twice)) %>%
  mutate(pRT_1 = 1*(response == 1 & fast_1 ==1) + 1*(response == 0 & fast_1 == 0),
         pRT_2 = 1*(response == 1 & fast_2 ==1) + 1*(response == 0 & fast_2 == 0),
         pRT_3 = 1*(response == 1 & fast_3 ==1) + 1*(response == 0 & fast_3 == 0))

saveRDS(dat,"empirical_data/total_emp_dat.rds")

#--------------------------------------------------------------------------------
# only the necessary columns for data analysis
dat_analysis = dat %>%
  select(student_id, item_id, response, pRT_1, pRT_2, pRT_3)

saveRDS(dat_analysis, "empirical_data/emp_dat_analysis.rds")

#-------------------------------------------------------------------------------


