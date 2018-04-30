# To create dataset of Trump tweets for subsequent analysis

setwd("/Users/philipp/Google Drive/Courses/Math 640 Bayesian Statistics/final project/")

library(jsonlite)
library(dplyr)
library(lubridate)
library(stringr)
library(ggplot2)

# https://github.com/bpb27/trump_tweet_data_archive

tweets <- fromJSON("condensed_2018.json") %>% 
	bind_rows(fromJSON("condensed_2017.json")) %>% 
	bind_rows(fromJSON("condensed_2016.json")) %>% 
	bind_rows(fromJSON("condensed_2015.json"))

tweets$date <- ymd_hms(paste(
	substr(tweets$created_at, 27, 30), 
	substr(tweets$created_at, 5, 10), 
	substr(tweets$created_at, 12, 19)), tz="EST")

saveRDS(tweets, "trump_tweets.rds")
