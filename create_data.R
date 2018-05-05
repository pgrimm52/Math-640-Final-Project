########################
# File: create_data.R
# Purpose: to create dataset of Trump tweets & tweetstorms for subsequent analysis
# Output: tweet_storms.rds
########################

# Load libraries
# setwd("/Users/philipp/Google Drive/Courses/Math 640 Bayesian Statistics/final project/")
library(jsonlite)
library(dplyr)
library(lubridate)
library(stringr)
library(ggplot2)
library(dbscan)
library(mcmcplots)
library(MCMCpack)

# Build dataset (JSON source: https://github.com/bpb27/trump_tweet_data_archive)
tweets <- fromJSON("condensed_2018.json") %>% 
	bind_rows(fromJSON("condensed_2017.json")) %>% 
	bind_rows(fromJSON("condensed_2016.json")) %>% 
	bind_rows(fromJSON("condensed_2015.json")) %>%
	tbl_df()

# Create proper date-time var 
tweets$date <- ymd_hms(paste(
	substr(tweets$created_at, 27, 30), 
	substr(tweets$created_at, 5, 10), 
	substr(tweets$created_at, 12, 19)), tz="EST")

# # Plotting overall tweet count by week
# tweets %>%
# 	ggplot(aes(date)) + 
# 	geom_histogram(binwidth = 86400*7)
# 
# # Subset data
# tweets %>%
# 	filter(source %in% c("Twitter Web Client", 
# 											 "Twitter for Android", 
# 											 "Twitter for iPhone")) %>%
# 	ggplot(aes(date)) +
# 	geom_histogram(binwidth = 86400*7) + 
# 	facet_grid(source ~ .)
# 
# # Tweets by date, source (last Android day: 2017-03-07)
# tweets %>%
# 	filter(source %in% c("Twitter for Android", 
# 											 "Twitter for iPhone"),
# 				 year(date) == 2017,
# 				 month(date)<5) %>%
# 	mutate(month_day = floor_date(date, "day")) %>%
# 	group_by(source, month_day) %>%
# 	summarize(n = n()) %>%
# 	ggplot(aes(x = month_day, y = n, color = source)) + geom_line()

# Build tweet subset (use android until Mar 8, 2017)
t1 <- tweets %>%
	filter(source == "Twitter for Android",
				 date > ymd(20151108),
				 date < ymd(20170308),
				 !is_retweet,
				 !str_detect(text, '^"'))

t2 <- tweets %>%
	filter(source == "Twitter for iPhone",
				 date > ymd(20170308),
				 date < ymd(20171108),
				 !is_retweet,
				 !str_detect(text, '^"'))

tweets_subset <- bind_rows(t1, t2) %>%
	mutate(post_election = date > ymd(20161109))

# Get tweetstorms via clustering, save full dataset for text analysis
tweets_subset$cluster <- dbscan(as.matrix(tweets_subset$date), 
																eps = 1000, 
																minPts = 3)$cluster

saveRDS(tweets_subset, "tweets_subset.rds")

# # Examine text of tweetstorms
# tweets_subset %>% 
# 	filter(cluster != 0) %>%
# 	dplyr::select(date, cluster, text) %>%
# 	View()

# # Date distribution of tweetstorms
# tweets_subset %>%
# 	filter(cluster != 0) %>%
# 	group_by(cluster) %>%
# 	summarize(date = mean(date)) %>%
# 	arrange(desc(date)) %>%
# 	mutate(previous = lead(date),
# 				 days_elapsed = as.duration(date - previous)/ddays(1)) %>% # time_elapsed converted to days
# 	slice(-n()) %>%
# 	ggplot(aes(date)) + geom_histogram(bins=24)

# Output simple tweet_storms dataset (measured in days elapsed), save dataset
tweet_storms <- tweets_subset %>%
	filter(cluster != 0) %>%
	group_by(cluster) %>%
	summarize(date = mean(date)) %>%
	arrange(desc(date)) %>%
	mutate(previous = lead(date),
				 days_elapsed = as.duration(date - previous)/ddays(1),
				 post_election = date > ymd(20161109)) %>% # time_elapsed converted to days
	slice(-n()) %>% # omit last value (= NA)
	dplyr::select(days_elapsed, post_election)

saveRDS(tweet_storms, "tweet_storms.rds")

# # Shape of tweetstorms is Gamma-like (magic!)
# qplot(days_elapsed, data=tweet_storms, color=post_election, geom="density")
# hist(tweet_storms$days_elapsed, xlim=c(0, 30))
# curve(5e2*dweibull(x, 0.91, 2.87), add=TRUE, col="red")
