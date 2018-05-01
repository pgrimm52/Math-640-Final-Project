# Preliminary exploration of dataset

# load libraries
library(jsonlite)
library(dplyr)
library(lubridate)
library(stringr)
library(ggplot2)
library(dbscan)

# load data
data <- readRDS("trump_tweets.rds") %>%
	tbl_df() %>%
	arrange(desc(date)) %>%
	filter(!is_retweet, 
				 source %in% c("Twitter for iPhone", "Twitter for Android"),
				 year(date) > 2015) %>%
	mutate(previous = lead(date),
				 hours_elapsed = (date - previous)/3600) %>%
	slice(-6285) # omit last entry with blank hours_elapsed
	
# seconds -> hours: div by 3600
# seconds -> days: div by 86400

summary(test$hours_elapsed)
quantile(test$hours_elapsed, c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975), na.rm=TRUE)

qplot(hours_elapsed, data=test, facets = year(date) ~ .)

summary(as.numeric(data$hours_elapsed)[year(data$date) > 2017])

# explore data
data %>% 
	count(source)

qplot(date, 
			data = data %>% 
				filter(source %in% c("Twitter Web Client", "Twitter for iPhone", "Twitter for Android")), 
			geom="histogram", 
			facets = source ~ . )

data %>%
	count(is_retweet)

# clustering tweetstorms
test_subset <- data %>% slice(1:1000)
test_subset <- data

test_subset$cluster <- dbscan(as.matrix(test_subset$date), eps = 1000, minPts = 3)$cluster

test_subset %>%
	count(cluster)

test_subset %>%
	filter(cluster != 0) %>%
	group_by(cluster) %>%
	summarize(date = mean(date)) %>%
	mutate(previous = lead(date),
				 time_elapsed = (date - previous)/3600) %>%
	qplot(time_elapsed, data=.)
