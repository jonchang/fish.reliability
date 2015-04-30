#!/usr/bin/env Rscript

# anonymizes and cleans up results from a sqlite file

library(stringr)
library(reshape2)
library(plyr)
library(dplyr)
library(rjson)
library(tidyr)
library(lubridate)
library(digest)
library(devtools)
library(readr)

secret <- read_lines("../../SECRET_KEY")
morphos <- read_lines("../../MORPHOS")

tb <- src_sqlite(path="reliability.sqlite")

log <- tbl(tb, "log")
plan <- tbl(tb, "plan")

res <- plan %>% left_join(log, by="assignmentId") %>% collect()  %>% group_by(workerId) %>% filter(type=="completed", max_sequence-1 == max(sequence)) %>% distinct(workerId, sequence)

alldata <- res %>% mutate(role=ifelse(workerId %in% morphos, "morphologist", ifelse(str_detect(workerId, "^A"), "turker", "student")), sequence=as.numeric(sequence))

get_marks <- function(blob, name) {
	jsoned <- fromJSON(blob)
	marks <- fromJSON(jsoned$marks) %>% lapply(unlist) %>% melt() %>%
	    rename(mark = L1) %>%	filter(str_detect(mark, "^[JEOPACD]")) %>%
			group_by(mark) %>% transmute(variable=c("x", "y"), value) %>%
			ungroup() %>% spread(variable, value) %>% cbind(name)
	marks
}

get_duration <- function(blob, name) {
	jsoned <- fromJSON(blob)
	logger <- fromJSON(jsoned$logger) %>% lapply(unlist) %>% lapply(data.frame) %>% lapply(t) %>%
		do.call(plyr::rbind.fill.matrix, .) %>% data.frame(stringsAsFactors=F) %>% tbl_df() %>%
		unite(payload, -X1, -X2, sep=",") %>% rename(timestamp=X1, action=X2) %>%
		mutate(timestamp = as.POSIXct(as.numeric(timestamp)/1000, origin="1970-01-01 00:00.00 UTC")) %>%
		filter(str_detect(action, "^data")) %>%
		summarise(duration=last(timestamp)-first(timestamp)) %>% cbind(name)
	logger
}

get_urls <- function(blob, name) {
	jsoned <- fromJSON(blob)
	res <- jsoned %>% melt() %>% transmute(sequence=L1-1, url=value)
	cbind(name, res)
}

allmarks <- rowwise(alldata) %>% do(get_marks(.$result, .[c("workerId","sequence")]))
alldurations <- rowwise(alldata) %>% do(get_duration(.$result, .[c("workerId","sequence")]))
allurls <- alldata %>% distinct(workerId, list) %>% rowwise() %>% do(get_urls(.$list, .["workerId"]))

hugedata <- ungroup(alldata) %>% transmute(workerId, sequence=as.integer(sequence), role) %>% inner_join(allmarks) %>% inner_join(alldurations) %>% inner_join(allurls)

fish_reliability <- hugedata %>% select(workerId) %>% distinct(workerId) %>% rowwise() %>% do(data.frame(workerId=.$workerId, sha1mac=hmac(secret, .$workerId, algo="sha1"))) %>% left_join(hugedata) %>% mutate(fn=str_match(url, "([-a-zA-Z0-9_])+\\.jpg$")[, 1]) %>% separate(fn, into="family", remove = T, extra = "drop") %>% select(-url, -workerId)

write_csv(fish_reliability, path = "fish_reliability.csv")

use_data(fish_reliability, overwrite = TRUE)
