#!/usr/bin/env Rscript

# anonymizes and cleans up results from the mturk results file

library(readr)
library(dplyr)
library(stringr)
library(rjson)
library(tidyr)
library(magrittr)
library(lubridate)
library(digest)
library(devtools)

secret <- read_lines("../../SECRET_KEY")

results <- read_tsv("fishturk.results") %>% filter(Answer.marks != "") %>% mutate(tip=str_match(annotation, ".*/(.*)/([A-Za-z]+_[A-Za-z]+)")[, 3], family=str_match(annotation, ".*/(.*)/([A-Za-z]+_[A-Za-z]+)")[, 2], Answer.marks=str_replace_all(Answer.marks, '""', '"'), start=parse_date_time(creationtime, "bdHMSY"), end=parse_date_time(assignmentsubmittime, "bdHMSY"), time_taken = end - start) %>% select(-start, -end)

get_marks <- function(blob, ...) {
  fromJSON(blob) %>% as.data.frame %>% t %>% as.data.frame %>% set_colnames(c("x", "y")) %>% mutate(mark=rownames(.)) %>% cbind(...)
}

fish_families <- rowwise(results) %>% do(get_marks(.$Answer.marks, worker=hmac(secret, .$workerId, algo="sha1"), family=.$family, tip=.$tip, time_taken=.$time_taken))

write_csv(fish_families, path = "fish_families.csv")

use_data(fish_families, overwrite = TRUE)
