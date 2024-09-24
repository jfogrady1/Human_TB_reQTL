#  Need to modify covariates as qtltools requires the first column to be called "id"
library(data.table)
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
data <- fread(args[1])

ids <- colnames(data)[5:52]

colnames(data)[1:4] <- c("#chr", "start", "end", "gid")
head(data)

data$pid <- NA
data$strand <- NA
data <- data %>% select(1,2,3,53,4,54,5:52)
head(data)
colnames(data)[1] <- "chr"

data$chr <- gsub("chr", "", data$chr)

rownames(data) <- 1:nrow(data)
head(data)
colnames(data)[1] <- "#chr"
write.table(data, file = args[2], row.names = F, col.names = T, sep = "\t", quote = F)
