library(ape)
library(optparse)
library(readr)
library(dplyr)

write_sample <- function(name, s, Tree, Log) {
  s.tree <- Tree[s]
  s.log <- Log[s, ]
  write.nexus(s.tree, file = paste0("sample/", data.set, ".", name, ".trees"))
  write_tsv(s.log, paste0("sample/", data.set, ".", name, ".log"))
}

args <- OptionParser() %>%
  add_option("--dataset") %>%
  parse_args()

data.set <- args$dataset

tree.file <- paste0(data.set, ".trees")
log.file <- paste0(data.set, ".log")
fine.tree.file <- paste0(data.set, ".fine.trees")
fine.log.file <- paste0(data.set, ".fine.log")

tree <- read.nexus(tree.file)
fine.tree <- read.nexus(fine.tree.file)
log <- read_tsv(log.file, comment = "#")
fine.log <- read_tsv(fine.log.file, comment = "#")

first.s <- 2:51
mid.s <- 103:152
last.s <- 204:253

random.s <- sample(seq(as.integer(length(tree) * .1), length(tree)), 250)

write_sample("first", first.s, fine.tree, fine.log)
write_sample("mid", mid.s, fine.tree, fine.log)
write_sample("last", last.s, fine.tree, fine.log)
write_sample("random", random.s, tree, log)

s.tree <- c(fine.tree[c(first.s, mid.s, last.s)], tree[random.s])
s.log <- bind_rows(fine.log[c(first.s, mid.s, last.s), ], log[random.s, ])
write.nexus(s.tree, file = paste0("sample/", data.set, ".all.trees"))
write_tsv(s.log, paste0("sample/", data.set, ".all.log"))
