library(tidyverse)
library(ape)
library(rwty) # https://github.com/r8roy/RWTY
library(parallel)
library(TreeDist)
library(phangorn)
library(phytools)
library(treespace)
library(Quartet)
library(distory)

set.seed(0)

alpha.trees <-read.nexus("data/Ireland_alpha.trees") %>% `[`(., -c(1:length(.) * .1)) # was _1000
delta.trees <- read.nexus("data/Ireland_delta.trees")  %>% `[`(., -c(1:length(.) * .1)) # was _1000
rsv.trees <- read.nexus("data/RSV2.trees")  %>% `[`(., -c(1:length(.) * .1))
hcv.trees <- read.nexus("data/hcv_coal.hcv.trees")  %>% `[`(., -c(1:length(.) * .1))

# parse hop moves
#hop.trees.raw <- readLines("data/random.trees.1000.hops.1.nwk")[-1]
#start.trees <- grep("^\\#", hop.trees.raw)
#all.trees <- hop.trees.raw[c(start.trees + 1, length(hop.trees.raw))] %>%
#  my_read_tree()

all.trees <- read.tree("data/random.trees.nwk")[1:1000]
#ntrees <- length(start.trees)
N <- length(all.trees)
#nhoptrees <- length(hop.trees.raw)

cat("copy\n")
# random copy number
#R <- rbinom(999, 5, 3 / 5)
R.geom2 <- rgeom(N, 1/2)
R.geom5 <- rgeom(N, 1/5)
R.geom9 <- rgeom(N, 0.9)
R.max <- unlist(lapply(1:N, function(i) max(R.geom2[i], R.geom5[i], R.geom9[i], 4)))

cat("NNI\n")
# NNI walks
#trees.walk <- lapply(1:N, function(i) do_walk(all.trees[[i]], R.max[i]))
#trees.walkbig <- lapply(1:N, function(i) do_walk(all.trees[[i]], R.max[i], step.size = 10))
trees.nni.1 <- do_walk(all.trees[[1]], 1000, step.size = 1)
trees.nni.10 <- do_walk(all.trees[[1]], 1000, step.size = 10)
trees.nni.100 <- do_walk(all.trees[[1]], 1000, step.size = 100)
trees.nni.1000 <- do_walk(all.trees[[1]], 1000, step.size = 1000)
trees.nni.10000 <- do_walk(all.trees[[1]], 1000, step.size = 10000)

cat("SPR\n")
#SPR walks
unrooted.tree <- unroot(all.trees[[1]])
trees.spr.1 <- do_walk(unrooted.tree, 1000, step.size = 1, foo = phangorn::rSPR) %>%
  lapply(phytools::midpoint_root)
trees.spr.10 <- do_walk(unrooted.tree, 1000, step.size = 10, foo = phangorn::rSPR) %>%
  lapply(phytools::midpoint_root)
trees.spr.100 <- do_walk(unrooted.tree, 1000, step.size = 100, foo = phangorn::rSPR) %>%
  lapply(phytools::midpoint_root)
trees.spr.1000 <- do_walk(unrooted.tree, 1000, step.size = 1000, foo = phangorn::rSPR) %>%
  lapply(phytools::midpoint_root)
trees.spr.10000 <- do_walk(unrooted.tree, 1000, step.size = 10000, foo = phangorn::rSPR) %>%
  lapply(phytools::midpoint_root)

cat("collect and write\n\n")
trees <- list(
  # Alpha
  alpha = alpha.trees,
  # Delta
  delta = delta.trees,
  # HCV
  hcv = hcv.trees,
  # RSV
  rsv = rsv.trees,
  # Cratopus
  `285` = read.nexus("data/285.trees"), 
  #`285.1000`="285.1000",
#  Cratopus = read.nexus("data/"), # already done
  # Toy
  random.1000 = all.trees,
  # Toy-copy2
  random.copy.1 = all.trees[c(unlist(lapply(1:N, rep, 2)))],
  # Toy-copy5
  random.copy.4 = all.trees[c(unlist(lapply(1:N, rep, 5)))],
  # Toy-copyary2
  random.copy.vary2 = all.trees[unlist(lapply(1:N, function(x) rep(x, R.geom2[x]+ 1)))],
  # Toy-copyary5
  random.copy.vary5 = all.trees[unlist(lapply(1:N, function(x) rep(x, R.geom5[x]+ 1)))],
  # Toy-copyary9
  random.copy.vary9 = all.trees[unlist(lapply(1:N, function(x) rep(x, R.geom9[x]+ 1)))],
  # Toy-walk2
  #  random.walk.1 = lapply(trees.walk, function(tree) tree[1:2]) %>% do.call(what = c),
  # Toy-walk5
  #  random.walk.4 = lapply(trees.walk, function(tree) tree[1:5]) %>% do.call(what = c),
  # Toy-walkvary2
  #  random.walk.vary2 = lapply(1:N, function(i) trees.walk[[i]][1:(R.geom2[i] + 1)]) %>% do.call(what = c),
  # Toy-walkvary5
  #  random.walk.vary5 = lapply(1:N, function(i) trees.walk[[i]][1:(R.geom5[i] + 1)]) %>% do.call(what = c),
  # Toy-walkvary2
  #  random.walk.vary9 = lapply(1:N, function(i) trees.walk[[i]][1:(R.geom9[i] + 1)]) %>% do.call(what = c),
  # Toy-walk2
  #  random.walkbig.1 = lapply(trees.walkbig, function(tree) tree[1:2]) %>% do.call(what = c),
  # Toy-walk5
  #  random.walkbig.4 = lapply(trees.walkbig, function(tree) tree[1:5]) %>% do.call(what = c),
  # Toy-walkvary2
  #  random.walkbig.vary2 = lapply(1:N, function(i) trees.walkbig[[i]][1:(R.geom2[i] + 1)]) %>% do.call(what = c),
  # Toy-walkvary5
  #  random.walkbig.vary5 = lapply(1:N, function(i) trees.walkbig[[i]][1:(R.geom5[i] + 1)]) %>% do.call(what = c),
  # Toy-walkvary2
  #  random.walkbig.vary9 = lapply(1:N, function(i) trees.walkbig[[i]][1:(R.geom9[i] + 1)]) %>% do.call(what = c),
  nni.1 = trees.nni.1,
  nni.10 = trees.nni.10,
  nni.100 = trees.nni.100,
  nni.1000 = trees.nni.1000,
  nni.10000 = trees.nni.10000,
  
  spr.1 = trees.spr.1,
  spr.10 = trees.spr.10,
  spr.100 = trees.spr.100,
  spr.1000 = trees.spr.1000,
  spr.10000 = trees.spr.10000
  # Toy-copyvary
  #  random.copy.vary = all.trees[c(unlist(lapply(1:ntrees, function(x) rep(x, R[x]+ 1))), ntrees + 1)],
  # Toy-hopfirst
  #  random.hopfirst.1 = hop.trees.raw[
  #    c(unlist(lapply(start.trees, function(x) x + 1:2)), nhoptrees)
  #  ] %>%
  #    my_read_tree(),
  # Toy-hopfirstvary
  # random.hopfirst.vary = hop.trees.raw[
  #   c(unlist(lapply(seq_along(start.trees), function(i) start.trees[i] + 1:(R[i] + 1))), nhoptrees)
  # ] %>%
  #   my_read_tree(),
  # # Toy-hoplast
  # random.hoplast.1 = hop.trees.raw[
  #   c(start.trees[1] + 1, unlist(lapply(start.trees[-1], function(x) x - 2:1)), nhoptrees + 1 - 2:1)
  # ] %>%
  #   my_read_tree(),
  # # Toy-hoplastvary
  # random.hoplast.vary = hop.trees.raw[
  #   c(start.trees[1] + 1, unlist(lapply(seq_along(start.trees)[-1], function(i) start.trees[i] - (R[i - 1] + 1):1)), nhoptrees + 1 - (R[999] + 1):1)
  # ] %>%
  #   my_read_tree(),
  # # Toy-hopmid
  #  random.hopmid.1 = hop.trees.raw[
  #    c(
  #      unlist(lapply(1:(ntrees - 1), function(i) c(start.trees[i] + 1, get_midpoints(start.trees[i], start.trees[i + 1], 1)))),
  #      start.trees[ntrees] + 1,
  #      get_midpoints(start.trees[ntrees], nhoptrees + 1, 1),
  #      nhoptrees
  #    )
  #  ] %>%
  #    my_read_tree()
  # # Toy-hopmidvary
  # random.hopmid.vary = hop.trees.raw[
  #   c(
  #     unlist(lapply(1:(ntrees - 1), function(i) c(start.trees[i] + 1, get_midpoints(start.trees[i], start.trees[i + 1], R[i])))),
  #     start.trees[ntrees] + 1,
  #     get_midpoints(start.trees[ntrees], nhoptrees + 1, R[999]),
  #     nhoptrees
  #   )
  # ] %>%
  #   my_read_tree()
)

prefix <- c(
#  all="all_8", shuffle="all_8_shuffle",
  alpha="Ireland_alpha", delta="Ireland_delta", hcv="hcv", rsv = 'RSV2',
  `285` = "285", #`285.1000` = "285.1000", Cratopus = "Cratopus",
  random="random_100", random.2="random_100_2", random.10="random_100_10", random.vary="random_100_vary",
  alpha.small="alpha_small", alpha.5 = "alpha_5", alpha.vary="alpha_vary",
  random.1000 = 'random',
  random.copy.1 = 'random_copy_1', random.copy.4 = 'random_copy_4',
  random.copy.vary2 = 'random_copy_vary2', random.copy.vary9 = 'random_copy_vary9', random.copy.vary5 = 'random_copy_vary5',
#  random.walk.1 = 'random_walk_1', random.walk.4 = 'random_walk_4',
#  random.walk.vary2 = 'random_walk_vary2',  random.walk.vary9 = 'random_walk_vary9', random.walk.vary5 = 'random_walk_vary5',
#  random.walkbig.1 = 'random_walkbig_1', random.walkbig.4 = 'random_walkbig_4',
#  random.walkbig.vary2 = 'random_walkbig_vary2',  random.walkbig.vary9 = 'random_walkbig_vary9', random.walkbig.vary5 = 'random_walkbig_vary5',
  nni.1 = 'nni_1', nni.10 = 'nni_10', nni.100 = 'nni_100', nni.1000 = 'nni_1000', nni.10000 = 'nni_10000',
  spr.1 = 'spr_1', spr.10 = 'spr_10', spr.100 = 'spr_100', spr.1000 = 'spr_1000', spr.10000 = 'spr_10000',
  
#  random.copy.vary = 'random_copy_vary',
#  random.hopfirst.1 = 'random_hopfirst_1', random.hopfirst.vary = 'random_hopfirst_vary',
#  random.hoplast.1 = 'random_hoplast_1', random.hoplast.vary = 'random_hoplast_vary',
#  random.hopmid.1 = 'random_hopmid_1', random.hopmid.vary = 'random_hopmid_vary'
)

lapply(seq_along(prefix), function(i) write.tree(trees[[i]], paste0("trees/", prefix[i], ".nwk"))) %>% invisible()