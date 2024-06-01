suppressPackageStartupMessages(
  {
    library(tidyverse)
    library(ape)
    library(rwty) # https://github.com/r8roy/RWTY
    library(parallel)
    library(TreeDist)
    library(phangorn)
    library(phytools)
    library(treespace)
 #   library(Kaphi)
    library(distory)
  }
)

# get evenly spaced points
get_midpoints <- function(start, end, n) {
  seq(start + 1, end - 1, length.out = n + 2)[-c(1, n + 2)] %>%
    as.integer()
}

# auxiliary function to read trees
my_read_tree <- function(x, tmp.file = tempfile(), remove.root = FALSE) {
  trees <- read.tree(text = x)
  lapply(
    trees,
    function (tree) {
      tree <- if (remove.root) tree else remove_root_edge(tree)
      tree$edge.length <- rep(1, nrow(tree$edge))
      tree$edge.length[tree$edge[,1] == (Ntip(tree) + 1)] <- 0.5
      tree
    }
  )
}

# remove root edge from tree
remove_root_edge <- function(tree) {
  tree$edge[tree$edge > Ntip(tree)] <- tree$edge[tree$edge > Ntip(tree)] - 1
  tree$edge <- tree$edge[-1, ]
  tree$Nnode <- tree$Nnode - 1
  tree
}

# discrete
discrete_metric <- function(x, y) ifelse(all.equal.phylo(x, y, use.edge.length = FALSE), 0, 1)

# NNI
NNI_lower <- function(x, y) NNIDist(x, y)[[2]]

# BHV
dist_BHV <- function(x, y) {distory::dist.multiPhylo(list(x, y))}

# wrapper function to compute ESS
compute_trees <- function(trees, burnin = 0, FOO = topological.pseudo.ess) {
  chain.trees <- mclapply(
    1:1,
    function(rep) {
      chain <- list(
        trees = trees,
        gens.per.tree = 1
      )
      class(chain) <- "rwty.chain"
      return(chain)
    }
  )
  
  # compute ESS
  message()
  message(Sys.time())
  message("PD")
  flush.console()
  res.PD <- FOO(chain.trees, treedist = "PD", burnin = burnin) %>%
    mutate(dist = "PD")
  message(Sys.time())
  message("RF")
  flush.console()
  res.RF <- FOO(chain.trees, treedist = "RF", burnin = burnin) %>%
    mutate(dist = "RF")
  message(Sys.time())
  message("JRF")
  flush.console()
  res.JRF <- FOO(chain.trees, treedist = "JRF", burnin = burnin) %>%
    mutate(dist = "JRF")
  message(Sys.time())
  message("WRF")
  flush.console()
  res.WRF <- FOO(chain.trees, treedist = "WRF", burnin = burnin) %>%
    mutate(dist = "WRF")
  message(Sys.time())
  message("KF")
  flush.console()
  res.KF <- FOO(chain.trees, treedist = "KF", burnin = burnin) %>%
    mutate(dist = "KF")
  message(Sys.time())
  message("KC")
  flush.console()
  res.KC <- FOO(chain.trees, treedist = treespace::treeDist, burnin = burnin) %>%
    mutate(dist = "KC")
  message(Sys.time())
  message("SPR")
  flush.console()
  res.SPR <- FOO(chain.trees, treedist = SPR.dist, burnin = burnin) %>%
    mutate(dist = "SPR")
  message(Sys.time())
  message("MSD")
  flush.console()
  res.MSD <- FOO(chain.trees, treedist = MatchingSplitDistance, burnin = burnin) %>%
    mutate(dist = "MSD")
  message(Sys.time())
  message("BHV")
  flush.console()
  res.BHV <- FOO(chain.trees, treedist = dist_BHV, burnin = burnin) %>%
    mutate(dist = "BHV")
  message(Sys.time())
  message("NNI")
  flush.console()
  res.NNI <- FOO(chain.trees, treedist = NNI_lower, burnin = burnin) %>%
    mutate(dist = "NNI")
  message(Sys.time())
  message("IRF")
  flush.console()
  res.IRF <- FOO(chain.trees, treedist = TreeDist::InfoRobinsonFoulds, burnin = burnin) %>%
    mutate(dist = "IRF")
  message(Sys.time())
  message("NS")
  flush.console()
  #res.NS <- FOO(chain.trees, treedist = TreeDist::NyeSimilarity, burnin = burnin) %>%
  #  mutate(dist = "NS")
  #message(Sys.time())
  #message("MSI")
  #flush.console()
  #res.MSI <- FOO(chain.trees, treedist = TreeDist::MatchingSplitInfoDistance, burnin = burnin) %>%
  #  mutate(dist = "MSI")
  #message(Sys.time())
#  message("Trip")
#  flush.console()
#  res.Trip <- FOO(chain.trees, treedist = Kaphi::Trip, burnin = burnin) %>%
#    mutate(dist = "Trip")
#  message(Sys.time())
#  message("discrete")
#  res.discrete <- FOO(chain.trees, treedist = discrete_metric, burnin = burnin) %>%
#    mutate(dist = "discrete")
#  message(Sys.time())
  
  df <- bind_rows(res.PD, res.RF, res.JRF, res.KF, res.KC, res.WRF, res.SPR, res.MSD, res.BHV, res.NNI, res.IRF, res.NS, res.MSI) #, res.discrete) #, res.Trip)
  
  return(df)
}

set.seed(0)
rwty.processors <<- 8
options(mc.cores = 8)

cat("prepare\n\n")

#all.trees <- allFurcTrees(8, to.plot = FALSE)
#all.trees <- all.trees[all.trees %>% lapply(is.binary) %>% unlist()] %>%
#  lapply(function(tree) {tree$edge.length <- rep(1, 13); tree %>% phangorn::midpoint()})

#random.trees <- rmtree(N = 1000, n = 100)

#write.tree(random.trees, "data/random.trees.nwk")

#alpha.trees <-read.nexus("data/Ireland_alpha_trim.trees") %>% `[`(., -c(1:length(.) * .1))
#delta.trees <- read.nexus("data/Ireland_delta_trim.trees")  %>% `[`(., -c(1:length(.) * .1))
alpha.trees <-read.nexus("data/Ireland_alpha_1000.trees") %>% `[`(., -c(1:length(.) * .1))
delta.trees <- read.nexus("data/Ireland_delta_1000.trees")  %>% `[`(., -c(1:length(.) * .1))
rsv.trees <- read.nexus("data/RSV2.trees")  %>% `[`(., -c(1:length(.) * .1))
hcv.trees <- read.nexus("data/hcv_coal.hcv.trees")  %>% `[`(., -c(1:length(.) * .1))

# parse hop moves
hop.trees.raw <- readLines("data/random.trees.1000.hops.1.nwk")[-1]
start.trees <- grep("^\\#", hop.trees.raw)
all.trees <- hop.trees.raw[c(start.trees + 1, length(hop.trees.raw))] %>%
  my_read_tree()
ntrees <- length(start.trees)
nhoptrees <- length(hop.trees.raw)

# random copy number
R <- rbinom(999, 5, 3 / 5)

trees <- list(
  # Alpha
  alpha = alpha.trees,
  # Delta
  delta = delta.trees,
  # RSV
  rsv = rsv.trees,
  # HCV
  hcv = hcv.trees,
  # Toy
  random.1000 = all.trees,
  # Toy-copy2
  random.copy.1 = all.trees[c(unlist(lapply(1:ntrees, rep, 2)), ntrees + 1)],
  # Toy-copyvary
  random.copy.vary = all.trees[c(unlist(lapply(1:ntrees, function(x) rep(x, R[x]+ 1))), ntrees + 1)],
  # Toy-hopfirst
  random.hopfirst.1 = hop.trees.raw[
    c(unlist(lapply(start.trees, function(x) x + 1:2)), nhoptrees)
  ] %>%
    my_read_tree(),
  # Toy-hopfirstvary
  random.hopfirst.vary = hop.trees.raw[
    c(unlist(lapply(seq_along(start.trees), function(i) start.trees[i] + 1:(R[i] + 1))), nhoptrees)
  ] %>%
    my_read_tree(),
  # Toy-hoplast
  random.hoplast.1 = hop.trees.raw[
    c(start.trees[1] + 1, unlist(lapply(start.trees[-1], function(x) x - 2:1)), nhoptrees + 1 - 2:1)
  ] %>%
    my_read_tree(),
  # Toy-hoplastvary
  random.hoplast.vary = hop.trees.raw[
    c(start.trees[1] + 1, unlist(lapply(seq_along(start.trees)[-1], function(i) start.trees[i] - (R[i - 1] + 1):1)), nhoptrees + 1 - (R[999] + 1):1)
  ] %>%
    my_read_tree(),
  # Toy-hopmid
  random.hopmid.1 = hop.trees.raw[
    c(
      unlist(lapply(1:(ntrees - 1), function(i) c(start.trees[i] + 1, get_midpoints(start.trees[i], start.trees[i + 1], 1)))),
      start.trees[ntrees] + 1,
      get_midpoints(start.trees[ntrees], nhoptrees + 1, 1),
      nhoptrees
    )
  ] %>%
    my_read_tree(),
  # Toy-hopmidvary
  random.hopmid.vary = hop.trees.raw[
    c(
      unlist(lapply(1:(ntrees - 1), function(i) c(start.trees[i] + 1, get_midpoints(start.trees[i], start.trees[i + 1], R[i])))),
      start.trees[ntrees] + 1,
      get_midpoints(start.trees[ntrees], nhoptrees + 1, R[999]),
      nhoptrees
    )
  ] %>%
    my_read_tree()
)

# filter
#trees <- trees[c("alpha", "delta")]

prefix <- c(
  all="all_8", shuffle="all_8_shuffle",
  alpha="Ireland_alpha", delta="Ireland_delta", hcv="hcv", rsv = 'RSV2',
  random="random_100", random.2="random_100_2", random.10="random_100_10", random.vary="random_100_vary",
  alpha.small="alpha_small", alpha.5 = "alpha_5", alpha.vary="alpha_vary",
  random.1000 = 'random',
  random.copy.1 = 'random_copy_1', random.copy.vary = 'random_copy_vary',
  random.hopfirst.1 = 'random_hopfirst_1', random.hopfirst.vary = 'random_hopfirst_vary',
  random.hoplast.1 = 'random_hoplast_1', random.hoplast.vary = 'random_hoplast_vary',
  random.hopmid.1 = 'random_hopmid_1', random.hopmid.vary = 'random_hopmid_vary'
)

# only run newest data sets
subset.to.run <- names(trees)
trees <- trees[subset.to.run]
prefix = prefix[subset.to.run]

lapply(seq_along(prefix), function(i) write.tree(trees[[i]], paste0("trees/", prefix[i], ".nwk"))) %>% invisible()

cat("pseudo\n\n")
pseudo.ess <- lapply(trees, compute_trees, FOO = function(x, ...) topological.pseudo.ess(x, ..., n = 100))
cat("approx\n\n")
approx.ess <- lapply(trees, compute_trees, FOO = topological.approx.ess)
cat("autocorr\n\n")
autocorr <- lapply(trees, compute_trees, FOO = topological.autocorr)

cat("write\n\n")
lapply(
  seq_along(prefix),
  function(i) {
    write_csv(pseudo.ess[[i]], paste0("stats/", prefix[i], ".psuedoess.csv"))
    write_csv(approx.ess[[i]], paste0("stats/", prefix[i], ".approxess.csv"))
    write_csv(autocorr[[i]], paste0("stats/", prefix[i], ".autocorr.csv"))
  }
) %>% invisible()