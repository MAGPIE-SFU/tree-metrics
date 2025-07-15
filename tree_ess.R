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
    library(Quartet)
    library(rrnni)
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

do_walk <- function(tree, steps, step.size = 1, foo = phangorn::rNNI) {
  path <- list(tree)
  
  for (i in 1:steps) {
    path[[i + 1]] <- foo(path[[i]], step.size, 1)
  }
  
  path
}

# discrete
discrete_metric <- function(x, y) ifelse(all.equal.phylo(x, y, use.edge.length = FALSE), 0, 1)

# NNI
NNI_lower <- function(x, y) NNIDist(x, y)[[2]]
NNI_upper <- function(x, y) NNIDist(x, y)[[4]]

# BHV
dist_BHV <- function(x, y) {distory::dist.multiPhylo(list(x, y))}

# Quartet and Triplet
Quart <- function(x, y, tempdir = "/Volumes/ram/", normalize = TRUE) {
  x.file <- tempfile(tmpdir = tempdir, fileext = as.character(Sys.getpid()))
  write.tree(x, x.file)
  
  y.file <- tempfile(tmpdir = tempdir, fileext = as.character(Sys.getpid()))
  write.tree(y, y.file)
  
  dist <- Quartet::QuartetDistance(x.file, y.file)
  
  if (normalize) {
    dist <- dist / choose(Ntip(x), 4)
  }
  
  try(
    {
      file.remove(x.file)
    },
    silent = TRUE
  )
  try(
    {
      file.remove(y.file)
    },
    silent = TRUE
  )
  
  dist
}

Trip <- function(x, y, tempdir = "/Volumes/ram/", normalize = TRUE) {
  x.file <- tempfile(tmpdir = tempdir, fileext = as.character(Sys.getpid()))
  write.tree(x, x.file)
  
  y.file <- tempfile(tmpdir = tempdir, fileext = as.character(Sys.getpid()))
  write.tree(y, y.file)
  
  dist <- Quartet::TripletDistance(x.file, y.file)
  
  if (normalize) {
    dist <- dist / choose(Ntip(x), 3)
  }
  
  try(
    {
      file.remove(x.file)
    },
    silent = TRUE
  )
  try(
    {
      file.remove(y.file)
    },
    silent = TRUE
  )
  
  dist
}

# ranked NNI
RNNI_dist <- function(x, y, tempdir = '/Volumes/ram/') {
  java.call <- "java -cp '../git/BEASTLabs/build/classes:../git/beast2/dist/beast2.jar:../git/BEASTLabs/lib/colt.jar:../git/BEASTLabs/lib/nashorn-core-15.3.jar:../git/BeastFX/dist/BeastFX.jar' hack.app.TreeComparer"
  
  x.file <- tempfile(tmpdir = tempdir, fileext = as.character(Sys.getpid()))
  write.tree(x, x.file)
  
  y.file <- tempfile(tmpdir = tempdir, fileext = as.character(Sys.getpid()))
  write.tree(y, y.file)
  
  out.file <- tempfile(tmpdir = tempdir, fileext = as.character(Sys.getpid()))
  
  system(paste(java.call, x.file, y.file, out.file, sep = " "))
  
  dist <- readLines(out.file) %>% as.numeric()

  try(
    {
      file.remove(x.file)
    },
    silent = TRUE
  )
  try(
    {
      file.remove(y.file)
    },
    silent = TRUE
  )
  try(
    {
      file.remove(out.file)
    },
    silent = TRUE
  )
  
  dist
}

# old ranked NNI
phylo_rnni <- function(x, y) {
  x.ranked <- phytools::force.ultrametric(x) %>%
    rrnni::as_ranked()
  y.ranked <- phytools::force.ultrametric(y) %>%
    rrnni::as_ranked()
  
  rrnni::rnni(x.ranked, y.ranked)
}

run_and_write <- function(trees, i, burnin, FOO, file.names, prefix) {
  message()
  message(prefix[i])

  val <- compute_trees(trees[[i]], burnin, FOO)
  
  write_csv(val, file.names[i])
}

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
  
 # rwty.processors <<- 1
  message(Sys.time())
  message("Trip")
  message(Sys.time())
  res.Trip <- FOO(chain.trees, treedist = Trip, burnin = burnin) %>%
    mutate(dist = "Trip")
  message(Sys.time())
  message("Quart")
  res.Quart <- FOO(chain.trees, treedist = Quart, burnin = burnin) %>%
    mutate(dist = "Quart")
  message(Sys.time())
  
 # rwty.processors <<- 8
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
  message("NNIup")
  flush.console()
  res.NNIup <- FOO(chain.trees, treedist = NNI_upper, burnin = burnin) %>%
    mutate(dist = "NNIup")
  message(Sys.time())
  message("IRF")
  flush.console()
  res.IRF <- FOO(chain.trees, treedist = TreeDist::InfoRobinsonFoulds, burnin = burnin) %>%
    mutate(dist = "IRF")
  message(Sys.time())
  message("RNNI")
  res.RNNI <- FOO(chain.trees, treedist = RNNI_dist, burnin = burnin) %>%
    mutate(dist = "RNNI")
  message(Sys.time())
#  message("NS")
#  flush.console()
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
  
  df <- bind_rows(res.PD, res.RF, res.JRF, res.KF, res.KC, res.WRF, res.SPR, res.MSD, res.BHV, res.NNI, res.NNIup, res.IRF, res.Trip, res.Quart, res.RNNI)#, res.MSI) #, res.discrete) #, res.Trip), res.NS)
  
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
cat("load data\n")
if (FALSE) {
  prefix <- c(
    alpha="Ireland_alpha", delta="Ireland_delta", hcv="hcv",rsv = 'RSV2',
    `285`="285", `285.1000`="285.1000", Cratopus="Cratopus",
    random.1000 = 'random',
    random.copy.1 = 'random_copy_1', random.copy.4 = 'random_copy_4',
    random.copy.vary2 = 'random_copy_vary2', random.copy.vary9 = 'random_copy_vary9', random.copy.vary5 = 'random_copy_vary5',
    nni.1 = 'nni_1', nni.10 = 'nni_10', nni.100 = 'nni_100', nni.1000 = 'nni_1000', nni.10000 = 'nni_10000',
    spr.1 = 'spr_1', spr.10 = 'spr_10', spr.100 = 'spr_100', spr.1000 = 'spr_1000', spr.10000 = 'spr_10000'
  )
  
  trees <- lapply(prefix, function(x) read.tree(paste0("trees/", x, ".nwk")))
} else {
  alpha.trees <-read.nexus("data/Ireland_alpha_1000.trees") %>% `[`(., -c(1:length(.) * .1))
  delta.trees <- read.nexus("data/Ireland_delta_1000.trees")  %>% `[`(., -c(1:length(.) * .1))
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
  # RSV
  rsv = rsv.trees,
  # HCV
  hcv = hcv.trees,
  `285`="285", `285.1000`="285.1000", Cratopus="Cratopus",
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
  all="all_8", shuffle="all_8_shuffle",
  alpha="Ireland_alpha", delta="Ireland_delta", hcv="hcv", rsv = 'RSV2',
  `285` = "285", `285.1000` = "285.1000", Cratopus = "Cratopus",
  random="random_100", random.2="random_100_2", random.10="random_100_10", random.vary="random_100_vary",
  alpha.small="alpha_small", alpha.5 = "alpha_5", alpha.vary="alpha_vary",
  random.1000 = 'random',
  random.copy.1 = 'random_copy_1', random.copy.4 = 'random_copy_4',
  random.copy.vary2 = 'random_copy_vary2', random.copy.vary9 = 'random_copy_vary9', random.copy.vary5 = 'random_copy_vary5',
  random.walk.1 = 'random_walk_1', random.walk.4 = 'random_walk_4',
  random.walk.vary2 = 'random_walk_vary2',  random.walk.vary9 = 'random_walk_vary9', random.walk.vary5 = 'random_walk_vary5',
  random.walkbig.1 = 'random_walkbig_1', random.walkbig.4 = 'random_walkbig_4',
  random.walkbig.vary2 = 'random_walkbig_vary2',  random.walkbig.vary9 = 'random_walkbig_vary9', random.walkbig.vary5 = 'random_walkbig_vary5',
  nni.1 = 'nni_1', nni.10 = 'nni_10', nni.100 = 'nni_100', nni.1000 = 'nni_1000', nni.10000 = 'nni_10000',
  spr.1 = 'spr_1', spr.10 = 'spr_10', spr.100 = 'spr_100', spr.1000 = 'spr_1000', spr.10000 = 'spr_10000',
  
random.copy.vary = 'random_copy_vary',
  random.hopfirst.1 = 'random_hopfirst_1', random.hopfirst.vary = 'random_hopfirst_vary',
  random.hoplast.1 = 'random_hoplast_1', random.hoplast.vary = 'random_hoplast_vary',
  random.hopmid.1 = 'random_hopmid_1', random.hopmid.vary = 'random_hopmid_vary'
)

lapply(seq_along(prefix), function(i) write.tree(trees[[i]], paste0("trees/", prefix[i], ".nwk"))) %>% invisible()
}

# filter
#trees <- list("Cratopus" = trees[[c("Cratopus")]])

# only run newest data sets
subset.to.run <- names(trees)
#trees <- trees[subset.to.run]
prefix = prefix[subset.to.run]



cat("pseudo\n\n")
lapply(
  seq_along(trees),
  run_and_write,
  trees = trees,
  burnin = 0,
  FOO = function(x, ...) topological.pseudo.ess(x, ..., n = 100),
  file.names = paste0("stats/", prefix, ".psuedoess.csv"),
  prefix = prefix
) %>% invisible()
cat("approx\n\n")
lapply(
  seq_along(trees),
  run_and_write,
  trees = trees,
  burnin = 0,
  FOO =topological.approx.ess,
  file.names = paste0("stats/", prefix, ".approxess.csv"),
  prefix = prefix
) %>% invisible()
if (FALSE) {
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
}
