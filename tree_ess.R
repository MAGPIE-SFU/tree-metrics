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
    library(distory)
    library(optparse)
    library(treess)
  }
)

#SCRATCH_DIR <- "/Volumes/ram/"
#JAVACP <- "../git/BEASTLabs/build/classes:../git/beast2/dist/beast2.jar:../git/BEASTLabs/lib/colt.jar:../git/BEASTLabs/lib/nashorn-core-15.3.jar:../git/BeastFX/dist/BeastFX.jar"
SCRATCH_DIR <- "./tmpdir"
JAVACP <- "~/.beast/2.7/BEASTLabs/lib/BEASTlabs.v2.0.3.jar:~/.beast/2.7/BEAST.base/lib/BEAST.base.jar:~/.beast/2.7/BEAST.app/lib/BEAST.app.jar:/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/beast/2.7.7/lib/launcher.jar:~/git/BEASTLabs/build"
JAVA_STACK_MEM <- "256m"
JAVA_MEM <- "1g"

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
Quart_all <- function(trees, tempdir = SCRATCH_DIR, normalize = TRUE) {
  trees.file <- tempfile(tmpdir = tempdir, fileext = as.character(Sys.getpid()))
  write.tree(trees, trees.file)
  
  dist <- Quartet::AllPairsQuartetDistance(trees.file)
  
  if (normalize) {
    dist <- dist / choose(Ntip(trees[[1]]), 4)
  }
  
  try(
    {
      file.remove(trees.file)
    },
    silent = TRUE
  )
  
  dist
}

Trip_all <- function(trees, tempdir = SCRATCH_DIR, normalize = TRUE) {
  trees.file <- tempfile(tmpdir = tempdir, fileext = as.character(Sys.getpid()))
  write.tree(trees, trees.file)
  
  dist <- Quartet::AllPairsTripletDistance(trees.file)
  
  if (normalize) {
    dist <- dist / choose(Ntip(trees[[1]]), 3)
  }
  
  try(
    {
      file.remove(trees.file)
    },
    silent = TRUE
  )
  
  dist
}

# Quartet and Triplet
Quart <- function(x, y, tempdir = SCRATCH_DIR, normalize = TRUE) {
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

Trip <- function(x, y, tempdir = SCRATCH_DIR, normalize = TRUE) {
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
RNNI_dist_all <- function(trees, tempdir = SCRATCH_DIR) {
  java.call <- paste0("java -Xss", JAVA_STACK_MEM, " -Xmx", JAVA_MEM, " -cp '", JAVACP,"' hack.app.TreeComparerNexus")
  
  trees.file <- tempfile(tmpdir = tempdir, fileext = as.character(Sys.getpid()))
  write.nexus(trees, file = trees.file)
  
  out.file <- tempfile(tmpdir = tempdir, fileext = as.character(Sys.getpid()))
  
  system(paste(java.call, trees.file, out.file, sep = " "))
  
  dist <- read_csv(out.file)
  
  if (any(is.na(dist))) {
    stop("RNNI failed:", trees.file, ", ", out.file, "\n")
  }
  
  try(
    {
      file.remove(trees.file)
    },
    silent = TRUE
  )
  try(
    {
      file.remove(out.file)
    },
    silent = TRUE
  )
  
  mat <- dist %>%
    pivot_wider(names_from=tree2, values_from=rNNI) %>%
    mutate(tree1 = 0) %>%
    setNames(., 1:nrow(.)) %>%
    bind_rows(., rep(0, nrow(.)) %>% setNames(., 1:length(.))) %>%
    as.matrix() %>%
    unname()
  mat[is.na(mat)] <- 0
  mat + t(mat)
}

RNNI_dist <- function(x, y, tempdir = SCRATCH_DIR) {
  java.call <- paste0("java -Xss256m -Xmx512m -cp '", JAVACP,"' hack.app.TreeComparer")
  
  x.file <- tempfile(tmpdir = tempdir, fileext = as.character(Sys.getpid()))
  write.tree(x, x.file)
  
  y.file <- tempfile(tmpdir = tempdir, fileext = as.character(Sys.getpid()))
  write.tree(y, y.file)
  
  out.file <- tempfile(tmpdir = tempdir, fileext = as.character(Sys.getpid()))
  
  system(paste(java.call, x.file, y.file, out.file, sep = " "))
  
  dist <- readLines(out.file) %>% as.numeric()
  
  if (is.na(dist)) {
    stop("RNNI failed:", x.file, ", ", y.file, ", ", out.file, "\n")
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

#run_and_write <- function(trees, i, burnin, FOO, file.names, prefix) {
#  message()
#  message(prefix[i])

#  val <- compute_trees(trees[[i]], burnin, FOO)

# write_csv(val, file.names[i])
run_and_write <- function(trees, burnin, FOO, file.name, use.all = FALSE, prefix, append = FALSE) {
  val <- compute_trees(trees, burnin, FOO, use.all = use.all, append = append, file.name = file.name)
  
  if (!append) {
  	write_csv(val, file.name)
  }
}

do_frechet_mean <- function(trees, treedist, burnin) {
  trees <- trees[[1]]$trees
  use.all <- FALSE
  
  if (burnin > 0) {
    trees <- trees[-c(1:(burnin * length(trees)))]
  }
  
  if (class(treedist) == "character") {
    if (treedist %in% c("RNNI", "Trip", "Quart")) {
      use.all <- TRUE
    }
    
    treedist = switch(
      treedist,
      "PD" = rwty:::path.distance,
      "RF" = rwty:::rf.distance,
      "JRF" = rwty:::jrf.distance,
      "WRF" = rwty:::wrf.distance,
      "KF" = rwty:::kf.distance,
      "RNNI" = RNNI_dist_all,
      "Trip" = Trip_all,
      "Quart" = Quart_all
    )
  }
  
  if (use.all) {
    dmat <- treedist(trees)
  } else {
  dmat <- mclapply(
    seq_along(trees),
    function(i) {
      res <- rep(0, i)
      
      if (i == length(trees)) {
        res
      } else {
        c(
          res, 
          lapply(
            (i+1):length(trees),
            function(j) treedist(trees[[i]], trees[[j]])
          ) %>%
            unlist()
        )
      }
    }
  ) %>%
    unlist() %>%
    matrix(nrow = length(trees))
  
  dmat <- dmat + t(dmat)
  }
  
  tibble(frechet = treess:::frechetCorrelationESS(dmat))
}

# wrapper function to compute ESS
compute_trees <- function(trees, burnin = 0, FOO = topological.pseudo.ess, use.all, append = FALSE, file.name = NA) {
  chain.trees <- lapply(
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
  message(Sys.time())
  message("PD")
  flush.console()
  res.PD <- FOO(chain.trees, treedist = "PD", burnin = burnin) %>%
    mutate(dist = "PD")
  if (append) {
    write_csv(res.PD, file.name)
  }
  message(Sys.time())
  message("RF")
  flush.console()
  res.RF <- FOO(chain.trees, treedist = "RF", burnin = burnin) %>%
    mutate(dist = "RF")
  if (append) {
    write_csv(res.RF, file.name, append = TRUE)
  }
  message(Sys.time())
  message("JRF")
  flush.console()
  res.JRF <- FOO(chain.trees, treedist = "JRF", burnin = burnin) %>%
    mutate(dist = "JRF")
  if (append) {
    write_csv(res.JRF, file.name, append = TRUE)
  }
  message(Sys.time())
  message("WRF")
  flush.console()
  res.WRF <- FOO(chain.trees, treedist = "WRF", burnin = burnin) %>%
    mutate(dist = "WRF")
  if (append) {
    write_csv(res.WRF, file.name, append = TRUE)
  }
  message(Sys.time())
  message("KF")
  flush.console()
  res.KF <- FOO(chain.trees, treedist = "KF", burnin = burnin) %>%
    mutate(dist = "KF")
  if (append) {
    write_csv(res.KF, file.name, append = TRUE)
  }
  message(Sys.time())
  message("KC")
  flush.console()
  res.KC <- FOO(chain.trees, treedist = treespace::treeDist, burnin = burnin) %>%
    mutate(dist = "KC")
  if (append) {
    write_csv(res.KC, file.name, append = TRUE)
  }
  message(Sys.time())
  message("SPR")
  flush.console()
  res.SPR <- FOO(chain.trees, treedist = SPR.dist, burnin = burnin) %>%
    mutate(dist = "SPR")
  if (append) {
    write_csv(res.SPR, file.name, append = TRUE)
  }
  message(Sys.time())
  message("MSD")
  flush.console()
  res.MSD <- FOO(chain.trees, treedist = MatchingSplitDistance, burnin = burnin) %>%
    mutate(dist = "MSD")
  if (append) {
    write_csv(res.MSD, file.name, append = TRUE)
  }
  message(Sys.time())
  message("BHV")
  flush.console()
  res.BHV <- FOO(chain.trees, treedist = dist_BHV, burnin = burnin) %>%
    mutate(dist = "BHV")
  if (append) {
    write_csv(res.BHV, file.name, append = TRUE)
  }
  message(Sys.time())
  message("NNI")
  flush.console()
  res.NNI <- FOO(chain.trees, treedist = NNI_lower, burnin = burnin) %>%
    mutate(dist = "NNI")
  if (append) {
    write_csv(res.NNI, file.name, append = TRUE)
  }
  message(Sys.time())
  message("NNIup")
  flush.console()
  res.NNIup <- FOO(chain.trees, treedist = NNI_upper, burnin = burnin) %>%
    mutate(dist = "NNIup")
  if (append) {
    write_csv(res.NNIup, file.name, append = TRUE)
  }
  message(Sys.time())
  message("IRF")
  flush.console()
  res.IRF <- FOO(chain.trees, treedist = TreeDist::InfoRobinsonFoulds, burnin = burnin) %>%
    mutate(dist = "IRF")
  if (append) {
    write_csv(res.IRF, file.name, append = TRUE)
  }
  message(Sys.time())
  message("RNNI")
  flush.console()
  res.RNNI <- FOO(chain.trees, treedist = if (use.all) 'RNNI' else RNNI_dist, burnin = burnin) %>%
    mutate(dist = "RNNI")
  if (append) {
    write_csv(res.RNNI, file.name, append = TRUE)
  }
  message(Sys.time())
  message("Trip")
  flush.console()
  res.Trip <- FOO(chain.trees, treedist = if (use.all) 'Trip' else Trip, burnin = burnin) %>%
    mutate(dist = "Trip")
  if (append) {
    write_csv(res.Trip, file.name, append = TRUE)
  }
  message(Sys.time())
  message("Quart")
  flush.console()
  res.Quart <- FOO(chain.trees, treedist = if (use.all) 'Quart' else Quart, burnin = burnin) %>%
    mutate(dist = "Quart")
  if (append) {
    write_csv(res.Quart, file.name, append = TRUE)
  }
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

args <- OptionParser() %>%
  add_option("--prefix") %>%
  add_option("--threads", type = 'integer', default = 2) %>%
  add_option("--tmpdir", default = SCRATCH_DIR)  %>%
  add_option("--javamem", default = JAVA_MEM) %>%
  add_option("--javastackmem", default = JAVA_STACK_MEM) %>%
  parse_args()

prefix <- args$prefix
threads <- args$threads
SCRATCH_DIR <- args$tmpdir
JAVA_MEM <- args$javamem
JAVA_STACK_MEM <- args$javastackmem

set.seed(0)
rwty.processors <<- threads
options(mc.cores = threads)

#cat("prepare\n\n")

#all.trees <- allFurcTrees(8, to.plot = FALSE)
#all.trees <- all.trees[all.trees %>% lapply(is.binary) %>% unlist()] %>%
#  lapply(function(tree) {tree$edge.length <- rep(1, 13); tree %>% phangorn::midpoint()})

#random.trees <- rmtree(N = 1000, n = 100)

#write.tree(random.trees, "data/random.trees.nwk")

#alpha.trees <-read.nexus("data/Ireland_alpha_trim.trees") %>% `[`(., -c(1:length(.) * .1))
#delta.trees <- read.nexus("data/Ireland_delta_trim.trees")  %>% `[`(., -c(1:length(.) * .1))
cat("load data\n")
  
#  trees <- lapply(prefix, function(x) read.tree(paste0("trees/", x, ".nwk")))
tree <- read.tree(paste0("trees/", prefix, ".nwk"))

# filter
#trees <- list("Cratopus" = trees[[c("Cratopus")]])

# only run newest data sets
#subset.to.run <- names(trees)
#trees <- trees[subset.to.run]
#prefix = prefix[subset.to.run]

cat("pseudo\n\n")
run_and_write (
  trees = tree,
  burnin = 0,
  FOO = function(x, ...) topological.pseudo.ess(x, ..., n = 100),
  file.name = paste0("stats/", prefix, ".psuedoess.csv"),
  prefix = prefix
) %>% invisible()

cat("approx\n\n")
run_and_write(
  trees = tree,
  burnin = 0,
  FOO =topological.approx.ess,
  file.name = paste0("stats/", prefix, ".approxess.csv"),
  prefix = prefix
) %>% invisible()

cat("frechet\n\n")
run_and_write (
  trees = tree,
  burnin = 0,
  FOO = do_frechet_mean,
  file.name = paste0("stats/", prefix, ".frechetess.csv"),
  use.all = TRUE,
  prefix = prefix
) %>% invisible()
