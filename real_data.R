library(tidyverse)
library(coda)

get_midpoints <- function(x, n) {
  seq(1, length(x), length.out = n + 2)[-c(1, n + 2)] %>%
    as.integer()
}

get_walk <- function(x, n, size = 0.01) {
  out <- rep(x, n + 1)
  rand <- runif(n, -size, size)
  
  for (i in 1:n) {
    out[i + 1] <- out[i] + rand[i]
  }
  
  out
}

#from https://github.com/beast-dev/beast-mcmc/blob/master/src/dr/inference/trace/TraceCorrelation.java
tracer_ESS <- function (values) {
  samples <- length(values)
  maxLag <- min(samples - 1, 2000)
  m <- mean(values)
  
  gammaStat <- rep(0, maxLag)
  varStat <- 0.0
  
  lag <- 0
  
  while (lag < maxLag) {
    del1 <- values[1 : (samples - lag)] - m;
    del2 <- values[1 : (samples - lag) + lag] - m;
    
    gammaStat[lag + 1] <- sum(del1 * del2) / (samples - lag)
    
    if (lag == 0) {
      varStat <- gammaStat[1];
    } else if (lag %% 2 == 0) {
      if (gammaStat[lag] + gammaStat[lag + 1] > 0) {
        varStat <- varStat + 2.0 * (gammaStat[lag] + gammaStat[lag + 1]);
      } else {
        maxLag <- lag;
      }
    }
    
    lag <- lag + 1
  }
  
  if (gammaStat[1] == 0) {
    1;
  } else {
    samples * gammaStat[1] / varStat;
  }
}

set.seed(0)

SIZE <- 100000
FUN <- "runif"

fun <- if (FUN == "rnorm") rnorm else runif

norm.nums <- fun(SIZE)
path.nums <- lapply(
  1:(SIZE-1),
  function(i) {
    c(
      seq(
        norm.nums[i],
        norm.nums[i + 1],
        by = sign(norm.nums[i + 1] - norm.nums[i]) * (if (FALSE) 0.01 else 0.001)
      ),
      norm.nums[i + 1]
    )
  }
)
R.temp <- rbinom(SIZE, 5, 3 / 5)
R.temp2 <- rgeom(SIZE, 1 / 2)
R.temp9 <- rgeom(SIZE, 0.9)
R.temp5 <- rgeom(SIZE, 1 / 5)
R <- unlist(lapply(seq_along(path.nums), function(i) min(R.temp[i], length(path.nums[[i]]) - 1)))
R2 <- unlist(lapply(seq_along(path.nums), function(i) min(R.temp2[i], length(path.nums[[i]]) - 1)))
R9 <- unlist(lapply(seq_along(path.nums), function(i) min(R.temp2[i], length(path.nums[[i]]) - 1)))
R5 <- unlist(lapply(seq_along(path.nums), function(i) min(R.temp2[i], length(path.nums[[i]]) - 1)))
R.max <- unlist(lapply(seq_along(R.temp), function(i) max(R.temp[i], R.temp2[i], R.temp9[i], R.temp5[i], 4)))

walks <- lapply(seq_along(R.max), function(i) get_walk(norm.nums[[i]], R.max[i]))

data <- list(
  random = norm.nums,
  random.copy.1 = norm.nums[c(unlist(lapply(1:SIZE, rep, 2)))],
  random.copy.4 = norm.nums[c(unlist(lapply(1:SIZE, rep, 4)))],
#  random.copy.vary = norm.nums[c(unlist(lapply(1:SIZE, function(i) rep(i, R.temp[i] + 1))))],
  random.copy.vary2 = norm.nums[c(unlist(lapply(1:SIZE, function(i) rep(i, R.temp2[i] + 1))))],
  random.copy.vary9 = norm.nums[c(unlist(lapply(1:SIZE, function(i) rep(i, R.temp9[i] + 1))))],
  random.copy.vary5 = norm.nums[c(unlist(lapply(1:SIZE, function(i) rep(i, R.temp5[i] + 1))))],
#  random.hopfirst.1 = c(unlist(lapply(path.nums, function(x) x[1:2])), norm.nums[SIZE]),
#  random.hopfirst.vary = c(unlist(lapply(seq_along(path.nums), function(i) path.nums[[i]][1:(R[i] + 1)])), norm.nums[SIZE]),
#  random.hopfirst.vary2 = c(unlist(lapply(seq_along(path.nums), function(i) path.nums[[i]][1:(R2[i] + 1)])), norm.nums[SIZE]),
#  random.hoplast.1 = c(unlist(lapply(path.nums, function(x) x[c(1,length(x) - 1)])), norm.nums[SIZE]),
#  random.hoplast.vary = c(norm.nums[1], unlist(lapply(seq_along(path.nums), function(i) path.nums[[i]][length(path.nums[[i]]) - (R[i]:0)]))),
#  random.hoplast.vary2 = c(norm.nums[1], unlist(lapply(seq_along(path.nums), function(i) path.nums[[i]][length(path.nums[[i]]) - (R2[i]:0)]))),
#  random.hopmid.1 = c(
#    unlist(lapply(path.nums, function(x) x[c(1, get_midpoints(x, 1))])),
#    norm.nums[SIZE]
#  ),
#  random.hopmid.vary = c(
#    unlist(lapply(seq_along(path.nums), function(i) path.nums[[i]][c(1, get_midpoints(path.nums[[i]], R[i]))])),
#    norm.nums[SIZE]
#  ),
#  random.hopmid.vary2 = c(
#    unlist(lapply(seq_along(path.nums), function(i) path.nums[[i]][c(1, get_midpoints(path.nums[[i]], R2[i]))])),
#    norm.nums[SIZE]
#  ),
  random.mid = c(
    unlist(lapply(1:(SIZE - 1), function(i) c(norm.nums[i], mean(norm.nums[i:(i + 1)])))),
    norm.nums[SIZE]
  ),
  random.walk.1 = c(
    unlist(lapply(walks, function(x) x[1:2]))
  ),
random.walk.4 = c(
  unlist(lapply(walks, function(x) x[1:5]))
),
#  random.walk.vary = c(
#    unlist(lapply(seq_along(walks), function(i) walks[[i]][1:(R.temp[i] + 1)]))
#  ),
  random.walk.vary2 = c(
    unlist(lapply(seq_along(walks), function(i) walks[[i]][1:(R.temp2[i] + 1)]))
  ),
  random.walk.vary9 = c(
    unlist(lapply(seq_along(walks), function(i) walks[[i]][1:(R.temp9[i] + 1)]))
  ),
  random.walk.vary5 = c(
    unlist(lapply(seq_along(walks), function(i) walks[[i]][1:(R.temp5[i] + 1)]))
  )
)

lapply(
  names(data),
  function(i) {
    write_tsv(
      tibble(Sample=seq_along(data[[i]]), x = data[[i]]),
      paste0("realdata.", FUN, ".", as.character(SIZE), "/", i, ".txt")
    )
  }
)

coda.ess <- lapply(data, . %>% mcmc(data = .) %>% effectiveSize) %>%
  unlist()
tracer.ess <- lapply(data, tracer_ESS) %>% unlist()

tibble(dataset = names(data) %>% gsub(".var1", "", .), coda.ESS = coda.ess, tracer.ESS = tracer.ess) %>%
  write_csv(paste0("realdata.", FUN, ".", as.character(SIZE), "/ESS.csv"))
