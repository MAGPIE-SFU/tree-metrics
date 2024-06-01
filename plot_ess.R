library(tidyverse)

my_ggsave <- function(p, file.name, width = 9, height = 8) {
  ggsave(
    plot = p,
    filename = file.name,
    'pdf',
    colormodel = 'rgb',
    useDingbats = FALSE,
    width = width,
    height = height
  )
}

# from https://stackoverflow.com/a/64011534/6490232
`-.gg` <- function(plot, layer) {
  if (missing(layer)) {
    stop("Cannot use `-.gg()` with a single argument. Did you accidentally put - on a new line?")
  }
  if (!is.ggplot(plot)) {
    stop('Need a plot on the left side')
  }
  plot$layers = c(layer, plot$layers)
  plot
}

plot_comb_ess <- function(pseudo.ess, approx.ess, data.labels, filter = NA, facet.scales = 'free') {
  if (!is.na(filter)[1]) {
    pseudo.ess <- filter(pseudo.ess, dataset %in% filter)
    approx.ess <- filter(approx.ess, dataset %in% filter)
    data.labels <- data.labels[names(data.labels) %in% filter]
  }
  
  ggplot(
    pseudo.ess,
    aes(x = dist, y = Chain.1, colour = dataset)
  ) +
    geom_boxplot(show.legend = FALSE, outliers = FALSE, alpha = 0) +
    geom_point(
      data = approx.ess,
      aes(x = dist, y = approx.ess),
      shape = 18,
      size = 3,
      show.legend = FALSE
    ) +
    facet_wrap(. ~ dataset, scales = facet.scales, labeller = labeller(dataset = data.labels)) +
    labs(x = "Metric", y = "Estimated ESS") +
    scale_x_discrete(
      breaks = c("RF", "WRF", "KF", "PD", "JRF", "SPR", "KC", "MSD", "IRF", "BHV", "NNI"),
      limits = c("RF", "WRF", "KF", "PD", "JRF", "SPR", "KC", "MSD", "IRF", "BHV", "NNI"),
      labels = c("RF", "wRF", "KF", "PD", "JRF", "SPR", "KC", "MS", "IRF", "BHV", "NNI")
    ) +
    theme_bw() +
    lims(y = c(0, NA)) +
    theme(axis.text.x = element_text(angle = 30, hjust = .75))
}

options("readr.show_col_types" = FALSE, "readr.show_progress" = FALSE)

prefix <- c(
  alpha="Ireland_alpha", delta="Ireland_delta", hcv="hcv", rsv = 'RSV2',
  random.1000 = 'random',
  random.copy.1 = 'random_copy_1', random.copy.vary = 'random_copy_vary',
  random.hopfirst.1 = 'random_hopfirst_1', random.hopfirst.vary = 'random_hopfirst_vary',
  random.hoplast.1 = 'random_hoplast_1', random.hoplast.vary = 'random_hoplast_vary',
  random.hopmid.1 = 'random_hopmid_1', random.hopmid.vary = 'random_hopmid_vary'
)
data.name <- c(
 "Alpha", "Delta", "HCV", "RSV",
  "Toy",
  "Toy-copy2", "Toy-copyvary",
  "Toy-hopfirst", "Toy-hopfirstvary",
  "Toy-hoplast", "Toy-hoplastvary",
  "Toy-hopmid", "Toy-hopmidvary"
) %>%
  setNames(prefix)

#top.limit <- c(rep(NA, 9), 600, 600)
top.limit <- rep(NA, length(prefix))

prefix.t <- tibble(dataset = factor(prefix, levels = prefix)) %>% mutate(id = as.character(1:n()))

pseudo.ess <- paste0("stats/", prefix, ".psuedoess.csv") %>%
  lapply(read_csv) %>%
  bind_rows(.id = 'id') %>%
  right_join(prefix.t, ., by = "id") %>%
  select(-id)

approx.ess <- paste0("stats/", prefix, ".approxess.csv") %>%
  lapply(read_csv) %>%
  bind_rows(.id = 'id') %>%
  right_join(prefix.t, ., by = "id") %>%
  select(-id)

auto.corr <- paste0("stats/", prefix, ".autocorr.csv") %>%
  lapply(read_csv) %>%
  bind_rows(.id = 'id') %>%
  right_join(prefix.t, ., by = "id") %>%
  select(-id) %>%
  group_by(dataset, dist) %>%
  mutate(scaled.distance = topo.distance / max(topo.distance))

approx.ess.p <- ggplot(
  approx.ess,
  aes(x = dist, y = approx.ess, colour = dist)
) +
  geom_point(show.legend = FALSE) +
  facet_wrap(. ~ dataset, scales = 'free', labeller = labeller(dataset = data.name)) +
  labs(x = "Metric", y = "Approximate ESS") +
  theme_bw() +
  lims(y = c(0, NA)) +
  theme(axis.text.x = element_text(angle = 30, hjust = .75))

pseudo.ess.p <- ggplot(
  pseudo.ess,
  aes(x = dist, y = Chain.1, colour = dist)
) +
  geom_boxplot(show.legend = FALSE) +
  facet_wrap(. ~ dataset, scales = 'free', labeller = labeller(dataset = data.name)) +
  labs(x = "Metric", y = "Pseudo ESS") +
  theme_bw() +
  lims(y = c(0, NA)) +
  theme(axis.text.x = element_text(angle = 30, hjust = .75))
  

ess.comb.p <-plot_comb_ess(pseudo.ess, approx.ess, data.labels = data.name)

ess.comb.emp.p <- plot_comb_ess(
  pseudo.ess,
  approx.ess,
  data.labels = data.name,
  filter = c("Ireland_alpha", "Ireland_delta", "hcv", "RSV2")
) +
  scale_colour_manual(
    name = "",
    breaks = c("Ireland_alpha", "Ireland_delta", "hcv", "RSV2"),
    limits = c("Ireland_alpha", "Ireland_delta", "hcv", "RSV2"),
    labels = c("Alpha", "Delta", "HCV", "RSV"),
    values = c("#E41A1C", "#7570B3", "#1B9E77", "#1F78B4"),
    guide = 'none'
  ) +
  geom_hline(
    aes(yintercept = y),
    data = tibble(y = c(4501, 4501, 901, 1801), dataset = factor(c("Ireland_alpha", "Ireland_delta", "hcv", "RSV2"))),
    linetype = 'dashed',
    colour = 'gray',
    size = 1
  )

ess.comb.toy.p <- plot_comb_ess(
  pseudo.ess,
  approx.ess,
  data.labels = data.name,
  filter =   c('random', 'random_copy_1', 'random_copy_vary', 'random_hopfirst_1', 'random_hopfirst_vary', 'random_hoplast_1', 'random_hoplast_vary', 'random_hopmid_1', 'random_hopmid_vary'),
  facet.scales = 'fixed'
) +
  scale_colour_manual(
    name = "",
    breaks = character(0),
    limits = character(0),
    values = character(0),
    na.value = 'black',
    guide = 'none'
  ) -
  geom_hline(aes(yintercept = y), data = tibble(y = 1000), linetype = 'dotted', colour = 'gray', size = 1)

autocorr.p <- ggplot(
  auto.corr,
  aes(x = sampling.interval, y = topo.distance, colour = dist, group = dist)
) + 
  geom_line() + 
  facet_wrap(. ~ dataset, scales = 'free', labeller = labeller(dataset = data.name)) + 
  labs(x = "Sampling Interval", y = "Average Pairwise Distance", colour = "Metric") + 
  lims(y = c(0, NA)) +
  theme_bw()

autocorr.scaled.p <- ggplot(
  auto.corr,
  aes(x = sampling.interval, y = scaled.distance, colour = dist, group = dist)
) + 
  geom_line(alpha = 0.75) + 
  facet_wrap(. ~ dataset, scales = 'free', labeller = labeller(dataset = data.name)) + 
  labs(x = "Sampling Interval", y = "Average Pairwise Distance (Scaled)", colour = "Metric") + 
  lims(y = c(0, NA)) +
  theme_bw()

my_ggsave(approx.ess.p, "plots/approx_ess.pdf")
my_ggsave(pseudo.ess.p, "plots/pseudo_ess.pdf")
my_ggsave(ess.comb.p, "plots/ess_comb.pdf")
my_ggsave(ess.comb.emp.p, "plots/ess_comb_emp.pdf", width = 6, height = 5.33)
my_ggsave(ess.comb.toy.p, "plots/ess_comb_toy.pdf")
my_ggsave(autocorr.p, "plots/dist_topo.pdf")
my_ggsave(autocorr.scaled.p, "plots/dist_scaled_topo.pdf")