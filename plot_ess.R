library(tidyverse)
library(scales)
library(optparse)

my_ggsave <- function(p, file.name, width = 9, height = 5.33) {
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

blank.colour.scale <- scale_colour_manual(
  name = "",
  breaks = character(0),
  limits = character(0),
  values = character(0),
  na.value = 'black',
  guide = 'none'
)

metrics <- c("RF", "WRF", "JRF", "IRF", "MSD", "KF", "PD", "KC", "BHV", "Trip", "Quart", "SPR", "NNI", "RNNI")
metric.labels <- c("RF", "wRF", "JRF", "IRF", "MS", "KF", "PD", "KC", "BHV", "Trip", "Quart", "SPR", "NNI", "RNNI")

plot_comb_ess <- function(pseudo.ess, approx.ess, data.labels, filter = NA, ymax = NA, facet.scales = 'free', ncol = NULL) {
  if (!is.na(filter)[1]) {
    pseudo.ess <- filter(pseudo.ess, dataset %in% filter) %>% mutate(dataset = factor(dataset, levels = filter))
    approx.ess <- filter(approx.ess, dataset %in% filter) %>% mutate(dataset = factor(dataset, levels = filter))
    data.labels <- data.labels[names(data.labels) %in% filter]
  }
  
  ggplot(
    pseudo.ess,
    aes(x = dist, y = Chain.1)#, colour = dataset)
  ) +
   geom_boxplot(show.legend = FALSE, outliers = FALSE, alpha = 0) +
    geom_point(
      data = approx.ess,
      aes(x = dist, y = approx.ess, colour = small),
      shape = 18,
      size = 3,
      show.legend = FALSE
    ) +
    geom_text(data = approx.ess, aes(x = dist, label = my.label), y = 100, colour = 'red') +
    facet_wrap(. ~ dataset, scales = facet.scales, labeller = labeller(dataset = data.labels), ncol = ncol) +
    labs(x = "Metric", y = "Estimated ESS") +
    scale_x_discrete(
      breaks = metrics,
      limits = metrics,
      labels = metric.labels,
      guide = guide_axis(angle = 45)
    ) +
    scale_colour_manual(breaks = c(TRUE, FALSE), limits = c(TRUE, FALSE), values = c("blue", "black")) +
    theme_bw() +
    lims(y = c(0, ymax)) +
    my.theme
}

options("readr.show_col_types" = FALSE, "readr.show_progress" = FALSE)

prefix <- c(
  alpha="Ireland_alpha", delta="Ireland_delta", hcv="hcv", rsv = 'RSV2',
  `285.1000` = "285.1000", `285` = "285", Cratopus = "Cratopus",
  random.1000 = 'random',
  random.copy.1 = 'random_copy_1', random.copy.4 = 'random_copy_4',
  random.copy.vary9 = 'random_copy_vary9', random.copy.vary2 = 'random_copy_vary2', random.copy.vary5 = 'random_copy_vary5',
#  random.walk.1 = 'random_walk_1', random.walk.4 = 'random_walk_4',
 # random.walk.vary9 = 'random_walk_vary9', random.walk.vary2 = 'random_walk_vary2', random.walk.vary5 = 'random_walk_vary5',
#  random.walkbig.1 = 'random_walkbig_1', random.walkbig.4 = 'random_walkbig_4',
#  random.walkbig.vary9 = 'random_walkbig_vary9', random.walkbig.vary2 = 'random_walkbig_vary2', random.walkbig.vary5 = 'random_walkbig_vary5',
#  random.copy.vary = 'random_copy_vary',
  nni.1 = 'nni_1', nni.10 = 'nni_10', nni.100 = 'nni_100', nni.1000 = 'nni_1000', nni.10000 = 'nni_10000',
spr.1 = 'spr_1', spr.10 = 'spr_10', spr.100 = 'spr_100', spr.1000 = 'spr_1000', sor.10000 = 'spr_10000'
#  random.hopfirst.1 = 'random_hopfirst_1', random.hopfirst.vary = 'random_hopfirst_vary',
#  random.hoplast.1 = 'random_hoplast_1', random.hoplast.vary = 'random_hoplast_vary',
#  random.hopmid.1 = 'random_hopmid_1', random.hopmid.vary = 'random_hopmid_vary'
)

theor.ess <- tibble(
  dataset = factor(
    c(
      'random',
      'random_copy_1', 'random_copy_4',
      'random_copy_vary9', 'random_copy_vary2', 'random_copy_vary5',
      'random_walk_1', 'random_walk_4',
      'random_walk_vary9', 'random_walk_vary2', 'random_walk_vary5',
      'random_walkbig_1', 'random_walkbig_4',
      'random_walkbig_vary9', 'random_walkbig_vary2', 'random_walkbig_vary5',
      'random_hopmid_1'
    ),
    levels = c(
      'random',
      'random_copy_1', 'random_copy_4',
      'random_copy_vary9', 'random_copy_vary2', 'random_copy_vary5',
      'random_walk_1', 'random_walk_4',
      'random_walk_vary9', 'random_walk_vary2', 'random_walk_vary5',
      'random_walkbig_1', 'random_walkbig_4',
      'random_walkbig_vary9', 'random_walkbig_vary2', 'random_walkbig_vary5',
      'random_hopmid_1'
    )
  ),
  y = c(1000, 1000, 1000, 909, 667, 556, 1000, 1000, 909, 667, 556, 1000, 1000, 909, 667, 556, NA)
)

#top.limit <- c(rep(NA, 9), 600, 600)
top.limit <- rep(NA, length(prefix))

prefix.t <- tibble(dataset = factor(prefix, levels = prefix)) %>% mutate(id = as.character(1:n()))

args <- OptionParser() %>%
  add_option("--pres", type = 'logical', default = FALSE, action = 'store_true') %>%
  parse_args()

pres <- args$pres

if (pres) {
  my.theme <- theme(
    text = element_text(size = 12)#,
    
#    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0, lineheight = 100, debug = TRUE)
  )
} else {
  my.theme <- theme()#axis.text.x = element_text(angle = 45, hjust = .75))
}

if (pres) {
  data.name <- c(
    "Alpha", "Delta", "HCV", "RSV",
    "Weevil (1000)", "Weevil (MrBayes)", "Weevil (BEAST)",
    "Random trees",
    "Random trees (2 copies)", "Random trees (5 copies)",
    "Random trees (Poisson p=0.9 copies)", "Random trees (Poisson p=0.5 copies)", "Random trees (Poisson p=0.2 copies)",
    # "Toy_walk-1", "Toy-walk-4",
    # "Toy-walkvary-0.9", "Toy-walkvary-0.5", "Toy-walkvary-0.2",
    # "Toy_walkbig1", "Toy-walkbig4",
    # "Toy-walkbigvary9", "Toy-walkbigvary2", "Toy-walkbigvary5",
    #  "Toy-copyvary",
    "1 NNI step", "10 NNI steps", "100 NNI steps", "1000 NNI steps", "Random-NNI-10000",
    "1 SPR step", "10 SPR steps", "100 SPR steps", "1000 SPR steps", "Random-SPR-10000"
    #  "Toy-hopfirst", "Toy-hopfirstvary",
    #  "Toy-hoplast", "Toy-hoplastvary",
    #  "Toy-hopmid", "Toy-hopmidvary"
  ) %>%
    setNames(prefix)
} else {
  data.name <- c(
    "Alpha", "Delta", "HCV", "RSV",
    "Weevil (1000)", "Weevil (MrBayes)", "Weevil (BEAST)",
    "Toy",
    "Toy-copy-2", "Toy-copy-5",
    "Toy-copyvary-0.9", "Toy-copyvary-0.5", "Toy-copyvary-0.2",
    # "Toy_walk-1", "Toy-walk-4",
    # "Toy-walkvary-0.9", "Toy-walkvary-0.5", "Toy-walkvary-0.2",
    # "Toy_walkbig1", "Toy-walkbig4",
    # "Toy-walkbigvary9", "Toy-walkbigvary2", "Toy-walkbigvary5",
    #  "Toy-copyvary",
    "Random-NNI-1", "Random-NNI-10", "Random-NNI-100", "Random-NNI-1000", "Random-NNI-10000",
    "Random-SPR-1", "Random-SPR-10", "Random-SPR-100", "Random-SPR-1000", "Random-SPR-10000"
    #  "Toy-hopfirst", "Toy-hopfirstvary",
    #  "Toy-hoplast", "Toy-hoplastvary",
    #  "Toy-hopmid", "Toy-hopmidvary"
  ) %>%
  setNames(prefix)
}

pseudo.ess <- paste0("stats/", prefix, ".psuedoess.csv") %>%
  lapply(read_csv) %>%
  bind_rows(.id = 'id') %>%
  right_join(prefix.t, ., by = "id") %>%
  select(-id) %>%
  filter(dist %in% metrics)

approx.ess <- paste0("stats/", prefix, ".approxess.csv") %>%
  lapply(read_csv) %>%
  bind_rows(.id = 'id') %>%
  right_join(prefix.t, ., by = "id") %>%
  select(-id) %>%
  filter(dist %in% metrics) %>%
  group_by(dataset) %>%
  mutate(small = approx.ess < median(approx.ess)) %>%
  ungroup() %>%
  mutate(
    my.label = ifelse(
      dist != "IRF" | ! dataset %in% c("random_copy_4", "random_copy_vary5", "spr_1", "nni_10", "nni_100"),
      NA,
      as.character(round(approx.ess, 0))
    ),
    approx.ess = ifelse(
      dist != "IRF" | ! dataset %in% c("random_copy_4", "random_copy_vary5", "spr_1", "nni_10", "nni_100"),
      approx.ess,
      NA
    )
  )

if (FALSE) {
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
  
}
ess.comb.p <-plot_comb_ess(pseudo.ess, approx.ess, data.labels = data.name)

ess.comb.weevil.p <- plot_comb_ess(
  pseudo.ess,
  approx.ess,
  data.labels = data.name,
  filter = c("285", "Cratopus")
) +
  geom_hline(
  aes(yintercept = y),
  data = tibble(y = c(30004, 2501), dataset = c("285", "Cratopus")),
  linetype = 'dashed',
  colour = 'gray',
  linewidth = 1
)

ess.comb.emp.p <- plot_comb_ess(
  pseudo.ess,
  approx.ess,
  data.labels = data.name,
  filter = c("Ireland_alpha", "Ireland_delta", "hcv", "RSV2") #, "285")
) +
#  scale_colour_manual(
#    name = "",
#    breaks = c("Ireland_alpha", "Ireland_delta", "hcv", "RSV2"),
#    limits = c("Ireland_alpha", "Ireland_delta", "hcv", "RSV2"),
#    labels = c("Alpha", "Delta", "HCV", "RSV"),
#    values = c("#E41A1C", "#7570B3", "#1B9E77", "#1F78B4"),
#    guide = 'none'
#  ) +
  geom_hline(
    aes(yintercept = y),
    data = tibble(y = c(4501, 4501, 901, 1801), dataset = factor(c("Ireland_alpha", "Ireland_delta", "hcv", "RSV2"))), #, 30004) //, "285"))),
    linetype = 'dashed',
    colour = 'gray',
    linewidth = 1
  )

toy.filter <- c(
  'random',
  'random_copy_1', 'random_copy_4',
  'random_copy_vary9', 'random_copy_vary2', 'random_copy_vary5'
  #   'random_walk_1', 'random_walk_4',
  #   'random_walk_vary9', 'random_walk_vary2', 'random_walk_vary5',
#  'random_walkbig_1', 'random_walkbig_4',
#  'random_walkbig_vary9', 'random_walkbig_vary2', 'random_walkbig_vary5'
  #    'random_hopmid_1'
)

ess.comb.toy.p <- plot_comb_ess(
  pseudo.ess,
  approx.ess,
  data.labels = data.name,
  filter = toy.filter,
  facet.scales = 'fixed'
) +
  blank.colour.scale -
  geom_hline(
    aes(yintercept = y),
    data = filter(theor.ess, dataset %in% toy.filter),
    linetype = 'solid',
    colour = 'gray',
    linewidth = 2
  )

ess.comb.nni.p <- plot_comb_ess(
  pseudo.ess,
  approx.ess,
  data.labels = data.name,
  filter =  if (pres)  {
    c('nni_1', 'nni_10', 'nni_100', 'nni_1000')
  } else {
    c('nni_1', 'nni_10', 'nni_100', 'nni_1000', 'nni_10000')
  },
  facet.scales = 'free',
  ncol = ifelse(pres, 1, 3)
) #+
#  blank.colour.scale

ess.comb.spr.p <- plot_comb_ess(
  pseudo.ess,
  approx.ess,
  data.labels = data.name,
  filter =  if (pres)  {
      c('spr_1', 'spr_10', 'spr_100', 'spr_1000')
    } else {
      c('spr_1', 'spr_10', 'spr_100', 'spr_1000', 'spr_10000')
    },
  facet.scales = 'free',
  ncol = ifelse(pres, 1, 3)
)# +
# blank.colour.scale 

if (FALSE) {
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
my_ggsave(autocorr.p, "plots/dist_topo.pdf")
my_ggsave(autocorr.scaled.p, "plots/dist_scaled_topo.pdf")
}

postfix = ifelse(pres, "_pres", "")

my_ggsave(ess.comb.weevil.p, paste0("plots/ess_comb_weevil", postfix, ".pdf"), width = 6, height = 3)
my_ggsave(ess.comb.emp.p, paste0("plots/ess_comb_emp", postfix, ".pdf"), width = 6)
my_ggsave(ess.comb.toy.p, paste0("plots/ess_comb_toy", postfix, ".pdf"))
my_ggsave(ess.comb.nni.p, paste0("plots/ess_comb_nni", postfix, ".pdf"), width = ifelse(pres, 3, 9), height = ifelse(pres, 10.66, 5.33))
my_ggsave(ess.comb.spr.p, paste0("plots/ess_comb_spr", postfix, ".pdf"), width = ifelse(pres, 3, 9), height = ifelse(pres, 10.66, 5.33))