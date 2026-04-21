# new likelihood vs distance plots feb 2026

# log files for middle 50
alpha.log.mid <- parse_beast_tracelog_file("treemetricsamples/Ireland_alpha_clean.mid.log")
hcv.log.mid   <- parse_beast_tracelog_file("treemetricsamples/hcv.mid.log")
rsv.log.mid   <- parse_beast_tracelog_file("treemetricsamples/RSV2.mid.log")
delta.log.mid <- parse_beast_tracelog_file("treemetricsamples/Ireland_delta_clean.mid.log")

# trees for middle 50
tree.sample.a.mid <- read.nexus("treemetricsamples/Ireland_alpha_clean.mid.trees")
tree.sample.d.mid <- read.nexus("treemetricsamples/Ireland_delta_clean.mid.trees")
tree.sample.h.mid <- read.nexus("treemetricsamples/hcv.mid.trees")
tree.sample.r.mid <- read.nexus("treemetricsamples/RSV2.mid.trees")

library(ape)

# ----- function to rename tip labels -----
rename_tips <- function(trees, prefix = "t") {
  if (!inherits(trees, "multiPhylo")) {
    stop("Input must be a multiPhylo object")
  }
  
  n_tips <- length(trees[[1]]$tip.label)
  new_labels <- paste0(prefix, 1:n_tips)
  
  trees <- lapply(trees, function(tr) {
    tr$tip.label <- new_labels
    tr
  })
  
  class(trees) <- "multiPhylo"
  return(trees)
}




## run metrics ##
# alpha 
alpha.diff.df.mid <- distmetrics(alpha.log.mid[,"likelihood"], tree.sample.a.mid)
# delta
delta.diff.df.mid <- distmetrics(delta.log.mid[,"likelihood"], tree.sample.d.mid)
# hcv
hcv.diff.df.mid   <- distmetrics(hcv.log.mid[,"likelihood"], tree.sample.h.mid)
# rsv
rsv.diff.df.mid   <- distmetrics(rsv.log.mid[,"likelihood"], tree.sample.r.mid)  

# --  add in BHV -- #
# rsv
r.bhv.mid <- distory::dist.multiPhylo(tree.sample.r.mid)
r.bhv.mid <- as.vector(r.bhv.mid)
r.bhv.df.mid <- data.frame(like=dplyr::filter(rsv.diff.df.mid, metric=="RF")$like,
                           metric=rep("BHV"),
                           diff=r.bhv.mid)

# hcv
h.bhv.mid <- distory::dist.multiPhylo(tree.sample.h.mid)
h.bhv.mid <- as.vector(h.bhv.mid)
h.bhv.df.mid <- data.frame(like=dplyr::filter(hcv.diff.df.mid, metric=="RF")$like,
                           metric=rep("BHV"),
                           diff=h.bhv.mid)

# alpha
a.bhv.mid <- distory::dist.multiPhylo(rename_tips(tree.sample.a.mid))
a.bhv.mid <- as.vector(a.bhv.mid)
a.bhv.df.mid <- data.frame(like=dplyr::filter(alpha.diff.df.mid, metric=="RF")$like,
                           metric=rep("BHV"),
                           diff=a.bhv.mid)

# delta
d.bhv.mid <- distory::dist.multiPhylo(rename_tips(tree.sample.d.mid))
d.bhv.mid <- as.vector(d.bhv.mid)
d.bhv.df.mid <- data.frame(like=dplyr::filter(delta.diff.df.mid, metric=="RF")$like,
                           metric=rep("BHV"),
                           diff=d.bhv.mid)                           

### add in ranked NNI ###
# middle 50 trees (51:100)

# rsv
rsv.rnni <- read.csv("treemetricsamples/RSV2.all.trees.rnni.csv")
rsv.rnni.mid <- dplyr::filter(rsv.rnni, tree1 %in% 51:100 & tree2 %in% 51:100)[,3]
rsv.rnni.mid.df <- data.frame(like=dplyr::filter(rsv.diff.df.mid, metric=="RF")$like,
                              metric=rep("RNNI"),
                              diff=rsv.rnni.mid)

# hcv                        
hcv.rnni <- read.csv("treemetricsamples/hcv.all.trees.rnni.csv")
hcv.rnni.mid <- dplyr::filter(hcv.rnni, tree1 %in% 51:100 & tree2 %in% 51:100)[,3]
hcv.rnni.mid.df <- data.frame(like=dplyr::filter(hcv.diff.df.mid, metric=="RF")$like,
                              metric=rep("RNNI"),
                              diff=hcv.rnni.mid)

# alpha
a.rnni <- read.csv("treemetricsamples/Ireland_alpha_clean.all.trees.rnni.csv")
a.rnni.mid <- dplyr::filter(a.rnni, tree1 %in% 51:100 & tree2 %in% 51:100)[,3]
a.rnni.mid.df <- data.frame(like=dplyr::filter(alpha.diff.df.mid, metric=="RF")$like,
                            metric=rep("RNNI"),
                            diff=a.rnni.mid)

# delta
d.rnni <- read.csv("treemetricsamples/Ireland_delta_clean.all.trees.rnni.csv")
d.rnni.mid <- dplyr::filter(d.rnni, tree1 %in% 51:100 & tree2 %in% 51:100)[,3]
d.rnni.mid.df <- data.frame(like=dplyr::filter(delta.diff.df.mid, metric=="RF")$like,
                            metric=rep("RNNI"),
                            diff=d.rnni.mid)

# combine
delta.mid <- rbind(delta.diff.df.mid, d.rnni.mid.df, d.bhv.df.mid)
alpha.mid <- rbind(alpha.diff.df.mid, a.rnni.mid.df, a.bhv.df.mid)
rsv.mid   <- rbind(rsv.diff.df.mid, rsv.rnni.mid.df, r.bhv.df.mid)
hcv.mid   <- rbind(hcv.diff.df.mid, hcv.rnni.mid.df, h.bhv.df.mid)

library(forcats)
metric.order <- c("RF", "wRF", "JRF", "IRF", "MS", "KF", "PD", "KC", "BHV", "Trip", "Quart", "SPR", "NNI", "RNNI")

rsv.mid   <- rsv.mid   %>% mutate(metric=fct_relevel(factor(metric), metric.order)) %>% arrange(metric)
hcv.mid   <- hcv.mid   %>% mutate(metric=fct_relevel(factor(metric), metric.order)) %>% arrange(metric)
delta.mid <- delta.mid %>% mutate(metric=fct_relevel(factor(metric), metric.order)) %>% arrange(metric)
alpha.mid <- alpha.mid %>% mutate(metric=fct_relevel(factor(metric), metric.order)) %>% arrange(metric)

saveRDS(list(alpha.mid,delta.mid,rsv.mid,hcv.mid),file="middle.rds")

## --- plot --- ##
library(latex2exp)
y.lab = expression(abs(L(T[i]) - L(T[j])))
x.lab = expression(d(T[i], T[j]))

# rsv 
rsv.p.mid <- ggplot(rsv.mid, aes(x=diff, y=like)) +
  geom_point(alpha=0.075, color="#1F78B4") +
  theme_bw(base_size=12) +
  labs(y=y.lab, x=x.lab, subtitle="d. RSV") +
  facet_wrap(metric~., scales="free_x", ncol=3) +
  theme(axis.text.x=element_text(angle=35, margin=margin(5,0,0,0)),
        panel.spacing=unit(0.5, "lines"))

hcv.p.mid <- ggplot(hcv.mid, aes(x=diff, y=like, group=metric)) +
  geom_point(alpha=0.075, color="#1B9E77") +
  theme_bw(base_size=12) +
  labs(y=y.lab, x=x.lab, subtitle="c. HCV") +
  facet_wrap(metric~., scales="free_x", ncol=3) +
  theme(axis.text.x=element_text(angle=35, margin=margin(5,0,0,0)),
        panel.spacing=unit(0.5, "lines"))

alpha.p.mid <- ggplot(alpha.mid, aes(x=diff, y=like, group=metric)) +
  geom_point(alpha=0.075, color="#E41A1C") +
  theme_bw(base_size=12) +
  labs(y=y.lab, x=x.lab, subtitle="a. Alpha") +
  facet_wrap(metric~., scales="free_x", ncol=3) +
  theme(axis.text.x=element_text(angle=35, margin=margin(5,0,0,0)),
        panel.spacing=unit(0.5, "lines"))

delta.p.mid <- ggplot(delta.mid, aes(x=diff, y=like, group=metric)) +
  geom_point(alpha=0.075, color="#7570B3") +
  theme_bw(base_size=12) +
  labs(y=y.lab, x=x.lab, subtitle="b. Delta") +
  facet_wrap(metric~., scales="free_x", ncol=3) +
  theme(axis.text.x=element_text(angle=35, margin=margin(5,0,0,0)),
        panel.spacing=unit(0.5, "lines"))

library(ggpubr)
mid.plots <- ggarrange(alpha.p.mid, delta.p.mid, hcv.p.mid, rsv.p.mid)
ggsave(filename = "likedist_plot_mid.svg",
       plot = mid.plots,
       width = 10, height = 15,   
       units = "in",
       device = "svg")
