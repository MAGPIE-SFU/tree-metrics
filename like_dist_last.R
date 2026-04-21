# new likelihood vs distance plots feb 2026

# log files for last 50
alpha.log.last <- parse_beast_tracelog_file("treemetricsamples/Ireland_alpha_clean.last.log")
hcv.log.last   <- parse_beast_tracelog_file("treemetricsamples/hcv.last.log")
rsv.log.last   <- parse_beast_tracelog_file("treemetricsamples/RSV2.last.log")
delta.log.last <- parse_beast_tracelog_file("treemetricsamples/Ireland_delta_clean.last.log")

# trees for last 50
tree.sample.a.last <- read.nexus("treemetricsamples/Ireland_alpha_clean.last.trees")
tree.sample.d.last <- read.nexus("treemetricsamples/Ireland_delta_clean.last.trees")
tree.sample.h.last <- read.nexus("treemetricsamples/hcv.last.trees")
tree.sample.r.last <- read.nexus("treemetricsamples/RSV2.last.trees")

## run metrics ##
# alpha 
alpha.diff.df.last <- distmetrics(alpha.log.last[,"likelihood"], tree.sample.a.last)
# delta
delta.diff.df.last <- distmetrics(delta.log.last[,"likelihood"], tree.sample.d.last)
# hcv
hcv.diff.df.last   <- distmetrics(hcv.log.last[,"likelihood"], tree.sample.h.last)
# rsv
rsv.diff.df.last   <- distmetrics(rsv.log.last[,"likelihood"], tree.sample.r.last)

# --  add in BHV -- #
# rsv
r.bhv.last <- distory::dist.multiPhylo(tree.sample.r.last)
r.bhv.last <- as.vector(r.bhv.last)
r.bhv.df.last <- data.frame(like=dplyr::filter(rsv.diff.df.last, metric=="RF")$like,
                            metric=rep("BHV"),
                            diff=r.bhv.last)

# hcv
h.bhv.last <- distory::dist.multiPhylo(tree.sample.h.last)
h.bhv.last <- as.vector(h.bhv.last)
h.bhv.df.last <- data.frame(like=dplyr::filter(hcv.diff.df.last, metric=="RF")$like,
                            metric=rep("BHV"),
                            diff=h.bhv.last)

# alpha
a.bhv.last <- distory::dist.multiPhylo(rename_tips(tree.sample.a.last))
a.bhv.last <- as.vector(a.bhv.last)
a.bhv.df.last <- data.frame(like=dplyr::filter(alpha.diff.df.last, metric=="RF")$like,
                            metric=rep("BHV"),
                            diff=a.bhv.last)
# delta
d.bhv.last <- distory::dist.multiPhylo(rename_tips(tree.sample.d.last))
d.bhv.last <- as.vector(d.bhv.last)
d.bhv.df.last <- data.frame(like=dplyr::filter(delta.diff.df.last, metric=="RF")$like,
                            metric=rep("BHV"),
                            diff=d.bhv.last)

### ranked NNI ###
# last 50 trees (151:200)  <-- adjust if needed

# rsv
rsv.rnni <- read.csv("treemetricsamples/RSV2.all.trees.rnni.csv")
rsv.rnni.last <- dplyr::filter(rsv.rnni, tree1 %in% 101:150 & tree2 %in% 101:150)[,3]
rsv.rnni.last.df <- data.frame(like=dplyr::filter(rsv.diff.df.last, metric=="RF")$like,
                               metric=rep("RNNI"),
                               diff=rsv.rnni.last)

# hcv
hcv.rnni <- read.csv("treemetricsamples/hcv.all.trees.rnni.csv")
hcv.rnni.last <- dplyr::filter(hcv.rnni, tree1 %in% 101:150 & tree2 %in% 101:150)[,3]
hcv.rnni.last.df <- data.frame(like=dplyr::filter(hcv.diff.df.last, metric=="RF")$like,
                               metric=rep("RNNI"),
                               diff=hcv.rnni.last)

# alpha
a.rnni <- read.csv("treemetricsamples/Ireland_alpha_clean.all.trees.rnni.csv")
a.rnni.last <- dplyr::filter(a.rnni, tree1 %in% 101:150 & tree2 %in% 101:150)[,3]
a.rnni.last.df <- data.frame(like=dplyr::filter(alpha.diff.df.last, metric=="RF")$like,
                             metric=rep("RNNI"),
                             diff=a.rnni.last)

# delta
d.rnni <- read.csv("treemetricsamples/Ireland_delta_clean.all.trees.rnni.csv")
d.rnni.last <- dplyr::filter(d.rnni, tree1 %in% 101:150 & tree2 %in% 101:150)[,3]
d.rnni.last.df <- data.frame(like=dplyr::filter(delta.diff.df.last, metric=="RF")$like,
                             metric=rep("RNNI"),
                             diff=d.rnni.last)

# combine
delta.last <- rbind(delta.diff.df.last, d.rnni.last.df, d.bhv.df.last)
alpha.last <- rbind(alpha.diff.df.last, a.rnni.last.df, a.bhv.df.last)
rsv.last   <- rbind(rsv.diff.df.last, rsv.rnni.last.df, r.bhv.df.last)
hcv.last   <- rbind(hcv.diff.df.last, hcv.rnni.last.df, h.bhv.df.last)

library(forcats)
metric.order <- c("RF", "wRF", "JRF", "IRF", "MS", "KF", "PD", "KC", "BHV", "Trip", "Quart", "SPR", "NNI", "RNNI")

rsv.last   <- rsv.last   %>% mutate(metric=fct_relevel(factor(metric), metric.order)) %>% arrange(metric)
hcv.last   <- hcv.last   %>% mutate(metric=fct_relevel(factor(metric), metric.order)) %>% arrange(metric)
delta.last <- delta.last %>% mutate(metric=fct_relevel(factor(metric), metric.order)) %>% arrange(metric)
alpha.last <- alpha.last %>% mutate(metric=fct_relevel(factor(metric), metric.order)) %>% arrange(metric)

saveRDS(list(alpha.last,delta.last,rsv.last,hcv.last),file="last.rds")

## --- plot --- ##
library(latex2exp)
y.lab = expression(abs(L(T[i]) - L(T[j])))
x.lab = expression(d(T[i], T[j]))

# rsv 
rsv.p.last <- ggplot(rsv.last, aes(x=diff, y=like)) +
  geom_point(alpha=0.075, color="#1F78B4") +
  theme_bw(base_size=12) +
  labs(y=y.lab, x=x.lab, subtitle="d. RSV") +
  facet_wrap(metric~., scales="free_x", ncol=3) +
  theme(axis.text.x=element_text(angle=35, margin=margin(5,0,0,0)),
        panel.spacing=unit(0.5, "lines"))

hcv.p.last <- ggplot(hcv.last, aes(x=diff, y=like, group=metric)) +
  geom_point(alpha=0.075, color="#1B9E77") +
  theme_bw(base_size=12) +
  labs(y=y.lab, x=x.lab, subtitle="c. HCV") +
  facet_wrap(metric~., scales="free_x", ncol=3) +
  theme(axis.text.x=element_text(angle=35, margin=margin(5,0,0,0)),
        panel.spacing=unit(0.5, "lines"))

alpha.p.last <- ggplot(alpha.last, aes(x=diff, y=like, group=metric)) +
  geom_point(alpha=0.075, color="#E41A1C") +
  theme_bw(base_size=12) +
  labs(y=y.lab, x=x.lab, subtitle="a. Alpha") +
  facet_wrap(metric~., scales="free_x", ncol=3) +
  theme(axis.text.x=element_text(angle=35, margin=margin(5,0,0,0)),
        panel.spacing=unit(0.5, "lines"))

delta.p.last <- ggplot(delta.last, aes(x=diff, y=like, group=metric)) +
  geom_point(alpha=0.075, color="#7570B3") +
  theme_bw(base_size=12) +
  labs(y=y.lab, x=x.lab, subtitle="b. Delta") +
  facet_wrap(metric~., scales="free_x", ncol=3) +
  theme(axis.text.x=element_text(angle=35, margin=margin(5,0,0,0)),
        panel.spacing=unit(0.5, "lines"))

library(ggpubr)
combined.plot <- ggarrange(alpha.p.last, delta.p.last, hcv.p.last, rsv.p.last,
                           ncol = 2, nrow = 2) 

# save as SVG
ggsave(filename = "likedist_plot_last.svg",
       plot = combined.plot,
       width = 10, height = 15,   
       units = "in",
       device = "svg")

