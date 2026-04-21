# new likelihood vs distance plots feb 2026
library(dplyr)
library(reshape2)
# log files for first 50
alpha.log.first <- parse_beast_tracelog_file("treemetricsamples/Ireland_alpha_clean.first.log")
hcv.log.first <- parse_beast_tracelog_file("treemetricsamples/hcv.first.log")
rsv.log.first <- parse_beast_tracelog_file("treemetricsamples/RSV2.first.log")
delta.log.first <- parse_beast_tracelog_file("treemetricsamples/Ireland_alpha_clean.first.log")

# trees for first 50
tree.sample.a.first <- read.nexus("treemetricsamples/Ireland_alpha_clean.first.trees")
tree.sample.d.first <- read.nexus("treemetricsamples/Ireland_delta_clean.first.trees")
tree.sample.h.first <- read.nexus("treemetricsamples/hcv.first.trees")
tree.sample.r.first <- read.nexus("treemetricsamples/RSV2.first.trees")

## run metrics ##
# alpha 
alpha.diff.df.first <- distmetrics(alpha.log.first[,"likelihood"],tree.sample.a.first)
# delta
delta.diff.df.first <- distmetrics(delta.log.first[,"likelihood"],tree.sample.d.first)
# hcv
hcv.diff.df.first <- distmetrics(hcv.log.first[,"likelihood"],tree.sample.h.first)
# rsv
rsv.diff.df.first <- distmetrics(rsv.log.first[,"likelihood"],tree.sample.r.first)


# --  add in BHV -- #
#rsv
r.bhv.first <- distory::dist.multiPhylo(tree.sample.r.first)
r.bhv.first <- as.vector(r.bhv.first)
r.bhv.df.first <- data.frame(like=dplyr::filter(rsv.diff.df.first, metric=="RF")$like,
                       metric=rep("BHV"),
                       diff=r.bhv.first)
#hcv
h.bhv.first <- distory::dist.multiPhylo(tree.sample.h.first)
h.bhv.first <- as.vector(h.bhv.first)
h.bhv.df.first <- data.frame(like=dplyr::filter(hcv.diff.df.first, metric=="RF")$like,
                       metric=rep("BHV"),
                       diff=h.bhv.first) 

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
#alpha
a.bhv.first <- distory::dist.multiPhylo(rename_tips(tree.sample.a.first))
a.bhv.first <- as.vector(a.bhv.first)
a.bhv.df.first <- data.frame(like=dplyr::filter(alpha.diff.df.first, metric=="RF")$like,
                       metric=rep("BHV"),
                       diff=a.bhv.first)  
#delta
d.bhv.first <- distory::dist.multiPhylo(rename_tips(tree.sample.d.first))
d.bhv.first <- as.vector(a.bhv.first)
d.bhv.df.first <- data.frame(like=dplyr::filter(delta.diff.df.first, metric=="RF")$like,
                       metric=rep("BHV"),
                       diff=d.bhv.first)  

### add in the ranked NNI ####
# first 50 trees
rsv.rnni <- read.csv("treemetricsamples/RSV2.all.trees.rnni.csv")
rsv.rnni.first <- dplyr::filter(rsv.rnni,tree1%in%c(1:50) & tree2%in%c(1:50))[,3]
rsv.rnni.first.df <- data.frame(like=dplyr::filter(rsv.diff.df.first, metric=="RF")$like,
                         metric=rep("RNNI"),
                         diff= rsv.rnni.first)
 #hcv                        
hcv.rnni <- read.csv("treemetricsamples/hcv.all.trees.rnni.csv")
hcv.rnni.first <- dplyr::filter(hcv.rnni,tree1%in%c(1:50) & tree2%in%c(1:50))[,3]
hcv.rnni.first.df <- data.frame(like=dplyr::filter(hcv.diff.df.first, metric=="RF")$like,
                         metric=rep("RNNI"),
                         diff= hcv.rnni.first)
 #alpha
a.rnni <- read.csv("treemetricsamples/Ireland_alpha_clean.all.trees.rnni.csv")
a.rnni.first <- dplyr::filter(a.rnni,tree1%in%c(1:50) & tree2%in%c(1:50))[,3]
a.rnni.first.df <- data.frame(like=dplyr::filter(alpha.diff.df.first, metric=="RF")$like,
                         metric=rep("RNNI"),
                         diff= a.rnni.first)   
#delta
d.rnni <- read.csv("treemetricsamples/Ireland_delta_clean.all.trees.rnni.csv")
d.rnni.first <- dplyr::filter(d.rnni,tree1%in%c(1:50) & tree2%in%c(1:50))[,3]
d.rnni.first.df <- data.frame(like=dplyr::filter(delta.diff.df.first, metric=="RF")$like,
                         metric=rep("RNNI"),
                         diff= d.rnni.first)                     

# combine
delta.first <- rbind(delta.diff.df.first,d.rnni.first.df,d.bhv.df.first)
alpha.first <- rbind(alpha.diff.df.first,a.rnni.first.df,a.bhv.df.first)
rsv.first <- rbind(rsv.diff.df.first,rsv.rnni.first.df,r.bhv.df.first)
hcv.first <- rbind(hcv.diff.df.first,hcv.rnni.first.df,h.bhv.df.first)

library(forcats)
metric.order <- c("RF", "wRF", "JRF", "IRF", "MS", "KF", "PD", "KC", "BHV", "Trip", "Quart", "SPR", "NNI", "RNNI")
rsv.first <- rsv.first %>% 
  mutate(metric=factor(metric)) %>% 
  mutate(metric=fct_relevel(metric,metric.order)) %>%
  arrange(metric)
hcv.first <- hcv.first %>% 
  mutate(metric=factor(metric)) %>% 
  mutate(metric=fct_relevel(metric,metric.order)) %>%
  arrange(metric)
delta.first <- delta.first %>% 
  mutate(metric=factor(metric)) %>% 
  mutate(metric=fct_relevel(metric,metric.order)) %>%
  arrange(metric)
alpha.first <- alpha.first %>% 
  mutate(metric=factor(metric)) %>% 
  mutate(metric=fct_relevel(metric,metric.order)) %>%
  arrange(metric)

saveRDS(list(alpha.first,delta.first,rsv.first,hcv.first),file="first.rds")

## --- plot --- ##

library(latex2exp)
y.lab =  expression(abs(L(T[i]) - L(T[j])))
x.lab =  expression(d(T[i],T[j]))
#rsv 
rsv.p.first <- ggplot(rsv.first ,aes(x=diff,y=like))+
  geom_point(alpha=0.075, color =  "#1F78B4") +
  theme_bw(base_size = 12) +
  labs(y=y.lab,x=x.lab,subtitle = "d. RSV")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35,margin = margin(5,0,0,0)),panel.spacing = unit(0.5, "lines")) 

hcv.p.first <- ggplot(hcv.first,aes(x=diff,y=like,group=metric))+
  geom_point(alpha=0.075, color = "#1B9E77") +
  theme_bw(base_size = 12) +
  labs(y=y.lab,x=x.lab,subtitle = "c. HCV")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35,margin = margin(5,0,0,0)),panel.spacing = unit(0.5, "lines")) 

alpha.p.first <- ggplot(alpha.first,aes(x=diff,y=like,group=metric))+
  geom_point(alpha=0.075, color = "#E41A1C") +
  theme_bw(base_size = 12) +
  labs(y=y.lab,x=x.lab,subtitle = "a. Alpha")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35,margin = margin(5,0,0,0)),panel.spacing = unit(0.5, "lines")) 

delta.p.first <- ggplot(delta.first,aes(x=diff,y=like,group=metric))+
  geom_point(alpha=0.075, color = "#7570B3") +
  theme_bw(base_size = 12) +
  labs(y=y.lab,x=x.lab,subtitle = "b. Delta")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35, margin = margin(5,0,0,0)),panel.spacing = unit(0.5, "lines")) 


library(ggpubr)
first.plots <- ggarrange(alpha.p.first,delta.p.first,hcv.p.first,rsv.p.first)
ggsave(filename = "likedist_plot_first.svg",
       plot = first.plots,
       width = 10, height = 15,   
       units = "in",
       device = "svg")
