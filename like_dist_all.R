# new likelihood vs distance plots feb 2026
library(ape)
library(Quartet)
### Functions 
#### Distances and Likelihood ####
##pairwise distances
distmetrics <- function(likeli,tree.sample){
RFpair=wRFpair=KFpair=PDpair=JRFpair=SPRpair=NNIpair=BHVpair=KCpair=MSpair=IRFpair=Qpair=Tpair=like=rep(NA,1225)
k <- 1
for (i in 1:(50 - 1)){
for (j in (i + 1):50) {
  RFpair[k] <- phangorn::RF.dist(tree.sample[[i]],tree.sample[[j]])
  wRFpair[k] <- phangorn::wRF.dist(tree.sample[[i]],tree.sample[[j]])
 KFpair[k] <- phangorn::KF.dist(tree.sample[[i]],tree.sample[[j]])
  PDpair[k] <- phangorn::path.dist(tree.sample[[i]],tree.sample[[j]])
 JRFpair[k] <- TreeDist::JaccardRobinsonFoulds(tree.sample[[i]],tree.sample[[j]])
 KCpair[k] <- TreeDist::KendallColijn(tree.sample[[i]],tree.sample[[j]])
 SPRpair[k] <- TreeDist::SPRDist(tree.sample[[i]],tree.sample[[j]])
  MSpair[k] <- TreeDist::MatchingSplitDistance(tree.sample[[i]],tree.sample[[j]])
 IRFpair[k] <- TreeDist::InfoRobinsonFoulds(tree.sample[[i]],tree.sample[[j]])
 # new ones
 NNIpair[k]  <- TreeDist::NNIDist(tree.sample[[i]],tree.sample[[j]])[4]
 Tpair[k] <- Quartet::TripletDistance(TQFile(tree.sample[[i]]),TQFile(tree.sample[[j]]))
 Qpair[k] <- Quartet::QuartetDistance(TQFile(tree.sample[[i]]),TQFile(tree.sample[[j]]))
 
#likelihood differences
  like[k] <- abs(likeli[i]-likeli[j])
  k <- k + 1
}}

#make a nice df
dist.data <- data.frame(RF=RFpair,
                        wRF=wRFpair,
                        JRF=JRFpair,
                        IRF=IRFpair,
                        MS=MSpair,
                        KF=KFpair,
                        PD=PDpair,
                        KC=KCpair,
                        Trip=Tpair,
                        Quart=Qpair,
                       SPR=SPRpair,
                        NNI=NNIpair)
dist.data <- melt(dist.data)
like.dist <- data.frame(rep(like,12),dist.data)
names(like.dist) <- c("like","metric","diff")
return(like.dist)
}

# ----- function to rename tip labels so BHV will work on alpha and delta -----
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

### FIRST 50 ####

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
rsv.diff.df.first <- disatmetrics(rsv.log.first[,"likelihood"],tree.sample.r.first)


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
ggsave(filename = "likedist_plot_first.png",
       plot = first.plots,
       width = 10, height = 15,   
       units = "in",
       device = "png")


### MIDDLE 50 ####

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
ggsave(filename = "likedist_plot_mid.png",
       plot = mid.plots,
       width = 10, height = 15,   
       units = "in",
       device = "png")


#### LAST 50 #####

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
# alpha
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
ggsave(filename = "likedist_plot_last.png",
       plot = combined.plot,
       width = 10, height = 15,   
       units = "in",
       device = "png")
