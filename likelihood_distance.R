alpha.log <- parse_beast_tracelog_file("Archive/Ireland_alpha.log")
hcv.log <- parse_beast_tracelog_file("TtB_fine_trees/hcv_coal.log")
rsv.log <- parse_beast_tracelog_file("TtB_fine_trees/RSV2.log")
delta.log <- parse_beast_tracelog_file("Archive/Ireland_delta.log")

#### Distances and Likelihood ####
##pairwise distances
library(Quartet)
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

# alpha 
alpha.diff.df <- distmetrics(alpha.log[c(99951:100000),"likelihood"],tree.sample.a)
# alpha 
delta.diff.df <- distmetrics(delta.log[c(99951:100000),"likelihood"],tree.sample.d)
# hcv
hcv.diff.df <- distmetrics(hcv.log[c(99951:100000),"likelihood"],tree.sample.h)
# rsv
rsv.diff.df <- distmetrics(rsv.log[c(99951:100000),"likelihood"],tree.sample.r)

### add in BHV 
r.bhv <- distory::dist.multiPhylo(tree.sample.r)
r.bhv <- as.vector(r.bhv)
r.bhv.df <- data.frame(like=filter(rsv.diff.df, metric=="RF")$like,
                       metric=rep("BHV"),
                       diff=r.bhv) 
h.bhv <- as.vector(distory::dist.multiPhylo(tree.sample.h))
h.bhv.df <- data.frame(like=filter(hcv.diff.df, metric=="RF")$like,
                       metric=rep("BHV"),
                       diff=h.bhv) 

d.bhv <- as.vector(distory::dist.multiPhylo(d.tree))

d.bhv <- treespace(d.tree, method="BHV")

d.bhv <- filter(delta.,metric=="BHV")
a.bhv <- filter(alpha.like.dist,metric=="BHV")


a.bhv <- as.vector(distory::dist.multiPhylo(tree.sample.a))



### add in the ranked NNI -- April 9, 2025 ####
# first 50 trees
rnnifiles <- list.files(path="sampled-trees-rnni",patt= "first50")
rnnifiles <- paste0("sampled-trees-rnni/",rnnifiles)
rnni.list.first <- lapply(rnnifiles, read.csv)

rnnifiles <- list.files(path="sampled-trees-rnni",patt= "mid50")
rnnifiles <- paste0("sampled-trees-rnni/",rnnifiles)
rnni.list.mid <- lapply(rnnifiles, read.csv)

rnnifiles <- list.files(path="sampled-trees-rnni",patt= "last50")
rnnifiles <- paste0("sampled-trees-rnni/",rnnifiles)
rnni.list.last <- lapply(rnnifiles, read.csv)

rnni.alpha <- data.frame(like=filter(alpha.diff.df, metric=="RF")$like,
                         metric=rep("RNNI"),
                         diff= rnni.list[[1]]$rNNI)
rnni.delta <- data.frame(like=filter(delta.mid, metric=="RF")$like,
                         metric=rep("RNNI"),
                         diff= rnni.list.mid[[2]]$rNNI)
rnni.hcv <- data.frame(like=filter(hcv.diff.df, metric=="RF")$like,
                         metric=rep("RNNI"),
                         diff= rnni.list[[3]]$rNNI)
rnni.rsv <- data.frame(like=filter(rsv.diff.df, metric=="RF")$like,
                         metric=rep("RNNI"),
                         diff= rnni.list[[4]]$rNNI)


delta.diff.df <- rbind(delta.diff.df,rnni.delta,d.bhv)
alpha.diff.df <- rbind(alpha.diff.df,rnni.alpha,a.bhv)
rsv.diff.df <- rbind(rsv.diff.df,rnni.rsv,r.bhv.df)
rsv.diff.df  <- filter(rsv.diff.df, metric != "rNNI")
hcv.diff.df <-rbind(hcv.diff.df,rnni.hcv,h.bhv.df)

rsv.diff.df <- filter(rsv.diff.df,metric != "BHV")
rsv.diff.df <- rbind(rsv.diff.df,r.bhv.df)

rsv.diff.df <- filter(rsv.diff.df,metric != "BHV")
rsv.diff.df <- rbind(rsv.diff.df,r.bhv.df)
#reorder the metrics 
library(forcats)
metric.order <- c("RF", "wRF", "JRF", "IRF", "MS", "KF", "PD", "KC", "BHV", "Trip", "Quart", "SPR", "NNI", "RNNI")
rsv.diff.df <- rsv.diff.df %>% 
  mutate(metric=factor(metric)) %>% 
  mutate(metric=fct_relevel(metric,metric.order)) %>%
  arrange(metric)
hcv.diff.df <- hcv.diff.df %>% 
  mutate(metric=factor(metric)) %>% 
  mutate(metric=fct_relevel(metric,metric.order)) %>%
  arrange(metric)
delta.diff.df <- delta.diff.df %>% 
  mutate(metric=factor(metric)) %>% 
  mutate(metric=fct_relevel(metric,metric.order)) %>%
  arrange(metric)

alpha.diff.df <- alpha.diff.df %>% 
  mutate(metric=factor(metric)) %>% 
  mutate(metric=fct_relevel(metric,metric.order)) %>%
  arrange(metric)




#delta 
delta.p <- ggplot(d.mid,aes(x=diff,y=like,group=metric))+
  geom_point(alpha=0.075, color = "#7570B3") +
  theme_bw(base_size = 12) +
 # labs(y=y.lab,x=x.lab,subtitle = "b. Delta")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35, margin = margin(5,0,0,0)),panel.spacing = unit(0.5, "lines")) 

#rsv 
rsv.p <- ggplot(rsv.diff.df,aes(x=diff,y=like))+
  geom_point(alpha=0.075, color =  "#1F78B4") +
  theme_bw(base_size = 12) +
  labs(y=y.lab,x=x.lab,subtitle = "d. RSV")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35,margin = margin(5,0,0,0)),panel.spacing = unit(0.5, "lines")) 

#alpha 
alpha.p <- ggplot(alpha.diff.df,aes(x=diff,y=like,group=metric))+
  geom_point(alpha=0.075, color = "#E41A1C") +
  theme_bw(base_size = 12) +
  labs(y=y.lab,x=x.lab,subtitle = "a. Alpha")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35,margin = margin(5,0,0,0)),panel.spacing = unit(0.5, "lines")) 

#hcv 
hcv.p <- ggplot(hcv.diff.df,aes(x=diff,y=like,group=metric))+
  geom_point(alpha=0.075, color = "#1B9E77") +
  theme_bw(base_size = 12) +
  labs(y=y.lab,x=x.lab,subtitle = "c. HCV")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35,margin = margin(5,0,0,0)),panel.spacing = unit(0.5, "lines")) 

ggarrange(alpha.p,delta.p,hcv.p,rsv.p)


