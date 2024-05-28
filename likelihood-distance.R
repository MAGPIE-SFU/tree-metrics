# likelihood v distance plots and correlation plot

library(treespace)
library(ape)
library(phangorn) #dist funs
library(TreeDist) #dist funs
library(distory) #bhv
library(ggplot2)
library(tracerer)
library(RColorBrewer)
library(dplyr)
library(ggpubr)
library(reshape2)
library(latex2exp) #for plot labels

# the first half of this code computes the data needed for the figures in the second half
# you can just load the rda files and make the plots (skip to line 100 -- figures section)

#### DATA ####
## RSV
r.tree <- read.nexus("TtB_fine_trees/RSV2.trees")
rsv.log <- parse_beast_tracelog_file("TtB_fine_trees/RSV2.log")
#first 50
tree.sample.50.rsv <- r.tree[c(49975:50024)]
likeli.50.rsv <- rsv.log[c(49975:50024),"likelihood"]

tree.sample.last.rsv <- r.tree[c(99951:100000)]
likeli.last.rsv <- rsv.log[c(99951:100000),"likelihood"]

### HCV
h.tree <- read.nexus("TtB_fine_trees/hcv_coal-hcv.trees")
hcv.log <- parse_beast_tracelog_file("TtB_fine_trees/hcv_coal.log")
tree.sample.50.hcv <- h.tree[c(49975:50024)]
likeli.50.hcv<- hcv.log[c(49975:50024),"likelihood"]

tree.sample.last.hcv <- h.tree[c(99951:100000)]
likeli.last.hcv<- hcv.log[c(99951:100000),"likelihood"]

### alpha
a.tree <- read.nexus("Archive/Ireland_alpha.trees")
alpha.log <- parse_beast_tracelog_file("Archive/Ireland_alpha.log")
tree.sample.50.alpha <- a.tree[c(49975:50024)]
likeli.50.alpha <- alpha.log[c(49975:50024),"likelihood"]

tree.sample.last.alpha <- a.tree[c(99951:100000)]
likeli.last.alpha <- alpha.log[c(99951:100000),"likelihood"]

### delta 
d.tree <- read.nexus("Archive/Ireland_delta.trees")
delta.log <- parse_beast_tracelog_file("Archive/Ireland_delta.log")
tree.sample <- delta.tree[c(1:50)]
likeli <- delta.log[c(1:50),"likelihood"]

tree.sample.50.delta <- d.tree[c(49975:50024)]
likeli.50.delta <- alpha.log[c(49975:50024),"likelihood"]

tree.sample.last.delta <- d.tree[c(99951:100000)]
likeli.last.delta <- delta.log[c(99951:100000),"likelihood"]

#### Distances and Likelihood ####
##pairwise distances
#RFpair=wRFpair=KFpair=PDpair=JRFpair=SPRpair= NNIpair=BHVpair=KCpair =MSpair=IRFpair=like=matrix(NA,ncol = 50,nrow=50)
#IRF.HCV1=IRF.RSV1=IRF.alpha1=IRF.delta1=IRF.HCV2=IRF.RSV2=IRF.alpha2=IRF.delta2=like.delta1=like.delta2=like.alpha1=like.alpha2=like.rsv1=like.rsv2=like.hcv1=like.hcv2=matrix(NA,ncol = 50,nrow=50)
#for (j in 1:50){
#for (i in 1:50) {
##  RFpair[i,j] <- phangorn::RF.dist(tree.sample[[i]],tree.sample[[j]])
##  wRFpair[i,j] <- phangorn::wRF.dist(tree.sample[[i]],tree.sample[[j]])
##  KFpair[i,j] <- phangorn::KF.dist(tree.sample[[i]],tree.sample[[j]])
##  PDpair[i,j] <- phangorn::path.dist(tree.sample[[i]],tree.sample[[j]])
##  JRFpair[i,j] <- TreeDist::JaccardRobinsonFoulds(tree.sample[[i]],tree.sample[[j]])
##  KCpair[i,j] <- TreeDist::KendallColijn(tree.sample[[i]],tree.sample[[j]])
##  SPRpair[i,j] <- TreeDist::SPRDist(tree.sample[[i]],tree.sample[[j]])
##  MSpair[i,j] <- TreeDist::MatchingSplitDistance(tree.sample[[i]],tree.sample[[j]])
##  NNIpair[i,j]  <- as.data.frame(TreeDist::NNIDist(tree.sample[[i]],tree.sample[[j]]))[2,]
#  #BHVpair[i,j] <- distory::dist.multiPhylo(c(tree.sample[[i]],tree.sample[[j]]))
#  IRFpair[i,j] <- TreeDist::InfoRobinsonFoulds(tree.sample[[i]],tree.sample[[j]])
##likelihood differences
#  like[i,j] <- abs(likeli[i]-likeli[j])
#}}
#
#
##make a nice df
#dist.data <- data.frame(RF=as.vector(RFpair[c(1:50),]),wRF=as.vector(wRFpair[c(1:50),]),
#                        KF=as.vector(KFpair[c(1:50),]),PD=as.vector(PDpair[c(1:50),]),
#                        JRF=as.vector(JRFpair[c(1:50),]),SPR=as.vector(SPRpair[c(1:50),]),
#                        KC=as.vector(KCpair[c(1:50),]),MS=as.vector(MSpair[c(1:50),]),
#                        IRF=as.vector(IRFpair[c(1:50),]),BHV=as.vector(MSpair[c(1:50),]), 
#                        NNI=as.vector(NNIpair[c(1:50),]))
#dist.data <- melt(dist.data)
#like.dist <- data.frame(rep(as.vector(like[c(1:50),]),11),dist.data)
#names(like.dist) <- c("like","metric","diff")
#
## save the data :) 
##rsv.like.dist <- like.dist
##save(rsv.like.dist,file="likelihooddata.rda")
#
#hcv.like.dist <- like.dist
#save(rsv.like.dist,hcv.like.dist,alpha.like.dist,delta.like.dist,file="likelihooddata-first50.rda")

### Figures ####
## load the differences/likelihoods for the first 50 trees 
load("likelihooddata.rda")

# reorder the metrics so theyre all in the same order 
alpha.like.dist$metric <- factor(alpha.like.dist$metric, levels=c("RF","wRF","KF","PD","JRF","SPR","KC","MS","IRF","BHV","NNI")) 
delta.like.dist$metric <- factor(delta.like.dist$metric, levels=c("RF","wRF","KF","PD","JRF","SPR","KC","MS","IRF","BHV","NNI")) 

x.lab <- TeX("$d(T_i,T_j)$")
y.lab <- TeX("$|L(T_i)-L(T_j)|$")

rsv.p <- ggplot(rsv.like.dist,aes(x=diff,y=like))+
  geom_point(alpha=0.025, color =  "#1F78B4") +
  theme_bw(base_size = 12) +
  labs(y=y.lab,x=x.lab,subtitle = "d. RSV")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35),panel.spacing = unit(0.2, "lines")) 

hcv.p <- ggplot(hcv.like.dist,aes(x=diff,y=like,group=metric))+
  geom_point(alpha=0.025, color = "#1B9E77") +
  theme_bw(base_size = 12) +
  labs(y=y.lab,x=x.lab,subtitle = "c. HCV")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35),panel.spacing = unit(0.2, "lines")) 

alpha.p <- ggplot(alpha.like.dist,aes(x=diff,y=like,group=metric))+
  geom_point(alpha=0.025, color = "#E41A1C") +
  theme_bw(base_size = 12) +
  labs(y=y.lab,x=x.lab,subtitle = "a. Alpha")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35),panel.spacing = unit(0.2, "lines")) 

delta.p <- ggplot(delta.like.dist,aes(x=diff,y=like,group=metric))+
  geom_point(alpha=0.025, color = "#7570B3") +
  theme_bw(base_size = 12) +
  labs(y=y.lab,x=x.lab,subtitle = "b. Delta")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35),panel.spacing = unit(0.2, "lines")) 

ggarrange(alpha.p,delta.p,hcv.p,rsv.p)

ggsave(,filename = "likedist.png",dpi=300,width=12,height = 12)



### Correlation with likelihood ####
### Read in the data for the middle and last 50
load("last-like-dist.rda")
load("middle-like-dist.rda")

### make one absolutely massive df for each dataset
rsv.all <- data.frame(rbind(rsv.like.dist,rsv.like.dist.middle,rsv.like.dist.last), Sample= rep(c("First 50","Middle 50","Last 50"),each = 27500))
hcv.all <- data.frame(rbind(hcv.like.dist,hcv.like.dist.middle,hcv.like.dist.last), Sample= rep(c("First 50","Middle 50","Last 50"),each = 27500))
alpha.all <- data.frame(rbind(alpha.like.dist,alpha.like.dist.middle,alpha.like.dist.last), Sample= rep(c("First 50","Middle 50","Last 50"),each = 27500))
delta.all <- data.frame(rbind(delta.like.dist,delta.like.dist.middle,delta.like.dist.last), Sample= rep(c("First 50","Middle 50","Last 50"),each = 27500))


rsv.cor <- data.frame(rsv.all %>% 
                        group_by(metric,Sample) %>% 
                        summarise(correlation = cor(like, diff,method="spearman"),.groups="keep"),data=rep("RSV"))
hcv.cor <- data.frame(hcv.all %>% 
                        group_by(metric,Sample) %>% 
                        summarise(correlation = cor(like, diff,method="spearman")),data=rep("HCV"))
alpha.cor <- data.frame(alpha.all %>% 
                        group_by(metric,Sample) %>% 
                        summarise(correlation = cor(like, diff,method="spearman")),data=rep("Alpha"))
delta.cor <- data.frame(delta.all %>% 
                        group_by(metric,Sample) %>% 
                        summarise(correlation = cor(like, diff,method="spearman")),data=rep("Delta"))

cordata <- rbind(rsv.cor,hcv.cor,alpha.cor,delta.cor)
levels(as.factor(cordata$Sample))
cordata$Sample <- factor(cordata$Sample, levels=c("First 50", "Middle 50","Last 50")) 
#create color spectrum
blues11 <- colorRampPalette(c("#E41A1C","pink","#7570B3","#1F78B4","#1B9E77"))(11)
# correlation plot 
ggplot(cordata,aes(x=data,y=correlation,fill=metric,color=metric))+
  geom_bar(stat = "identity",position = "dodge",alpha=0.7, width=0.75) +
  scale_fill_manual(values = blues11) +
  scale_color_manual(values = blues11) +
  theme_bw(base_size = 12) +
  theme(legend.title = element_blank()) +
  labs(y="Correlation",x="Data") +
  facet_grid(Sample~. )

ggsave(filename = "like-corr-all.png",dpi=300,width = 6,height = 6)

#### Scatterplots for middle and last samples ####

rsv.middle <- ggplot(rsv.like.dist.middle,aes(x=diff,y=like))+
  geom_point(alpha=0.025, color =  "#1F78B4") +
  theme_bw(base_size = 12) +
  labs(y=y.lab,x=x.lab,subtitle = "d. RSV")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35),panel.spacing = unit(0.2, "lines")) 

hcv.middle <- ggplot(hcv.like.dist.middle,aes(x=diff,y=like,group=metric))+
  geom_point(alpha=0.025, color = "#1B9E77") +
  theme_bw(base_size = 12) +
  labs(y=y.lab,x=x.lab,subtitle = "c. HCV")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35),panel.spacing = unit(0.2, "lines")) 

alpha.middle <- ggplot(alpha.like.dist.middle,aes(x=diff,y=like,group=metric))+
  geom_point(alpha=0.025, color = "#E41A1C") +
  theme_bw(base_size = 12) +
  labs(y=y.lab,x=x.lab,subtitle = "a. Alpha")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35),panel.spacing = unit(0.2, "lines")) 

delta.middle <- ggplot(delta.like.dist.middle,aes(x=diff,y=like,group=metric))+
  geom_point(alpha=0.025, color = "#7570B3") +
  theme_bw(base_size = 12) +
  labs(y=y.lab,x=x.lab,subtitle = "b. Delta")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35),panel.spacing = unit(0.2, "lines")) 

ggarrange(alpha.middle,delta.middle,hcv.middle,rsv.middle)

ggsave(,filename = "likedist-middle.png",dpi=300,width=12,height = 12)

rsv.last <- ggplot(rsv.like.dist.last,aes(x=diff,y=like))+
  geom_point(alpha=0.025, color =  "#1F78B4") +
  theme_bw(base_size = 12) +
  labs(y=y.lab,x=x.lab,subtitle = "d. RSV")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35),panel.spacing = unit(0.2, "lines")) 

hcv.last <- ggplot(hcv.like.dist.last,aes(x=diff,y=like,group=metric))+
  geom_point(alpha=0.025, color = "#1B9E77") +
  theme_bw(base_size = 12) +
  labs(y=y.lab,x=x.lab,subtitle = "c. HCV")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35),panel.spacing = unit(0.2, "lines")) 

alpha.last <- ggplot(alpha.like.dist.last,aes(x=diff,y=like,group=metric))+
  geom_point(alpha=0.025, color = "#E41A1C") +
  theme_bw(base_size = 12) +
  labs(y=y.lab,x=x.lab,subtitle = "a. Alpha")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35),panel.spacing = unit(0.2, "lines")) 

delta.last <- ggplot(delta.like.dist.last,aes(x=diff,y=like,group=metric))+
  geom_point(alpha=0.025, color = "#7570B3") +
  theme_bw(base_size = 12) +
  labs(y=y.lab,x=x.lab,subtitle = "b. Delta")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35),panel.spacing = unit(0.2, "lines")) 

ggarrange(alpha.last,delta.last,hcv.last,rsv.last)

ggsave(,filename = "likedist-last.png",dpi=300,width=12,height = 12)


#### Histograms #######
rsv.last.hist <- ggplot(rsv.like.dist.last,aes(x=diff))+
  geom_histogram(alpha=0.4, color =  "#1F78B4", fill = "#1F78B4") +
  theme_bw(base_size = 12) +
  labs(y="Count",x=x.lab,subtitle = "d. RSV")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35),panel.spacing = unit(0.2, "lines")) 

hcv.last.hist <- ggplot(hcv.like.dist.last,aes(x=diff,group=metric))+
  geom_histogram(alpha=0.4, color = "#1B9E77", fill= "#1B9E77") +
  theme_bw(base_size = 12) +
  labs(y="Count",x=x.lab,subtitle = "c. HCV")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35),panel.spacing = unit(0.2, "lines")) 

alpha.last.hist <- ggplot(alpha.like.dist.last,aes(x=diff,group=metric))+
  geom_histogram(alpha=0.4, color = "#E41A1C", fill="#E41A1C") +
  theme_bw(base_size = 12) +
  labs(y="Count",x=x.lab,subtitle = "a. Alpha")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35),panel.spacing = unit(0.2, "lines")) 

delta.last.hist <- ggplot(delta.like.dist.last,aes(x=diff,group=metric))+
  geom_histogram(alpha=0.4, color = "#7570B3", fill="#7570B3") +
  theme_bw(base_size = 12) +
  labs(y=y.lab,x=x.lab,subtitle = "b. Delta")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35),panel.spacing = unit(0.2, "lines")) 

ggarrange(alpha.last.hist,delta.last.hist,hcv.last.hist,rsv.last.hist)

ggsave(,filename = "dist-last-hist.png",dpi=300,width=12,height = 12)


rsv.middle.hist <- ggplot(rsv.like.dist.middle,aes(x=diff))+
  geom_histogram(alpha=0.4, color =  "#1F78B4", fill = "#1F78B4") +
  theme_bw(base_size = 12) +
  labs(y="Count",x=x.lab,subtitle = "d. RSV")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35),panel.spacing = unit(0.2, "lines")) 

hcv.middle.hist <- ggplot(hcv.like.dist.middle,aes(x=diff,group=metric))+
  geom_histogram(alpha=0.4, color = "#1B9E77", fill= "#1B9E77") +
  theme_bw(base_size = 12) +
  labs(y="Count",x=x.lab,subtitle = "c. HCV")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35),panel.spacing = unit(0.2, "lines")) 

alpha.middle.hist <- ggplot(alpha.like.dist.middle,aes(x=diff,group=metric))+
  geom_histogram(alpha=0.4, color = "#E41A1C", fill="#E41A1C") +
  theme_bw(base_size = 12) +
  labs(y="Count",x=x.lab,subtitle = "a. Alpha")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35),panel.spacing = unit(0.2, "lines")) 

delta.middle.hist <- ggplot(delta.like.dist.middle,aes(x=diff,group=metric))+
  geom_histogram(alpha=0.4, color = "#7570B3", fill="#7570B3") +
  theme_bw(base_size = 12) +
  labs(y=y.lab,x=x.lab,subtitle = "b. Delta")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35),panel.spacing = unit(0.2, "lines")) 

ggarrange(alpha.middle.hist,delta.middle.hist,hcv.middle.hist,rsv.middle.hist)

ggsave(,filename = "dist-middle-hist.png",dpi=300,width=12,height = 12)

rsv.hist <- ggplot(rsv.like.dist,aes(x=diff))+
  geom_histogram(alpha=0.4, color =  "#1F78B4", fill = "#1F78B4") +
  theme_bw(base_size = 12) +
  labs(y="Count",x=x.lab,subtitle = "d. RSV")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35),panel.spacing = unit(0.2, "lines")) 

hcv.hist <- ggplot(hcv.like.dist,aes(x=diff,group=metric))+
  geom_histogram(alpha=0.4, color = "#1B9E77", fill= "#1B9E77") +
  theme_bw(base_size = 12) +
  labs(y="Count",x=x.lab,subtitle = "c. HCV")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35),panel.spacing = unit(0.2, "lines")) 

alpha.hist <- ggplot(alpha.like.dist,aes(x=diff,group=metric))+
  geom_histogram(alpha=0.4, color = "#E41A1C", fill="#E41A1C") +
  theme_bw(base_size = 12) +
  labs(y="Count",x=x.lab,subtitle = "a. Alpha")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35),panel.spacing = unit(0.2, "lines")) 

delta.hist <- ggplot(delta.like.dist,aes(x=diff,group=metric))+
  geom_histogram(alpha=0.4, color = "#7570B3", fill="#7570B3") +
  theme_bw(base_size = 12) +
  labs(y=y.lab,x=x.lab,subtitle = "b. Delta")+
  facet_wrap(metric~.,scales = "free_x",ncol = 3) +
  theme(axis.text.x = element_text(angle = 35),panel.spacing = unit(0.2, "lines")) 

ggarrange(alpha.hist,delta.hist,hcv.hist,rsv.hist)

ggsave(,filename = "dist-first-hist.png",dpi=300,width=12,height = 12)
