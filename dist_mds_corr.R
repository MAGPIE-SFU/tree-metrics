# MDS and corr plots for tree metrics paper
# May 9 2024

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

#### Data ####
rsv.tree <- read.nexus("TtB_fine_trees/RSV2.trees")
rsv.log <- parse_beast_tracelog_file("TtB_fine_trees/RSV2.log")

#### Look at the MCMC ####
# plot posterior 
plot(rsv.log$posterior,xlim = c(0,1000))
plot(rsv.log$likelihood, xlim = c(0,2000))


#### Sample trees ####
# fist 50, middle 50, and last 50 ?
set.seed(666)
random.trees <- floor(runif(250, min = 50, max = 100000))
tree.sample <- rsv.tree[c(1:49,49975:50025,99951:100000,random.trees)]

# sample first 400
#tree.sample <- rsv.tree[c(1:400)]

#### Distance matrices ####
# distance matrices
RFmat <- phangorn::RF.dist(tree.sample)
wRFmat <- phangorn::wRF.dist(tree.sample)
KFmat <- phangorn::KF.dist(tree.sample)
#SPRmat <- phangorn::SPR.dist(tree.sample)
PDmat <- phangorn::path.dist(tree.sample)
JRFmat <- TreeDist::JaccardRobinsonFoulds(tree.sample)
KCmat <- TreeDist::KendallColijn(tree.sample)
SPRmat <- TreeDist::SPRDist(tree.sample)
MSmat <- TreeDist::MatchingSplitDistance(tree.sample)
NNImat <- TreeDist::NNIDist(tree.sample)
BHVmat <- distory::dist.multiPhylo(tree.sample)

# I saved the above since they took a while to run
#load("treemetrics-mds.Rdata")

#### MDS plots ####
calc.mds <- function(D,name){
  # fun to compute MDS and return dataframe with name of method as a col
  # D = dist matrix
  # name = name of method (character string)
  dist.mds <- cmdscale(D)
  colnames(dist.mds) <- c("V1","V2")
  dist.mds <- data.frame(Metric = rep(name),Group = c(rep("First 50",50),rep("Mid 50",50),rep("Last 50",50),rep("Random",250)),dist.mds)
  dist.mds <- data.frame(Metric = rep(name),dist.mds)
  
  return(dist.mds)
}

KF.mds <- calc.mds(KFmat,"KF")
RF.mds <- calc.mds(RFmat,"RF")
wRF.mds <- calc.mds(wRFmat,"wRF")
JRF.mds <- calc.mds(JRFmat,"JRF")
MS.mds <- calc.mds(MSmat,"MS")
KC.mds <- calc.mds(KCmat,"KC")
PD.mds <- calc.mds(PDmat,"PD")
NNI.mds <- calc.mds(NNImat$best_lower,"NNI")
SPR.mds <- calc.mds(SPRmat,"SPR")
BHV.mds <- calc.mds(BHVmat,"BHV")

mds.data <- rbind(KF.mds,RF.mds,wRF.mds,JRF.mds,KC.mds,MS.mds,PD.mds,NNI.mds,SPR.mds,BHV.mds)
mds.data <- rbind(KF.mds,RF.mds,wRF.mds,JRF.mds,KC.mds,MS.mds,PD.mds)

mds.data$Group <- factor(mds.data$Group, levels = c("Random","First 50","Last 50","Mid 50"))

mds.data %>%
  arrange(Group) %>%
ggplot( aes(y = V1, x = V2, color=Group)) +  
  geom_point(show.legend = T, stroke = 1, alpha=0.4) +
  theme_bw(base_size=12)+
  theme(legend.position = "none")+
  scale_color_manual(values = c("lightgrey","#7FC97F", "#BEAED4","#FDC086")) +
  facet_wrap(.~Metric, scales = "free") +
  xlab("MDS2") +
  ylab("MDS1")

ggsave(file="mdsplots.png",width=9,height=8,dpi=300,bg="white")


#### Corr plots ####
#RFpair=wRFpair=KFpair=PDpair=JRFpair=KCpair=SPRpair=MSpair= NNIpair=BHVpair=rep(NA,398)
for (i in 1:399) {
#  RFpair[i] <- phangorn::RF.dist(tree.sample[[i]],tree.sample[[i+1]])
#  wRFpair[i] <- phangorn::wRF.dist(tree.sample[[i]],tree.sample[[i+1]])
#  KFpair[i] <- phangorn::KF.dist(tree.sample[[i]],tree.sample[[i+1]])
#  PDpair[i] <- phangorn::path.dist(tree.sample[[i]],tree.sample[[i+1]])
#  JRFpair[i] <- TreeDist::JaccardRobinsonFoulds(tree.sample[[i]],tree.sample[[i+1]])
#  KCpair[i] <- TreeDist::KendallColijn(tree.sample[[i]],tree.sample[[i+1]])
 # SPRpair[i] <- TreeDist::SPRDist(tree.sample[[i]],tree.sample[[i+1]])
# MSpair[i] <- TreeDist::MatchingSplitDistance(tree.sample[[i]],tree.sample[[i+1]])
 # NNIpair[i] <- TreeDist::NNIDist(tree.sample[[i]],tree.sample[[i+1]])
 # BHVpair[i] <- distory::dist.multiPhylo(tree.sample[[i]],tree.sample[[i+1]])
}



library(GGally)


metrics.df <- cbind(RFpair,wRFpair,KFpair,PDpair,JRFpair,KCpair, MSpair,NNIpair,SPRpair)
colnames(metrics.df) <- c("RF","wRF","KF","PD","JRF","KC","MS","NNI","SPR")



pairp <- ggpairs(metrics.df, aes(alpha = 0.4),
        upper = list(na = "na")) + theme_bw(legend.title=element_blank()) + scale_color_manual(values = c("lightgrey"))





