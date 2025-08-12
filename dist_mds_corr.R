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
write.tree(tree.sample, "rsv-randomtrees-mds.trees")

# sample first 400
#tree.sample <- rsv.tree[c(1:400)]

#### Distance matrices ####
# distance matrices
RFmat <- phangorn::RF.dist(tree.sample)
wRFmat <- phangorn::wRF.dist(tree.sample)
KFmat <- phangorn::KF.dist(tree.sample)
SPRmat <- phangorn::SPR.dist(tree.sample)
PDmat <- phangorn::path.dist(tree.sample)
JRFmat <- TreeDist::JaccardRobinsonFoulds(tree.sample)
KCmat <- TreeDist::KendallColijn(tree.sample)
SPRmat <- TreeDist::SPRDist(tree.sample)
MSmat <- TreeDist::MatchingSplitDistance(tree.sample)
NNImat <- TreeDist::NNIDist(tree.sample)
IRFmat <- TreeDist::InfoRobinsonFoulds(tree.sample)
BHVmat <- distory::dist.multiPhylo(tree.sample)


rnnimat <- read.csv("rsv-randomtrees-mds.csv")

# turn rnni into a matrix
mat <- matrix(0, nrow = 400, ncol = 400)
mat[lower.tri(mat, diag = FALSE)] <- rnnimat$rNNI
RNNImat <- as.dist(mat)

### TO DO May 24 2025 ####
### --> add trip and quart dist and change nni use [4]
### then generate the mds plots 

library(Quartet)

NNImat <- TreeDist::NNIDist(tree.sample)[4]
NNImat <- as.matrix(NNImat$best_upper)

Tmat <- Quartet::TripletDistance(TQFile(tree.sample),TQFile(tree.sample))
Qmat <- Quartet::QuartetDistance(TQFile(tree.sample),TQFile(tree.sample))

Tmat <- matrix(nrow=400,ncol = 400)
for (j in 1:400){
  for (i in 1:400) {
Tmat[i,j] <- Quartet::TripletDistance(TQFile(tree.sample[[i]]),TQFile(tree.sample[[j]]))
  }}

Qmat <- matrix(nrow=400,ncol = 400)
for (j in 1:400){
  for (i in 1:400) {
    Qmat[i,j] <- Quartet::QuartetDistance(TQFile(tree.sample[[i]]),TQFile(tree.sample[[j]]))
  }}






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
IRF.mds <- calc.mds(IRFmat,"IRF")
NNI.mds <- calc.mds(NNImat,"NNI")
SPR.mds <- calc.mds(SPRmat,"SPR")
BHV.mds <- calc.mds(BHVmat,"BHV")
Trip.mds <- calc.mds(Tmat,"Trip")
Quart.mds <-  calc.mds(Qmat,"Quart")
RNNI.mds <- calc.mds(RNNImat,"RNNI")


mds.data <- rbind(RF.mds,wRF.mds,JRF.mds,IRF.mds,
                  MS.mds,KF.mds,PD.mds,KC.mds,
                  BHV.mds,Trip.mds,Quart.mds,SPR.mds,NNI.mds,RNNI.mds)

#mds.data <- rbind(KF.mds,RF.mds,wRF.mds,JRF.mds,KC.mds,MS.mds,PD.mds)

mds.data$Group <- factor(mds.data$Group, levels = c("Random","First 50","Last 50","Mid 50"))
mds.data$Metric <- factor(mds.data$Metric, levels = c("RF", "wRF", "JRF", "IRF", "MS", "KF", "PD", "KC", "BHV", "Trip", "Quart", "SPR", "NNI", "RNNI"))
levels(as.factor(mds.data$Metric))

library(dplyr)
library(ggplot2)
mds.data %>%
  arrange(Group) %>%
ggplot( aes(y = V1, x = V2, color=Group)) +  
  geom_point(show.legend = T, stroke = 1, alpha=0.4) +
  theme_bw(base_size=12)+
  scale_color_manual(values = c("lightgrey","#7FC97F", "#BEAED4","#FDC086")) +
  facet_wrap(.~Metric, scales = "free") +
  xlab("MDS2") +
  ylab("MDS1") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(.75, .1),
    legend.title = element_blank()
   # legend.justification = c("bottom")
   # legend.box.just = "left",
   # legend.margin = margin(6, 6, 6, 6)
  )

ggsave(file="mdsplots-2025.png",width=9,height=9,dpi=800,bg="white")


#### Corr plots ####
#RFpair=wRFpair=KFpair=PDpair=JRFpair=KCpair=SPRpair=MSpair= NNIpair=BHVpair=rep(NA,398)
for (i in 1:399) {
  RFpair[i] <- phangorn::RF.dist(tree.sample[[i]],tree.sample[[i+1]])
  wRFpair[i] <- phangorn::wRF.dist(tree.sample[[i]],tree.sample[[i+1]])
  KFpair[i] <- phangorn::KF.dist(tree.sample[[i]],tree.sample[[i+1]])
  PDpair[i] <- phangorn::path.dist(tree.sample[[i]],tree.sample[[i+1]])
  JRFpair[i] <- TreeDist::JaccardRobinsonFoulds(tree.sample[[i]],tree.sample[[i+1]])
  KCpair[i] <- TreeDist::KendallColijn(tree.sample[[i]],tree.sample[[i+1]])
  SPRpair[i] <- TreeDist::SPRDist(tree.sample[[i]],tree.sample[[i+1]])
 MSpair[i] <- TreeDist::MatchingSplitDistance(tree.sample[[i]],tree.sample[[i+1]])
  NNIpair[i] <- TreeDist::NNIDist(tree.sample[[i]],tree.sample[[i+1]])
  BHVpair[i] <- distory::dist.multiPhylo(tree.sample[[i]],tree.sample[[i+1]])
}



library(GGally)


metrics.df <- cbind(RFpair,wRFpair,KFpair,PDpair,JRFpair,KCpair, MSpair,NNIpair,SPRpair)
colnames(metrics.df) <- c("RF","wRF","KF","PD","JRF","KC","MS","NNI","SPR")

pairp <- ggpairs(metrics.df, aes(alpha = 0.4),
        upper = list(na = "na")) + theme_bw(legend.title=element_blank()) + scale_color_manual(values = c("lightgrey"))





