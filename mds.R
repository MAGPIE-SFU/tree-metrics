# new trees 2026

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
library(Quartet)



# load in all the RSV trees -- the extra 250 are already sampled
tree.sample.rsv <- read.nexus("treemetricsamples/RSV2.all.trees")
tree.sample <- tree.sample.rsv

## Distance matrices ##
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

# read in rnni values
rnnimat <- read.csv("treemetricsamples/RSV2.all.trees.rnni.csv")

# turn rnni into a matrix
mat <- matrix(0, nrow = 400, ncol = 400)
mat[lower.tri(mat, diag = FALSE)] <- rnnimat$rNNI
RNNImat <- as.dist(mat)


# NNI 
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

ggsave(file="mdsplots-2026.png",width=9,height=9,dpi=800,bg="white")


## correlation bar plot 
RFpair=wRFpair=KFpair=PDpair=JRFpair=KCpair=SPRpair=MSpair= NNIpair=BHVpair=rep(NA,398)
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
