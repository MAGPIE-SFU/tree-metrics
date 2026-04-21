# First run the like_dist scripts for first, mid, last 
# to get all of the likelihood and distance values. 
# I saved them as .rds on 21 April 2026, so load them


### make one absolutely massive df for each dataset
rsv.all <- data.frame(rbind(rsv.first,rsv.mid,rsv.last), Sample= rep(c("First 50","Middle 50","Last 50"),each = 17150))
hcv.all <- data.frame(rbind(hcv.first,hcv.mid,hcv.last), Sample= rep(c("First 50","Middle 50","Last 50"),each = 17150))
alpha.all <- data.frame(rbind(alpha.first,alpha.mid,alpha.last), Sample= rep(c("First 50","Middle 50","Last 50"),each = 17150))
delta.all <- data.frame(rbind(delta.first,delta.mid,delta.last), Sample= rep(c("First 50","Middle 50","Last 50"),each = 17150))

# corrlation for each dataset
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

head(cordata)
unique(cordata$metric)

#the metrics should already be in the correct order
metric.order <- c("RF", "wRF", "JRF", "IRF", "MS", "KF", "PD", "KC", "BHV", "Trip", "Quart", "SPR", "NNI", "RNNI")



#create color spectrum
blues13 <- colorRampPalette(c("#E41A1C","#FFB7C3","#7570B3",
                            "#F6D8AE","#1F78B4","#28112B",
                            "#1B9E77","#D95D39","#E8C547",
                            "#0FA3B1","#AC3931","#bc6c25",
                            "#6a994e","#A53860"), bias=1)(14)

# correlation plot 
ggplot(cordata,aes(x=data,y=correlation,fill=metric,color=metric))+
  geom_bar(stat = "identity",position = "dodge",alpha=0.7, width=0.75) +
  scale_fill_manual(values = blues13) +
  scale_color_manual(values = blues13) +
  theme_bw(base_size = 12) +
  theme(legend.title = element_blank()) +
  labs(y="Correlation",x="Data") +
  facet_grid(Sample~. )

ggsave(filename = "like-corr-all-2026.png",dpi=300,width = 6,height = 6)
