####Rarefaction analyses####
library(iNEXT)
library(gridExtra)
library(grid)
library(ggplot2)

polTPint <- read.table("inex_polxindvTP.txt" , h=T)

polTP1i <- polTPint[,1]
out1i <- iNEXT(polTP1i, q=0, datatype="abundance" , endpoint = 291)
ChaoRichness(polTP1i)
p1i<-ggiNEXT(out1i , type = 2)+
  labs(title="Network 1" , x="Number of interactions")+
  theme(legend.position="none")
grid.arrange(p1i, ncol=1)

polTP5i <- polTPint[,2]
out5i <- iNEXT(polTP5i, q=0, datatype="abundance" , endpoint = 534)
ChaoRichness(polTP5i)
p5i<-ggiNEXT(out5i , type = 2)+
  labs(title="Network 2",x="Number of interactions")+
  theme(legend.position="none")
grid.arrange(p5i, ncol=1)

polTP6i <- polTPint[,3]
out6i <- iNEXT(polTP6i, q=0, datatype="abundance" , endpoint = 1752)
ChaoRichness(polTP6i)
p6i<-ggiNEXT(out6i , type = 2)+
  labs(title="Network 3",x="Number of interactions")+
  theme(legend.position="none")
grid.arrange(p6i, ncol=1)

polTP8i <- polTPint[,4]
out8i <- iNEXT(polTP8i, q=0, datatype="abundance" , endpoint = 461)
ChaoRichness(polTP8i)
p8i<-ggiNEXT(out8i , type = 2)+
  labs(title="Network 4",x="Number of interactions")+
  theme(legend.position="none")
grid.arrange(p8i, ncol=1)

polTP10i <- polTPint[,5]
out10i <- iNEXT(polTP10i, q=0, datatype="abundance" , endpoint = 503)
ChaoRichness(polTP10i)
p10i<-ggiNEXT(out10i , type = 2)+
  labs(title="Network 5",x="Number of interactions")+
  theme(legend.position="none")
grid.arrange(p10i, ncol=1)

polTP11i <- polTPint[,6]
out11i <- iNEXT(polTP11i, q=0, datatype="abundance"  , endpoint = 124)
ChaoRichness(polTP11i)
p11i<-ggiNEXT(out11i , type = 2)+
  labs(title="Network 6",x="Number of interactions")+
  theme(legend.position="none")
grid.arrange(p11i, ncol=1)

polTP14i <- polTPint[,7]
out14i <- iNEXT(polTP14i, q=0, datatype="abundance" , endpoint = 78)
ChaoRichness(polTP14i)
p14i<-ggiNEXT(out14i , type = 2)+
  labs(title="Network 7",x="Number of interactions")+
  theme(legend.position="none")
grid.arrange(p14i, ncol=1)

polTP15i <- polTPint[,8]
out15i <- iNEXT(polTP15i, q=0, datatype="abundance" , endpoint = 549)
ChaoRichness(polTP15i)
p15i<-ggiNEXT(out15i , type = 2)+
  labs(title="Network 8",x="Number of interactions")+
  theme(legend.position="none")
grid.arrange(p15i, ncol=1)

polTP16i <- polTPint[,9]
out16i <- iNEXT(polTP16i, q=0, datatype="abundance" , endpoint = 480)
ChaoRichness(polTP16i)
p16i<-ggiNEXT(out16i , type = 2)+
  labs(title="Network 9",x="Number of interactions")+
  theme(legend.position="none")
grid.arrange(p16i, ncol=1)

polTP17i <- polTPint[,10]
out17i <- iNEXT(polTP17i, q=0, datatype="abundance" , endpoint = 489)
ChaoRichness(polTP17i)
p17i<-ggiNEXT(out17i , type = 2)+
  labs(title="Network 10",x="Number of interactions")+
  theme(legend.position="none")
grid.arrange(p17i, ncol=1)

polTP21i <- polTPint[,11]
out21i <- iNEXT(polTP21i, q=0, datatype="abundance" , endpoint = 385)
ChaoRichness(polTP21i)
p21i<-ggiNEXT(out21i , type = 2)+
  labs(title="Network 11",x="Number of interactions")+
  theme(legend.position="none")
grid.arrange(p21i, ncol=1)

polTP23i <- polTPint[,12]
out23i <- iNEXT(polTP23i, q=0, datatype="abundance" , endpoint = 78)
ChaoRichness(polTP23i)
p23i<-ggiNEXT(out23i , type = 2)+
  labs(title="Network 12",x="Number of interactions")+
  theme(legend.position="none")
grid.arrange(p23i, ncol=1)

polTP25i <- polTPint[,13]
out25i <- iNEXT(polTP25i, q=0, datatype="abundance" , endpoint = 492)
ChaoRichness(polTP25i)
p25i<-ggiNEXT(out25i , type = 2)+
  labs(title="Network 13",x="Number of interactions")+
  theme(legend.position="none")
grid.arrange(p25i, ncol=1)

polTP28i <- polTPint[,14]
out28i <- iNEXT(polTP28i, q=0, datatype="abundance" , endpoint = 189)
ChaoRichness(polTP28i)
p28i<-ggiNEXT(out28i , type = 2)+
  labs(title="Network 14",x="Number of interactions")+
  theme(legend.position="none")
grid.arrange(p28i, ncol=1)

polTP30i <- polTPint[,15]
out30i <- iNEXT(polTP30i, q=0, datatype="abundance" , endpoint = 193)
ChaoRichness(polTP30i)
p30i<-ggiNEXT(out30i , type = 2)+
  labs(title="Network 15" , x="Number of interactions")+
  theme(legend.position="none")
grid.arrange(p30i, ncol=1)

png(file="rareintnetworks.png", width=15, height=7, units="in", res=300)
grid.arrange(p1i, p5i, p6i, p8i, p10i, p11i, p14i, p15i, p16i, p17i, p21i, p23i, p25i, p28i, p30i, nrow=3, ncol=5)
dev.off()

####Network metrics####
require(bipartite)

PolC <- read.table("PolC.txt" , h=T)
int_polc <- frame2webs(PolC)
attach(int_polc)
names(int_polc)
####Network plots####
plotweb(int_polc$network1, labsize = 2, text.rot = 90, y.lim = c(-0.1,2.5),
        col.high = 'chartreuse3', col.low = 'tan1',
        bor.col.high = 'white', bor.col.low = 'white',
        col.interaction = 'khaki1', bor.col.interaction = 'khaki1',
        arrow = 'no',)

plotweb(int_polc$network2, labsize = 2, text.rot = 90, y.lim = c(-0.1,2.5),
        col.high = 'chartreuse3', col.low = 'tan1',
        bor.col.high = 'white', bor.col.low = 'white',
        col.interaction = 'khaki1', bor.col.interaction = 'khaki1',
        arrow = 'no',)

plotweb(int_polc$network3, labsize = 2, text.rot = 90, y.lim = c(-0.1,2.5),
        col.high = 'chartreuse3', col.low = 'tan1',
        bor.col.high = 'white', bor.col.low = 'white',
        col.interaction = 'khaki1', bor.col.interaction = 'khaki1',
        arrow = 'no',)

plotweb(int_polc$network4, labsize = 2, text.rot = 90, y.lim = c(-0.1,2.5),
        col.high = 'chartreuse3', col.low = 'tan1',
        bor.col.high = 'white', bor.col.low = 'white',
        col.interaction = 'khaki1', bor.col.interaction = 'khaki1',
        arrow = 'no',)

plotweb(int_polc$network5, labsize = 2, text.rot = 90, y.lim = c(-0.1,2.5),
        col.high = 'chartreuse3', col.low = 'tan1',
        bor.col.high = 'white', bor.col.low = 'white',
        col.interaction = 'khaki1', bor.col.interaction = 'khaki1',
        arrow = 'no',)

plotweb(int_polc$network6, labsize = 2, text.rot = 90, y.lim = c(-0.1,2.5),
        col.high = 'chartreuse3', col.low = 'tan1',
        bor.col.high = 'white', bor.col.low = 'white',
        col.interaction = 'khaki1', bor.col.interaction = 'khaki1',
        arrow = 'no',)

plotweb(int_polc$network7, labsize = 2, text.rot = 90, y.lim = c(-0.1,2.5),
        col.high = 'chartreuse3', col.low = 'tan1',
        bor.col.high = 'white', bor.col.low = 'white',
        col.interaction = 'khaki1', bor.col.interaction = 'khaki1',
        arrow = 'no',)

plotweb(int_polc$network8, labsize = 2, text.rot = 90, y.lim = c(-0.1,2.5),
        col.high = 'chartreuse3', col.low = 'tan1',
        bor.col.high = 'white', bor.col.low = 'white',
        col.interaction = 'khaki1', bor.col.interaction = 'khaki1',
        arrow = 'no',)

plotweb(int_polc$network9, labsize = 2, text.rot = 90, y.lim = c(-0.1,2.5),
        col.high = 'chartreuse3', col.low = 'tan1',
        bor.col.high = 'white', bor.col.low = 'white',
        col.interaction = 'khaki1', bor.col.interaction = 'khaki1',
        arrow = 'no',)

plotweb(int_polc$network10, labsize = 2, text.rot = 90, y.lim = c(-0.1,2.5),
        col.high = 'chartreuse3', col.low = 'tan1',
        bor.col.high = 'white', bor.col.low = 'white',
        col.interaction = 'khaki1', bor.col.interaction = 'khaki1',
        arrow = 'no',)

plotweb(int_polc$network11, labsize = 2, text.rot = 90, y.lim = c(-0.1,2.5),
        col.high = 'chartreuse3', col.low = 'tan1',
        bor.col.high = 'white', bor.col.low = 'white',
        col.interaction = 'khaki1', bor.col.interaction = 'khaki1',
        arrow = 'no',)

plotweb(int_polc$network12, labsize = 2, text.rot = 90, y.lim = c(-0.1,2.5),
        col.high = 'chartreuse3', col.low = 'tan1',
        bor.col.high = 'white', bor.col.low = 'white',
        col.interaction = 'khaki1', bor.col.interaction = 'khaki1',
        arrow = 'no',)

plotweb(int_polc$network13, labsize = 2, text.rot = 90, y.lim = c(-0.1,2.5),
        col.high = 'chartreuse3', col.low = 'tan1',
        bor.col.high = 'white', bor.col.low = 'white',
        col.interaction = 'khaki1', bor.col.interaction = 'khaki1',
        arrow = 'no',)

plotweb(int_polc$network14, labsize = 2, text.rot = 90, y.lim = c(-0.1,2.5),
        col.high = 'chartreuse3', col.low = 'tan1',
        bor.col.high = 'white', bor.col.low = 'white',
        col.interaction = 'khaki1', bor.col.interaction = 'khaki1',
        arrow = 'no',)

plotweb(int_polc$network15, labsize = 2, text.rot = 90, y.lim = c(-0.1,2.5),
        col.high = 'chartreuse3', col.low = 'tan1',
        bor.col.high = 'white', bor.col.low = 'white',
        col.interaction = 'khaki1', bor.col.interaction = 'khaki1',
        arrow = 'no',)

####H2####
# Network1
obs <- unlist(networklevel(network1, index="H2"))
nulls <- nullmodel(network1, N=1000, method=1)

null <- unlist(sapply(nulls, networklevel, index="H2"))
mean.null<-mean(null); min.null<-min(null); max.null<-max(null); sd.null<-sd(null)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(null>obs) / length(null)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(null), xlim=c(min(obs, min(null)), max(obs, max(null))),
     main="network1"); abline(v=obs, col="blue", lwd=2)  

res.network1 <- cbind(obs, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)

#result1
res.network1

# Network2
obs <- unlist(networklevel(network5, index="H2"))
nulls <- nullmodel(network5, N=1000, method=1)

null <- unlist(sapply(nulls, networklevel, index="H2"))
mean.null<-mean(null); min.null<-min(null); max.null<-max(null); sd.null<-sd(null)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(null>obs) / length(null)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(null), xlim=c(min(obs, min(null)), max(obs, max(null))),
     main="network5"); abline(v=obs, col="blue", lwd=2)  

res.network5 <- cbind(obs, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)

#result2
res.network5

# network3
obs <- unlist(networklevel(network6, index="H2"))
nulls <- nullmodel(network6, N=1000, method=1)

null <- unlist(sapply(nulls, networklevel, index="H2"))
mean.null<-mean(null); min.null<-min(null); max.null<-max(null); sd.null<-sd(null)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(null>obs) / length(null)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(null), xlim=c(min(obs, min(null)), max(obs, max(null))),
     main="network6"); abline(v=obs, col="blue", lwd=2)  

res.network6 <- cbind(obs, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)

#result3
res.network6

# network4
obs <- unlist(networklevel(network8, index="H2"))
nulls <- nullmodel(network8, N=1000, method=1)

null <- unlist(sapply(nulls, networklevel, index="H2"))
mean.null<-mean(null); min.null<-min(null); max.null<-max(null); sd.null<-sd(null)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(null>obs) / length(null)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(null), xlim=c(min(obs, min(null)), max(obs, max(null))),
     main="network8"); abline(v=obs, col="blue", lwd=2)  

res.network8 <- cbind(obs, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)

#result4
res.network8

# network5
obs <- unlist(networklevel(network10, index="H2"))
nulls <- nullmodel(network10, N=1000, method=1)

null <- unlist(sapply(nulls, networklevel, index="H2"))
mean.null<-mean(null); min.null<-min(null); max.null<-max(null); sd.null<-sd(null)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(null>obs) / length(null)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(null), xlim=c(min(obs, min(null)), max(obs, max(null))),
     main="network10"); abline(v=obs, col="blue", lwd=2)  

res.network10 <- cbind(obs, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)

#result5
res.network10

# network6
obs <- unlist(networklevel(network11, index="H2"))
nulls <- nullmodel(network11, N=1000, method=1)

null <- unlist(sapply(nulls, networklevel, index="H2"))
mean.null<-mean(null); min.null<-min(null); max.null<-max(null); sd.null<-sd(null)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(null>obs) / length(null)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(null), xlim=c(min(obs, min(null)), max(obs, max(null))),
     main="network11"); abline(v=obs, col="blue", lwd=2)  

res.network11 <- cbind(obs, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)

#result6
res.network11

# network7
obs <- unlist(networklevel(network14, index="H2"))
nulls <- nullmodel(network14, N=1000, method=1)

null <- unlist(sapply(nulls, networklevel, index="H2"))
mean.null<-mean(null); min.null<-min(null); max.null<-max(null); sd.null<-sd(null)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(null>obs) / length(null)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(null), xlim=c(min(obs, min(null)), max(obs, max(null))),
     main="network14"); abline(v=obs, col="blue", lwd=2)  

res.network14 <- cbind(obs, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)

#result7
res.network14

# network8
obs <- unlist(networklevel(network15, index="H2"))
nulls <- nullmodel(network15, N=1000, method=1)

null <- unlist(sapply(nulls, networklevel, index="H2"))
mean.null<-mean(null); min.null<-min(null); max.null<-max(null); sd.null<-sd(null)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(null>obs) / length(null)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(null), xlim=c(min(obs, min(null)), max(obs, max(null))),
     main="network15"); abline(v=obs, col="blue", lwd=2)  

res.network15 <- cbind(obs, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)

#result8
res.network15

# network9
obs <- unlist(networklevel(network16, index="H2"))
nulls <- nullmodel(network16, N=1000, method=1)

null <- unlist(sapply(nulls, networklevel, index="H2"))
mean.null<-mean(null); min.null<-min(null); max.null<-max(null); sd.null<-sd(null)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(null>obs) / length(null)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(null), xlim=c(min(obs, min(null)), max(obs, max(null))),
     main="network16"); abline(v=obs, col="blue", lwd=2)  

res.network16 <- cbind(obs, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)

#result9
res.network16

# network10
obs <- unlist(networklevel(network17, index="H2"))
nulls <- nullmodel(network17, N=1000, method=1)

null <- unlist(sapply(nulls, networklevel, index="H2"))
mean.null<-mean(null); min.null<-min(null); max.null<-max(null); sd.null<-sd(null)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(null>obs) / length(null)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(null), xlim=c(min(obs, min(null)), max(obs, max(null))),
     main="network17"); abline(v=obs, col="blue", lwd=2)  

res.network17 <- cbind(obs, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)

#result10
res.network17

# network11
obs <- unlist(networklevel(network21, index="H2"))
nulls <- nullmodel(network21, N=1000, method=1)

null <- unlist(sapply(nulls, networklevel, index="H2"))
mean.null<-mean(null); min.null<-min(null); max.null<-max(null); sd.null<-sd(null)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(null>obs) / length(null)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(null), xlim=c(min(obs, min(null)), max(obs, max(null))),
     main="network21"); abline(v=obs, col="blue", lwd=2)  

res.network21 <- cbind(obs, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)

#result11
res.network21

# network12
obs <- unlist(networklevel(network23, index="H2"))
nulls <- nullmodel(network23, N=1000, method=1)

null <- unlist(sapply(nulls, networklevel, index="H2"))
mean.null<-mean(null); min.null<-min(null); max.null<-max(null); sd.null<-sd(null)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(null>obs) / length(null)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(null), xlim=c(min(obs, min(null)), max(obs, max(null))),
     main="network23"); abline(v=obs, col="blue", lwd=2)  

res.network23 <- cbind(obs, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)

#result12
res.network23

# network13
obs <- unlist(networklevel(network25, index="H2"))
nulls <- nullmodel(network25, N=1000, method=1)

null <- unlist(sapply(nulls, networklevel, index="H2"))
mean.null<-mean(null); min.null<-min(null); max.null<-max(null); sd.null<-sd(null)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(null>obs) / length(null)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(null), xlim=c(min(obs, min(null)), max(obs, max(null))),
     main="network25"); abline(v=obs, col="blue", lwd=2)  

res.network25 <- cbind(obs, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)

#result13
res.network25

# network14
obs <- unlist(networklevel(network28, index="H2"))
nulls <- nullmodel(network28, N=1000, method=1)

null <- unlist(sapply(nulls, networklevel, index="H2"))
mean.null<-mean(null); min.null<-min(null); max.null<-max(null); sd.null<-sd(null)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(null>obs) / length(null)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(null), xlim=c(min(obs, min(null)), max(obs, max(null))),
     main="network28"); abline(v=obs, col="blue", lwd=2)  

res.network28 <- cbind(obs, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)

#result14
res.network28

# network15
obs <- unlist(networklevel(network30, index="H2"))
nulls <- nullmodel(network30, N=1000, method=1)

null <- unlist(sapply(nulls, networklevel, index="H2"))
mean.null<-mean(null); min.null<-min(null); max.null<-max(null); sd.null<-sd(null)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(null>obs) / length(null)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(null), xlim=c(min(obs, min(null)), max(obs, max(null))),
     main="network30"); abline(v=obs, col="blue", lwd=2)  

res.network30 <- cbind(obs, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)

#result15
res.network30

####Modularity####

#network1
mod1.1<-computeModules(network1);mod1.2<-computeModules(network1);mod1.3<-computeModules(network1);mod1.4<-computeModules(network1);mod1.5<-computeModules(network1);mod1.6<-computeModules(network1);mod1.7<-computeModules(network1);mod1.8<-computeModules(network1);mod1.9<-computeModules(network1);mod1.10<-computeModules(network1)
mod.obs<-c(mod1.1@likelihood,mod1.2@likelihood,mod1.3@likelihood,mod1.4@likelihood,mod1.5@likelihood,mod1.6@likelihood,mod1.7@likelihood,mod1.8@likelihood,mod1.9@likelihood,mod1.10@likelihood)
obs<-max(mod.obs);sd.obs<-sd(mod.obs)

plotModuleWeb(mod1.1)

nulls <- nullmodel(network1, N=1000, method=1)
modules.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
zscore <- (obs - mean(like.nulls))/sd(like.nulls)

mean.null<-mean(like.nulls); min.null<-min(like.nulls); max.null<-max(like.nulls); sd.null<-sd(like.nulls)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(like.nulls>obs) / length(like.nulls)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(like.nulls), xlim=c(min((obs), min(like.nulls)), max((obs), max(like.nulls))),
     main="network1"); abline(v=(obs), col="blue", lwd=2)

mod.network1 <- cbind(obs, sd.obs,zscore, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)
mod.network1

# network2

mod1.1<-computeModules(network5);mod1.2<-computeModules(network5);mod1.3<-computeModules(network5);mod1.4<-computeModules(network5);mod1.5<-computeModules(network5);mod1.6<-computeModules(network5);mod1.7<-computeModules(network5);mod1.8<-computeModules(network5);mod1.9<-computeModules(network5);mod1.10<-computeModules(network5)
mod.obs<-c(mod1.1@likelihood,mod1.2@likelihood,mod1.3@likelihood,mod1.4@likelihood,mod1.5@likelihood,mod1.6@likelihood,mod1.7@likelihood,mod1.8@likelihood,mod1.9@likelihood,mod1.10@likelihood)
obs<-max(mod.obs);sd.obs<-sd(mod.obs)

plotModuleWeb(mod1.1)

nulls <- nullmodel(network5, N=1000, method=1)
modules.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
zscore <- (obs - mean(like.nulls))/sd(like.nulls)

mean.null<-mean(like.nulls); min.null<-min(like.nulls); max.null<-max(like.nulls); sd.null<-sd(like.nulls)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(like.nulls>obs) / length(like.nulls)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(like.nulls), xlim=c(min((obs), min(like.nulls)), max((obs), max(like.nulls))),
     main="network5"); abline(v=(obs), col="blue", lwd=2)

mod.network2 <- cbind(obs, sd.obs,zscore, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)
mod.network2

# network3

mod1.1<-computeModules(network6);mod1.2<-computeModules(network6);mod1.3<-computeModules(network6);mod1.4<-computeModules(network6);mod1.5<-computeModules(network6);mod1.6<-computeModules(network6);mod1.7<-computeModules(network6);mod1.8<-computeModules(network6);mod1.9<-computeModules(network6);mod1.10<-computeModules(network6)
mod.obs<-c(mod1.1@likelihood,mod1.2@likelihood,mod1.3@likelihood,mod1.4@likelihood,mod1.5@likelihood,mod1.6@likelihood,mod1.7@likelihood,mod1.8@likelihood,mod1.9@likelihood,mod1.10@likelihood)
obs<-max(mod.obs);sd.obs<-sd(mod.obs)

plotModuleWeb(mod1.1)

nulls <- nullmodel(network5, N=1000, method=1)
modules.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
zscore <- (obs - mean(like.nulls))/sd(like.nulls)

mean.null<-mean(like.nulls); min.null<-min(like.nulls); max.null<-max(like.nulls); sd.null<-sd(like.nulls)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(like.nulls>obs) / length(like.nulls)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(like.nulls), xlim=c(min((obs), min(like.nulls)), max((obs), max(like.nulls))),
     main="network6"); abline(v=(obs), col="blue", lwd=2)

mod.network3 <- cbind(obs, sd.obs,zscore, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)
mod.network3

# network4

mod1.1<-computeModules(network8);mod1.2<-computeModules(network8);mod1.3<-computeModules(network8);mod1.4<-computeModules(network8);mod1.5<-computeModules(network8);mod1.6<-computeModules(network8);mod1.7<-computeModules(network8);mod1.8<-computeModules(network8);mod1.9<-computeModules(network8);mod1.10<-computeModules(network8)
mod.obs<-c(mod1.1@likelihood,mod1.2@likelihood,mod1.3@likelihood,mod1.4@likelihood,mod1.5@likelihood,mod1.6@likelihood,mod1.7@likelihood,mod1.8@likelihood,mod1.9@likelihood,mod1.10@likelihood)
obs<-max(mod.obs);sd.obs<-sd(mod.obs)

plotModuleWeb(mod1.1)

nulls <- nullmodel(network8, N=1000, method=1)
modules.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
zscore <- (obs - mean(like.nulls))/sd(like.nulls)

mean.null<-mean(like.nulls); min.null<-min(like.nulls); max.null<-max(like.nulls); sd.null<-sd(like.nulls)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(like.nulls>obs) / length(like.nulls)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(like.nulls), xlim=c(min((obs), min(like.nulls)), max((obs), max(like.nulls))),
     main="network8"); abline(v=(obs), col="blue", lwd=2)

mod.network4 <- cbind(obs, sd.obs,zscore, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)
mod.network4

# network5

mod1.1<-computeModules(network10);mod1.2<-computeModules(network10);mod1.3<-computeModules(network10);mod1.4<-computeModules(network10);mod1.5<-computeModules(network10);mod1.6<-computeModules(network10);mod1.7<-computeModules(network10);mod1.8<-computeModules(network10);mod1.9<-computeModules(network10);mod1.10<-computeModules(network10)
mod.obs<-c(mod1.1@likelihood,mod1.2@likelihood,mod1.3@likelihood,mod1.4@likelihood,mod1.5@likelihood,mod1.6@likelihood,mod1.7@likelihood,mod1.8@likelihood,mod1.9@likelihood,mod1.10@likelihood)
obs<-max(mod.obs);sd.obs<-sd(mod.obs)

plotModuleWeb(mod1.1)

nulls <- nullmodel(network10, N=1000, method=1)
modules.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
zscore <- (obs - mean(like.nulls))/sd(like.nulls)

mean.null<-mean(like.nulls); min.null<-min(like.nulls); max.null<-max(like.nulls); sd.null<-sd(like.nulls)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(like.nulls>obs) / length(like.nulls)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(like.nulls), xlim=c(min((obs), min(like.nulls)), max((obs), max(like.nulls))),
     main="network10"); abline(v=(obs), col="blue", lwd=2)

mod.network5 <- cbind(obs, sd.obs,zscore, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)
mod.network5

# network6

mod1.1<-computeModules(network11);mod1.2<-computeModules(network11);mod1.3<-computeModules(network11);mod1.4<-computeModules(network11);mod1.5<-computeModules(network11);mod1.6<-computeModules(network11);mod1.7<-computeModules(network11);mod1.8<-computeModules(network11);mod1.9<-computeModules(network11);mod1.10<-computeModules(network11)
mod.obs<-c(mod1.1@likelihood,mod1.2@likelihood,mod1.3@likelihood,mod1.4@likelihood,mod1.5@likelihood,mod1.6@likelihood,mod1.7@likelihood,mod1.8@likelihood,mod1.9@likelihood,mod1.10@likelihood)
obs<-max(mod.obs);sd.obs<-sd(mod.obs)

plotModuleWeb(mod1.1)

nulls <- nullmodel(network11, N=1000, method=1)
modules.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
zscore <- (obs - mean(like.nulls))/sd(like.nulls)

mean.null<-mean(like.nulls); min.null<-min(like.nulls); max.null<-max(like.nulls); sd.null<-sd(like.nulls)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(like.nulls>obs) / length(like.nulls)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(like.nulls), xlim=c(min((obs), min(like.nulls)), max((obs), max(like.nulls))),
     main="network11"); abline(v=(obs), col="blue", lwd=2)

mod.network6 <- cbind(obs, sd.obs,zscore, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)
mod.network6

# network7

mod1.1<-computeModules(network14);mod1.2<-computeModules(network14);mod1.3<-computeModules(network14);mod1.4<-computeModules(network14);mod1.5<-computeModules(network14);mod1.6<-computeModules(network14);mod1.7<-computeModules(network14);mod1.8<-computeModules(network14);mod1.9<-computeModules(network14);mod1.10<-computeModules(network14)
mod.obs<-c(mod1.1@likelihood,mod1.2@likelihood,mod1.3@likelihood,mod1.4@likelihood,mod1.5@likelihood,mod1.6@likelihood,mod1.7@likelihood,mod1.8@likelihood,mod1.9@likelihood,mod1.10@likelihood)
obs<-max(mod.obs);sd.obs<-sd(mod.obs)

plotModuleWeb(mod1.5)

nulls <- nullmodel(network14, N=1000, method=1)
modules.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
zscore <- (obs - mean(like.nulls))/sd(like.nulls)

mean.null<-mean(like.nulls); min.null<-min(like.nulls); max.null<-max(like.nulls); sd.null<-sd(like.nulls)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(like.nulls>obs) / length(like.nulls)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(like.nulls), xlim=c(min((obs), min(like.nulls)), max((obs), max(like.nulls))),
     main="network14"); abline(v=(obs), col="blue", lwd=2)

mod.network7 <- cbind(obs, sd.obs,zscore, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)
mod.network7

# network8

mod1.1<-computeModules(network15);mod1.2<-computeModules(network15);mod1.3<-computeModules(network15);mod1.4<-computeModules(network15);mod1.5<-computeModules(network15);mod1.6<-computeModules(network15);mod1.7<-computeModules(network15);mod1.8<-computeModules(network15);mod1.9<-computeModules(network15);mod1.10<-computeModules(network15)
mod.obs<-c(mod1.1@likelihood,mod1.2@likelihood,mod1.3@likelihood,mod1.4@likelihood,mod1.5@likelihood,mod1.6@likelihood,mod1.7@likelihood,mod1.8@likelihood,mod1.9@likelihood,mod1.10@likelihood)
obs<-max(mod.obs);sd.obs<-sd(mod.obs)

plotModuleWeb(mod1.1)

nulls <- nullmodel(network15, N=1000, method=1)
modules.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
zscore <- (obs - mean(like.nulls))/sd(like.nulls)

mean.null<-mean(like.nulls); min.null<-min(like.nulls); max.null<-max(like.nulls); sd.null<-sd(like.nulls)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(like.nulls>obs) / length(like.nulls)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(like.nulls), xlim=c(min((obs), min(like.nulls)), max((obs), max(like.nulls))),
     main="network15"); abline(v=(obs), col="blue", lwd=2)

mod.network8 <- cbind(obs, sd.obs,zscore, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)
mod.network8

# network9

mod1.1<-computeModules(network16);mod1.2<-computeModules(network16);mod1.3<-computeModules(network16);mod1.4<-computeModules(network16);mod1.5<-computeModules(network16);mod1.6<-computeModules(network16);mod1.7<-computeModules(network16);mod1.8<-computeModules(network16);mod1.9<-computeModules(network16);mod1.10<-computeModules(network16)
mod.obs<-c(mod1.1@likelihood,mod1.2@likelihood,mod1.3@likelihood,mod1.4@likelihood,mod1.5@likelihood,mod1.6@likelihood,mod1.7@likelihood,mod1.8@likelihood,mod1.9@likelihood,mod1.10@likelihood)
obs<-max(mod.obs);sd.obs<-sd(mod.obs)

plotModuleWeb(mod1.1)

nulls <- nullmodel(network16, N=1000, method=1)
modules.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
zscore <- (obs - mean(like.nulls))/sd(like.nulls)

mean.null<-mean(like.nulls); min.null<-min(like.nulls); max.null<-max(like.nulls); sd.null<-sd(like.nulls)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(like.nulls>obs) / length(like.nulls)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(like.nulls), xlim=c(min((obs), min(like.nulls)), max((obs), max(like.nulls))),
     main="network16"); abline(v=(obs), col="blue", lwd=2)

mod.network9 <- cbind(obs, sd.obs,zscore, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)
mod.network9

# network10

mod1.1<-computeModules(network17);mod1.2<-computeModules(network17);mod1.3<-computeModules(network17);mod1.4<-computeModules(network17);mod1.5<-computeModules(network17);mod1.6<-computeModules(network17);mod1.7<-computeModules(network17);mod1.8<-computeModules(network17);mod1.9<-computeModules(network17);mod1.10<-computeModules(network17)
mod.obs<-c(mod1.1@likelihood,mod1.2@likelihood,mod1.3@likelihood,mod1.4@likelihood,mod1.5@likelihood,mod1.6@likelihood,mod1.7@likelihood,mod1.8@likelihood,mod1.9@likelihood,mod1.10@likelihood)
obs<-max(mod.obs);sd.obs<-sd(mod.obs)

plotModuleWeb(mod1.1)

nulls <- nullmodel(network17, N=1000, method=1)
modules.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
zscore <- (obs - mean(like.nulls))/sd(like.nulls)

mean.null<-mean(like.nulls); min.null<-min(like.nulls); max.null<-max(like.nulls); sd.null<-sd(like.nulls)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(like.nulls>obs) / length(like.nulls)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(like.nulls), xlim=c(min((obs), min(like.nulls)), max((obs), max(like.nulls))),
     main="network17"); abline(v=(obs), col="blue", lwd=2)

mod.network10 <- cbind(obs, sd.obs,zscore, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)
mod.network10

# network11

mod1.1<-computeModules(network21);mod1.2<-computeModules(network21);mod1.3<-computeModules(network21);mod1.4<-computeModules(network21);mod1.5<-computeModules(network21);mod1.6<-computeModules(network21);mod1.7<-computeModules(network21);mod1.8<-computeModules(network21);mod1.9<-computeModules(network21);mod1.10<-computeModules(network21)
mod.obs<-c(mod1.1@likelihood,mod1.2@likelihood,mod1.3@likelihood,mod1.4@likelihood,mod1.5@likelihood,mod1.6@likelihood,mod1.7@likelihood,mod1.8@likelihood,mod1.9@likelihood,mod1.10@likelihood)
obs<-max(mod.obs);sd.obs<-sd(mod.obs)

plotModuleWeb(mod1.4)

nulls <- nullmodel(network21, N=1000, method=1)
modules.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
zscore <- (obs - mean(like.nulls))/sd(like.nulls)

mean.null<-mean(like.nulls); min.null<-min(like.nulls); max.null<-max(like.nulls); sd.null<-sd(like.nulls)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(like.nulls>obs) / length(like.nulls)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(like.nulls), xlim=c(min((obs), min(like.nulls)), max((obs), max(like.nulls))),
     main="network21"); abline(v=(obs), col="blue", lwd=2)

mod.network11 <- cbind(obs, sd.obs,zscore, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)
mod.network11

# network12

mod1.1<-computeModules(network23);mod1.2<-computeModules(network23);mod1.3<-computeModules(network23);mod1.4<-computeModules(network23);mod1.5<-computeModules(network23);mod1.6<-computeModules(network23);mod1.7<-computeModules(network23);mod1.8<-computeModules(network23);mod1.9<-computeModules(network23);mod1.10<-computeModules(network23)
mod.obs<-c(mod1.1@likelihood,mod1.2@likelihood,mod1.3@likelihood,mod1.4@likelihood,mod1.5@likelihood,mod1.6@likelihood,mod1.7@likelihood,mod1.8@likelihood,mod1.9@likelihood,mod1.10@likelihood)
obs<-max(mod.obs);sd.obs<-sd(mod.obs)

plotModuleWeb(mod1.6)

nulls <- nullmodel(network23, N=1000, method=1)
modules.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
zscore <- (obs - mean(like.nulls))/sd(like.nulls)

mean.null<-mean(like.nulls); min.null<-min(like.nulls); max.null<-max(like.nulls); sd.null<-sd(like.nulls)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(like.nulls>obs) / length(like.nulls)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(like.nulls), xlim=c(min((obs), min(like.nulls)), max((obs), max(like.nulls))),
     main="network23"); abline(v=(obs), col="blue", lwd=2)

mod.network12 <- cbind(obs, sd.obs,zscore, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)
mod.network12

# network13

mod1.1<-computeModules(network25);mod1.2<-computeModules(network25);mod1.3<-computeModules(network25);mod1.4<-computeModules(network25);mod1.5<-computeModules(network25);mod1.6<-computeModules(network25);mod1.7<-computeModules(network25);mod1.8<-computeModules(network25);mod1.9<-computeModules(network25);mod1.10<-computeModules(network25)
mod.obs<-c(mod1.1@likelihood,mod1.2@likelihood,mod1.3@likelihood,mod1.4@likelihood,mod1.5@likelihood,mod1.6@likelihood,mod1.7@likelihood,mod1.8@likelihood,mod1.9@likelihood,mod1.10@likelihood)
obs<-max(mod.obs);sd.obs<-sd(mod.obs)

plotModuleWeb(mod1.1)

nulls <- nullmodel(network25, N=1000, method=1)
modules.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
zscore <- (obs - mean(like.nulls))/sd(like.nulls)

mean.null<-mean(like.nulls); min.null<-min(like.nulls); max.null<-max(like.nulls); sd.null<-sd(like.nulls)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(like.nulls>obs) / length(like.nulls)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(like.nulls), xlim=c(min((obs), min(like.nulls)), max((obs), max(like.nulls))),
     main="network25"); abline(v=(obs), col="blue", lwd=2)

mod.network13 <- cbind(obs, sd.obs,zscore, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)
mod.network13

# network14

mod1.1<-computeModules(network28);mod1.2<-computeModules(network28);mod1.3<-computeModules(network28);mod1.4<-computeModules(network28);mod1.5<-computeModules(network28);mod1.6<-computeModules(network28);mod1.7<-computeModules(network28);mod1.8<-computeModules(network28);mod1.9<-computeModules(network28);mod1.10<-computeModules(network28)
mod.obs<-c(mod1.1@likelihood,mod1.2@likelihood,mod1.3@likelihood,mod1.4@likelihood,mod1.5@likelihood,mod1.6@likelihood,mod1.7@likelihood,mod1.8@likelihood,mod1.9@likelihood,mod1.10@likelihood)
obs<-max(mod.obs);sd.obs<-sd(mod.obs)

plotModuleWeb(mod1.7)

nulls <- nullmodel(network28, N=1000, method=1)
modules.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
zscore <- (obs - mean(like.nulls))/sd(like.nulls)

mean.null<-mean(like.nulls); min.null<-min(like.nulls); max.null<-max(like.nulls); sd.null<-sd(like.nulls)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(like.nulls>obs) / length(like.nulls)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(like.nulls), xlim=c(min((obs), min(like.nulls)), max((obs), max(like.nulls))),
     main="network28"); abline(v=(obs), col="blue", lwd=2)

mod.network14 <- cbind(obs, sd.obs,zscore, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)
mod.network14

# network15

mod1.1<-computeModules(network30);mod1.2<-computeModules(network30);mod1.3<-computeModules(network30);mod1.4<-computeModules(network30);mod1.5<-computeModules(network30);mod1.6<-computeModules(network30);mod1.7<-computeModules(network30);mod1.8<-computeModules(network30);mod1.9<-computeModules(network30);mod1.10<-computeModules(network30)
mod.obs<-c(mod1.1@likelihood,mod1.2@likelihood,mod1.3@likelihood,mod1.4@likelihood,mod1.5@likelihood,mod1.6@likelihood,mod1.7@likelihood,mod1.8@likelihood,mod1.9@likelihood,mod1.10@likelihood)
obs<-max(mod.obs);sd.obs<-sd(mod.obs)

plotModuleWeb(mod1.7)

nulls <- nullmodel(network30, N=1000, method=1)
modules.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
zscore <- (obs - mean(like.nulls))/sd(like.nulls)

mean.null<-mean(like.nulls); min.null<-min(like.nulls); max.null<-max(like.nulls); sd.null<-sd(like.nulls)
upperCI<-mean.null+(qnorm(0.975)*sd.null/sqrt(1000))
lowerCI<-mean.null-(qnorm(0.975)*sd.null/sqrt(1000))
praw <- sum(like.nulls>obs) / length(like.nulls)
pvalue<- ifelse(praw > 0.5, 1-praw, praw); pvalue
plot(density(like.nulls), xlim=c(min((obs), min(like.nulls)), max((obs), max(like.nulls))),
     main="network30"); abline(v=(obs), col="blue", lwd=2)

mod.network15 <- cbind(obs, sd.obs,zscore, mean.null, min.null, max.null, sd.null, upperCI, lowerCI, pvalue)
mod.network15

####Network metrics results####

#results H2
res.network1
res.network2
res.network3
res.network4
res.network5
res.network6
res.network7
res.network8
res.network9
res.network10
res.network11
res.network12
res.network13
res.network14
res.network15


#results Modularity

mod.network1
mod.network2
mod.network3
mod.network4
mod.network5
mod.network6
mod.network7
mod.network8
mod.network9
mod.network10
mod.network11
mod.network12
mod.network13
mod.network14
mod.network15

####Structural Equation Models####
require(piecewiseSEM)

dados<-read.table("Models.txt", h=T)
names(dados)

shapiro.test(dados$DistMedBorda)
shapiro.test(dados$nButHa)
shapiro.test(dados$DistVMed)
shapiro.test(dados$CArboPlog)

shapiro.test(dados$H2)
shapiro.test(dados$Mod_zs)
shapiro.test(dados$abdPM)
shapiro.test(dados$abdPh)
shapiro.test(dados$abdCr)
shapiro.test(dados$riqM)
shapiro.test(dados$riqPh)
shapiro.test(dados$riqCr)


#H2'
modC<- psem(
  lm(H2 ~ DistVMed + CArboPlog + nButHa + DistMedBorda, data = dados),
  lm(CArboPlog ~ DistMedBorda + nButHa, data = dados),
  DistVMed %~~% nButHa,
  DistMedBorda %~~% DistVMed
); summary(modC)

#Q_zs
modC<- psem(
  lm(Mod_zs ~ DistVMed + CArboPlog + nButHa + DistMedBorda, data = dados),
  lm(CArboPlog ~ DistMedBorda + nButHa, data = dados),
  DistVMed %~~% nButHa
); summary(modC)

#Pollinators mean abundance (log)
abdPM<- psem(
  lm(log(abdPM) ~ DistVMed + CArboPlog + nButHa + DistMedBorda, data = dados),
  lm(CArboPlog ~ DistMedBorda + nButHa, data = dados),
  DistVMed %~~% nButHa
); summary(abdPM)

#Peripheral pollinators mean abundance
abdPh<- psem(
  lm(abdPh ~ DistVMed + CArboPlog + nButHa + DistMedBorda, data = dados),
  lm(CArboPlog ~ DistMedBorda + nButHa, data = dados),
  DistVMed %~~% nButHa
); summary(abdPh)

#Core pollinators mean abundance (log)
abdCr<- psem(
  lm(log(abdCr) ~ DistVMed + CArboPlog + nButHa + DistMedBorda, data = dados),
  lm(CArboPlog ~ DistMedBorda + nButHa, data = dados),
  DistVMed %~~% nButHa
); summary(abdCr)

#Pollinators mean richness
riqM<- psem(
  lm(riqM ~ DistVMed + CArboPlog + nButHa + DistMedBorda, data = dados),
  lm(CArboPlog ~ DistMedBorda + nButHa, data = dados),
  DistVMed %~~% nButHa
); summary(riqM)

#Peripheral pollinators mean richness (log)
riqPh<- psem(
  lm(log(riqPh) ~ DistVMed + CArboPlog + nButHa + DistMedBorda, data = dados),
  lm(CArboPlog ~ DistMedBorda + nButHa, data = dados),
  DistVMed %~~% nButHa
); summary(riqPh)

#Core pollinators mean richness
riqCr<- psem(
  lm(riqCr ~ DistVMed + CArboPlog + nButHa + DistMedBorda, data = dados),
  lm(CArboPlog ~ DistMedBorda + nButHa, data = dados),
  DistVMed %~~% nButHa
); summary(riqCr)

