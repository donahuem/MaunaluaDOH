#Analysis of Maunalua DOH data from August 2019

library(tidyverse)
library(vegan)
library(BiodiversityR)
library(lubridate)
library(GGally)
library(ggplot2)
library(MASS)
library(ellipse)
library(viridis)

# Get water sample data -----------------------------------------------------

source("MBayDOH_curation.R")

#ADD IN NEWER DATASET FROM KIM - V14

# Name useful subsets of variables ----------------------------------------

c.bgc <- c("depth","temp","salinity", "pcDO","pH","TA","turbidity","TSS")
c.nuts <- c("logN","logtotN","logNH4","logtotP","logSi","chla")
c.fdom <- c("logBIX","logHIX","logMC", "logFI")
c.ben <- c("silt","sand","rock","pavement","coral","algae.nc","algae.calc","cca","other")
c.phrm <-  c("IBU","SMX","CBZ","GLY")
#data collected "everywhere"
c.ew <- c("temp","salinity", "pcDO","pH","TA","logturb","logN","logtotN",
          "logNH4","logtotP","logSi","logchla","logBIX","logHIX","logMC", "logFI")
#all response vars
c.all <- c(c.bgc,c.nuts,c.fdom,c.ben,c.phrm)
#driver vars
c.drv <- c("depth","offshore","longshore","otpN","sgd","streams")
c.sgd <- c("temp","salinity","logSi","HIX")
i.xtr <- d$site==267|d$site==277|d$site==263


#Clustering
library(ade4)
library(vegan)
library(gclus)
library(cluster)
library(RColorBrewer)
library(labdsv)
#library(mvpart)
#library(MVPARTwrap)

#first, cluster with all sites using c.ew (vars measured everywhere)
r.ew <- complete.cases(d[,c.ew])
d.ew <- d[r.ew,c("site",c.ew)]
ew.norm <- decostand(d.ew,"normalize")
ew.stand <- decostand(d.ew,"standardize")
ew.ch <- vegdist(ew.norm,"euc")

#hierarchical clustering(single linkage agglomerative)
ew.ch.single <- hclust(ew.ch, method="single")
plot(ew.ch.single)
ew.ch.single.cophen <- cophenetic(ew.ch.single)
cor(ew.ch, ew.ch.single.cophen)  #.755
gow.dist.single <- sum((ew.ch-ew.ch.single.cophen)^2) #15.14

#complete linkage clustering
ew.ch.complete <- hclust(ew.ch,method="complete")
plot(ew.ch.complete)
ew.ch.complete.cophen <- cophenetic(ew.ch.complete)
cor(ew.ch, ew.ch.complete.cophen) #0.781
gow.dist.complete <- sum((ew.ch-ew.ch.complete.cophen)^2) #24.5
ew.complete.g6 <- cutree(ew.ch.complete,6)

#UPGMA clustering (avg agglomerative)  -THIS LOOKS LIKE THE BEST FOR HIERARCHICAL
ew.ch.UPGMA <- hclust(ew.ch, method="average")
plot(ew.ch.UPGMA)
ew.ch.UPGMA.cophen <- cophenetic(ew.ch.UPGMA)
cor(ew.ch, ew.ch.UPGMA.cophen) #0.786
gow.dist.UPGMA <- sum((ew.ch-ew.ch.UPGMA.cophen)^2)  #2.66
#identify best number of groups for this method
asw=numeric(nrow(d.ew))
for(k in 2:(nrow(d.ew)-1)) {
  sil <- silhouette(cutree(ew.ch.UPGMA,k=k),ew.ch)
  asw[k] <- summary(sil)$avg.width
}
plot(1:nrow(d.ew),asw,type="h")
ew.chwo <- reorder.hclust(ew.ch.UPGMA,ew.ch)
ew.dend <- as.dendrogram(ew.chwo)
ew.colorbydepth <- ifelse(d.ew$depth<=2,"cyan","blue")
heatmap(as.matrix(ew.ch),Rowv=ew.dend, symm=TRUE,margin=c(3,3),RowSideColors=ew.colorbydepth)

#min variance clustering (Ward)
ew.ch.ward <- hclust(ew.ch, method="ward.D")
plot(ew.ch.ward)
ew.ch.ward.cophen <- cophenetic(ew.ch.ward)
cor(ew.ch, ew.ch.ward.cophen)  #0.757
sum((ew.ch-ew.ch.ward.cophen)^2)

#k-means partitioning - normalized:  12 (or 14) groups
ew.km.cascade <- cascadeKM(ew.norm,inf.gr=2,sup.gr=20,iter=100,criterion="ssi")
plot(ew.km.cascade,sortg=TRUE)
ew.kmeans <- kmeans(ew.norm, centers=12,nstart=100)
d.ew[order(ew.kmeans$cluster),]

#k-means partitioning - standardized
ew.s.km.cascade <- cascadeKM(ew.stand,inf.gr=2,sup.gr=20,iter=100,criterion="ssi")
plot(ew.s.km.cascade,sortg=TRUE)
ew.s.kmeans <- kmeans(ew.stand, centers=13,nstart=100)
d.ew[order(ew.kmeans$cluster),]



#PAM clustering
asw.pam <- numeric(nrow(d.ew))
for (k in 2:(nrow(d.ew)/2)){
  asw.pam[k] <- pam(ew.ch,k,diss=TRUE)$silinfo$avg.width
  k.best <- which.max(asw.pam)
  plot(1:nrow(d.ew),asw.pam,type="h")
}

#Compare Clusters from kmeans partitioning (12 groups)
d.ew$km <- as.factor(ew.kmeans$cluster)
d.ew.loc <- left_join(d.ew,d[,c("site","lat","long")])
write.csv(d.ew.loc,"Data_EWloc_20210909.csv")

d.ew.long <-
  pivot_longer(d.ew,depth:logFI,names_to = "env_parm",values_to="env_data")

d.ew.long %>% 
  ggplot(aes(y=env_data)) +
    facet_wrap(vars(env_parm),scales="free") +
    geom_boxplot(aes(x=km, group=km,color=km))
  
d.ew.loc %>%  ggplot(aes(kw)) +
  #scale_colour_viridis(option = "D",direction=-1,begin =0.15,end=0.85) +
  geom_point(mapping=aes(x=long,y=lat,color=km)) 


# Correlation Plots -------------------------------------------------------
pairs(d[(d$benthic==0 & !i.xtr),c.ew])
e<- ggplot(data=d)
e + geom_text(mapping=aes(x=logN,y=log(BIX),label=site))
e + geom_text(mapping=aes(x=BIX,y=FI,label=site))
e + geom_text(mapping=aes(x=logtotP,y=logtotN,label=site))
e + geom_text(mapping=aes(x=logN,y=logtotN,label=site))
e + geom_text(mapping=aes(x=logN,y=logtotP,label=site))
e + geom_text(mapping=aes(x=logtotP,y=logSi,label=site))
e + geom_text(mapping=aes(x=logtotN,y=NH4,label=site))
e + geom_text(mapping=aes(x=logturb,y=log(TSS),label=site))
e + geom_text(mapping=aes(x=logturb,y=silt,label=site))
e + geom_text(mapping=aes(x=DO,y=pcDO,label=site))
e + geom_text(mapping=aes(x=TA,y=pH,label=site))
e + geom_text(mapping=aes(x=TA,y=pH,label=site))

ei <- ggplot(data=d[!i.xtr,])
ei + geom_text(mapping=aes(x=logtotN,y=NH4,label=site))

# Histograms and maps of vars from DOH points ------------------------------------------------------

#DOH map for nitrogen
#plot log10(N) transform and histogram
filter(d, (doh>0)) %>%
  ggplot() + 
  geom_histogram(aes(logN))
#exclude logN>1 from color scale, plot high values as x's
filter(d, (doh>0 & logN<1)) %>%
  ggplot() +
  scale_colour_viridis(option = "D",direction=-1,begin =0.15,end=.8) +
  geom_point(mapping=aes(x=long,y=lat,color=logN))+ 
  filter(d, (doh>0 & logN>=1)) %>%
  geom_text(mapping=aes(x=long,y=lat,label="x"))

#DOH map for totP
#log transform and histogram
filter(d, (doh>0)) %>%
  ggplot() + 
  geom_histogram(aes(logtotP,color=doh))
#no exclusions
filter(d, (doh>0)) %>%
  ggplot() +
  scale_colour_viridis(option = "D",direction=-1,begin =0.15,end=0.85) +
  geom_point(mapping=aes(x=long,y=lat,color=logtotP))

#DOH map for NH4
#log transform and histogram
filter(d, (doh>0)) %>%
  ggplot() + 
  geom_histogram(aes(NH4))
#no exclusions
filter(d, (doh>0)) %>%
  ggplot() +
  scale_colour_viridis(option = "D",direction=-1,begin =0.15,end=0.85) +
  geom_point(mapping=aes(x=long,y=lat,color=NH4))

#DOH map for silicate
#log transform and histogram
filter(d, (doh>0)) %>%
  ggplot() + 
  geom_histogram(aes(logSi,color=doh))
#exclude logSi >3 (?)
filter(d, (doh>0)) %>%  # && logSi > 2.3)) %>%
  ggplot() +
  scale_colour_viridis(option = "D",direction=-1,begin =0.15,end=0.85) +
  geom_point(mapping=aes(x=long,y=lat,color=logSi))

#DOH map for turbidity
#log transform and histogram
filter(d, (doh>0)) %>%
  ggplot() + 
  geom_histogram(aes(logturb)) 

filter(d, (doh>0)) %>%
  ggplot() +
  scale_colour_viridis(option = "D",direction=-1,begin =0.15,end=0.85) +
  geom_point(mapping=aes(x=long,y=lat,color=logturb))

#DOH map for TSS
#log transform and histogram
filter(d, (doh>0 & TSS>5)) %>%
  ggplot() + 
  geom_histogram(aes(log(TSS))) 

filter(d, (doh>0 & TSS>5)) %>%
  ggplot() +
  scale_colour_viridis(option = "D",direction=-1,begin =0.15,end=0.85) +
  geom_point(mapping=aes(x=long,y=lat,color=log(TSS)))

#DOH map for chla
#log transform and histogram
filter(d, (doh>0)) %>%
  ggplot() + 
  geom_histogram(aes(logchla)) 

filter(d, (doh>0)) %>%
  ggplot() +
  scale_colour_viridis(option = "D",direction=-1,begin =0.15,end=0.85) +
  geom_point(mapping=aes(x=long,y=lat,color=logchla))

#DOH map for HIX
#log transform and histogram
filter(d, (doh>0)) %>%
  ggplot() + 
  geom_histogram(aes(log(HIX))) 

filter(d, (doh>0)) %>%
  ggplot() +
  scale_colour_viridis(option = "D",direction=-1,begin =0.15,end=0.85) +
  geom_point(mapping=aes(x=long,y=lat,color=log(HIX)))

#map for BIX,HIX, FI
filter(d, doh>0) %>% 
  ggplot() + 
  scale_color


#DOH map for pcDO
#log transform and histogram
filter(d, (doh>0)) %>%
  ggplot() + 
  geom_histogram(aes(pcDO)) 

filter(d, (doh>0)) %>%
  ggplot() +
  scale_colour_viridis(option = "D",direction=-1,begin =0.15,end=0.85) +
  geom_point(mapping=aes(x=long,y=lat,color=pcDO))

#DOH map for DO
#log transform and histogram
filter(d, (doh>0)) %>%
  ggplot() + 
  geom_histogram(aes(DO)) 

filter(d, (doh>0)) %>%
  ggplot() +
  scale_colour_viridis(option = "D",direction=-1,begin =0.15,end=0.85) +
  geom_point(mapping=aes(x=long,y=lat,color=DO))

#DOH map for salinity
#log transform and histogram
filter(d, (doh>0)) %>%
  ggplot() + 
  geom_histogram(aes(salinity)) 
#exclude salinities below 30
filter(d, (doh>0 & salinity>30)) %>%
  ggplot() +
  scale_colour_viridis(option = "D",direction=-1,begin =0.15,end=0.85) +
  geom_point(mapping=aes(x=long,y=lat,color=salinity)) +
  filter(d, (doh>0 & salinity<30)) %>%
  geom_text(mapping=aes(x=long,y=lat,label="x"))

#Map for silt

filter(d, (doh>0)) %>%
  ggplot() + 
  geom_histogram(aes(silt)) 

filter(d, (doh>0)) %>%
  ggplot() +
  scale_colour_viridis(option = "D",direction=-1,begin =0.15,end=0.85) +
  geom_point(mapping=aes(x=long,y=lat,color=silt))

filter(d, (doh>0 & benthic==0)) %>%
  ggplot() +
  scale_colour_viridis(option = "D",direction=-1,begin =0.15,end=0.85) +
  geom_point(mapping=aes(x=d$long,y=d$lat,color=(!is.na(d$GLY))*1))

#Map for depth

filter(d, (doh>0)) %>%
  ggplot() + 
  geom_histogram(aes(depth)) 

filter(d, (doh>0 & benthic==0)) %>%
  ggplot() +
  scale_colour_viridis(option = "D",direction=-1,begin =0.15,end=0.85) +
  geom_point(mapping=aes(x=long,y=lat,color=depth))

# Calculate distance matrices ---------------------------------------------

#Create distance matrices& clusters for columns 9:36, excluding pharms
ind<-complete.cases(d[,c.ew])

d.euclid <- vegdist(d[ind,c.ew],method="euclidean")
clust.euclid <-hclust(d.euclid)
plot(clust.euclid,labels=d$site[ind],cex=0.5)

d.scl.eu <- vegdist(scale(d[ind,c.ew]),method="euclidean")
clust.scl.eu <- hclust(d.scl.eu)
plot(clust.scl.eu,labels=d$site[ind],cex=0.5)

d.bray <- vegdist(d[ind,c.ew],method="bray",na.rm=TRUE)
clust.bray <- hclust(d.bray)
plot(clust.bray,labels=d$site[ind],cex=0.5)

#PCA including benthic
d.pca <- rda(scale(d[ind,c.ew]),scale=TRUE)
summary(d.pca)

biplot(d.pca,scaling=1,display=c("sites","species"),type=c("text","points"))


#Exclude Bottom Samples
ind_nb<-complete.cases(d[,c.ew])&(d$benthic==0)
d.nb.euclid <- vegdist(d[ind_nb,c.ew],method="euclidean")
clust.nb.euclid <-hclust(d.nb.euclid)
plot(clust.nb.euclid,labels=d$site[ind_nb],cex=0.5)

d.nb.scl.eu <- vegdist(scale(d[ind_nb,c.ew]),method="euclidean")
clust.nb.scl.eu <- hclust(d.nb.scl.eu)
plot(clust.nb.scl.eu,labels=d$site[ind_nb],cex=0.5)

d.nb.bray <- vegdist(d[ind_nb,c.ew],method="bray",na.rm=TRUE)
clust.nb.bray <- hclust(d.nb.bray)
plot(clust.nb.bray,labels=d$site[ind_nb],cex=0.5)

#PCA excluding benthic
d.nb.pca <- rda(scale(d[ind_nb,c.ew]),scale=TRUE)
summary(d.nb.pca)

biplot(d.nb.pca,scaling=1,display=c("sites","species"),type=c("text","points"))
#note - very little difference with and without bottom samples

#PCA on only fdom
d.fdom.pca <- rda(scale(d[ind_nb,c.fdom]),scale=TRUE)
summary(d.fdom.pca)
biplot(d.fdom.pca,scaling=1,display=c("sites","species"),type=c("text","points"))
#MC and HIX are pretty inverse

d.fdom.pca <- rda(scale(d[ind_nb,c("logHIX","logBIX","logFI")]),scale=TRUE)
summary(d.fdom.pca)
biplot(d.fdom.pca,scaling=1,display=c("sites","species"),type=c("text","points"))


#par(mfrow = c(3,2))
for (i in seq(1,4)){
  for (j in seq(2,4)){
    print(i)
    print(j)
    if (j>i){
    biplot(d.nb.pca,scaling=1,display=c("sites","species"),type=c("text","points"),choices=c(i,j))
    }
  }
}

#NMDS
indexR <- ind_nb
d.nb.nmds <- metaMDS(d[indexR,c.ew],distance = "bray")
par(mfrow=c(1,1))
stressplot(d.nb.nmds,main="Shepard plot")
gof=goodness(d.nb.nmds)
plot(d.nb.nmds,type="n", main="NMDS")
points(d.nb.nmds,display="sites",col=d$clust[indexR],pch=16)
text(d.nb.nmds, display="species",cex=.8)

#add Ward cluster
d.nb.nmds.ward <- hclust(d.nb.bray,"ward.D")
d.nb.nmds.ward.groups <- cutree(d.nb.nmds.ward,k=5)
grp.lev <-levels(factor(d.nb.nmds.ward.groups))
sit.sc <- scores(d.nb.nmds)
p<- ordiplot(sit.sc, type="n", main="NMDS/Bray - clusters Ward/Bray")
for (i in 1:length(grp.lev)){
  points(sit.sc[d.nb.nmds.ward.groups==i,],pch=(14-i),cex=2,col=(i+1))
}

#RDA - explanatory vars
#we consider 4 explanatory variables:  dist2coast, distHK, sgd, otpN, streams
#we consider all RVs except pharmaceuticals (cols_nph),then physical vars (cols_phys),
#  and biochem vars (bc1 = nuts+carbonate chem, bc2=cells + fDOM), and pharma
rda.nph <- rda(scale(d[ind_nb,c.ew]) ~ d$sgdx[ind_nb]+d$otpN[ind_nb] + d$streamsx[ind_nb]
               +d$offshore[ind_nb] + d$longshore[ind_nb])
summary(rda.nph)
RsquareAdj(rda.nph)$r.squared
plot(rda.nph,scaling=1,main="Triplot RDA rda.nph")

###RDA FIGURES FOR OCEAN SCIENCES (toggle indices for nearshore and pharm plots)
#index for nearshore sites (exclude Kim's clusters 2,11)
ind_near<-complete.cases(d[,c.ew])&(d$benthic==0) #&(d$clust!=2)&(d$clust!=11)
#index for pharm sites
ind_phrm <-complete.cases(d[,c.phrm])&(d$benthic==0)

indexR= ind_near #ind_phrm #
indexC= c.ew #c.phrm#
yRDA <- data.frame(scale(d[indexR,indexC]))
xRDA <- data.frame(SGD=d$sgdx[indexR],OSDS=d$otpN[indexR],STREAM=d$streamsx[indexR],
                   OFFSHORE=d$offshore[indexR], LONGSHORE=d$longshore[indexR])
attach(xRDA)
mod.rda <- rda(yRDA ~ SGD + OSDS + STREAM + OFFSHORE) 
summary(mod.rda)
RsquareAdj(mod.rda)$r.squared
mod.rda.sc <- scores(mod.rda,choices=1:2,display=c("sp","cn","bp"),scaling=2)
plot(mod.rda,choices=1:2,display=c("sp","cn","bp"),scaling=2,type="none",lwd=2)
xarrow=sign(mod.rda.sc$species[,1])*(abs(mod.rda.sc$species[,1])-0.02)
yarrow=sign(mod.rda.sc$species[,2])*(abs(mod.rda.sc$species[,2])-0.02)
arrows(0,0,xarrow,yarrow,length=0,lty=1,lwd=2,col="red")
text(mod.rda,dis="cn",scaling=2,lwd=2)
text(mod.rda,"species",col="red",cex=0.9,scaling=2)
detach(xRDA)

##Linear Discriminant Analysis using Kim's clusters
#Back to using sites that exclude bottom (ind_nb) and all non-pharm sites
indexR <- ind_near
indexC <- c.ew

xLDA <- data.frame(d[indexR,c.ew])
xLDA$logN <- log(xLDA$N)
xLDA$logNH4 <- log(xLDA$NH4)
xLDA$logtotN <- log(xLDA$totN)
xLDA$logtotP <- log(xLDA$totP)
xLDA$logsi <- log(xLDA$silicate)
xLDA <-data.frame(scale(xLDA),na.rm=TRUE)
mod.lda <- lda(d$clust[indexR]~depth+temp+salinity+pcDO+pH+TA+turbidity+TSS+
                              logN+logtotN+logNH4+logtotP+logSi+chla+
                              FCM+BIX+HIX+MC+FI,data=xLDA)
Cs <- mod.lda$scaling #normalized eigenvectors
mod.lda$svd^2 #canonical eigenvalues
(Fp <- predict(mod.lda)$x)
mod.class <- predict(mod.lda)$class
mod.post <- predict(mod.lda)$posterior
mod.table <- table(d$clust[indexR],mod.class)
diag(prop.table(mod.table,1))

plot(Fp[,1],Fp[,2],type="n")
text(Fp[,1],Fp[,2],d$site[indexR],col=as.numeric(mod.class))
abline(v=0, lty="dotted")
abline(h=0, lty="dotted")
for(i in 1:length(levels(mod.class))){
  cov <-cov(Fp[mod.class==i,])
  centre<- apply(Fp[mod.class==i,],2,mean)
  lines(ellipse(cov,center=centre,level=0.95))
}

##This
