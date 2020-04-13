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

data <- read_csv("Data_MASTER_20200331_v9.csv")
data <- mutate_at(data,vars(TSS:Chla,EnteroCFU100ml),as.numeric)

#subset to data columns of use; reformat times and dates
d <- dplyr::select(data, site=SiteNum,lat=LatGPS,long=LongGPS,date=Date,time=Time,team=Team,
           benthic=Benthic, benthicpair=BenthicPair, clust = clusterfit,depth=Depth_FINAL,
           temp=TempFinal,salinity=Salinity_Silbiger,DO=ODO,pcDO=ODOPerFinal,pH=pH_insitu_Silbiger,TA=TA_Raw_Silbiger,
           chla=Chla,turbidity=Turb_Final,TSS=TSS,NH4=Ammonia,N=NNN,totN=TotalN,totP=TotalP,silicate=Silica,
           FCM=`FCM_Concentration (cells/uL)`,BIX=fDOM_BIX,HIX=fDOM_HIX,FI=fDOM_FI,MC=`fDOM_M:C`,
           IBU=IBU_ngmL,SMX=SMX_ngmL,CBZ=CBZ_ngmL,GLY=GLY_ngmL,
           offshore=dist2coast,longshore=DistanceHK_m,otpN=otpN,sgd=gw_wells,streams=streams)
d <- mutate(d,doh=0+1*(site>=1 & site<=160)+2*(site>=201 & site<=260))  #DOH offshore=1, DOH nearshore=2
d <- mutate(d,date = ymd(date))
d <- mutate(d,datetime = as.POSIXct(paste(d$date,d$time), format="%Y-%m-%d %H:%M:%S"))
d <- mutate(d,sgdx = exp(-0.001*d$sgd))
d <- mutate(d,streamsx = exp(-0.001*streams))
d <- mutate(d,logN = log10(d$N))
d <- mutate(d,logNH4 = log10(d$NH4))
d <- mutate(d,logtotN = log10(d$totN))
d <- mutate(d,logtotP = log10(d$totP))
d <- mutate(d,logSi = log10(d$silicate))
d <- mutate(d,logturb = log10(turbidity))
d <- mutate(d,logchla = log10(chla))

#Outliers
#high nutrient sites (N,P,Si):  267, 277,263 - likely to be real, close to SGD
#high BIX: 252
#high FI: 248
#THE FOLLOWING WERE FIXED IN V9 OF DATASET
#both FI(site=248) and BIX(site=252) are 10^4 larger than median have the exact value 65535; replace with NAs
#FI(site=251) has an exact value of 0, which I think is impossible and no other fDOM measurement is zero; replace with NA
#d$FI[d$site==248]<-NA
#d$BIX[d$site==252]<-NA
#d$FI[d$site==251] <-NA
#d$salinity[d$site==115]=35.16

#Benthic Cover Data
benthic <- read_csv("PercentCovers_PM.csv")
b <- dplyr::select(benthic,
                   quad = 'Quadrad number',
                   site = 'Site Number',
                   coral = Coral,
                   sand = Sand,
                   silt = Silt,
                   pavement = Pav,
                   rock = Rock,
                   unknown = Unk,
                   a.spic = ASPIC,
                   other.calc = CA_O,
                   cca = CCA,
                   algae.gfa = GFA,
                   g.sal = GSAL,
                   halimeda = Hal,
                   turf = Turf,
                   other.algae = UnkMacroa)
#simplify so the variables of interest are:  
#silt, sand, rock, pavement, coral, non-calcifying algae, calcifying non-crustose algae, cca, other
#non-calcifying algae = Gracilaria salicornia, Acanthophora spicifera, turf, and other macroalgae
#calcifying non-crustose = halimeda and padina (most of other.calc)
b <- mutate(b, algae.nc = a.spic+algae.gfa+turf+other.algae,
            algae.calc = halimeda+other.calc,
            other = unknown)
b.sum <- b %>% 
         group_by(site) %>% 
         summarise(silt=mean(silt),sand=mean(sand),rock=mean(rock),pavement=mean(pavement),
                   coral=mean(coral),algae.nc=mean(algae.nc),algae.calc=mean(algae.calc),
                   cca=mean(cca), other=mean(other))
d <- left_join(d,b.sum, by="site")

save(d,file="Datav9_plusbenthic.Rdata")

cols_phys <- c("depth","temp","salinity", "turbidity","TSS", "DO")
#cols_bc1 <- c("pH","TA","logN","logtotN","logNH4","logtotP","logSi")
cols_bc1 <- c("pH","TA","N","totN","NH4","totP","silicate")
cols_bc2 <- c("chla","FCM","BIX","HIX","MC", "FI")
cols_ben <- c("silt","sand","rock","pavement","coral","algae.nc","algae.calc","cca","other")
cols_phrm <-  c("IBU","SMX","CBZ","GLY")
cols_nph <- c(cols_phys,cols_bc1,cols_bc2)
cols_drv <- c("dist2coast","distHK","otpN","sgd","streams")
cols_sgd <- c("salinity","silicate","HIX")
ind_xtr <- d$site==267|d$site==277

#correlation plots
# ggpairs(d,columns=cols_phys)
# ggpairs(d,columns=cols_bc1)
# ggpairs(d,columns=cols_bc2)

e<- ggplot(data=d)
e + geom_text(mapping=aes(x=log(N),y=log(BIX),label=site))
e + geom_text(mapping=aes(x=BIX,y=FI,label=site))
e + geom_text(mapping=aes(x=log(totP),y=log(totN),label=site))
e + geom_text(mapping=aes(x=log(N),y=log(totP),label=site))
e + geom_text(mapping=aes(x=log(totP),y=log(silicate),label=site))
e + geom_text(mapping=aes(x=log(N),y=log(NH4),label=site))
e + geom_text(mapping=aes(x=log(turbidity),y=log(TSS),label=site))

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
  geom_point(mapping=aes(x=long,y=lat,color=logN))

#DOH map for silicate
#log transform and histogram
filter(d, (doh>0)) %>%
  ggplot() + 
  geom_histogram(aes(logSi,color=doh))
#exclude logSi >3 (?)
filter(d, (doh>0 & logSi < 2.3)) %>%
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
  geom_histogram(aes(HIX)) 

filter(d, (doh>0 & HIX<3)) %>%
  ggplot() +
  scale_colour_viridis(option = "D",direction=-1,begin =0.15,end=0.85) +
  geom_point(mapping=aes(x=long,y=lat,color=HIX))


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

#Create distance matrices& clusters for columns 9:36, excluding pharms
ind<-complete.cases(d[,cols_nph])
d.euclid <- vegdist(d[ind,cols_nph],method="euclidean")
clust.euclid <-hclust(d.euclid)
plot(clust.euclid,labels=d$site[ind],cex=0.5)

d.scl.eu <- vegdist(scale(d[ind,cols_nph]),method="euclidean")
clust.scl.eu <- hclust(d.scl.eu)
plot(clust.scl.eu,labels=d$site[ind],cex=0.5)

d.bray <- vegdist(d[ind,cols_nph],method="bray",na.rm=TRUE)
clust.bray <- hclust(d.bray)
plot(clust.bray,labels=d$site[ind],cex=0.5)

#Exclude Benthic
ind_nb<-complete.cases(d[,cols_nph])&(d$benthic==0)
d.nb.euclid <- vegdist(d[ind_nb,cols_nph],method="euclidean")
clust.nb.euclid <-hclust(d.nb.euclid)
plot(clust.nb.euclid,labels=d$site[ind_nb],cex=0.5)

d.nb.scl.eu <- vegdist(scale(d[ind_nb,cols_nph]),method="euclidean")
clust.nb.scl.eu <- hclust(d.nb.scl.eu)
plot(clust.nb.scl.eu,labels=d$site[ind_nb],cex=0.5)

d.nb.bray <- vegdist(d[ind_nb,cols_nph],method="bray",na.rm=TRUE)
clust.nb.bray <- hclust(d.nb.bray)
plot(clust.nb.bray,labels=d$site[ind_nb],cex=0.5)

#PCA excluding benthic
d.nb.pca <- rda(scale(d[ind_nb,cols_nph]),scale=TRUE)
summary(d.nb.pca)

biplot(d.nb.pca,scaling=1,display=c("sites","species"),type=c("text","points"))

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
d.nb.nmds <- metaMDS(d[indexR,cols_nph],distance = "bray")
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
rda.nph <- rda(scale(d[ind_nb,cols_nph]) ~ d$sgdx[ind_nb]+d$otpN[ind_nb] + d$streamsx[ind_nb]
               +d$offshore[ind_nb] + d$longshore[ind_nb])
summary(rda.nph)
RsquareAdj(rda.nph)$r.squared
plot(rda.nph,scaling=1,main="Triplot RDA rda.nph")

###RDA FIGURES FOR OCEAN SCIENCES (toggle indices for nearshore and pharm plots)
#index for nearshore sites (exclude Kim's clusters 2,11)
ind_near<-complete.cases(d[,cols_nph])&(d$benthic==0)&(d$clust!=2)&(d$clust!=11)
#index for pharm sites
ind_phrm <-complete.cases(d[,cols_phrm])&(d$benthic==0)

indexR=ind_phrm #ind_near
indexC=cols_phrm#cols_nph
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
indexC <- cols_nph

xLDA <- data.frame(d[indexR,cols_nph])
xLDA$logN <- log(xLDA$N)
xLDA$logNH4 <- log(xLDA$NH4)
xLDA$logtotN <- log(xLDA$totN)
xLDA$logtotP <- log(xLDA$totP)
xLDA$logsi <- log(xLDA$silicate)
xLDA <-data.frame(scale(xLDA),na.rm=TRUE)
mod.lda <- lda(d$clust[indexR]~depth+temp+salinity+turbidity+TSS+DO+
                              pH+TA+logN+logtotN+logNH4+logtotP+logsi+
                              FCM+BIX+HIX+MC,data=xLDA)
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