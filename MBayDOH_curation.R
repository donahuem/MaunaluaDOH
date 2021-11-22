#Data Curation

data <- read_csv("Data_MASTER_20200331_v9.csv")

#found mistake in TA data.  Site 235 (random inshore site) was measured immediately after site 181 (Paiko lagoon site)
#  and had the exact same, anomalous TA value as the Paiko site.  Paiko probably does have weird TA.  Exclude TA for site 235

data$TA_Raw_Silbiger[data$SiteNum==235] <- NA

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
d <- mutate(d,logHIX = log10(HIX))
d <- mutate(d,logBIX = log10(BIX))
d <- mutate(d,logMC = log10(MC))
d <- mutate(d,logFI = log10(FI))


# Get benthic data, curate, and merge -------------------------------------------------------------
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

#remove old geo vars from Kim
d <- subset(d,select=-c(longshore, offshore,otpN,sgd,streams,sgdx,streamsx))

#add geo data from Kim
d2.kim <- read_csv("ml_data_v14_addedflowlines2.csv",col_names=TRUE)  #updated from Kim Sept 2021; includes flow lines
d2.kim.geo <- subset(d2.kim, select =c(SiteNum, clusterfit:urbanperc))
names(d2.kim.geo)[1] <- 'site'
d <- left_join(d,d2.kim.geo,by='site')

rm(b,b.sum,benthic,d2.kim,d2.kim.geo,data)
