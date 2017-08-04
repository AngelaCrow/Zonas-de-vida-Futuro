#Ecoregion ID and names 
#9 = Grandes Planicies (GP)
#10 = Desiertos de America de Norte (DAN)
#11 = California Mediterranea (CM)
#12 = Elevaciones Semiaridas Meridionales (ESM)
#13 = Sierras Templadas (ST)
#14 = Selvas Calido-secas (SCS)
#15 = Selvas Calido-humedas (SCH)

estables45<-read.csv("Resample/HAD_ESM_LR_rcp45_2015_2039_bio/estables45.csv")
estables85<-read.csv("Resample/HAD_ESM_LR_rcp85_2015_2039_bio/estables85.csv")

cambio45<-read.csv("Resample/HAD_ESM_LR_rcp45_2015_2039_bio/cambio45.csv")
cambio85<-read.csv("Resample/HAD_ESM_LR_rcp85_2015_2039_bio/cambio85.csv")

mexbio.ecoanp<-read.csv("Mexbio_eco_stats.csv", header = T)
mexbio.ecoanp.sel<-c("X","Median_10")
mexbio.median <-mexbio.ecoanp[mexbio.ecoanp.sel]
head(mexbio.median)

dir.create("Resample/HAD_ESM_LR_rcp45_2015_2039_bio/median_mexbio")

GP_estables45 <- subset(estables45, ecoregiones == 9)
GP_estables45.mexbio <- GP_estables45[which(GP_estables45$MEXBIO_2010_gw_pr >= 0.1702225),]
GP_cambio45 <- subset(cambio45, ecoregiones == 9)
GP_cambio45.mexbio <- GP_cambio45[which(GP_cambio45$MEXBIO_2010_gw_pr >= 0.1702225),]

DAN_estables45 <- subset(estables45, ecoregiones == 10)
DAN_estables45.mexbio <- DAN_estables45[which(DAN_estables45$MEXBIO_2010_gw_pr >= 0.6300000),]
DAN_cambio45 <- subset(cambio45, ecoregiones == 10)
DAN_cambio45.mexbio <- DAN_cambio45[which(DAN_cambio45$MEXBIO_2010_gw_pr >= 0.6300000),]

CM_estables45 <- subset(estables45, ecoregiones == 11)
CM_estables45.mexbio <- CM_estables45[which(CM_estables45$MEXBIO_2010_gw_pr >= 0.5061778),]
CM_cambio45 <- subset(cambio45, ecoregiones == 11)
CM_cambio45.mexbio <- CM_cambio45[which(CM_cambio45$MEXBIO_2010_gw_pr >= 0.5061778),]

ESM_estables45 <- subset(estables45, ecoregiones == 12)
ESM_estables45.mexbio <- ESM_estables45[which(ESM_estables45$MEXBIO_2010_gw_pr >= 0.2754000),]
ESM_cambio45 <- subset(cambio45, ecoregiones == 12)
ESM_cambio45.mexbio <- ESM_cambio45[which(ESM_cambio45$MEXBIO_2010_gw_pr >= 0.2754000),]

ST_estables45 <- subset(estables45, ecoregiones == 13)
ST_estables45.mexbio <- ST_estables45[which(ST_estables45$MEXBIO_2010_gw_pr >= 0.3000000),]
ST_cambio45 <- subset(cambio45, ecoregiones == 13)
ST_cambio45.mexbio <- ST_cambio45[which(ST_cambio45$MEXBIO_2010_gw_pr >= 0.3000000),]

SCS_estables45 <- subset(estables45, ecoregiones == 14)
SCS_estables45.mexbio <- SCS_estables45[which(SCS_estables45$MEXBIO_2010_gw_pr >= 0.1056500),]
SCS_cambio45 <- subset(cambio45, ecoregiones == 14)
SCS_cambio45.mexbio <- SCS_cambio45[which(SCS_cambio45$MEXBIO_2010_gw_pr >= 0.1056500),]

SCH_estables45 <- subset(estables45, ecoregiones == 15)
SCH_estables45.mexbio <- SCH_estables45[which(SCH_estables45$MEXBIO_2010_gw_pr >= 0.1000000),]
SCH_cambio45 <- subset(cambio45, ecoregiones == 15)
SCH_cambio45.mexbio <- SCH_cambio45[which(SCH_cambio45$MEXBIO_2010_gw_pr >= 0.1000000),]


bajo.mexbio.estables45<-rbind(GP_estables45.mexbio,DAN_estables45.mexbio,
                         CM_estables45.mexbio,ESM_estables45.mexbio,
                         ST_estables45.mexbio,SCS_estables45.mexbio,
                         SCH_estables45.mexbio)
write.csv(bajo.mexbio.estables45, "Resample/HAD_ESM_LR_rcp45_2015_2039_bio/median_mexbio/bajo_mexbio_estables45.csv")

bajo.mexbio.cambio45<-rbind(GP_cambio45.mexbio,DAN_cambio45.mexbio,
                       CM_cambio45.mexbio,ESM_cambio45.mexbio,
                       ST_cambio45.mexbio,SCS_cambio45.mexbio,
                       SCH_cambio45.mexbio)
write.csv(bajo.mexbio.cambio45, "Resample/HAD_ESM_LR_rcp45_2015_2039_bio/median_mexbio/bajo_mexbio_cambio45.csv")

####RCP85####

dir.create("Resample/HAD_ESM_LR_rcp85_2015_2039_bio/median_mexbio")

GP_estables85 <- subset(estables85, ecoregiones == 9)
GP_estables85.mexbio <- GP_estables85[which(GP_estables85$MEXBIO_2010_gw_pr >= 0.1702225),]
GP_cambio85 <- subset(cambio85, ecoregiones == 9)
GP_cambio85.mexbio <- GP_cambio85[which(GP_cambio85$MEXBIO_2010_gw_pr >= 0.1702225),]

DAN_estables85 <- subset(estables85, ecoregiones == 10)
DAN_estables85.mexbio <- DAN_estables85[which(DAN_estables85$MEXBIO_2010_gw_pr >= 0.6300000),]
DAN_cambio85 <- subset(cambio85, ecoregiones == 10)
DAN_cambio85.mexbio <- DAN_cambio85[which(DAN_cambio85$MEXBIO_2010_gw_pr >= 0.6300000),]

CM_estables85 <- subset(estables85, ecoregiones == 11)
CM_estables85.mexbio <- CM_estables85[which(CM_estables85$MEXBIO_2010_gw_pr >= 0.5061778),]
CM_cambio85 <- subset(cambio85, ecoregiones == 11)
CM_cambio85.mexbio <- CM_cambio85[which(CM_cambio85$MEXBIO_2010_gw_pr >= 0.5061778),]

ESM_estables85 <- subset(estables85, ecoregiones == 12)
ESM_estables85.mexbio <- ESM_estables85[which(ESM_estables85$MEXBIO_2010_gw_pr >= 0.2754000),]
ESM_cambio85 <- subset(cambio85, ecoregiones == 12)
ESM_cambio85.mexbio <- ESM_cambio85[which(ESM_cambio85$MEXBIO_2010_gw_pr >= 0.2754000),]

ST_estables85 <- subset(estables85, ecoregiones == 13)
ST_estables85.mexbio <- ST_estables85[which(ST_estables85$MEXBIO_2010_gw_pr >= 0.3000000),]
ST_cambio85 <- subset(cambio85, ecoregiones == 13)
ST_cambio85.mexbio <- ST_cambio85[which(ST_cambio85$MEXBIO_2010_gw_pr >= 0.3000000),]

SCS_estables85 <- subset(estables85, ecoregiones == 14)
SCS_estables85.mexbio <- SCS_estables85[which(SCS_estables85$MEXBIO_2010_gw_pr >= 0.1056500),]
SCS_cambio85 <- subset(cambio85, ecoregiones == 14)
SCS_cambio85.mexbio <- SCS_cambio85[which(SCS_cambio85$MEXBIO_2010_gw_pr >= 0.1056500),]

SCH_estables85 <- subset(estables85, ecoregiones == 15)
SCH_estables85.mexbio <- SCH_estables85[which(SCH_estables85$MEXBIO_2010_gw_pr >= 0.1000000),]
SCH_cambio85 <- subset(cambio85, ecoregiones == 15)
SCH_cambio85.mexbio <- SCH_cambio85[which(SCH_cambio85$MEXBIO_2010_gw_pr >= 0.1000000),]


bajo.mexbio.estables85<-rbind(GP_estables85.mexbio,DAN_estables85.mexbio,
                              CM_estables85.mexbio,ESM_estables85.mexbio,
                              ST_estables85.mexbio,SCS_estables85.mexbio,
                              SCH_estables85.mexbio)
write.csv(bajo.mexbio.estables85, "Resample/HAD_ESM_LR_rcp85_2015_2039_bio/median_mexbio/bajo_mexbio_estables85.csv")

bajo.mexbio.cambio85<-rbind(GP_cambio85.mexbio,DAN_cambio85.mexbio,
                            CM_cambio85.mexbio,ESM_cambio85.mexbio,
                            ST_cambio85.mexbio,SCS_cambio85.mexbio,
                            SCH_cambio85.mexbio)
write.csv(bajo.mexbio.cambio85, "Resample/HAD_ESM_LR_rcp85_2015_2039_bio/median_mexbio/bajo_mexbio_cambio85.csv")