#Este codigo se debe usar para calcular la magnitud de cambio de las zonas de vida.
#calcula el cambio para cada variable de las que depende el diagrama de holdridge
#y un valor de vector total. Esto lo hace para cada pixel

library(raster)
library(reshape)
setwd("C:/CONABIO/CLIMA/LifeZones/LZ_R_uniatmos/")


#### Preparacion de matrices ####
#Importat raster de las zonas de vida de cada tiempo
#### Preparacion de matrices ####
#llamar raster de las zonas de vida de cada tiempo y mexbio
zvrcp45 <- raster("C:/CONABIO/CLIMA/LifeZones/LZ_R_uniatmos/Resample/CNRMCM5_rcp45_2015_2039_bio/zvhrcp45.tif")
zvrcp85 <- raster("C:/CONABIO/CLIMA/LifeZones/LZ_R_uniatmos/Resample/HAD_ESM_LR_rcp85_2015_2039_bio/zvhrcp85.tif")
zv <- raster("C:/CONABIO/CLIMA/LifeZones/LZ_R_uniatmos/Resample/Presente/zvh.tif")

mexbio2010 <- raster("F:/COBERTURAS/MexBio/MEXBIO_corregido/MEXBIO_2010_gw_pr.tif") #Global Human Footprint 
ecoregiones <- raster("C:/CONABIO/CLIMA/LifeZones/LZ_R_uniatmos/ecoregiones.tif")
compareRaster(zv, zvrcp45,mexbio2010,zvrcp85,ecoregiones)

matriz<-rasterToPoints(zv) #area de estudio
m<-matriz[,1:2]
zonas<-stack(zv,zvrcp45,zvrcp85, ecoregiones,mexbio2010) 
m.zonas<-extract(zonas,m) #matriz con zv del presente y futuro por pixel
coor.m.zonas<-as.data.frame(cbind(m,m.zonas))

#para separar los numeros del ID en columnas; LatReg (biotemperatura), 
#Ppt_P (Precipitación) y Phum(Piso de Humedad):

zvh<-as.factor(coor.m.zonas$zvh)
zvh <- colsplit(zvh, names=c("LatReg_P","Ppt_P","Phum_P"))

zvh45<-as.factor(coor.m.zonas$zvhrcp45)
zvh45 <- colsplit(zvh45, names=c("LatReg_45","Ppt_45","Phum_45"))

zvh85<-as.factor(coor.m.zonas$zvhrcp85)
zvh85 <- colsplit(zvh85, names=c("LatReg_85","Ppt_85","Phum_85"))

#Matriz con coordenadas (x y), codigos de las zonas de vida en el presente y futuro (Rcp 4.5 y 8.5)
# más las las columnas de las tres variables para cada tiempo

m_final <-as.data.frame(cbind(coor.m.zonas,zvh,zvh45,zvh85))
names(m_final)

#Rcp 4.5 - Calcular la diferencia entre los vectores de cada variable entre tiempos.
Lat.vector45<-abs(data.frame(m_final$LatReg_P - m_final$LatReg_45))
Ppt.vector45<-abs(data.frame(m_final$Ppt_P - m_final$Ppt_45))
Phum.vector45<-abs(data.frame(m_final$Phum_P - m_final$Phum_45))

#Rcp 4.5 -Se adiciona una columna: suma.vec45, con la cual es posible identificar 
#las zonas que cambian (suma.vec45 = número diferente a 0)
#o que se mantienen estables (suma.vec45 = 0)
suma.vec45<-(Lat.vector45+Ppt.vector45+Phum.vector45)
names(suma.vec45)
suma.vec45<-rename(suma.vec45,c(m_final.LatReg_P...m_final.LatReg_45 ="suma.vec45"))

#Rcp 8.5 - Calcular la diferencia entre los vectores de cada variable entre tiempos.
Lat.vector85<-abs(data.frame(m_final$LatReg_P - m_final$LatReg_85))
Ppt.vector85<-abs(data.frame(m_final$Ppt_P - m_final$Ppt_85))
Phum.vector85<-abs(data.frame(m_final$Phum_P - m_final$Phum_85))

#Rcp 8.5 -Se adiciona una columna: suma.vec85, con la cual es posible identificar 
#las zonas que cambian (suma.vec85 = número diferente a 0)
#o que se mantienen estables (suma.vec85 = 0)
suma.vec85<-(Lat.vector85+Ppt.vector85+Phum.vector85)
names(suma.vec85)
suma.vec85<-rename(suma.vec85,c(m_final.LatReg_P...m_final.LatReg_85 ="suma.vec85"))

m_final45<-as.data.frame(cbind(m_final,
                               Lat.vector45,
                               Ppt.vector45,
                               Phum.vector45,
                               suma.vec45 ))
names(m_final45)
m_final45<-rename(m_final45, c(m_final.LatReg_P...m_final.LatReg_45="Lat.vec45",
                             m_final.Ppt_P...m_final.Ppt_45="Ppt.vec45",
                             m_final.Phum_P...m_final.Phum_45="Phum45.vec"
                           ))

m_final85<-as.data.frame(cbind(m_final,
                               Lat.vector85,
                               Ppt.vector85,
                               Phum.vector85,
                               suma.vec85 ))
names(m_final85)
m_final85<-rename(m_final85, c(m_final.LatReg_P...m_final.LatReg_85="Lat.vec85",
                             m_final.Ppt_P...m_final.Ppt_85="Ppt.vec85",
                             m_final.Phum_P...m_final.Phum_85="Phum85.vec"
))

#Identificar zonas de clima estable (suma.vec45 = 0)
estables45 <-subset(m_final45, suma.vec45 == 0)
estables45 <-estables45[,1:7]
write.csv(estables45,"Resample/CNRMCM5_rcp45_2015_2039_bio/estables45.csv")
estables45 <-estables45[,1:3]
coor.ref<-crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
estables<-rasterFromXYZ(estables45, crs = coor.ref, digits=0)
plot(estables)
writeRaster(estables,filename="estables45.tif", format="GTiff", overwrite = TRUE)

estables85 <-subset(m_final85, suma.vec85 == 0) 
estables85 <-estables85[,1:7]
write.csv(estables45,"Resample/HAD_ESM_LR_rcp85_2015_2039_bio/estables85.csv")
estables85 <-estables45[,1:3]
estables<-rasterFromXYZ(estables85, crs = coor.ref, digits=0)
writeRaster(estables,filename="estables85.tif", format="GTiff", overwrite = TRUE)

#Identificar zonas de clima que cambian (suma.vec45 != 0)
cambio45 <-subset(m_final45, suma.vec45 != 0)
cambio.45select<-c("x","y","MEXBIO_2010_gw_pr","zvh","zvhrcp45",
                  "ecoregiones","Lat.vec45","Ppt.vec45","Phum45.vec","suma.vec45")
cambio45<-cambio45[cambio.45select]
write.csv(cambio45,"Resample/CNRMCM5_rcp45_2015_2039_bio/cambio45.csv")
cambio45 <-cambio45[,1:2]
cambio<-rasterFromXYZ(cambio45, crs = coor.ref, digits=0)
writeRaster(cambio,filename="cambio45.tif", format="GTiff", overwrite = TRUE)

cambio85 <-subset(m_final85, suma.vec85 != 0) 
cambio.85select<-c("x","y","MEXBIO_2010_gw_pr","zvh","zvhrcp85",
                   "ecoregiones","Lat.vec85","Ppt.vec85","Phum85.vec","suma.vec85")
cambio85<-cambio85[cambio.85select]
write.csv(cambio85,"Resample/HAD_ESM_LR_rcp85_2015_2039_bio/cambio85.csv")
cambio85 <-cambio85[,1:2]
cambio<-rasterFromXYZ(cambio85, crs = coor.ref)
writeRaster(cambio,filename="cambio85.tif", format="GTiff", overwrite = TRUE)



