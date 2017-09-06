# Script for calculating basic statistics for ecoregions and protected areas inside ecoregions
# Author: D. Ramírez-Mejía
# Date: 09.2017

rm(list=ls(all=TRUE))


#--- Packages ---

library(rgeos)
library(raster)
library(maptools)
library(sp)
library(plyr)
library(tidyr)
library(dplyr)


# --- Read raster files ---

projLambert_WSG84 <- "+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +datum=WGS84 +units=m +no_defs" #Mexico
edos <- readShapePoly("Z:/Coberturas/DEST_2012/DEst_2012.shp", IDvar=NULL, proj4string=CRS(projLambert_WSG84))
malla_raster <-  raster("Z:/CoberturasRestringidas/Malla_Mexico/Malla_pr_1km/malla_1kmPR_ccl1.tif")      
impact <- raster("W:/O_godinez/CORREDORES/RESISTENCIA/RESIST_ELEGIDOS/usv_roads_ca_caEucDist_accumulated_cost_surface_INVERTIDO.tif")


# --- Read ECOREGIONS and ANP shape files ---

ecoreg <- readShapePoly("Z:/Coberturas/Ecorre_terr_2008/ecorregiones_2008_c/ecort08cw.shp", IDvar=NULL, proj4string=CRS(projLambert_WSG84))
anp <- readShapePoly("Z:/Coberturas/ANP_2017_federales/181ANP_WGS1984_c.shp", IDvar=NULL, proj4string=CRS(projLambert_WSG84))

### Extract protected areas by year 
anp_t <- crop(anp, edos)
anp_t@data$YEAR <- gsub("-.*","",anp_t@data$PRIM_DEC)
anp_t@data <- (arrange(anp_t@data, YEAR))
anp_t@data$ID <- seq.int(nrow(anp_t@data))

anp12 <- anp_t[anp_t@data$YEAR <= "2012",]

# --- Dissolve to ecoregion level N1 ---
ecoreg_n1 <- gUnaryUnion(ecoreg, id = ecoreg@data$CVEECON1)
levels(ecoreg@data$DESECON1)

calif_med <- gUnaryUnion(ecoreg[ecoreg@data$DESECON1=="California Mediterranea",])
des_am_norte <- gUnaryUnion(ecoreg[ecoreg@data$DESECON1=="Desiertos de America del Norte",])
elev_semar_meri <- gUnaryUnion(ecoreg[ecoreg@data$DESECON1=="Elevaciones Semiaridas Meridionales",])
grandes_plan <- gUnaryUnion(ecoreg[ecoreg@data$DESECON1=="Grandes Planicies",])
selvas_cal_hum <- gUnaryUnion(ecoreg[ecoreg@data$DESECON1=="Selvas Calido-Humedas",])
selvas_cal_sec <- gUnaryUnion(ecoreg[ecoreg@data$DESECON1=="Selvas Calido-Secas",])
sierras_tem <- gUnaryUnion(ecoreg[ecoreg@data$DESECON1=="Sierras Templadas",])

ecoN1 <- c(calif_med, des_am_norte, elev_semar_meri, grandes_plan, selvas_cal_hum, selvas_cal_sec, sierras_tem)

# --- Erase ANP from Ecoregion shape ---

eco_noANP <- function (eco_stack, anp_file) {
  
  eco_buffer <- list()
  eco_diff <- list()
  anp_buff <- gBuffer(anp_file, byid=TRUE, width=0)
  
  for (i in 1:length(eco_stack)){
    eco_buffer[[i]] <- gBuffer(eco_stack[[i]], byid=TRUE, width=0)
    eco_diff[[i]] <- gDifference(eco_buffer[[i]], anp_buff)
  }
  
  return(eco_diff)
  
}

ecoN1_noANP_2012 <- eco_noANP(ecoN1, anp12)

# --- Extract ANP by Ecoregion ---

eco_ANP <- function (eco_stack, anp_file) {
  
  eco_buffer <- list()
  eco_anp <- list()
  anp_buff <- gBuffer(anp_file, byid=TRUE, width=0)
  
  for (i in 1:length(eco_stack)){
    eco_buffer[[i]] <- gBuffer(eco_stack[[i]], byid=TRUE, width=0)
    eco_anp[[i]] <- gIntersection(eco_buffer[[i]], anp_buff)
  }
  
  return(eco_anp)
  
}

ecoN1_ANP_2012 <- eco_ANP(ecoN1, anp12)


# --- IMPACT: Ecoregion (without Protected Areas) --- 

set_extent_eco <- function(raster, shps_stack){
  
  crop_x <- list()
  set_x <- list()
  for (i in 1:length(shps_stack)) {
    crop_x[[i]] <- crop(raster, extent(shps_stack[[i]]))
    set_x[[i]] <- mask(crop_x[[i]], shps_stack[[i]])
  }
  
  return(set_x)
  
}

impact_eco <- set_extent_eco(impact, ecoN1_noANP_2012)


# --- IMPACT: Protected Areas by Ecoregion --- 

impact_eco_anp <- set_extent_eco(impact, ecoN1_ANP_2012)


# --- STATS ---
# Set ecoregion name

names_eco <- c("calif_med","des_am_norte","elev_semar_meri","grandes_plan","selvas_cal_hum","selvas_cal_sec","sierras_tem")

eco_name <-  function(raster_stack){
  for (i in 1:length(raster_stack)){
    names(raster_stack[[i]]) <- names_eco[i]
  }
  return(raster_stack)
}

impact12_eco <- eco_name(impact_eco)
impact12_eco_anp <- eco_name(impact_eco_anp)


### Function to calculate stats

impact_eco_anp_stats <- function (raster_stack) {
  
  df_anp_mean <- c()
  df_anp_median <- c()
  df_anp_quantile <- c()
  df_anp_sd <- c()
  eco_anp_stats <- c()
  
  for (i in 1:length(raster_stack)){
    
    df_anp_mean[[i]] <- data.frame(cellStats(raster_stack[[i]], 'mean'))
    colnames(df_anp_mean[[i]]) <- "mean"
    df_anp_median[[i]] <-data.frame(median(raster_stack[[i]], na.rm = TRUE))
    colnames(df_anp_median[[i]]) <- "median"
    df_anp_quantile[[i]] <- t(data.frame(quantile(raster_stack[[i]])))
    colnames(df_anp_quantile[[i]]) <- c("Q0","Q25","Q50","Q75","Q100")
    rownames(df_anp_quantile[[i]]) <- i
    df_anp_sd[[i]] <- data.frame(cellStats(raster_stack[[i]], 'sd'))
    colnames(df_anp_sd[[i]]) <- "sd"
    eco_anp_stats[[i]] <- cbind(df_anp_mean[[i]],df_anp_median[[i]],df_anp_quantile[[i]],df_anp_sd[[i]])
  }
  
  return(ldply(eco_anp_stats,data.frame))
  
}

Impact_eco_stats <- impact_eco_anp_stats(impact12_eco)
Impact_eco_stats$Status <- "Ecoregion"
Impact_eco_stats$Ecoregion <- names_eco
Impact_eco_anp_stats <- impact_eco_anp_stats(impact12_eco_anp)
Impact_eco_anp_stats$Status <- "ANP"
Impact_eco_anp_stats$Ecoregion <- names_eco

Mexico_impact_eco_anp_stats <- rbind(Impact_eco_stats,Impact_eco_anp_stats)

write.table(Mexico_impact_eco_anp_stats, 
            "J:/USUARIOS/ANALISIS/SubCoord_Evaluacion_Ecosistemas/Biblioteca_codigosR/Impacto_Ecorregiones_ANP/Mexico_impact_eco_anp_stats.CSV", 
            sep = ",",  col.names = NA, row.names = TRUE)

