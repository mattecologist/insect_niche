#Matt Hill
#Background preparation
#Builds background data based on Clim1 and Clim2 data from file_prep.R

library (raster)
library (maptools)
library (rgeos)
library (dismo)

#load shapefiles of continent/country and koppen zones
worldshp <- readShapePoly ("/media/matt/OS/layers/Natural Earth/10m_admin_countries/countries.shp")
koppen <- readShapePoly ("/media/matt/OS/layers/world_koppen/world_koppen.shp")

#load predictors
setwd ("/media/matt/OS/layers/CRU/bioclim")
files <- list.files(pattern='asc', full.names=TRUE )
predictors <- stack(files)
predictors <- subset(predictors, c("bio_2","bio_3","bio_5","bio_6","bio_7","bio_13","bio_14","bio_15"))

ref_rast <- predictors[[1]]

#setwd for project
wd = "/media/matt/OS/output/"; setwd(wd) 
occur.files = list.files('occurences')

occur = NULL

#switch for using buffers

buffers = TRUE

for (tfile in occur.files) { cat(tfile,'\n'); occur = rbind(occur,read.csv(paste('occurences/',tfile,'/occur.csv',sep=''),as.is=TRUE)[1:4]) } 
colnames (occur) <- c("spp", "X", "Y", "status")
species = unique(occur$spp)

for (spp in species) {  cat(spp,'\n') 
  tt = which(occur$spp == spp); pnts = cbind(occur$X[tt],occur$Y[tt])  
  
  all_clim <- extract (predictors, pnts)
  clim <- data.frame (cbind (occur[tt,], all_clim))
  clim$spp <- NULL
  
  clim1 <- data.frame (clim[clim$status==0,])
  clim2 <- data.frame (clim[clim$status==1,])
  
  #dump this output for other analyses
  clim1$status <- NULL
  clim2$status <- NULL
  write.csv(clim1, file=(paste('occurences/',spp,'/clim1.csv',sep='')), row.names=F)
  write.csv(clim2, file=(paste('occurences/',spp,'/clim2.csv',sep='')), row.names=F)
  colnames (pnts) <- c("X", "Y")
  write.csv(pnts, file=(paste('occurences/',spp,'/pnts.csv',sep='')), row.names=F)
  
  #begin geoprocessing of points
  native <- clim1
  invasive <- clim2
  
  
  #determine unique koppen zones
  coordinates(native) <- ~X+Y
  coordinates(invasive) <- ~X+Y
  natK <- na.exclude (over (native, koppen))
  invK <- na.exclude (over (invasive, koppen))
  natK <- unique(natK$DN)
  invK <- unique(invK$DN)
  natKpoly <- koppen[koppen$DN %in% natK, ]
  invKpoly <- koppen[koppen$DN %in% invK, ]
  
  #determine unique CONTINENTS
  natC <- na.exclude (over (native, worldshp))
  invC <- na.exclude (over (invasive, worldshp))
  natC <- unique(natC$continent)
  invC <- unique(invC$continent)
  natCpoly <- worldshp[worldshp$continent %in% natC, ]
  natCpoly <- gUnaryUnion(natCpoly)
  projection(natCpoly) <- CRS('+proj=longlat')
  invCpoly <- worldshp[worldshp$continent %in% invC, ]
  invCpoly <- gUnaryUnion(invCpoly)
  projection(invCpoly) <- CRS('+proj=longlat')
  
  
  #break down Koppen polygons
  natKpoly <- gUnaryUnion(natKpoly)
  projection(natKpoly) <- CRS('+proj=longlat')
  
  invKpoly <- gUnaryUnion(invKpoly)
  projection(invKpoly) <- CRS('+proj=longlat')
  
  
  natpoly <- gIntersection (natKpoly, natCpoly, byid=TRUE)
  invpoly <- gIntersection (invKpoly, invCpoly, byid=TRUE)
  
  if (buffers == TRUE){
    #create buffer around points [d= radius of circle in metres (default = 50000)]
    #set some parameters; width is width of circle radius in metres; crd_lim is boundary of longitudes
    #boundary is required otherwise circles go off extents of world projection
    #these two parameters were balanced
    width = 1000000
    crd_lim = 169
    
    #select points that are above -crd_lim and belmow crd_lim
    native.1 <- native[native$X>-crd_lim & native$X<crd_lim, ]
    native.2 <- native[!native$X>-crd_lim | !native$X<crd_lim, ]
    
    invasive.1 <- invasive[invasive$X>-crd_lim & invasive$X<crd_lim, ]
    invasive.2 <- invasive[!invasive$X>-crd_lim | !invasive$X<crd_lim, ]
    
    natB <- circles(native.1, d=width, lonlat=TRUE)
    natBpoly <- gUnaryUnion(natB@polygons)
    invB <- circles(invasive.1, d=width, lonlat=TRUE)
    invBpoly <- gUnaryUnion(invB@polygons)
    
    if (length(native.2) >0 ){
      natB.2 <- circles(native.2, d=(width/2)*0.185, lonlat=TRUE) 
      natBpoly.2 <- gUnaryUnion(natB.2@polygons)
      natBpoly <- gUnion (natBpoly, natBpoly.2, byid=TRUE)
    }
    
    if (length(invasive.2) >0 ){
      invB.2 <- circles(invasive.2, d=(width/2)*0.185, lonlat=TRUE) 
      invBpoly.2 <- gUnaryUnion(invB.2@polygons)
      invBpoly <- gUnion (invBpoly, invBpoly.2, byid=TRUE)
    }
  
  
  #create single shape based on geographic intersections
  #use following if using the buffers
    natpoly <- gIntersection (natpoly, natBpoly, byid=TRUE)
    invpoly <- gIntersection (invpoly, invBpoly, byid=TRUE)
  
  }
  
  
  #write these out as polygons for later processing
  natID <- sapply(slot(natpoly, "polygons"), function(x) slot(x, "ID"))
  natdf <- data.frame(rep(0, length(natID)), row.names=natID)
  natSPDF <- SpatialPolygonsDataFrame(natpoly, natdf)
  
  invID <- sapply(slot(invpoly, "polygons"), function(x) slot(x, "ID"))
  invdf <- data.frame(rep(0, length(invID)), row.names=invID)
  invSPDF <- SpatialPolygonsDataFrame(invpoly, invdf)
  
  tmpf1 <- paste0('occurences/',spp,'/natpoly')
  tmpf2 <- paste0('occurences/',spp,'/invpoly')
  
  writePolyShape(natSPDF, tmpf1)
  writePolyShape(invSPDF, tmpf2)
  
  
  #sample 10000 random cells across backgrounds & write to .csv
  natpts <- data.frame (spsample (natpoly, 10000, type="random"))
  invpts <- data.frame (spsample (invpoly, 10000, type="random"))
  write.csv(natpts, file=(paste('occurences/',spp,'/natpts.csv',sep='')), row.names=F)
  write.csv(invpts, file=(paste('occurences/',spp,'/invpts.csv',sep='')), row.names=F)
  
  ##gives 1 for occupied cells and 0 for unoccupied cells across background
  #native background
  clim1$cell <- cellFromXY(ref_rast, clim1[,1:2])
  clim1$pa <- 1
  nat_df <- clim1[,c("X", "Y", "cell", "pa")]
  natpts$cell <- cellFromXY(ref_rast, natpts)
  colnames (natpts) <- c("X", "Y", "cell")
  nat_df <- merge (nat_df, natpts, all=T)
  nat_df$pa[is.na(nat_df$pa)] <- 0
  #reverse order and remove duplicates (0s in pa here)
  nat_df <- nat_df[order(nat_df$cell, nat_df$pa, decreasing=TRUE),]
  nat_df <- nat_df[!duplicated(nat_df$cell),]
  write.csv(nat_df, file=(paste('occurences/',spp,'/nat_df.csv',sep='')), row.names=F)
  
  #invasive background
  clim2$cell <- cellFromXY(ref_rast, clim2[,1:2])
  clim2$pa <- 1
  inv_df <- clim2[,c("X", "Y", "cell", "pa")]
  invpts$cell <- cellFromXY(ref_rast, invpts)
  colnames (invpts) <- c("X", "Y", "cell")
  inv_df <- merge (inv_df, invpts, all=T)
  inv_df$pa[is.na(inv_df$pa)] <- 0
  #reverse order and remove duplicates (0s in pa here)
  inv_df <- inv_df[order(inv_df$cell, inv_df$pa, decreasing=TRUE),]
  inv_df <- inv_df[!duplicated(inv_df$cell),]
  
  write.csv(inv_df, file=(paste('occurences/',spp,'/inv_df.csv',sep='')), row.names=F)
  
}

