######
# PCA niche overlap analysis with species loop.
# Original code is from Broennimann et al. 2012, now in the package 'ecospat'
# the following contains my own modifications to run over a species loop, 
# perform MESS analysis and output a table with all the scores from the niche tests

library(ecospat)
library(biomod2)
library(ade4)
library(adehabitatHS)
library(MASS)
library (maptools)
library (dismo)
library (fields)

#######################################################################################
###function fix for arrows
## I was having some display difficulties, so I quickly edited the function here.
# Draw arrows linking the centroid of the native and exotic (non-native) distribution (black continuous line)
# and between native and invaded extent (red dashed line).

matt.arrows <- function (sp1, sp2, clim1, clim2) 
{
  arrows(median(sp1[, 1]), median(sp1[, 2]), median(sp2[, 1]), 
         median(sp2[, 2]), col = "black", lwd = 2, length = 0.1)
  arrows(median(clim1[, 1]), median(clim1[, 2]), median(clim2[,1]), 
         median(clim2[, 2]), lty = 11, col = "red", lwd = 2, 
         length = 0.1)
}

#######################################################################################
# Load predictors 
# These ascii grids are from the CRU 3.0 database and at 10' resolution.

setwd ("/Volumes/Matt2015/layers/CRU/bioclim")
files <- list.files(pattern='asc', full.names=TRUE )
cruclim <- stack(files)

names (cruclim) <- c("bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", 
                     "bio01", "bio02","bio03","bio04", "bio05","bio06","bio07","bio08","bio09")

cruclim <- subset(cruclim, c("bio02","bio03","bio05","bio06","bio07","bio13","bio14","bio15"))


## this is the threshold at which to develop the niche metrics. This selection allows of running across
## a variety of thresholds. Use pcat = 0.25 (75th percentile) as default when first optimising the analysis.
pcathresh <- c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.25)


for (pcat in pcathresh) { cat(pcat,'\n') 

#  Folder for the species occurences. Each species is in its own folder, 
#  this lists those species folders to create the species list, and then the loop cycles through each species folder
#  with the associated distribution files

setwd ("/Volumes/Matt2015/Postdoc backup/output")
species = list.files('occurences_nobuffer')

## outdata here is the table that will be added to with all the species niche metrics
outdata = data.frame(species=species, overlap=NA, expansion=NA, stability=NA, unfilling=NA, euc_dist=NA, euc_ext=NA) 

ident = data.frame(species=species, a=NA, ap=NA, b=NA, bp=NA, b2=NA, b2p=NA, D=NA)

for (spp in species[1:22]) { cat(spp,'\n') 

# load climate variable for all site of the study area 1 = native area, 2 = invasive area (these data frames have been resampled to grid cell level)
#here the spp1 file is p/a dataset, just select the 0s as used in SDMs
    spp1 <- read.csv (paste0('./occurences_nobuffer/',spp,'/',"nat_df.csv"))
    spp1$cell <- NULL
    spp2 <- read.csv (paste0('./occurences_nobuffer/',spp,'/',"inv_df.csv"))
    spp2$cell <- NULL

#extract climate information across these data
spp1_vals <- extract (cruclim, spp1[,1:2])
spp1 <-cbind (spp1, spp1_vals)
spp1 <- na.omit (spp1)
rm (spp1_vals)

spp2_vals <- extract (cruclim, spp2[,1:2])
spp2 <-cbind (spp2, spp2_vals)
spp2 <- na.omit (spp2)
rm (spp2_vals)

########################### MESS ANALYSIS ###############################
#########################################################################
#load up the polygons
# natpoly <- readShapePoly (paste0("./occurences/", spp,"/natpoly.shp"))
# invpoly <- readShapePoly (paste0("./occurences/", spp, "/invpoly.shp"))
# 
# #convert to raster
# natmask <- rasterize (natpoly, cruclim)
# invmask <- rasterize (invpoly, cruclim)
# 
# e <- c(invpoly@bbox[1,1], invpoly@bbox[1,2], invpoly@bbox[2,1], invpoly@bbox[2,2])
# invmask <- crop (invmask, e)
# inv_cru <- crop (cruclim, invmask)
# inv_cru <- mask (inv_cru, invmask)
# 
# mmm <- mess (inv_cru, spp1[,4:11], full=FALSE)
# 
# ## The value here (0 or -10 or -20) reflects the MESS threshold to remove non-analog climates at
# mmc <- reclassify (mmm, c (-Inf, -10, NA))
# 
# ## I was having issues with this species - it doesn't have enoguh presence data inside analog climate space
# ## so had to be skipped for the MESS tests
# if (spp == "a_tessellatus"){
#   mmc <- mmm
# }
# 
# par (mfrow=c(2, 1))
# plot (mmm)
# plot (mmc)
# 
# spp_mmc <- extract (mmc, spp2[,1:2])
# spp_mmc[spp_mmc==Inf] <- NA
# spp2$mess <- spp_mmc
# spp2 <- spp2[!is.na(spp2$mess),]
# spp2$mess <- NULL
#########################################################################

#setup dataframes for analysis (clim1 = all cells P&A, occ.sp1 = presences across these cells)
clim1 <- spp1
occ.sp1 <- spp1[spp1$pa==1,]

clim2 <- spp2
occ.sp2 <- spp2[spp2$pa==1,]

#global species occurences
occ.sp.aggr <- rbind (spp1[spp1$pa==1,], spp2[spp2$pa==1,])
#global climate backgrounds
clim12<-rbind(clim1,clim2)

###### presence/absence datasets already made, just formatting here.
pa1 <- data.frame(clim1[,"pa"])
names (pa1) <-"pa"
pa2 <- data.frame(clim2[,"pa"])
names (pa2) <-"pa"
pa1 <- cbind (clim1, pa1)
pa2 <- cbind (clim2, pa2)

##get rid of "pa" elsewhere
clim1$pa <- NULL
clim2$pa <- NULL
occ.sp1$pa <- NULL
occ.sp2$pa <- NULL

#################################################################################################
############################## ANALYSIS - selection of parameters ###############################
#################################################################################################

# selection of the type of analysis.
# If PROJ =F, the models are calibrated on both ranges.
# If PROJ =T, the models are calibrated on species 1 range only and projected to range 2. 
# Analyses where both ranges are needed (ex: LDA) are not done
PROJ = F

# selection of variables to include in the analyses
names(clim12)
Xvar<-c(3:10)
nvar<-length(Xvar)

#number of interation for the tests of equivalency and similarity
iterations<-100

#resolution of the gridding of the climate space
R=100


#################################################################################################
################### row weigthing and grouping factors for ade4 functions  ######################
#################################################################################################

# if PROJ = F
row.w.1.occ<-1-(nrow(occ.sp1)/nrow(rbind(occ.sp1,occ.sp2))) # prevalence of occ1
row.w.2.occ<-1-(nrow(occ.sp2)/nrow(rbind(occ.sp1,occ.sp2))) # prevalence of occ2
row.w.occ<-c(rep(0, nrow(clim1)),rep(0, nrow(clim2)),rep(row.w.1.occ, nrow(occ.sp1)),rep(row.w.2.occ, nrow(occ.sp2)))

row.w.1.env<-1-(nrow(clim1)/nrow(clim12))  # prevalence of clim1
row.w.2.env<-1-(nrow(clim2)/nrow(clim12))  # prevalence of clim2
row.w.env<-c(rep(row.w.1.env, nrow(clim1)),rep(row.w.2.env, nrow(clim2)),rep(0, nrow(occ.sp1)),rep(0, nrow(occ.sp2)))

fac<-as.factor(c(rep(1, nrow(clim1)),rep(2, nrow(clim2)),rep(1, nrow(occ.sp1)),rep(2, nrow(occ.sp2))))

# if PROJ = T

#row.w.occ.PROJT<-c(rep(0, nrow(clim1)),rep(0, nrow(clim2)),rep(1, nrow(occ.sp1)),rep(0, nrow(occ.sp2)))
#row.w.env.PROJT<-c(rep(1, nrow(clim1)),rep(0, nrow(clim2)),rep(0, nrow(occ.sp1)),rep(0, nrow(occ.sp2)))

# global dataset for the analysis and rows for each sub dataset
data.env.occ<-rbind(clim1,clim2,occ.sp1,occ.sp2)[Xvar]
row.clim1<-1:nrow(clim1)
row.clim2<-(nrow(clim1)+1):(nrow(clim1)+nrow(clim2))
row.clim12<-1:(nrow(clim1)+nrow(clim2))
row.sp1<-(nrow(clim1)+nrow(clim2)+1):(nrow(clim1)+nrow(clim2)+nrow(occ.sp1))
row.sp2<-(nrow(clim1)+nrow(clim2)+nrow(occ.sp1)+1):(nrow(clim1)+nrow(clim2)+nrow(occ.sp1)+nrow(occ.sp2))
row.sp2 <- na.omit (row.sp2)



#################################################################################################
#################################### PCA-ENV ####################################################
#################################################################################################

# measures niche overlap along the two first axes of a PCA calibrated on all the pixels of the study areas

if(PROJ == F){  #fit of the analyse using occurences from both ranges		
  pca.cal <-dudi.pca(data.env.occ,row.w = row.w.env, center = T, scale = T, scannf = F, nf = 2)
}
if(PROJ == T){	#fit of the analyse using occurences from range 1		
  pca.cal <-dudi.pca(data.env.occ,row.w = row.w.env.PROJT, center = T, scale = T, scannf = F, nf = 2)
}
# predict the scores on the axes
scores.clim12<- pca.cal$li[row.clim12,]
scores.clim1<- pca.cal$li[row.clim1,]
scores.clim2<- pca.cal$li[row.clim2,]
scores.sp1<- pca.cal$li[row.sp1,]
scores.sp2<- pca.cal$li[row.sp2,]
scores.sp2 <- na.omit (scores.sp2)

# calculation of occurence density and test of niche equivalency and similarity 
## order for the data into ecospat.grid.clim.dyn = global climate, native climate, native climate at presence points
#th.sp = 0.1 will take out 10% of records (related to PCA centroids) 
z1<- ecospat.grid.clim.dyn(scores.clim12,scores.clim1,scores.sp1,R, th.sp=0)
z2<- ecospat.grid.clim.dyn(scores.clim12,scores.clim2,scores.sp2,R, th.sp=0)
a<-ecospat.niche.equivalency.test(z1,z2,rep=100)# test of niche equivalency and similarity according to Warren et al. 2008
b<-ecospat.niche.similarity.test(z1,z2,rep=100)
b2<-ecospat.niche.similarity.test(z2,z1,rep=100)

###add to table
ident_i = which(ident$species==spp)

#species name
ident$species[ident_i] = spp
#niche equivalency
ident$a[ident_i] = round(a$obs$D,3)
#niche equivalency p value
ident$ap[ident_i] = round(a$p.D,3)

#niche equivalency and P value (native -> invasive)
ident$b[ident_i] = round(b$obs$D,3)
ident$bp[ident_i] = round(b$p.D,3)

ident$b2[ident_i]= round(b2$obs$D,3)
ident$b2p[ident_i]= round(b2$p.D,3)

ident$D[ident_i] = round(as.numeric(ecospat.niche.overlap(z1,z2,cor=T)[1]),3)

#plot			
postscript(file=paste0(spp,"_nichechanges.eps"),horizontal=FALSE,paper="special",height=8, width=8)
#x11(); 
layout(matrix(c(1,1,2,2,1,1,2,2,3,3,4,5,3,3,6,7), 4, 4, byrow = TRUE))
ecospat.plot.niche(z1,title="PCA-env - native niche",name.axis1="PC1",name.axis2="PC2")
ecospat.plot.niche(z2,title="PCA-env - invasive niche",name.axis1="PC1",name.axis2="PC2")
ecospat.plot.contrib(pca.cal$co,pca.cal$eig)
#plot.new();
#test <- text(0.5,0.5,paste("niche overlap:","\n","D=",round(as.numeric(ecospat.niche.overlap(z1,z2,cor=T)[1]),3)))
ecospat.plot.overlap.test(a,"D","Equivalency")
ecospat.plot.overlap.test(b,"D","Similarity 2->1")
ecospat.plot.overlap.test(b2,"D","Similarity 1->2")
dev.off()

################################################################################################
### Niche metric outputs
################################################################################################
# z1.dyn<-ecospat.grid.clim.dyn (scores.clim12, scores.clim1, scores.sp1, R=100)
# z2.dyn<-ecospat.grid.clim.dyn (scores.clim12, scores.clim2, scores.sp2, R=100)
 a.dyn<-ecospat.niche.equivalency.test(z1=z1 , z2=z2, rep=100)

#turn on plotting device
postscript(file=paste0(spp,"_niche_equiv.eps"),horizontal=FALSE,paper="special",height=8, width=8)
ecospat.plot.niche.dyn(z1, z2, title=paste0(spp," Niche Categories and Species Density"),quant=0.75, 
                       colz1="#66C2A5", colz2="#FC8D62", colinter="#8DA0CB", colZ1="#1B9E77", colZ2="#D95F02")
matt.arrows (scores.sp1, scores.sp2, scores.clim1, scores.clim2)
#turn off plotting device
dev.off()

#turn on plotting device
postscript(file=paste0(spp,"_niche_equiv2.eps"),horizontal=FALSE,paper="special",height=8, width=8)
ecospat.plot.overlap.test(a.dyn,"D","Equivalency")
#turn off plotting device
dev.off()

### Here is where the gradients are set to examine niche metrics across different percentiles

#measure euclidean distance between niche centroud (points of arrows)
#between native and invasive distributions
x1 <- cbind(median(scores.sp1[, 1]), median(scores.sp1[, 2]))
x2 <- cbind(median(scores.sp2[, 1]), median(scores.sp2[, 2]))
euc.dist <- rdist (x1, x2)

#between native and invasive extents
y1 <- cbind(median(scores.clim1[, 1]), median(scores.clim1[, 2]))
y2 <- cbind(median(scores.clim2[, 1]), median(scores.clim2[, 2]))
euc.ext <- rdist (y1, y2)

#Niche scores calculation
niche_scores<-ecospat.niche.dyn.index (z1, z2, intersection=pcat)

kk <- (as.numeric(ecospat.niche.overlap(z1,z2,cor=T)[1]))
ii = which(outdata$species==spp)

jj <- niche_scores$dynamic.index.w[1]
ll <- niche_scores$dynamic.index.w[2]
mm <- niche_scores$dynamic.index.w[3]
outdata$species[ii] = spp
outdata$overlap[ii] = kk
outdata$expansion[ii] = jj
outdata$stability[ii] = ll
outdata$unfilling[ii] = mm
outdata$euc_dist[ii] = euc.dist
outdata$euc_ext[ii] = euc.ext

}
write.csv (outdata, file=paste0("~/Documents/niche/pca_nobuff/pca_env",pcat,".csv"))
write.csv (ident, file=paste0("~/Documents/niche/pca_nobuff/pca_env_ident",pcat,".csv"))

}






       

