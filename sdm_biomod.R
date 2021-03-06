## Species distribution models/ ENMs to test the ability of native range data to project across to invasive range
## This will take the native range data and build the models just on that, and project to the world surfaces.
## As we are only interested in making the best model and assessing it post-hoc with the Boyce index, 
## all data is used in training and no AUC or TSS thresholds are employed to determine model fit before projections.
## Backgrounds are the same as used for the PCA-ENV analysis
## 
## Scripts largely follow the BIOMOD2 tutorial with some slight modifications
 
library (biomod2)
library (maptools)

#setup CRU predictors
setwd ("/Volumes/Matt2015/layers/CRU/bioclim")
files <- list.files(pattern='asc', full.names=TRUE )
predictors <- stack(files)

#take out relevant predictors and rename to SWD names
cruclim <- subset(predictors, c("bio_2","bio_3","bio_5","bio_6","bio_7","bio_13","bio_14","bio_15"))
names (cruclim) <- c("bio02","bio03","bio05","bio06","bio07","bio13","bio14","bio15")


###Loop begins here###
setwd ("/Volumes/Matt2015/Postdoc backup/output")
species = list.files('occurences_nobuffer')

for (spp in species[20:22]) { cat(spp,'\n') 

setwd ("/Volumes/Matt2015/Postdoc backup/output")
#read in the native and invasive occurences + background
#native occurences
spp1 <- read.csv (paste0('./occurences_nobuffer/',spp,'/',"nat_df.csv"))
#invasive occurences
#spp2 <- read.csv (paste0('./occurences_nobuffer/',spp,'/',"inv_df.csv"))

#remove the redundant cell column
spp1$cell <- NULL
#spp2$cell <- NULL

##Make pa 0 = NA for pseudoabsences
spp1$pa[spp1$pa==0] <- NA
#spp2$pa[spp2$pa==0] <- NA

#spp2 becomes entire distribution
#spp2 <- rbind(spp1, spp2)

## background of the native range
bg1_vals <- extract (cruclim, spp1[,1:2], df=TRUE)
bg1_vals$ID <- NULL

## background of the native + invasive range
#bg2_vals <- extract (cruclim, spp2[,1:2], df=TRUE)
#bg2_vals$ID <- NULL

#create Spatial Points Data Frames
coordinates (spp1) <- ~X+Y
#coordinates (spp2) <- ~X+Y
spp1_df <- SpatialPointsDataFrame(spp1, data=bg1_vals)
#spp2_df <- SpatialPointsDataFrame(spp2, data=bg2_vals)

#setup BIOMOD data
#following tutorial conventions for simplicity

setwd("/Volumes/Matt2015/niche_nobuff")

###################################################################################################
##GLOBAL MODELLING OPTIONS
###################################################################################################
myBiomodOption <- BIOMOD_ModelingOptions(
  MAXENT = list( path_to_maxent.jar = "/Volumes/Matt2015/maxent/maxent.jar",
                 maximumiterations = 200,
                 visible = FALSE,
                 linear = FALSE,
                 quadratic = FALSE,
                 product = FALSE,
                 threshold = FALSE,
                 hinge = TRUE,
                 lq2lqptthreshold = 80,
                 l2lqthreshold = 10,
                 hingethreshold = 15,
                 #beta_threshold = -1, #numeric (default -1.0), regularization parameter to be applied to all linear, quadratic and product features; negative value enables automatic setting
                 beta_categorical = -1,
                 beta_lqp = -1,
                 beta_hinge = 2,
                 defaultprevalence = 0.5))


###################################################################################################
##NATIVE RANGE MODELLING
###################################################################################################


myRespName <-paste0(spp)
myResp <- spp1
myExpl <- spp1_df
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.name = myRespName)
                              

myBiomodModelOut <- BIOMOD_Modeling(
                                        myBiomodData,
                                        models = c('GLM', 'GBM', 'RF', 'MAXENT'), #four models used
                                        models.options = myBiomodOption,
                                        NbRunEval=10, # how many runs to do = normally 10, but here is 1
                                        DataSplit=70, # data split =100. normally 70 for training, 30 for testing
                                        Yweights=NULL,
                                        Prevalence=0.5,
                                        VarImport=10,
                                        #models.eval.meth = c('TSS','ROC'),
                                        SaveObj = TRUE,
                                        rescal.all.models = TRUE)


#write out evaluation metrics
#eval<- get_evaluations(myBiomodModelOut)
#write.csv (eval, file=paste0(spp, "_nat_eval.csv"))

#projections

myBiomodProj <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                    new.env = cruclim,
                                    proj.name = 'curr_nat',
                                    #xy.new.env = cruclim,
                                    selected.models = 'all',
                                    #binary.meth = c('TSS', 'ROC'),
                                    filtered.meth = NULL,
                                    compress = 'gzip',
                                    clamping.mask = T,
                                    do.stack=T)

myEnsemble <- BIOMOD_EnsembleModeling (modeling.output =myBiomodModelOut, 
                                       chosen.models='all', 
                                       em.by = 'all',
                                       eval.metric=c('TSS','ROC'),
                                       eval.metric.quality.threshold=c(0.5, 0.7), ##this TSS score is low, but sometimes shit is broke
                                       models.eval.meth = c('TSS', 'ROC'), 
                                       prob.mean = FALSE,
                                       prob.ci.alpha = 0.05,
                                       committee.averaging = FALSE,
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')

ensembleBiomodProj <- BIOMOD_EnsembleForecasting(EM.output=myEnsemble,
                                    projection.output=myBiomodProj)
                                   
####play with the outdata now

proj_stack <- get_predictions(ensembleBiomodProj)

#nat_proj_TSS <- proj_stack[[1]]
#nat_proj_ROC <- proj_stack[[2]]

#1 = TSS by mean, 2 = TSS by wmean, 3 = ROC by mean, 4 = ROC by wmean

writeRaster (proj_stack[[1]], file=paste0(spp,"_nat_TSS.asc"), overwrite=TRUE)
writeRaster (proj_stack[[2]], file=paste0(spp,"_nat_ROC.asc"), overwrite=TRUE)

##clear memory & temp files
removeTmpFiles(h = 0.01)
gc()
}

library (raster); library (ecospat);library (maptools) 

# set up a blank data frame to write boyce scores to
rm(outdata)
outdata = data.frame(species=species, boyce=NA) 

setwd ("/Volumes/Matt2015/Postdoc backup/output")


# begin species loop
for (spp in species) { cat(spp,'\n') 
  
  ##read in the shape files for the distributions  
  invpoly <- readShapePoly (paste0("./occurences/", spp, "/invpoly.shp"))     
  natpoly <- readShapePoly (paste0("./occurences/", spp, "/natpoly.shp"))     
  
  #read in the occurences
  occur1 <- read.csv (paste0("./occurences/", spp, "/clim1.csv"))
  occur2 <- read.csv (paste0("./occurences/", spp, "/clim2.csv"))
  occur <- rbind (occur1, occur2)
  
  #single data frame of lats and longs for observations
  occur_pts <- occur2[,1:2]
  
  #read in the raster of the consensus model for the species 
  #(_nat_ROC.asc = model built on native range only and projected globally)                    
  r1 <- raster (paste0("/Volumes/Matt2015/niche_nobuff/",spp,"_nat_ROC.asc"))
  
  #rasterize the polygon of the invasive range
  invR <- rasterize (invpoly, r1)
  
  #natR <- rasterize (natpoly, r1)
  #allR <- merge (natR, invR)
  
  # mask the raster of the consensus model 
  r2 <- mask (r1, invR)
  
  #take the predictions within the mask, create a data.frame and remove the NAs
  r3<- as.data.frame (r2)
  fit <- na.omit (r3)
  
  # extract the modelled values at the presence points
  obs <- extract (r2, occur_pts)
  obs <- na.omit (obs)
  
  
  #
  boyce_val <- ecospat.boyce (fit, obs, nclass=0, 
                              window.w="default",
                              res=100, PEplot=FALSE)
  
  jj <- boyce_val$Spearman.cor
  jj
  
  ii = which(outdata$species==spp)
  
  outdata$species[ii] = spp
  outdata$boyce[ii] = jj
}

write.csv (outdata, file="~/Documents/niche/boyce.csv")
