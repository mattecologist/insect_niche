library (biomod2)
library (maptools)

#setup CRU predictors
setwd ("/media/matt/OS/layers/CRU/bioclim")
files <- list.files(pattern='asc', full.names=TRUE )
predictors <- stack(files)

#take out relevant predictors and rename to SWD names
cruclim <- subset(predictors, c("bio_2","bio_3","bio_5","bio_6","bio_7","bio_13","bio_14","bio_15"))
names (cruclim) <- c("bio02","bio03","bio05","bio06","bio07","bio13","bio14","bio15")

#cruclim <- subset(cruclim, c("bio_5","bio_6","bio_13","bio_14"))
#names (cruclim) <- c("bio_5","bio_6","bio_13","bio_14")



###Loop begins here###
setwd ("/media/matt/OS/output/")
species = list.files('occurences')



for (spp in species) { cat(spp,'\n') 

setwd ("/media/matt/OS/output/")

#read in the native and invasive occurences + background
#native occurences
spp1 <- read.csv (paste0('./occurences/',spp,'/',"nat_df.csv"))
#invasive occurences
spp2 <- read.csv (paste0('./occurences/',spp,'/',"inv_df.csv"))

#remove the redundant cell column
spp1$cell <- NULL
spp2$cell <- NULL

##Make pa 0 = NA for pseudoabsences
spp1$pa[spp1$pa==0] <- NA
spp2$pa[spp2$pa==0] <- NA

#spp2 becomes entire distribution
spp2 <- rbind(spp1, spp2)

## background of the native range
bg1_vals <- extract (cruclim, spp1[,1:2], df=TRUE)
bg1_vals$ID <- NULL

## background of the native + invasive range
bg2_vals <- extract (cruclim, spp2[,1:2], df=TRUE)
bg2_vals$ID <- NULL

#create Spatial Points Data Frames
coordinates (spp1) <- ~X+Y
coordinates (spp2) <- ~X+Y
spp1_df <- SpatialPointsDataFrame(spp1, data=bg1_vals)
spp2_df <- SpatialPointsDataFrame(spp2, data=bg2_vals)

#setup BIOMOD data
#following tutorial conventions for simplicity

setwd("/media/matt/OS/output/sept")

###################################################################################################
##GLOBAL MODELLING OPTIONS
###################################################################################################
myBiomodOption <- BIOMOD_ModelingOptions(
  MAXENT = list( path_to_maxent.jar = "/media/matt/OS/maxent/maxent.jar",
                 maximumiterations = 200))
                 #visible = FALSE,
                 #linear = FALSE,
                 #quadratic = FALSE,
                 #product = FALSE,
                 #threshold = FALSE,
                 #hinge = TRUE,
                 #lq2lqptthreshold = 80,
                 #l2lqthreshold = 10,
                 #hingethreshold = 15,
                 #beta_threshold = -1, #numeric (default -1.0), regularization parameter to be applied to all linear, quadratic and product features; negative value enables automatic setting
                 #beta_categorical = -1,
                 #beta_lqp = -1,
                 #beta_hinge = 2,
                 #defaultprevalence = 0.5))


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
                                        models = c('GLM', 'GBM', 'RF', 'MAXENT'),
                                        models.options = myBiomodOption,
                                        NbRunEval=1,
                                        DataSplit=100,
                                        Yweights=NULL,
                                        Prevalence=0.5,
                                        VarImport=10,
                                        models.eval.meth = c('TSS','ROC'),
                                        SaveObj = TRUE,
                                        rescal.all.models = TRUE)


#write out evaluation metrics
eval<- get_evaluations(myBiomodModelOut)
write.csv (eval, file=paste0(spp, "_nat_eval.csv"))

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
                                       #eval.metric=c('TSS','ROC'),
                                       #eval.metric.quality.threshold=c(0.5, 0.8), ##this TSS score is low, but sometimes shit is broke
                                       #models.eval.meth = c('TSS', 'ROC'), 
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

###################################################################################################
##NATIVE + INVASIVE RANGE MODELLING
###################################################################################################

myRespName <-paste0(spp)
myResp <- spp2
myExpl <- spp2_df

                myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.name = myRespName)

                myBiomodModelOut <- BIOMOD_Modeling(
                                     myBiomodData,
                                     models = c('GLM', 'GBM', 'RF', 'MAXENT'),
                                     models.options = myBiomodOption,
                                     NbRunEval=1,
                                     DataSplit=70,
                                     Yweights=NULL,
                                     Prevalence=0.5,
                                     VarImport=10,
                                     models.eval.meth = c('TSS','ROC'),
                                     SaveObj = TRUE,
                                     rescal.all.models = TRUE)

#write out evaluation metrics
    eval <- get_evaluations(myBiomodModelOut)
    write.csv (eval, file=paste0(spp, "_natinv_eval.csv"))

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
                                       #eval.metric=c('TSS','ROC'),
                                       #eval.metric.quality.threshold=c(0.5, 0.8), ##this TSS score is low, but sometimes shit is broke
                                       #models.eval.meth = c('TSS', 'ROC'), 
                                       prob.mean = TRUE,
                                       
                                       prob.ci.alpha = 0.05,
                                       
                                       committee.averaging = FALSE,
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')

ensembleBiomodProj <- BIOMOD_EnsembleForecasting(EM.output=myEnsemble,
                                                 projection.output=myBiomodProj)

####play with the outdata now

proj_stack <- get_predictions(ensembleBiomodProj)

writeRaster (proj_stack[[1]], file=paste0(spp,"_natinv_TSS.asc"),overwrite=TRUE)
writeRaster (proj_stack[[2]], file=paste0(spp,"_natinv_ROC.asc"), overwrite=TRUE)

## h is hours, 0.01 removes 
removeTmpFiles(h = 0.01)
gc()
}