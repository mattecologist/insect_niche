library (maptools)
library (raster)

## reference raster which can be used to build up rasters to the same extent as the project.
refrast <- raster ("/Volumes/Matt2015/layers/CRU/bioclim/bio_1.asc")

setwd ("/Volumes/Matt2015/Postdoc backup/output")
species = list.files('occurences_nobuffer')

out_size = data.frame(species=species, cells=NA, enm=NA, poleward=NA, equator=NA) 

for (spp in species) { cat(spp,'\n') 

nat <- readShapePoly (paste0('./occurences_nobuffer/',spp,'/',"natpoly.shp"))
natpts <- read.csv (paste0('./occurences_nobuffer/',spp,'/',"clim1.csv"))
natpts$status <- 1

invpts <- read.csv (paste0('./occurences_nobuffer/',spp,'/',"clim2.csv"))

##simple calculation of cells with a presence added up in native range.
natR <- rasterize (natpts[,1:2], refrast, field=natpts$status)
jj <- cellStats (natR, 'sum', na.rm=TRUE)

## calculate area based on ENM prediction
## need to add in thresholding here too. 

r1 <- raster (paste0("/Volumes/Matt2015/niche_nobuff/",spp,"_nat_ROC.asc"))
r2 <- crop (r1, nat)
kk <- cellStats (r2, 'sum', na.rm=TRUE)


############# latitude shifts
lat_nat.1 <- natpts$Y[natpts$Y<=0]
lat_nat.1 <- lat_nat.1*-1
lat_nat.2 <- natpts$Y[natpts$Y>0]
lat_nat <- (c(lat_nat.1, lat_nat.2))

lat_inv.1 <- invpts$Y[invpts$Y<=0]
lat_inv.1 <- lat_inv.1*-1
lat_inv.2 <- invpts$Y[invpts$Y>0]
lat_inv <- (c(lat_inv.1, lat_inv.2))


ii = which(out_size$species==spp)

out_size$species[ii] = spp
out_size$cells[ii] = jj
out_size$enm[ii] = kk
out_size$poleward[ii] = max(lat_inv)-max(lat_nat)
out_size$equator[ii] = min(lat_inv)-min(lat_nat)
}

out_size$poleward[out_size$poleward<=0] <- 0
out_size$equator[out_size$equator<=0] <- 0

#out_size$enm <- log10(out_size$enm)

plot (enm~cells, data=out_size)
text(enm~cells, data=out_size, labels=species)

boyce_val <- read.csv ("~/Documents/niche/boyce.csv")
pca_val <- read.csv ("~/Documents/niche/pca_nobuff/pca_env0.25.csv")

### need to find this file and rebuild...
#nSDM <- read.csv("/media/matt/OS/output/out/nichshiftSDM.csv")

#nSDM<- nSDM[c("species", "nat_ROC", "nat_TSS")]

boyce_val <- boyce_val[,c("species", "boyce")]
pca_val <- pca_val[,c("species", "overlap", "expansion", "stability", "unfilling")]

#test <- merge (boyce_val, nSDM, by="species")
test <- merge (boyce_val, pca_val, by="species")

dataset <- merge (test, out_size, by="species")

## need to find this as well.. but perhaps can be excluded....
# enfa_nat <- read.csv ("/media/matt/OS/output/out/ENFA.native.csv")
# enfa_inv <- read.csv ("/media/matt/OS/output/out/ENFA.invasive.csv")
# 
# dataset$nat_mar <- enfa_nat$marginality/1.96
# dataset$inv_mar <- enfa_inv$marginality/1.96
# dataset$nat_spec <- enfa_nat$specialization
# dataset$inv_spec <- enfa_inv$specialization
# 
# dataset$delta_mar <- dataset$inv_mar-dataset$nat_mar  
# dataset$delta_spec <- dataset$inv_spec-dataset$nat_spec 

dataset <- cbind(dataset[1], round (dataset[,2:length(dataset),],3))

#boxplot native Y ~ boxplot invasive Y
#qnatify expansion towards equator and towards poles.

dataset$exp_bin <- 1
dataset$exp_bin [dataset$expansion <=0.1] <- 0


dataset$unf_bin <- 1
dataset$unf_bin [dataset$unfilling <=0.1] <- 0

write.csv (dataset, file="~/Documents/niche/summary_stats.csv")

dataset <- read.csv ("~/Documents/niche/summary_stats.csv")

mod1 = glm(exp_bin ~ cells + enm + poleward + equator + HII_mean_native +
          HII_mean_invasive + POP_mean_native + POP_mean_invasive + PORT_mean_native +
            PORT_mean_invasive + ROAD_mean_native + ROAD_mean_invasive + HANPP_mean_native +HANPP_mean_invasive, 
          family = binomial(link = "logit"), data = dataset)
summary(mod1)

mod2 <- nlme (expansion ~ unfilling (cells, poleward, equator, delta_mar, delta_spec), data=dataset,
              fixed= cells + poleward + equator + delta_mar + delta_spec)

mod3 <- gls (unfilling ~ cells + poleward + equator, data=dataset)
summary (mod3)

mod4 <- lme (boyce ~ Family, data=dataset, random=~1)

model1=gam(dataset$expansion~s(dataset$HII_mean_native, 4))
summary(model1) # Dev= 57% P<0.001 edf= 2.64 (BUT dev= 9.3% without outlier)
plot(model1,residuals=T,pch=20,shade=T,xlab="HII native range",ylab="Expansion",col=Order)
legend(17,0.8,c("Acari","Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera"),col=c(1:6),pch=20)


library(MuMIn)
#data(Cement)
fm1 <- glm(dataset$boyce ~ ., data = dataset[,7:10])
dd <- dredge(fm1)



