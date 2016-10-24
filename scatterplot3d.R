### The following is code for plots and calculting the amount of niche shift across the different insect orders
### in the study.
### Base code is from Diedrick Strubbe - cheers! (see references in main manuscript for publication details)

## Produces Figure 1 and Figure 2 in main manuscript.

library (scatterplot3d)

### !!! this boyden table is out of date, rebuild!!!!
boyce_val <- read.csv ("~/Documents/niche/boyce.csv")
#at the moment this is just taking in the 75% overlap one...
pca_val <- read.csv ("~/Documents/niche/pca_nobuff/pca_env0.25.csv")
pca_val$unfilling <- pca_val$unfilling * -1

pca_df <- pca_val

#to add colour between groups.
pca_df$type <- "A"

pca_df$type[pca_df$species=="a_taeniatulus"] <- "B"
pca_df$type[pca_df$species=="a_tessellatus"] <- "B"
pca_df$type[pca_df$species=="l_decemlineata"] <- "B"
pca_df$type[pca_df$species=="h_axyridis"] <- "B"
pca_df$type[pca_df$species=="a_glabripennis"] <- "B"

pca_df$type[pca_df$species=="h_destructor"] <- "C"

pca_df$type[pca_df$species=="t_oleracea"] <- "D"
pca_df$type[pca_df$species=="t_paludosa"] <- "D"
pca_df$type[pca_df$species=="a_albopictus"] <- "D"
pca_df$type[pca_df$species=="c_capitata"] <- "D"
pca_df$type[pca_df$species=="b_dorsalis"] <- "D"
pca_df$type[pca_df$species=="z_indianus"] <- "D"

pca_df$type[pca_df$species=="l_humile"] <- "E"
pca_df$type[pca_df$species=="p_chinensis"] <- "E"
pca_df$type[pca_df$species=="p_megacephala"] <- "E"
pca_df$type[pca_df$species=="s_invicta"] <- "E"
pca_df$type[pca_df$species=="t_speciesE"] <- "E"
pca_df$type[pca_df$species=="t_tsushimae"] <- "E"
pca_df$type[pca_df$species=="w_auropunctata"] <- "E"

test <- merge (boyce_val, pca_df, by="species")

test$points <- 0


test$X.y <- NULL

colnames (test) <- c("species", "ID", "boyce", "overlap", "expansion", "stability", "unfilling", "euc_dist", "euc_ext", "type", "points") 


#Add in the colours.

test$col <- "dark grey"

test$col[test$type=="B"] <- "#fb0009"
test$col[test$type=="C"] <- "#21b305"
test$col[test$type=="D"] <- "#e6c800"
test$col[test$type=="E"] <- "#007aa6"




#test$unfilling <- -1*test$unfilling

niche.data.temp1 <- test[,c("species", "boyce", "overlap", "expansion", "ID", "points")]
colnames(niche.data.temp1)[4] <- "change"
niche.data.temp2 <- test[,c("species", "boyce", "overlap", "unfilling", "ID", "points")]
colnames(niche.data.temp2)[4] <- "change"


pdf (file = "~/Documents/niche/figs/Fig1.pdf", paper="special", width=10, height=8)


plot.niche <-scatterplot3d(test$overlap,test$boyce,test$points,
                           pch=20,
                           cex.symbols=1.5,
                           color=test$col,
                           box=FALSE,grid=FALSE, lab = c(10,10,10),zlim=c(-1,1),
                           xlab = "niche overlap D", ylab= "Boyce index", zlab="unfilling           expansion",
                           main="Niche Metrics",
                           angle=45)



forexp <- plot.niche$xyz.convert(niche.data.temp1$overlap,niche.data.temp1$boyce,niche.data.temp1$change)
forunf <- plot.niche$xyz.convert(niche.data.temp2$overlap,niche.data.temp2$boyce,niche.data.temp2$change)
segments(forexp$x,forexp$y,forunf$x,forunf$y, col=test$col, lwd=3) 

plot.niche$plane3d(0,0,0, col="grey",lty="dashed")



text(forexp$x,forexp$y, 
     labels= pca_val$X,
     #outdata$lab,
    adj=c(0, 2),
     cex=.75)

dev.off()

########
#Summary plots
test$unfilling <- -1*test$unfilling

#expansion
exp.A <- test[test$type=="A",]
exp.B <- test[test$type=="B",]
exp.C <- test[test$type=="C",]
exp.D <- test[test$type=="D",]
exp.E <- test[test$type=="E",]

x1 <- ((length(which(test$expansion > 0.1))/length(test[,1]))*100)
x2 <- ((length(which(exp.B$expansion > 0.1))/length(exp.B[,1]))*100)
x3 <- ((length(which(exp.D$expansion > 0.1))/length(exp.D[,1]))*100)
x4 <- ((length(which(exp.E$expansion > 0.1))/length(exp.E[,1]))*100)

exp <- data.frame (rbind (x1, x2, x3, x4))
row.names (exp) <- c("Total", "Coleoptera", "Diptera", "Formicidae")



y1 <- ((length(which(test$unfilling > 0.1))/length(test[,1]))*100)
y2 <- ((length(which(exp.B$unfilling > 0.1))/length(exp.B[,1]))*100)
y3 <- ((length(which(exp.D$unfilling > 0.1))/length(exp.D[,1]))*100)
y4 <- ((length(which(exp.E$unfilling > 0.1))/length(exp.E[,1]))*100)

unf <- data.frame (rbind (y1, y2, y3, y4))
row.names (unf) <- c("Total", "Coleoptera", "Diptera", "Formicidae")


pdf (file = "~/Documents/niche/figs/Fig2.pdf", paper="special", width=10, height=6)

par(mfrow=c(1, 2))
barplot (exp[,1], col=c("light grey", "#fb0009", "#e6c800","#007aa6" ),names.arg=rownames(exp), main="Expansion", ylim=c(0,100))
barplot (unf[,1], col=c("light grey", "#fb0009", "#e6c800","#007aa6"),names.arg=rownames(unf), main="Unfilling", ylim=c(0,100))

dev.off()


#######################################################
## BEYOND THIS POINT IS JIBBERISH 


##640 400


test$col[test$type=="B"] <- "#fb0009"
test$col[test$type=="C"] <- "#21b305"
test$col[test$type=="D"] <- "#e6c800"
test$col[test$type=="E"] <- "#007aa6"


##################################
#Unfilling

test$col <- "black"
test$col[test$unfilling<=-0.1] <- "#66C2A5"

plot.niche <-scatterplot3d(test$overlap,test$boyce,test$points,pch = 20,
                           color=test$col,
                           cex.symbols=1.5,
                           box=FALSE,grid=FALSE, lab = c(10,10,10),zlim=c(-1,1),
                           xlab = "niche overlap D", ylab= "Boyce index", zlab="unfilling           expansion",
                           main="Unfilling")



forexp <- plot.niche$xyz.convert(niche.data.temp1$overlap,niche.data.temp1$boyce,niche.data.temp1$change)
forunf <- plot.niche$xyz.convert(niche.data.temp2$overlap,niche.data.temp2$boyce,niche.data.temp2$change)
segments(forexp$x,forexp$y,forunf$x,forunf$y, col=test$col, lwd=3) 

plot.niche$plane3d(0,0,0, col="grey",lty="dashed")

text(forexp$x,forexp$y, labels=niche.data.temp1$ID,cex=.75,pos=3)



scatterplot3d (test$overlap, test$boyceglob, test$expansion,  
               color=test$col, 
               #highlight.3d=TRUE,
               pch=20,
               box=F,
               type="h",
               
               main="Niche Metrics",
               xlab="Niche Overlap D",
               ylab="Boyce Index",
               zlab="Expansion",
               cex.axis=1, 
               cex.lab=1,
               label.tick.marks=T
)

library (ggplot2)

exp <- pca_val[with(pca_val, order(expansion, species)),]
exp$X <- seq(from=1, to=length(exp[,1]), by=1)

exp.barplot <- qplot(x=X, y=expansion, fill=type,
                       data=exp, geom="bar", stat="identity",
                       position="dodge") + coord_flip()


unf <- pca_val[with(pca_val, order(unfilling, species)),]
unf$X <- seq(from=1, to=length(exp[,1]), by=1)

unf.barplot <- qplot(x=type, y=unfilling, fill=type,
                     data=unf, geom="bar", stat="identity",
                     position="dodge") + coord_flip()

