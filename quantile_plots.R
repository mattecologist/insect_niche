##############################################################
library (ggplot2)

mytheme<- theme(axis.title.y = element_blank(),
                axis.title.x = element_blank(),
                axis.text.x  = element_text(angle=0, size=15, hjust=1, colour="#000000"),
                axis.text.y  = element_text(angle=0, size=15, colour="#000000"),
                title  = element_text(angle=0, size=20, colour="#000000"),
                axis.line = element_line (colour="black"),
                #panel.grid.minor=element_blank(), 
                #panel.grid.major=element_blank(),
                axis.text = element_text (colour="black", size=20))

#with MESS
#setwd("/media/matt/OS/output/out")
#withoutMESS
setwd("~/Documents/niche/pca_nobuff")

pca100 <- read.csv ("pca_env0.csv")
pca95 <- read.csv ("pca_env0.05.csv")
pca90 <- read.csv ("pca_env0.1.csv")
pca85 <- read.csv ("pca_env0.15.csv")
pca80 <- read.csv ("pca_env0.2.csv")
pca75 <- read.csv ("pca_env0.25.csv")
pca70 <- read.csv ("pca_env0.3.csv")

pca100$Quantile <- 100
pca95$Quantile <- 95
pca90$Quantile <- 90
pca85$Quantile <- 85
pca80$Quantile <- 80
pca75$Quantile <- 75
pca70$Quantile <- 70

pca_df <- rbind (pca100, pca95, pca90, pca85, pca80,pca75, pca70)

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

q <- ggplot (pca_df, aes(x=Quantile, y=expansion), group=species) +
  geom_point(aes(size=2, colour=type)) +
  mytheme +
  geom_line(aes(group=species, colour=type)) +
  ggtitle ("Expansion")+
  scale_color_manual(values=c("dark grey", "#fb0009","#21b305", "#e6c800",  "#007aa6")) +
  geom_abline (intercept=0.10, slope=0, alpha=0.9, colour="black", size=2, linetype=3)

p <- ggplot (pca_df, aes(x=Quantile, y=unfilling), group=species) +
  geom_point(aes(size=2, colour=type)) +
  mytheme +
  geom_line(aes(group=species, colour=type)) +
  ggtitle ("Unfilling")+
  scale_color_manual(values=c("dark grey", "#fb0009","#21b305", "#e6c800",  "#007aa6")) +
  geom_abline (intercept=0.10, slope=0, alpha=0.9, colour="black", size=2, linetype=3)


library (grid)
pdf (file = "~/Documents/niche/figs/Fig3.pdf", paper="special", width=16, height=9)
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2, widths=c(1,1))))



print(q, vp = vplayout(1, 1))
print(p, vp = vplayout(1, 2))
#grid.text("(a)", vp = vplayout(1,1), y=unit(0.97, "npc"), x=unit(0, "npc"), gp=gpar(cex=3))
#grid.text("(b)", vp = vplayout(1,2), y=unit(0.97, "npc"), x=unit(0, "npc"), gp=gpar(cex=3))


dev.off()