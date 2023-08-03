#==============================================================================
# Supplementary vegetation analysis for bachelor's thesis
# Institution: University of Leipzig
# Title: Variability in plant species composition in the floodplain meadows 'Papitzer Lachen'
# Author: Lukas Erzfeld
# Date: 26.01.2023
#===============================================================================

# required packages may have to be installed primarily
# install.packages(c("ggplot2", "ggpubr", "gridExtra", "readxl", "tidyverse", "vegan", "VennDiagram))
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(readxl)
library(tidyverse)
library(vegan)
library(VennDiagram)
#-------------------------------------------------------------------------------
# set current directory to working directory
path = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)
#-------------------------------------------------------------------------------

# Read vegetation data
veg_data <- read_excel("data/PL_vegetation_Braun_Blanquet_dots.xlsx")
veg_data[is.na(veg_data)] <- 0
dat <- as.data.frame(veg_data)
rownames(dat) <- as.character(dat[,1])
dat <- dat[,-1]

# Read elevation data
elevation <- read_csv("data/PL_veg_plots_elevation.csv")
elevation <- as.data.frame(elevation)
elevation <- elevation[order(elevation$PlotNr),]
rownames(elevation) <- as.character(elevation[,1])
elevation <- elevation[,-1]

# Read least cost path data
least_cost <- read_excel("data/PL_least_cost_paths.xlsx")

# Read groundwater-surface distance data
gw <- read_csv("data/PL_groundwater_surface_dist_mean.csv")
gw <- as.data.frame(gw)
gw <- gw[order(gw$PlotNr), ]
rownames(gw) <- as.character(gw[, 2])
gw <- gw[, -2]

# Read alluvial clay thickness data
clay <- read_csv("data/PL_alluvial_clay_thickness_mean.csv")
clay <- as.data.frame(clay)
clay <- clay[order(clay$PlotNr), ]
rownames(clay) <- as.character(clay[, 2])
clay <- clay[, -2]

# Ellenberg indicator values
ellenberg <- read.table("data/PL_veg_ellenberg.txt", header = TRUE, sep = "\t", na.string = "NA")
rownames(ellenberg) <- ellenberg[, 1]
ellenberg <- ellenberg[, -1]

ellenberg_mean <- function(deckung, eb){
  eb_non_na <- is.na(eb) == F
  deckung_tot <- sum(deckung[eb_non_na])
  sum(eb * deckung, na.rm = T) / deckung_tot
}

eb_plot <- matrix(0, 20, 5)
rownames(eb_plot) <- rownames(dat)
colnames(eb_plot) <- colnames(ellenberg)
for(i in 1:20){
  for(j in 1:5){
    eb_plot[i, j] <- ellenberg_mean(dat[i, ], ellenberg[, j])
  }
}

#-------------------------------------------------------------------------------

# Calculate species intersections and draw Venn Diagram
venn_count <- as.data.frame(t(dat))

# Scutellario-Veronicetum: 153001, 153002, 153003
# Arrhenatheretum elatioris: 143001, 163001, 163002, 173001, 183001, 193001, 
#                            203001, 213001, 223001, 233001, 243001
# Silaum silaus: 143002, 143003, 173002, 173003, 253001, 263001

# sum up abundances for each association
venn_count$scut_vero <- rowSums(venn_count[, c("153001", "153002", "153003")])
venn_count$arrh_elat <- rowSums(venn_count[, c("143001", "163001", "163002", "173001", 
                                               "183001", "193001", "203001", "213001", 
                                               "223001", "233001", "243001")])
venn_count$sila_sila <- rowSums(venn_count[, c("143002", "143003", "173002", "173003", 
                                               "253001", "263001")])
venn_count <- venn_count[, -c(1:20)]
# presence = 1, absence = NA
venn_count[venn_count > 0] <- 1
venn_count[venn_count == 0] <- NA

# count number of species per association
cat("Number of different species per classified association:\n")
colSums(venn_count, na.rm = T)
nspecies <- colSums(venn_count, na.rm = T)
sv <- nspecies[[1]]
ae <- nspecies[[2]]
ss <- nspecies[[3]]

# count number of species only occurring in one association
sv_excl <- venn_count[is.na(venn_count$arrh_elat) & is.na(venn_count$sila_sila), ]
cat("Species only occurring in Scutellario-Veronicetum: \n", rownames(sv_excl))
ae_excl <- venn_count[is.na(venn_count$scut_vero) & is.na(venn_count$sila_sila), ]
cat("Species only occurring in Arrhenatheretum elatius: \n", rownames(ae_excl))
ss_excl <- venn_count[is.na(venn_count$arrh_elat) & is.na(venn_count$scut_vero), ]
cat("Species only occurring in Silaum silaus: \n", rownames(ss_excl))

# col calculating number of intersecting species for all three associations
venn_count$all_intersect = rowSums(venn_count[, ])
venn_count[venn_count > 0] <- 1

# intersection Scutellario-Veronicetum, Arrhenatheretum elatioris and Silaum silaus
sv_ae_ss <- sum(venn_count$all_intersect, na.rm = T)
# intersection Scutellario-Veronicetum and Arrhenatheretum elatioris
sv_ae <- sum(venn_count$scut_vero == venn_count$arrh_elat, na.rm = T)
sv_ae <- sv_ae - sv_ae_ss
# intersection Scutellario-Veronicetum and Silaum silaus
sv_ss <- sum(venn_count$scut_vero == venn_count$sila_sila, na.rm = T)
sv_ss <- sv_ss - sv_ae_ss
# intersection Arrhenatheretum elatioris and Silaum silaus
ae_ss <- sum(venn_count$arrh_elat == venn_count$sila_sila, na.rm = T)
ae_ss <- ae_ss - sv_ae_ss

# draw Venn Diagram
grid.newpage()
png(filename="data/plots/Papitzer_Lachen_venn_species.png",
    width=3.25, height=3.25,
    units="in",
    res=600,
    pointsize=4)
venn.plot <- draw.triple.venn(area.vector = c((sv - sv_ae - sv_ss - sv_ae_ss), 
                                              sv_ss, 
                                              (ss - sv_ss - ae_ss - sv_ae_ss), 
                                              sv_ae, 
                                              sv_ae_ss, 
                                              ae_ss, 
                                              (ae - sv_ae - ae_ss - sv_ae_ss)), 
                              category = c("Scutellario-\nVeronicetum", "Silaum\nsilaus", 
                                           "Arrhenatheretum elatioris"), 
                              fill = c("orange", "#00a52a", "#9a12b2"),
                              lty = "blank", cex = 2, cat.cex = 2,
                              cat.col = c("orange", "#00a52a", "#9a12b2"), 
                              direct.area = T)
dev.off()
#-------------------------------------------------------------------------------

# Detrended Correspondence Analysis
dca <- decorana(dat)
dca
dca_sc <- scores(dca)

# Plot DCA and save file
grid.newpage()
png(filename="data/plots/Papitzer_Lachen_DCA.png",
    width=3.25, height=3.25,
    units="in",
    res=600,
    pointsize=4)

plot(dca_sc, xlab="DCA axis 1", ylab="DCA axis 2")
# TEXT FORMATTING
dca_sc_text <- cbind(dca_sc)
# all labels higher
dca_sc_text[, 2] <- dca_sc_text[, 2] + 0.05
# separate overlapping labels
# 163002
dca_sc_text[8, 1] <- dca_sc_text[8, 1] + 0.15
dca_sc_text[8, 2] <- dca_sc_text[8, 2] + -0.01
# 173001
dca_sc_text[11, 1] <- dca_sc_text[11, 1] + 0.2
dca_sc_text[11, 2] <- dca_sc_text[11, 2] + 0.05
# 243001
dca_sc_text[18, 1] <- dca_sc_text[18, 1] - 0.15
text(dca_sc_text[,1], dca_sc_text[,2] - 0.127, cex = 0.85, rownames(dat))

dev.off()

# Plot abundance of plant species (circle size) of each plot in ordination space
# uncomment lines 185-189 to plot
# for(i in 1:116){
#   # ggf. sqrt ändern
#   plot(dca_sc, cex = sqrt(dat[, i]) + 0.1, main = colnames(dat)[i])
#   Sys.sleep(2)
# }

# Plot species and plots in ordination space
png(filename="data/plots/Papitzer_Lachen_DCA_spec_eb.png",
    width=3.25, height=3.25,
    units="in",
    res=600,
    pointsize=4)
plot(dca, type = "n", xlim = (c(-2, 5)), ylim = c(-2, 2.5),
     xlab = "DCA axis 1", ylab = "DCA axis 2")
points(dca, display = "sites", cex = 0.8, pch = 15, col = "black", bg = "white")
text(dca, display = "spec", cex = 0.7, col = "black", select = c(2, 7, 9, 16, 20, 
                                                                24, 31, 32, 33,
                                                                35, 36, 39, 42, 
                                                                46, 50, 63, 65,
                                                                71, 75, 78, 83, 
                                                                89, 92, 110))
plot(envfit(dca_sc, eb_plot), col = "dodgerblue3")
dev.off()
#-------------------------------------------------------------------------------

# Create ordination plot with phytosociological associations and env-vars
dca_sc_gg = as.data.frame(dca_sc)

assoc <- c("Arrhenatheretum elatioris", "Silaum silaus", "Silaum silaus", "Scutellario-Veronicetum", 
           "Scutellario-Veronicetum", "Scutellario-Veronicetum", "Arrhenatheretum elatioris", "Arrhenatheretum elatioris", 
           "Arrhenatheretum elatioris", "Silaum silaus", "Silaum silaus", "Arrhenatheretum elatioris", 
           "Arrhenatheretum elatioris", "Arrhenatheretum elatioris", "Arrhenatheretum elatioris", "Arrhenatheretum elatioris", 
           "Arrhenatheretum elatioris", "Arrhenatheretum elatioris", "Silaum silaus", "Silaum silaus")
dca_sc_gg$assoc = assoc

# create envfit arrows
en_elev = envfit(dca, elevation$mean, permutations = 999, na.rm = TRUE)
en_elev_coord <- as.data.frame(scores(en_elev, "vectors")) * ordiArrowMul(en_elev)
en_gw = envfit(dca, gw$mean, permutations = 999, na.rm = TRUE)
en_gw_coord <- as.data.frame(scores(en_gw, "vectors")) * ordiArrowMul(en_gw)
en_lc = envfit(dca, least_cost$cost, permutations = 999, na.rm = TRUE)
en_lc_coord <- as.data.frame(scores(en_lc, "vectors")) * ordiArrowMul(en_lc)
en_clay = envfit(dca, clay$mean, permutations = 999, na.rm = TRUE)
en_clay_coord <- as.data.frame(scores(en_clay, "vectors")) * ordiArrowMul(en_clay)

gg = ggplot(data = dca_sc_gg, aes(x = DCA1, y = DCA2)) +
  scale_y_continuous(breaks=seq(-2, 2.0, 0.5)) + 
  scale_x_continuous(limits=c(-2.5, 4), breaks=seq(-2, 5.0, 1)) + 
  geom_point(data = dca_sc_gg, aes(shape = assoc), size = 3, alpha = 1) +
  scale_shape_manual(values = c(1, 2, 3)) +
  geom_segment(aes(x = 0, y = 0, xend = DCA1, yend = DCA2), 
               data = en_elev_coord, size =1, alpha = 0.5, colour = "blue") +
  geom_text(data = en_elev_coord, aes(x = DCA1, y = DCA2), colour = "blue", 
            fontface = "bold", label = "elevation") + 
  geom_segment(aes(x = 0, y = 0, xend = DCA1, yend = DCA2), 
               data = en_gw_coord, size =1, alpha = 0.5, colour = "blue") +
  geom_text(data = en_gw_coord, aes(x = DCA1 + 0.3, y = DCA2), colour = "blue", 
            fontface = "bold", label = "groundwater-surface dist") + 
  geom_segment(aes(x = 0, y = 0, xend = DCA1, yend = DCA2), 
               data = en_lc_coord, size =1, alpha = 0.5, colour = "blue") +
  geom_text(data = en_lc_coord, aes(x = DCA1, y = DCA2), colour = "blue", 
            fontface = "bold", label = "least cost") + 
  geom_segment(aes(x = 0, y = 0, xend = DCA1, yend = DCA2), 
               data = en_clay_coord, size =1, alpha = 0.5, colour = "blue") +
  geom_text(data = en_clay_coord, aes(x = DCA1, y = DCA2 - 0.04), colour = "blue", 
            fontface = "bold", label = "clay thickness") + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"),
        axis.text = element_text(color="black", size=10),
        legend.key = element_blank(),
        legend.title = element_text(size = 10, face = "bold", colour = "black"),
        legend.text = element_text(size = 9, colour = "black")) +
  labs(shape = "classified plant community")
gg

ggsave(filename = "Papitzer_Lachen_dca_env_var.png",
       plot = last_plot(),
       device = "png",
       path = "data/plots")
#-------------------------------------------------------------------------------

# correlations between DCA1 axis and env-variables
dca_sc_scatter <- as.data.frame(dca_sc)
dca_sc_scatter$gw_mean <- gw$mean
dca_sc_scatter$elev <- elevation$mean
dca_sc_scatter$lc <- least_cost$cost
dca_sc_scatter$clay_mean <- clay$mean

### DCA axis 1
# gw_mean
png(filename="data/plots/DCA1_vs_gw_mean.png",
    width=600, height=350)
scatter <- ggscatter(dca_sc_scatter, x = "DCA1", y = "gw_mean", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "DCA axis 1", ylab = "distance to groundwater table [m]",
          add.params = list(color = "black", fill = "lightgray"),
          cor.coef.size = 5)
scatter + font("xylab", size = 12)
dev.off()

# clay_mean
png(filename="data/plots/DCA1_vs_clay_mean.png",
    width=600, height=350)
scatter <- ggscatter(dca_sc_scatter, x = "DCA1", y = "clay_mean", 
                     add = "reg.line", conf.int = TRUE, 
                     cor.coef = TRUE, cor.method = "pearson",
                     xlab = "DCA axis 1", ylab = "thickness of alluvial clay layer [m]",
                     add.params = list(color = "black", fill = "lightgray"),
                     cor.coef.size = 5, cor.coef.coord = c(-1.15, 2.4))
scatter + font("xylab", size = 12)
dev.off()

# elevation a.s.l.
png(filename="data/plots/DCA1_vs_elev.png",
    width=600, height=350)
scatter <- ggscatter(dca_sc_scatter, x = "DCA1", y = "elev", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "DCA axis 1", ylab = "ground elevation a.s.l. [m]",
          add.params = list(color = "black", fill = "lightgray"),
          cor.coef.size = 5)
scatter + font("xylab", size = 12)
dev.off()

# least cost path
png(filename="data/plots/DCA1_vs_leastcost.png",
    width=600, height=350)
scatter <- ggscatter(dca_sc_scatter, x = "DCA1", y = "lc", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "DCA axis 1", ylab = "least cost to next water body",
          add.params = list(color = "black", fill = "lightgray"),
          cor.coef.size = 5, cor.coef.coord = c(-1.25, 20000))
scatter + font("xylab", size = 12)
dev.off()

### DCA axis 2
# gw_mean
scatter <- ggscatter(dca_sc_scatter, x = "DCA2", y = "gw_mean", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "DCA axis 2", ylab = "distance to groundwater table [m]",
          add.params = list(color = "black", fill = "lightgray"),
          cor.coef.size = 5)
scatter + font("xylab", size = 12)

# clay_mean
scatter <- ggscatter(dca_sc_scatter, x = "DCA2", y = "clay_mean", 
                     add = "reg.line", conf.int = TRUE, 
                     cor.coef = TRUE, cor.method = "pearson",
                     xlab = "DCA axis 2", ylab = "thickness of alluvial clay layer [m]",
                     add.params = list(color = "black", fill = "lightgray"),
                     cor.coef.size = 5, cor.coef.coord = c(-1.15, 2.5))
scatter + font("xylab", size = 12)

# elevation a.s.l.
scatter <- ggscatter(dca_sc_scatter, x = "DCA2", y = "elev", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "DCA axis 2", ylab = "elevation a.s.l. [m]",
          add.params = list(color = "black", fill = "lightgray"),
          cor.coef.size = 5)
scatter + font("xylab", size = 12)

# least cost path
scatter <- ggscatter(dca_sc_scatter, x = "DCA2", y = "lc", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "DCA axis 2", ylab = "least cost",
          add.params = list(color = "black", fill = "lightgray"),
          cor.coef.size = 5)
scatter + font("xylab", size = 12)


# correlations between env-variables
# elevation vs. groundwater surface distance
scatter <- ggscatter(dca_sc_scatter, x = "gw_mean", y = "elev", size = 1,
                     add = "reg.line", conf.int = TRUE, 
                     cor.coef = TRUE, cor.method = "pearson",
                     xlab = "distance to groundwater table [m]", ylab = "ground elevation [m]",
                     add.params = list(color = "black", fill = "lightgray"),
                     cor.coef.size = 2, cor.coef.coord = c(0.0, 98.0))
plot1 <- scatter + theme(text = element_text(size = 6)) 

# elevation vs. clay
scatter <- ggscatter(dca_sc_scatter, x = "clay_mean", y = "elev", size = 1,
                     add = "reg.line", conf.int = TRUE, 
                     cor.coef = TRUE, cor.method = "pearson",
                     xlab = "thickness of alluvial clay layer [m]", ylab = "ground elevation [m]",
                     add.params = list(color = "black", fill = "lightgray"),
                     cor.coef.size = 2, cor.coef.coord = c(1.0, 98.0))
plot2 <- scatter + theme(text = element_text(size = 6)) 

# elevation vs. least cost path
scatter <- ggscatter(dca_sc_scatter, x = "lc", y = "elev", size = 1,
                     add = "reg.line", conf.int = TRUE, 
                     cor.coef = TRUE, cor.method = "pearson",
                     xlab = "least cost to next water body", ylab = "ground elevation [m]",
                     add.params = list(color = "black", fill = "lightgray"),
                     cor.coef.size = 2, cor.coef.coord = c(8000, 97.8))
plot3 <- scatter + theme(text = element_text(size = 6)) + ylim(96.5, 98.0)

grid.arrange(plot1, plot2, plot3, ncol=2, nrow=2)

# save file
g <- arrangeGrob(plot1, plot2, plot3, nrow=2, ncol=2)
ggsave(g,
       file = "elev_vs_env_var.png",
       path = "data/plots")


# clay vs. ground water surface distance
scatter <- ggscatter(dca_sc_scatter, x = "gw_mean", y = "clay_mean", 
                     add = "reg.line", conf.int = TRUE, 
                     cor.coef = TRUE, cor.method = "pearson",
                     xlab = "distance to groundwater table [m]", ylab = "thickness of alluvial clay layer [m]",
                     add.params = list(color = "black", fill = "lightgray"),
                     cor.coef.size = 5)
scatter + font("xylab", size = 12)

# clay vs. least cost path
scatter <- ggscatter(dca_sc_scatter, x = "lc", y = "clay_mean", 
                     add = "reg.line", conf.int = TRUE, 
                     cor.coef = TRUE, cor.method = "pearson",
                     xlab = "least cost to next water body", ylab = "thickness of alluvial clay layer [m]",
                     add.params = list(color = "black", fill = "lightgray"),
                     cor.coef.size = 5)
scatter + font("xylab", size = 12)

# groundwater surface distance vs. least cost path
scatter <- ggscatter(dca_sc_scatter, x = "lc", y = "gw_mean", 
                     add = "reg.line", conf.int = TRUE, 
                     cor.coef = TRUE, cor.method = "pearson",
                     xlab = "least cost to next water body", ylab = "distance to groundwater table [m]",
                     add.params = list(color = "black", fill = "lightgray"),
                     cor.coef.size = 5)
scatter + font("xylab", size = 12)
