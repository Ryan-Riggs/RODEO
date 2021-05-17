##Figure 6. 
library(sf)
library(geosphere)
library(tidyverse)
library(sf)
library(sp)
library(ggplot2)
library(rgeos)
library(dplyr)
library(gmt)
library(vroom)
library(BAMMtools)
library(hydroGOF)
library(data.table)
library(parallel)
library(foreach)
library(tidync)
library(dataRetrieval)
library(ncdf4)
library(vioplot)

if (!"foreign" %in% rownames(installed.packages())){
  install.packages("foreign")}; require(foreign)
if (!"rgdal" %in% rownames(installed.packages())){
  install.packages("rgdal")}; require(rgdal)
if (!"shapefiles" %in% rownames(installed.packages())){
  install.packages("shapefiles")}; require(shapefiles)
if (!"RColorBrewer" %in% rownames(installed.packages())){
  install.packages("RColorBrewer")}; require(RColorBrewer)
if (!"zyp" %in% rownames(installed.packages())){
  install.packages("zyp")}; require(zyp)

png("E:\\research\\RODEO\\Figures\\boxplots.png",
    units = "in",
    width = 5.5,
    height = 5.5,
    res = 500,
    pointsize = 10)

dev.off()














##All results. 
width = 400
n = boxplot(gage_stats$NRMSE)
rB = boxplot(gage_stats$rBias)
rR = boxplot(gage_stats$RRMSE)
R = boxplot(gage_stats$R)
K = boxplot(gage_stats$KGE)
par(mai = c(1,1,1,1))
par(mfrow = c(1,1), pty = "s")
myCol = c("white", "lightskyblue2")
myCol1 = myCol
limits = c(-100, 300)
vColor = "grey90"
vOutline = "grey90"
vWdth = 1
vTy = 1
thickness = 2.65
b_Thickness = 0.6
vioplot(gage_stats$NRMSE[!(gage_stats$NRMSE%in%n$out)],gage_stats$rBias[!(gage_stats$rBias%in%rB$out)], gage_stats$RRMSE[!is.infinite(gage_stats$RRMSE)&!(gage_stats$RRMSE%in%rR$out)],lty = vTy,
        col = vColor,border = vOutline, lwd = vWdth,xlim = c(0,16), at = c(1,4,7), ylim = limits,wex = thickness, drawRect = FALSE, side = "left", yaxt = "n")
par(new = TRUE)
vioplot(gage_stats$NRMSE[gage_stats$GRWL_width_m>=width&!(gage_stats$NRMSE%in%n$out)],gage_stats$rBias[gage_stats$GRWL_width_m>=width&!(gage_stats$rBias%in%rB$out)], gage_stats$RRMSE[!is.infinite(gage_stats$RRMSE)&gage_stats$GRWL_width_m>=width&!(gage_stats$RRMSE%in%rR$out)],
        lty = vTy,col = vColor, lwd = vWdth,border = vOutline,xlim = c(0,16), at = c(2,5,8), ylim = limits,wex = thickness, drawRect = FALSE, side = "right", yaxt = "n")
par(new = TRUE)
vioplot(gage_stats$KGE[!(gage_stats$KGE%in%K$out)],gage_stats$R[!(gage_stats$R%in%R$out)],
        lty =vTy,col = vColor, lwd = vWdth,border = vOutline,xlim = c(0,16), at = c(10,13), ylim =c(-1,1), limits,wex = thickness, drawRect = FALSE, side = "left", yaxt = "n")
par(new = TRUE)
vioplot(gage_stats$KGE[gage_stats$GRWL_width_m>=width&!(gage_stats$KGE%in%K$out)],gage_stats$R[gage_stats$GRWL_width_m>=width&!(gage_stats$KGE%in%K$out)],
        lty = vTy,col = vColor, lwd = vWdth,border = vOutline,xlim = c(0,16), at = c(11,14), ylim =c(-1,1), limits,wex = thickness, drawRect = FALSE, side = "right", yaxt = "n")
par(new = TRUE)
boxplot(outline = FALSE,gage_stats$NRMSE,gage_stats$NRMSE[gage_stats$GRWL_width_m>=width], gage_stats$rBias,
        gage_stats$rBias[gage_stats$GRWL_width_m>=width], gage_stats$RRMSE,gage_stats$RRMSE[gage_stats$GRWL_width_m>=width],
        NA,NA,NA,NA, ylim = limits,xlim = c(0,16),
        ylab = "", las=2,boxwex = b_Thickness,
        at = c(1,2,4,5,7,8,10,11,13,14), col = rep(myCol, 4), medcol = "black", staplelty = 0, lwd =2, whisklty = 1,medlwd =2,bty = "n",xaxt = "n")
box(lwd = 2)
mtext("Percent (%)", side = 2,line = 3, col="black", cex = 1.25)
medList = list(gage_stats$NRMSE, gage_stats$NRMSE[gage_stats$GRWL_width_m>=width], gage_stats$rBias, gage_stats$rBias[gage_stats$GRWL_width_m>=width],
               gage_stats$RRMSE, gage_stats$RRMSE[gage_stats$GRWL_width_m>=width])
meds = lapply(medList, median, na.rm = TRUE)
#text(c(1,2,4,5,7,8), unlist(meds)+10, round(unlist(meds)))
numb = paste0("n=", nrow(gage_stats))
numb1 = paste0("n=", nrow(gage_stats[gage_stats$GRWL_width_m>=width,]))
numbCombined = rep(c(numb, numb1),4)
#text(c(1,2,4,5,7,8), unlist(meds)-5, numbCombined[1:6])
par(new = TRUE)
boxplot(outline = FALSE,NA, NA, NA,NA,NA,NA, gage_stats$KGE,gage_stats$KGE[gage_stats$GRWL_width_m>=width],gage_stats$R,gage_stats$R[gage_stats$GRWL_width_m>=width], ylim = c(-1,1), axes = FALSE, col = rep(myCol, 4),
        boxwex = b_Thickness,
        medcol = "red", staplelty = 0, border = "red", lwd = 2,whisklty = 1,medlwd = 2, at = c(1,2,4,5,7,8,10,11, 13,14), xlim = c(0,16))
medList = list(gage_stats$KGE,gage_stats$KGE[gage_stats$GRWL_width_m>=width],gage_stats$R, gage_stats$R[gage_stats$GRWL_width_m>=width])
meds = lapply(medList, median, na.rm = TRUE)
##text(c(10,11,13,14), unlist(meds)+.05, signif(unlist(meds),1))
#text(10:11, unlist(meds)-.05, numbCombined[7:8])

axis(4, col = "red",col.ticks = "red",col.axis = "red", las = TRUE, lwd = 2)
axis(4, col = "red",col.ticks = "red",col.axis = "red", las = TRUE, lwd = 2, at =seq(-1.5, 1.5))
axis(1, at = c(1.5,4.5,7.5), labels= c("NRMSE","rBias","RRMSE"))
axis(1, at = c(10.5), labels = c("KGE"), col = "red", col.ticks = "red", col.axis = "red")
axis(1, at = c(13.5), labels = c("R"), col = "red", col.ticks = "red", col.axis = "red")
mtext("R / KGE", side = 4,line = 3, col="red", cex = 1.25)
n = nrow(gage_stats)
text = paste("N = ", n)
top = expression(Width >= paste(100, " m"), 
                 atop(paste("N = 456 Gauges")))
bottom =expression(Width >= paste(400, " m"), 
                   atop(paste("N = 53 Gauges")))
legend("topleft", legend = c(top, bottom), xpd = TRUE, fill = c(myCol1[1], NA, myCol1[2], NA),
       horiz = FALSE, bty = "n", inset = c(0,0),
       border =c("black", NA, "black", NA), cex = 0.8, y.intersp = .5)
