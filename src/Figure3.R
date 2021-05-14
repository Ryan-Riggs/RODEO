##Follow GRDC gauge validation steps in Qvalidation.R first. 
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
source("E:/research/2019_08_30_rivObs/git/src/Error_stats_functions.R")


##################################################################################################################################################################################################################################
##Figure 3.
##################################################################################################################################################################################################################################

##GRDC
i = 18
print(i)
xSecIDcol = Site_number_xsections[i,]
xSecID=xSecIDcol[2:ncol(xSecIDcol)]
mInd = match(xSecID, data$ID) #changed from tab
xSecw = as.data.frame(seq(1:100))
xSecq = as.data.frame(seq(1:100))
mInd_paired = which(data_val_8415$ID %in% xSecID)
paired_df = data_val_8415[mInd_paired,]
###Pair largest GRADES flowline with same day Landsat widths. 
pdf_comid = grep("COMID", names(paired_df))
pdf_filtering = unique(paired_df[,pdf_comid])
list = na.omit(unlist(pdf_filtering))
testing_df = as.vector(as.numeric(length(list)))
for(l in 1:length(list)){
  paired_df$COMID = list[l]
  Gauge_comid_filt = Gauge_comid[.(paired_df$COMID), allow.cartesian = TRUE]
  joining = left_join(paired_df, Gauge_comid_filt, by = c("COMID" = "COMID", "Date" = "date"))
  testing_df[l] = mean(joining$Q, na.rm = TRUE)
}
place = as.data.frame(cbind(testing_df, list))
placement = place$list[place$testing_df==max(place$testing_df, na.rm = TRUE)]
paired_df$COMID = placement[1]
Gauge_comid_filter1  = Gauge_comid[.(paired_df$COMID), allow.cartesian = TRUE]
paired_df = left_join(paired_df, Gauge_comid_filter1, by = c("COMID" = "COMID", "Date" = "date"))
paired_df = distinct(paired_df, Date, .keep_all = TRUE)
p_q = quantile(paired_df$Q, probs = seq(percentiles[1],percentiles[2],.01), na.rm = TRUE)
p_w = quantile(paired_df$calc_mean, probs = seq(percentiles[1],percentiles[2],.01), na.rm = TRUE)
l_q = length(p_q)
l_w = length(p_w)
gage_quants_q[i,1:l_q] = p_q
gage_quants_w[i,1:l_w] = p_w
###Determine widths post 2014 and aggregate to same day in case of overlapping observations. 
data_val_subset = subset(data_val, xSecID[,1] ==data_val$ID)
data_val_sub_agg = aggregate(data_val_subset$calc_mean, by=list(data_val_subset$Date), FUN=max)#} else{next}
w_rc= p_w
q_rc = p_q
y = q_rc
x = w_rc
paired_df = paired_df[paired_df$Q>=p_q[[1]]&paired_df$Q<=rev(p_q)[[1]],]
paired_df = paired_df[!is.na(paired_df$`system:index`),]
##Create function to estimate Q from a Landsat width based on quantile pair up. 
spl = approxfun(x, y)
###Read in either USGS or GRDC gauges and process them. 
usgs_q = try(usgs_q_processing(rawDailyData <- readNWISdv(Site_number_xsections$`gageinfo$SITE_NUM`[i],"00060","1984-01-01"))) ##Usgs Make a switch for USGS vs Canadian. 
####################################################################################  
usgs_q = usgs_q[usgs_q$q>=0,]
usgs_q = usgs_q[usgs_q$Date<as.Date("2021-01-01"),]
usgs_q$Year = format(usgs_q$Date, "%Y")
AnnualAvg = aggregate(usgs_q$q, by=list(usgs_q$Year), mean, na.rm = TRUE)
gage_stats$Q_50[i]=mean(AnnualAvg$x, na.rm = TRUE)
data_val_sub_agg$Date = data_val_sub_agg$Group.1
joined = inner_join(data_val_sub_agg, usgs_q)
joined$model = spl(joined$x)
plot.new()
par(pty = "s", mar = c(8, 4.1,4.1,2.1))
plot(c(min(x, na.rm = TRUE), max(x, na.rm = TRUE)),c(min(spl(x), na.rm = TRUE), max(spl(x), na.rm = TRUE)), col = "white", 
     log = "", type="p", ylab="",
     xlab="",lwd = 2, bty = "n", yaxt = "n", xaxt = "n")
box(lwd = 2)
axis(side = 1, labels = TRUE, at = c(260,300,340,380))
axis(side = 2, las = 2, at = c(2000,3500,5000,6500))
mtext("GRADES Discharge (cms)",side = 2, cex = 1.25, line = 3.25)
mtext("RODEO Width (m)", side = 1, cex = 1.25, line = 3)
points(paired_df$calc_mean, paired_df$Q, col = "red")
colorPal = rev(colorRampPalette(brewer.pal(n = 8, name = "Spectral"))(95))
segments(x[-length(x)],spl(x)[-length(x)],x[-1L],spl(x)[-1L], col=  colorPal, lwd=4)
med = c(x[41], x[51])
loc = paired_df[paired_df$calc_mean>=med[1] & paired_df$calc_mean<=med[2],]
loc$dist = loc$Q - spl(loc$calc_mean)
loc = loc[loc$dist==max(loc$dist, na.rm = TRUE),]
arrows(loc$calc_mean, loc$Q, loc$calc_mean, spl(loc$calc_mean), code = 3, angle = 90, length = 0.10, lwd = 3)
#points(loc$calc_mean, loc$Q, col = "red", pch = 16)
legend <- bquote("Q"[i] - hat(Q)[i])  
legend(x = (loc$calc_mean*.85),y = (loc$Q* 1.05), legend = legend,
       col = "Red",bty = "n", cex= 1.5)
#segments(x[-length(x)],125,x[-1L],125, col="black",lty = 2, lwd=1, xpd = TRUE)
rng = seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), 2)
#lines(rng, rep(100, length(rng)), lty = 3, xpd = TRUE)
axis(side = 1, xpd = TRUE, line =6.5, hadj = 1, labels = FALSE, tick = TRUE)
segments(x[-length(x)],-800,x[-1L],-800, col=  colorPal, lwd=10, xpd = TRUE)
mtext("Percentile", side = 1, cex =1.25, line = 4.75)
axis(side = 2, pos = 500, xpd = TRUE)
mtext("5", side = 1, line = 7, at= c(min(x, na.rm = TRUE), 100))
mtext("95", side = 1, line = 7, at= c(max(x, na.rm = TRUE), 100))
mid = (min(x, na.rm = TRUE)+max(x, na.rm= TRUE))/2
mtext("50", side = 1, line = 7, at = c(mid, 100))
