
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

##Error functions. 
source("E:/research/2019_08_30_rivObs/git/src/Error_stats_functions.R")

Eff_widths = map_df(list.files("E:\\research\\RODEO\\VariableBuffers\\6xLength3xWidth\\USGS\\VariableBuffers_3xLength_1.5Width_USGS", full.names = TRUE), ~vroom(.x))
gageinfo = read.csv("E:\\research\\GRWL\\GRWL_2015_present\\width_val\\input\\gaugeData\\USGS\\gaugeTable.csv")

####################################################################################################################
##Assign closest gauge location to data points. 
####################################################################################################################
gauging = "USGS"
Eff_widths = Eff_widths[Eff_widths$Difference==0&Eff_widths$cloud<=10,]
Eff_widths$Date = as.Date(as.POSIXct(Eff_widths$`system:time_start`/1000, origin = "1970-01-01"))
test = Eff_widths[Eff_widths$Date> as.Date("2014-12-31"),]
data = distinct(test, ID, .keep_all = TRUE)
tab = data
Lat_lon_df = read.csv("E:\\research\\GRWL\\GRWL_2015_present\\East_west_together\\NA_quantiles_combined.csv")
Lat_lon_df = read.dbf("E:\\research\\GRWL\\Subset_GRWL_1spc\\NA_100m_min\\na_sj_using_R_min_100.dbf")
Lat_lon_df = as.data.frame(Lat_lon_df)
data$GRWL_width_m = Lat_lon_df$dbf.width_m[match(data$ID, Lat_lon_df$dbf.ID_2)]
data$lon_dd = Lat_lon_df$dbf.lon_dd[match(data$ID, Lat_lon_df$dbf.ID_2)]
data$lat_dd= Lat_lon_df$dbf.lat_dd[match(data$ID, Lat_lon_df$dbf.ID_2)]
tab = data
rm(Lat_lon_df)
gageinfo = gageinfo[gageinfo$GRWL_width>99,]
gageinfo = gageinfo[gageinfo$lakeFlag==0,]
nGRWL = as.numeric(1)
gageinfo_coords = cbind(gageinfo$LONG, gageinfo$LAT)
data_coords = cbind(data$lon_dd, data$lat_dd)
#Determine nearest xsections
closestDF=as.data.frame(array(NA, c(nrow(as.vector(gageinfo)), nGRWL)))
distanceDF=as.data.frame(array(NA, c(nrow(as.vector(gageinfo)), nGRWL)))
for(i in 1:nrow(gageinfo_coords)){
  dist = distGeo(gageinfo_coords[i, 1:2], data_coords[,1:2])##### seems to be working. 
  close=order(dist, decreasing=FALSE)[1:nGRWL]
  closeID=data$ID[close]
  distance_to = sort(dist)[1:nGRWL]
  closestDF[i,]=closeID
  distanceDF[i,]=distance_to 
}
#determine closest xsections to each gage
Site_number_xsections=cbind(gageinfo$SITE_NUM, closestDF)
Site_number_distances=cbind(gageinfo$SITE_NUM, distanceDF)
#filter out gages with no xsections within 500m. 
gage_filter = Site_number_distances[,(nGRWL+1)] < 500
Site_number_xsections = Site_number_xsections[gage_filter,]
##################################################################################################
##Functions.
##################################################################################################
usgs_processing = function(usgs_w){
  usgs_w$measurement_dt <- substr(usgs_w$measurement_dt, 0, 10)
  usgs_w_filter = subset(usgs_w, usgs_w$measured_rating_diff!= "Poor")# & usgs_w$chan_loc_dist <= 500)
  ##Just using the dishcharge values found in the width dataset - slightly different from discharge values.
  w = as.vector(usgs_w_filter$chan_width)
  q = as.vector(usgs_w_filter$discharge_va)
  ##convert to proper units for plotting. 
  w_m = as.numeric(w)*0.3048
  q_cms = as.numeric(q)*0.02832
  q_w = cbind(q_cms, w_m, usgs_w_filter$measurement_dt, usgs_w_filter$measured_rating_diff)
  return(q_w)
}

auto.legend.pos <- function(x,y,xlim=NULL,ylim=NULL) {
  if (dev.cur() > 1) {
    p <- par('usr')
    if (is.null(xlim)) xlim <- p[1:2]
    if (is.null(ylim)) ylim <- p[3:4]
  } else {
    if (is.null(xlim)) xlim <- range(x, na.rm = TRUE)
    if (is.null(ylim)) ylim <- range(y, na.rm = TRUE)
  }
  countIt <- function(a) {
    tl <- sum(x <= xlim[1]*(1-a)+xlim[2]*a & y >= ylim[1]*a+ylim[2]*(1-a))
    tr <- sum(x >= xlim[1]*a+xlim[2]*(1-a) & y >= ylim[1]*a+ylim[2]*(1-a))
    bl <- sum(x <= xlim[1]*(1-a)+xlim[2]*a & y <= ylim[1]*(1-a)+ylim[2]*a)
    br <- sum(x >= xlim[1]*a+xlim[2]*(1-a) & y <= ylim[1]*(1-a)+ylim[2]*a)
    c(topleft=tl,topright=tr,bottomleft=bl,bottomright=br)
  }
  for (k in seq(0.5,0.1,by=-0.05)) {
    a <- countIt(k)
    if (sum(a==0)>0) break
    #if (!is.na(sum(a))) break
    
  }
  names(a)[which(a==0)][1]   # may delete "[1]"
}

is.error <- function(
  expr,
  tell=FALSE,
  force=FALSE
)
{
  expr_name <- deparse(substitute(expr))
  test <- try(expr, silent=TRUE)
  iserror <- inherits(test, "try-error")
  if(tell) if(iserror) message("Note in is.error: ", test)
  if(force) if(!iserror) stop(expr_name, " is not returning an error.", call.=FALSE)
  # output:
  iserror
}

#################################################################################################
##Set up empty dataframes for collecting observations. 
#################################################################################################
start = as.Date("1979-01-01")
data_val = Eff_widths
RC_End = as.Date("2014-12-31") ###WAs 2014
RC_year = format(RC_End, "%Y")
RC_year_1 = format(RC_End+1, "%Y")
data_val$ID = data_val$ID
data_val$ID_2 = data_val$ID
data_val$calc_mean = data_val$Effective_width
xSecq=as.data.frame(matrix(numeric(), nrow =5, ncol = 11))
xSecw=as.data.frame(matrix(numeric(), nrow =5, ncol = 11))
xSecIDcol=grep("V", names(Site_number_xsections))
mInd = array(5, dimnames = NULL)
gage_stats = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 22))
l_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
u_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
width_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
colnames(gage_stats)= c("Site_number", "GRWL_width_m","n_Landsat_obs","R_2", "R", "RMSE", "p_val","Bias", "RRMSE","avg_std", "change", "RRMSE_median", "std_Q","STDE", "KGE", "NSE", "rBias",
                        "SDRR", "MRR", "NRMSE", "Q_50", "W_50")
usgs_all = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
as.data.frame(gage_stats)
gage_stats_col1 = as.vector(1)
gage_stats_col2 = as.vector(1)
gage_stats_GRADES_col1 = as.vector(1)
gage_stats_GRADES_col2 = as.vector(1)
paired_df_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
rmse = 1
width_grouping = 30
percentiles = c(0.05, 0.95)
training = 0.8

##Make sure the proper number of characters are in each gauge id. 
for(i in 1:nrow(Site_number_xsections)){
  if(nchar(Site_number_xsections$`gageinfo$SITE_NUM`[i])==7){
    Site_number_xsections$`gageinfo$SITE_NUM`[i] = paste0("0", Site_number_xsections$`gageinfo$SITE_NUM`[i])
  }
}
################################################################################################################
##processing each gauge location.
################################################################################################################
url = paste0("http://waterdata.usgs.gov/nwis/measurements?site_no=", 
             Site_number_xsections$`gageinfo$SITE_NUM`, "&agency_cd=USGS&format=rdb_expanded")
outPath = paste(dailyWdir, '/', gauges_ryan, ".csv", sep="")
for(i in 1:nrow(Site_number_xsections)){
  ###Filter to widths near gauge location.   
  print(i)
  xSecIDcol = Site_number_xsections[i,]
  xSecID=xSecIDcol[2:ncol(xSecIDcol)]
  mInd = match(xSecID, data$ID) #changed from tab
  xSecw = as.data.frame(seq(1:100))
  xSecq = as.data.frame(seq(1:100))
  mInd_paired = which(data_val$ID %in% xSecID)
  paired_df = data_val[mInd_paired,]
  usgs_q = try(usgs_processing(try(read.table(url[i], sep="\t", fill=T, skip=14, header=T), silent=T))) ##Usgs Make a switch for USGS vs Canadian. 
  if(is.error(usgs_q)){next}else{}
  usgs_q = as.data.frame(usgs_q)
  usgs_q$w_m = as.numeric(as.character(usgs_q$w_m))
  usgs_q$Date = as.Date(usgs_q$V3, format = "%Y-%m-%d")
  usgs_q = usgs_q[usgs_q$Date<as.Date("2021-01-01"),]
  usgs_q_1984 = usgs_q[usgs_q$Date>as.Date("1984-03-01", format = "%Y-%m-%d"),]
  joined = inner_join(paired_df, usgs_q)
  if(nrow(joined)>0){
    usgs_all[i,1:nrow(usgs_q_1984)] = usgs_q_1984$w_m
    mn = min(c(joined$calc_mean, joined$w_m), na.rm = TRUE)
    mx = max(c(joined$calc_mean, joined$w_m), na.rm = TRUE)
    plot(joined$calc_mean, joined$w_m, xlim = c(mn, mx), ylim = c(mn, mx))
    abline(0,1)
    gage_stats$Site_number[i]=Site_number_xsections[i,1]
    gage_stats$n_Landsat_obs[i] = nrow(joined)
    l_vals[i,1:nrow(joined)] = joined$calc_mean
    u_vals[i,1:nrow(joined)] = joined$w_m
  }
}
#######################################################################################################
##Filter results and plot Figure 4. 
#######################################################################################################
landsat = unlist(l_vals)
usgs = unlist(u_vals)
landsat = landsat[usgs>0]
usgs = usgs[usgs>0]
mn = min(c(landsat, usgs), na.rm = TRUE)
mx = max(c(landsat, usgs), na.rm = TRUE)
mn = 0
mx = 1300
par(pty = "s")
plot(usgs, landsat, yaxt = "n", xaxt = "n", pch = 19, col = rgb(red = 0, green = 0, blue = 0, alpha = 0.1),
     xlim = c(mn, mx), ylim = c(mn, mx), log = "", asp = 1, lwd = 2, bty = "n", xlab = "", ylab = "")
box(lwd = 2)
mtext("RODEO Width (m)", side = 2, line = 3, cex = 1.25)
mtext("In Situ Width (m)", side = 1, line = 3, cex = 1.25)
spacing = c(0,300,600,900,1200)
axis(side = 1, labels = TRUE, at = spacing, cex.axis = 1)#, at = c(0, 300, 600, 900))
axis(side = 2, las = 2, at = spacing, cex.axis = 1)
#abline(0,1, col = "black", lwd = 2, lty = 2)
regression = lm(landsat~usgs)
reg_sen = mblm::mblm(landsat~usgs)
regression = reg_sen
text(x = (420), y = (mx*.95), paste("Y =", signif(regression$coefficients[[2]],3),"X", signif(regression$coefficients[[1]],3)), col = "red",
     cex = 1.25)

###Stats
error = landsat - usgs
error = na.omit(error)
mae = mean(abs(error))
bias = mean(error)
rmse = sqrt(mean(error^2))









