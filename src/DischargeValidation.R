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

#########################################################################################################################
##USGS or Canada validation gauges. 
#########################################################################################################################
gauging = "Canada"
###Read in effective widths csv files. 
if(gauging == "USGS"){
  Eff_widths = map_df(list.files("USGSvalidation", full.names = TRUE), ~vroom(.x))
  gageinfo = read.csv("gaugeTable.csv")
  } else{
  Eff_widths = map_df(list.files("GRDCvalidation", full.names = TRUE), ~vroom(.x))
  gageinfo = read.csv("GRDCgaugeTable.csv")
  gageinfo$SITE_NUM = gageinfo$Station_Num
}

PathToGRADES7 = "PathToGRADES7"
PathToGRADES8 = "PathToGRADES8"

##########################################################################################################################
##Filter out Landsat widths and determine which COMIDs are in dataset. 
##########################################################################################################################
Eff_widths$Date = as.Date(as.POSIXct(Eff_widths$`system:time_start`/1000, origin = "1970-01-01"))
Eff_widths = Eff_widths[Eff_widths$cloud<=10&Eff_widths$Difference==0,]
Eff_widths = data.table(Eff_widths)
comidCols = grep("COMID", colnames(Eff_widths))
Eff_widths = as.data.frame(Eff_widths)
uniqueComid = unique(na.omit(unlist(Eff_widths[,comidCols])))
uniqueComid7 = uniqueComid[grep("^7", uniqueComid)]
uniqueComid8 = uniqueComid[grep("^8", uniqueComid)]
uniqueComid = uniqueComid7

#########################################################################################################################
##Read in GRADES_7 flowlines. 
#########################################################################################################################
# read in netCDF file:
ncIn = nc_open(PathToGRADES7)
allComid = ncvar_get(ncIn, "COMID")

roi = allComid%in%uniqueComid

start = 1
end = ncIn$var$Q$varsize[2]

j = 1
interval  = 50 
year = interval # assuming no leap years

output=list()
while (start < end){
  print(paste("chunk", j))
  # for last year: 
  if ((start+interval) > ncIn$var$Q$varsize[2]){
    interval = ncIn$var$Q$varsize[2] - start
  }
  nc = ncvar_get(ncIn, "Q", 
                 start=c(1, start), # starting index of netCDF to read in 
                 count=c(ncIn$var$Q$varsize[1], interval)) # ending index of netCDF to read in 
  nc = cbind(allComid, nc)
  timeSeries = seq(start, start +interval-1,1)
  nc_filt = nc[roi,]
  nc_filt = as.data.frame(nc_filt)
  colnames(nc_filt) = c("comid", timeSeries)
  tall = melt(nc_filt, id.var ="comid")
  output[[j]] = tall
  start = start+interval
  j = j+1
}
out7 = rbindlist(output)

uniqueComid = uniqueComid8
#########################################################################################################################
##Read in GRADES_8 flowlines and join with GRADES_7
########################################################################################################################## read in netCDF file:
ncIn = nc_open(PathToGRADES8)
allComid = ncvar_get(ncIn, "COMID")
roi = allComid%in%uniqueComid
start = 1
end = ncIn$var$Q$varsize[2]
j = 1
interval  = 50 
year = interval # assuming no leap years
output=list()
while (start < end){
  print(paste("chunk", j))
  if ((start+interval) > ncIn$var$Q$varsize[2]){
    interval = ncIn$var$Q$varsize[2] - start
  }
  nc = ncvar_get(ncIn, "Q", 
                 start=c(1, start), # starting index of netCDF to read in 
                 count=c(ncIn$var$Q$varsize[1], interval)) # ending index of netCDF to read in 
  nc = cbind(allComid, nc)
  timeSeries = seq(start, start +interval-1,1)
  nc_filt = nc[roi,]
  nc_filt = as.data.frame(nc_filt)
  colnames(nc_filt) = c("comid", timeSeries)
  tall = melt(nc_filt, id.var ="comid")
  output[[j]] = tall
  start = start+interval
  j = j+1
}
out8 = rbindlist(output)
out = bind_rows(out7, out8)
Gauge_comid = out

Gauge_comid = as.data.frame(Gauge_comid)
Gauge_comid$time = as.numeric(as.character(Gauge_comid$variable))
Gauge_comid$time = Gauge_comid$time-1
start = as.Date("1979-01-01")
Gauge_comid$date = start+Gauge_comid$time
Gauge_comid$Q = as.numeric(Gauge_comid$value)
Gauge_comid$COMID = Gauge_comid$comid

#########################################################################################################################
##Determine closest data points for each gauge location.  
#########################################################################################################################
test = Eff_widths#[Eff_widths$Date> as.Date("2014-12-31"),]
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
plot(gageinfo$e, gageinfo$n, type="p")
closestDF=as.data.frame(array(NA, c(nrow(as.vector(gageinfo)), nGRWL)))
distanceDF=as.data.frame(array(NA, c(nrow(as.vector(gageinfo)), nGRWL)))
for(i in 1:nrow(gageinfo_coords)){
  dist = distGeo(gageinfo_coords[i, 1:2], data_coords[,1:2])##### seems to be working. 
  close=order(dist, decreasing=FALSE)[1:nGRWL]
  closeID=data$ID[close]
  distance_to = sort(dist)[1:nGRWL]
  closestDF[i,]=closeID
  distanceDF[i,]=distance_to #was close
}
#determine closest xsections to each gage
Site_number_xsections=cbind(gageinfo$SITE_NUM, closestDF)
Site_number_distances=cbind(gageinfo$SITE_NUM, distanceDF)
#filter out gages with no xsections within 500m. 
gage_filter = Site_number_distances[,(nGRWL+1)] < 500
Site_number_xsections = Site_number_xsections[gage_filter,]
Gages_char = as.character(Site_number_xsections$`gageinfo$SITE_NUM`)
wd_q = "E:/research/GRDC/DailyQ/"
Gages_char = as.character(Site_number_xsections$`gageinfo$SITE_NUM`)
usgs_q_list = paste(wd_q, Gages_char, "_Q_Day.Cmd", ".txt", sep="")
files = list.files(wd_q)
gauges = gsub("_Q_Day.Cmd.txt", "",files)

#Make sure there are proper number of characters in USGS gauge numbers. 
if(gauging == "USGS"){
for(i in 1:nrow(Site_number_xsections)){
  if(nchar(Site_number_xsections$`gageinfo$SITE_NUM`[i])==7){
    Site_number_xsections$`gageinfo$SITE_NUM`[i] = paste0("0", Site_number_xsections$`gageinfo$SITE_NUM`[i])
  }
}
} else{}
##################################################################################################
##Functions. 
##################################################################################################
usgs_q_processing = function(usgs_q){
  q_v = as.vector(usgs_q[,4])
  q_c = as.character(usgs_q[4])
  q_n = as.numeric(q_v)
  q= q_n *0.02832
  usgs_q = cbind(usgs_q, q)
  as.character(usgs_q$datetime)
  return(usgs_q)
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

#########################################################################################################################
##Setup empty dataframes and get data in correct divisions for validation. 
#########################################################################################################################
start = as.Date("1979-01-01")
data_val = Eff_widths
RC_End = as.Date("2013-12-31") ###WAs 2014
RC_year = format(RC_End, "%Y")
RC_year_1 = format(RC_End+1, "%Y")
data_val$ID = data_val$ID
data_val$ID_2 = data_val$ID
data_val$calc_mean = data_val$Effective_width
data_val_8415 = data_val[data_val$Date< (RC_End +1),] ##Changed from 2015-01-01
data_val_1520 = data_val[data_val$Date> RC_End,] ##Changed from 2014-12-31
data_val = data_val_1520
tab$ID = tab$id
tab$ID_2 = tab$id
tab$width_m = data_val$width_m[match(tab$ID, data_val$ID)]
tab$width_m = data$width_m[match(tab$ID_2, data$ID_2)]
data_val_8415$Date = as.Date(data_val_8415$Date)
xSecq=as.data.frame(matrix(numeric(), nrow =5, ncol = 11))
xSecw=as.data.frame(matrix(numeric(), nrow =5, ncol = 11))
xSecIDcol=grep("V", names(Site_number_xsections))
mInd = array(5, dimnames = NULL)
rangedf_1 = as.data.frame(matrix(numeric(), nrow = 1, ncol = 4))
gage_stats = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 22))
gage_stats_GRADES = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 14))
colnames(gage_stats_GRADES)= c("Site_number", "GRWL_width_m","n_Landsat_obs","R_2", "R", "RMSE", "p_val","Bias", "RRMSE", "avg_std", "change", 'RRMSE_median', "std_Q", "STDE")
as.data.frame(gage_stats_GRADES)
l_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
Mean_grades = as.vector(nrow(Site_number_xsections))
grades_sd = as.vector(nrow(Site_number_xsections))
u_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
cloud_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
sd_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
sd_vals_1 = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
corrected_sd = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
width_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
q_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
gage_quants_q = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 100))
gage_quants_w = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 100))
colnames(gage_stats)= c("Site_number", "GRWL_width_m","n_Landsat_obs","R_2", "R", "RMSE", "p_val","Bias", "RRMSE","avg_std", "change", "RRMSE_median", "std_Q","STDE", "KGE", "NSE", "rBias",
                        "SDRR", "MRR", "NRMSE", "Q_50", "W_50")
all_q_df= as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
all_dt_df = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
as.data.frame(gage_stats)
gage_stats_col1 = as.vector(1)
gage_stats_col2 = as.vector(1)
gage_stats_GRADES_col1 = as.vector(1)
gage_stats_GRADES_col2 = as.vector(1)
paired_df_vals = as.data.frame(matrix(numeric(), nrow =nrow(Site_number_xsections), ncol = 20))
rmse = 1
width_grouping = 30
percentiles = c(0.05, 0.95)
Gauge_comid = as.data.table(Gauge_comid)
setkey(Gauge_comid, COMID)

#########################################################################################################################
##Process each gauge location. 
#########################################################################################################################
for(i in 1:nrow(Site_number_xsections)){
  ###Filter to widths near gauge location.   
  print(i)
  xSecIDcol = Site_number_xsections[i,]
  xSecID=xSecIDcol[2:ncol(xSecIDcol)]
  mInd = match(xSecID, data$ID) #changed from tab
  xSecw = as.data.frame(seq(1:100))
  xSecq = as.data.frame(seq(1:100))
  mInd_paired = which(data_val_8415$ID %in% xSecID)
  paired_df = data_val_8415[mInd_paired,]
  all_q = Eff_widths[which(Eff_widths$ID%in%xSecID),]
  ###Pair largest GRADES flowline with same day Landsat widths. 
  ################################################################################################    
  if(nrow(paired_df)>1){
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
    if(all(is.na(place$testing_df))){next} else{
      placement = place$list[place$testing_df==max(place$testing_df, na.rm = TRUE)]
      paired_df$COMID = placement[1]
    }
    Gauge_comid_filter1  = Gauge_comid[.(paired_df$COMID), allow.cartesian = TRUE]
    paired_df = left_join(paired_df, Gauge_comid_filter1, by = c("COMID" = "COMID", "Date" = "date"))
    paired_df = distinct(paired_df, Date, .keep_all = TRUE)
    if(all(is.na(paired_df$Q))){next}
    #########################################################################################################    
    p_q = quantile(paired_df$Q, probs = seq(percentiles[1],percentiles[2],.01), na.rm = TRUE)
    p_w = quantile(paired_df$calc_mean, probs = seq(percentiles[1],percentiles[2],.01), na.rm = TRUE)
    l_q = length(p_q)
    l_w = length(p_w)
    gage_quants_q[i,1:l_q] = p_q
    gage_quants_w[i,1:l_w] = p_w
    ###Determine widths post 2014 and aggregate to same day in case of overlapping observations. 
    data_val_subset = subset(data_val, xSecID[,1] ==data_val$ID)
    if(nrow(data_val_subset)>0){
      data_val_sub_agg = aggregate(data_val_subset$calc_mean, by=list(data_val_subset$Date), FUN=max)} else{next}
    w_rc= p_w
    q_rc = p_q
    y = q_rc
    x = w_rc
    paired_df = paired_df[paired_df$Q>=p_q[[1]]&paired_df$Q<=rev(p_q)[[1]],]
    paired_df = paired_df[!is.na(paired_df$`system:index`),]
    if(length(unique(x))>1){
      ###Create function to estimate Q from a Landsat width based on quantile pair up. 
      spl = approxfun(x, y)
    } else{next}
    ###Read in either USGS or GRDC gauges and process them. 
    if(gauging == "USGS"){
      usgs_q = try(usgs_q_processing(rawDailyData <- readNWISdv(Site_number_xsections$`gageinfo$SITE_NUM`[i],"00060","1984-01-01"))) ##Usgs Make a switch for USGS vs Canadian. 
    }else{
      ##Read in validation data. 
      ################################################################################GRDC
      usgs_q = try(read.table(usgs_q_list[i], stringsAsFactors = FALSE))##GRDC
      if(is.error(usgs_q)){next}else{
      usgs_q$V1 = substr(usgs_q$V1, 0, 10)
      usgs_q$datetime = usgs_q$V1
      usgs_q$q = as.numeric(usgs_q$V2)
      usgs_q$Date = as.Date(usgs_q$datetime, format = "%Y-%m-%d")
      }
    }
    ####################################################################################  
    if(is.error(usgs_q)){next}else{}
    usgs_q = usgs_q[usgs_q$q>=0,]
    usgs_q = usgs_q[usgs_q$Date<as.Date("2021-01-01"),]
    usgs_q$Year = format(usgs_q$Date, "%Y")
    if(nrow(usgs_q)==0)next
    AnnualAvg = aggregate(usgs_q$q, by=list(usgs_q$Year), mean, na.rm = TRUE)
    gage_stats$Q_50[i]=mean(AnnualAvg$x, na.rm = TRUE)
    data_val_sub_agg$Date = data_val_sub_agg$Group.1
    joined = inner_join(data_val_sub_agg, usgs_q)
    joined$model = spl(joined$x)
    if(length(na.omit(joined$model))>2&length(na.omit(joined$q))>2){} else{next}
    ###Plot the rating curve (black line) and GRADES Q (red circles), width bins (gray vertical lines), and estimated Landsat Q (blue circles).             
    par(mfrow =c(2, 2), mar = c(4, 4, 5, 4))
    layout(matrix(c(3,3,3,3,3,3,3,3,3,3,3,1,1,1,1,1,0,2,2,2,2,2), 2, 11, byrow = TRUE))
    plot(c(min(x, na.rm = TRUE), max(x, na.rm = TRUE)),c(min(spl(x), na.rm = TRUE), max(spl(x), na.rm = TRUE)), col = "white", 
         main=paste("USGS Gage:", Site_number_xsections$`gageinfo$SITE_NUM`[i]),
         log = "", type="p", ylab="Discharge (cms)",
         xlab="W (m)",)
    points(paired_df$calc_mean, paired_df$Q, col = "red")
    lines(x ,spl(x), col = "black", lwd = 2)
    points(joined$x,joined$model,pch = 1, col = "blue")
    modelQ = spl(paired_df$calc_mean)
    all_q_df[i,1:nrow(all_q)] = spl(all_q$Effective_width)
    all_dt_df[i,1:nrow(all_q)] = as.character(all_q$Date)
    difference = spl(paired_df$calc_mean)-paired_df$Q
    q_vals[i,1:length(modelQ)] = modelQ
    Mean_grades[i] = mean(paired_df$Q, na.rm = TRUE)
    
      ###Create a dataframe of Landsat width and the SD from GRADES within that width
      landsat_df = as.data.frame(joined$x)

      ###calculate group mean from rsq.
      paired_df$rsq = spl(paired_df$calc_mean)
     
      estimated_error = spl(paired_df$calc_mean) - paired_df$Q
      sd_vals_1[i,1:nrow(paired_df)] = estimated_error
      
      estimated_rmse = sqrt(mean((estimated_error^2), na.rm= TRUE))
      estimated_bias = mean(estimated_error, na.rm = TRUE)
      estimated_stde = sd.p(var(estimated_error))
      
      corrected_bias = (estimated_bias)*(1/.045)
      corrected_stde = estimated_stde*(1/1.29)
      corrected_rmse = sqrt((corrected_stde^2) + abs(corrected_bias^2))

      ###Replace the width value with the estimated Q from Landsat.             
      landsat_df[,1] = joined$model
      landsat_df$minimum = landsat_df[,1] - estimated_rmse
      landsat_df$maximum = landsat_df[,1] + estimated_rmse
      
      arrows(y0 = landsat_df$minimum, x0= joined$x, 
             y1 = landsat_df$maximum, x1 = joined$x, col = "indianred" , 
             code=3, angle = 180, length = 0, lwd = 0.5)
      ###Create a dataframe for error estimates on future plots based on date.             
      landsat_df = cbind(joined$model, landsat_df)
      ###Add in legend.             
      location = auto.legend.pos(na.omit(x), na.omit(y))
      legend(location,
             legend = c(paste("Rating curve (1984-", RC_year, ")", sep = ""),paste("Landsat widths (", RC_year_1, "-2020)", sep = ""),
                        "GRADES"), 
             col = c("black", "blue",
                     "red"),
             lty=c(1, NA, NA),
             lwd = c(2, NA, NA),
             pch = c(NA, 01, 01), 
             bty = "n", 
             text.col = "black", 
             horiz = FALSE, cex = 0.6)
      u_l_mn = min(spl(x), na.rm = TRUE)
      u_l_mx = max(spl(x), na.rm = TRUE)
      
      
      plot(c(u_l_mn, u_l_mx),c(u_l_mn, u_l_mx), type = "n", xlab = "In situ Discharge (cms)", ylab = "Landsat Discharge (cms)")
      points(joined$q, joined$model, col = "blue")
      abline(0,1)
      arrows(x0 = joined$q, y0=  landsat_df$minimum, 
             x1 = joined$q, y1 = landsat_df$maximum, col = "indianred",
             code=3, angle = 180, length = 0, lwd = 0.5)
      ###Calculate various statistics and add them to the gage_stats dataframe.             
      gage_stats$n_Landsat_obs[i] = length(na.omit(joined$model))
      mtext(side = 3, line = 1.5, adj = 0, text = paste0("R=",signif(Rvalue(joined$model, joined$q), 3)), cex = 0.75)
      mtext(side = 3, line = 1.5, adj = 1, text = paste0("RRMSE=",signif(RRMSE(joined$model, joined$q), 3)), cex = 0.75)
      mtext(side = 3, line = 4.5, adj = 1, text = paste0("rBias=",signif(rBias(joined$model, joined$q), 3)), cex = 0.75)
      mtext(side = 3, line = 4.5, adj = 0, text = paste0("KGE=",signif(KGE(joined$model, joined$q), 3)), cex = 0.75)
      mtext(side = 3, line = 3, adj = 1, text = paste0("NRSME=",signif(NRMSE(joined$model, joined$q), 3)), cex = 0.75)
      error = joined$model-joined$q
      rmse = sqrt(mean(error^2, na.rm = TRUE))
      gage_stats$RMSE[i] = rmse
      bias = mean(error, na.rm = TRUE)
      gage_stats$Bias[i] = bias
      stde = sqrt(var(error, na.rm = TRUE))
      stde = sd.p(error)
      gage_stats$STDE[i] = stde
      gage_stats$RRMSE[i] = RRMSE(joined$model, joined$q)
      gage_stats$KGE[i] = KGE(joined$model, joined$q)
      gage_stats$rBias[i] = rBias(joined$model, joined$q)
      gage_stats$NRMSE[i] = NRMSE(joined$model, joined$q)
      gage_stats$R[i]=Rvalue(joined$model, joined$q)
      gage_stats$GRWL_width_m[i] = data$GRWL_width_m[mInd]
      gage_stats$Site_number[i]=Site_number_xsections[i,1]
      l_vals[i,1:nrow(joined)] = joined$model
      u_vals[i,1:nrow(joined)] = joined$q
      width_vals[i,1:nrow(joined)] = joined$x
      
      ###add in hydrograph
      usgs_q$date = usgs_q$Date
      usgs_q_2015 = usgs_q %>%
        filter(usgs_q$date > (RC_End+1) & usgs_q$date < as.Date('2020-12-31'))
      plot(usgs_q_2015$date, usgs_q_2015$q,ylim = c(min(spl(x), na.rm = TRUE),
                                                    max(spl(x), na.rm = TRUE)), xlab = "", ylab = "Discharge (cms)", type = "l", col = "lightgray")
      points(joined$Date, joined$model, col = "blue")
      landsat_df$Date = joined$Date
      arrows(x0 = landsat_df$Date, y0 =  landsat_df$minimum, 
             x1 = landsat_df$Date, y1 = landsat_df$maximum, col = "indianred",
             code=3, angle = 180, length = 0, lwd = 0.5)
      landsat_diff = landsat_df$maximum - landsat_df$minimum
      title(paste("Hydrograph ", RC_year_1, "-2020", sep = ""), line = 0.5)
      legend("bottom", inset = c(0, -.2), legend = c("Rating curve", "In situ"),
             col = c("blue", "lightgray"),lty = c(NA, 1), pch = c(1, NA), bty = "n", xpd = TRUE, horiz = TRUE)
  } else{next} ### this one for min df. 
}
#########################################################################################################################################################
##combine USGS and GRDC outputs and filter any NA values. 
#########################################################################################################################################################
##Store first run values of USGS gauges. 
gage_stats= lapply(gage_stats, as.numeric)
gage_stats_can = gage_stats
u_vals_can = u_vals
l_vals_can = l_vals
sd_vals = sd_vals_1
sd_vals_can = sd_vals
grades_can = Mean_grades
width_can = width_vals
q_can = q_vals
all_q_can = all_q_df
all_dt_can = all_dt_df

##Combined with second run values of GRDC gauges. 
grades_usa = Mean_grades
sd_vals = sd_vals_1
sd_vals = bind_rows(sd_vals, sd_vals_can)
u_vals = bind_rows(u_vals, u_vals_can)
l_vals = bind_rows(l_vals, l_vals_can)
gage_stats = bind_rows(gage_stats, gage_stats_can)
grades = c(grades_usa, grades_can)
width_vals = bind_rows(width_vals, width_can)
q_vals = bind_rows(q_vals, q_can)
all_q_df = bind_rows(all_q_df, all_q_can)
all_dt_df = bind_rows(all_dt_df, all_dt_can)

##Filter out any NA values. 
filtering = !is.na(gage_stats$Site_number)
gage_stats= lapply(gage_stats, as.numeric)
gage_stats = as.data.frame(gage_stats)
nrow(gage_stats[filtering,])
apply(gage_stats[filtering,], 2, median, na.rm = TRUE)
