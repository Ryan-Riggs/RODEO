##Follow guidelines for validating USGS and GRDC gauges, in that order and joining them. 
library(caret)
library(sfsmisc)
library(scales)
library(ggpmisc)
library(scales)
library(gridExtra)
library(ggpubr)
library(grid)
library(cowplot)

#############################################################################################################################################
##Combined USGS and GRDC gauges and determine cal/val gauges. 
###############################################################################################################################
filtering = !is.na(gage_stats$Site_number)
gage_stats= lapply(gage_stats, as.numeric)
gage_stats = as.data.frame(gage_stats)
nrow(gage_stats[filtering,])
apply(gage_stats[filtering,], 2, median, na.rm = TRUE)
test = gage_stats[filtering,]
gage_stats_all = gage_stats
gage_stats_index = gage_stats_all==filtering
gage_stats = gage_stats[filtering, ]

l_vals = l_vals[filtering,]
u_vals = u_vals[filtering,]
sd_vals = sd_vals[filtering,]
#######################################################################################
##Subset into training/validation gauges. 
#######################################################################################
set.seed(1)
inTrain = createDataPartition(
  y = gage_stats$GRWL_width_m,
  p = 0.70,
  list = FALSE
)
# 
training = gage_stats[inTrain,]
testing = gage_stats[-inTrain,]
l_vals_testing = l_vals[-inTrain,]
l_vals_training = l_vals[inTrain,]
u_vals_testing = u_vals[-inTrain,]
u_vals_training = u_vals[inTrain,]
sd_vals_testing = sd_vals[-inTrain,]
sd_vals_training = sd_vals[inTrain,]

############################################################################################
##Error figures and calibration. Make sure to reselect testing data so it is with the calibrated errors. 
############################################################################################
###Total estimated error = dataframe where each row represents a gauge location and contains error bar lengths for each landsat observation. 
Total_estimated_error = sd_vals_training

###Estimated bias is the mean of the error bar lengths at each gauge. 
Bias_estimate = apply(Total_estimated_error, 1, FUN = mean, na.rm = TRUE)

###Estimated STDE is the sqrt of the variance of the error bar lengths at each gauge. 
sd.p=function(x){sd(na.omit(x))*sqrt((length(na.omit(x))-1)/length(na.omit(x)))}
STDE_estimate = apply(Total_estimated_error,1, FUN = sd.p)
###Estimated rmse is the rmse of the error bar lengths at each gauge. 
rmse_estimate = Total_estimated_error^2
rmse_estimate = sqrt(apply(rmse_estimate, 1, FUN = mean, na.rm = TRUE))

STDE_estimate[STDE_estimate==0] = NA
Bias_estimate[Bias_estimate==0] = NA
STDE_estimate[STDE_estimate==0] = NA
rmse_estimate[rmse_estimate==0] = NA

###Converts infinite values to NA. 
Bias_estimate[is.infinite(Bias_estimate)] = NA
STDE_estimate[is.infinite(STDE_estimate)] = NA
rmse_estimate[is.infinite(rmse_estimate)] = NA

####For later figure. 
Bias_estimate_train = Bias_estimate
STDE_estimate_train = STDE_estimate
rmse_estimate_train = rmse_estimate
###plot rmse vs abs(Bias_estimate)^2+(STDE_estimate^2) to double check that they line up. 
plot((rmse_estimate), sqrt((STDE_estimate^2) + abs(Bias_estimate^2)), xlim = c(min(rmse_estimate, na.rm = TRUE), max(rmse_estimate, na.rm = TRUE)),
     ylim = c(min(rmse_estimate, na.rm = TRUE), max(rmse_estimate, na.rm = TRUE)))
abline(0,1)

###Set up estimated error metrics for comparison.  
error_metrics = c("Bias", "STDE", "RMSE")
error_colors = c("gold", "blue", "forestgreen")
error_colors = c("firebrick", "dodgerblue4", "darkmagenta")
error_colors = brewer.pal(3,"Dark2")
estimated_df = cbind(abs(Bias_estimate), STDE_estimate, rmse_estimate)
estimated_df = as.data.frame(estimated_df)
colnames(estimated_df) = error_metrics

###Set up actual error metrics for comparison. 
actual_df = cbind(abs(training$Bias), training$STDE, training$RMSE)
actual_df = as.data.frame(actual_df)
colnames(actual_df) = error_metrics


########################################################################################################################
##Determine calibration values. 
########################################################################################################################
combined_values = c(abs(Bias_estimate), abs(training$Bias), STDE_estimate, training$STDE, rmse_estimate, training$RMSE)
mn = min(combined_values, na.rm = TRUE)
mx = max(combined_values, na.rm = TRUE)
output = as.vector(3)
error_plot = function(f,v){
  x = actual_df[,f]
  y = estimated_df[,f]
  lm = lm(y~x+0)
  output = lm$coefficients[1]
  r2 = summary(lm)$r.squared
  print(r2)
  return(output)
}

calibration = c(1:3)
out = lapply(calibration, error_plot)

#########################################################################################################################
##Determine validation gauges and apply calibration. 
#########################################################################################################################
Total_estimated_error = sd_vals_testing

###Estimated bias is the mean of the error bar lengths at each gauge. 
Bias_estimate = apply(Total_estimated_error, 1, FUN = mean, na.rm = TRUE)

###Estimated STDE is the sqrt of the variance of the error bar lengths at each gauge. 
STDE_estimate = apply(Total_estimated_error,1, FUN = sd.p)

###Estimated rmse is the rmse of the error bar lengths at each gauge. 
rmse_estimate = Total_estimated_error^2
rmse_estimate = sqrt(apply(rmse_estimate, 1, FUN = mean, na.rm = TRUE))

STDE_estimate[STDE_estimate==0] = NA
Bias_estimate[Bias_estimate==0] = NA
STDE_estimate[STDE_estimate==0] = NA
rmse_estimate[rmse_estimate==0] = NA

###Converts infinite values to NA. 
Bias_estimate[is.infinite(Bias_estimate)] = NA
STDE_estimate[is.infinite(STDE_estimate)] = NA
rmse_estimate[is.infinite(rmse_estimate)] = NA

corr_bias = (Bias_estimate)*(1/out[[1]])
corr_bias = abs(corr_bias)
corr_stde = STDE_estimate*(1/out[[2]])
corr_rmse = sqrt((corr_stde^2) + abs(corr_bias^2))

Total_estimated_error_all = sd_vals

###Estimated bias is the mean of the error bar lengths at each gauge. 
Bias_estimate_all = apply(Total_estimated_error_all, 1, FUN = mean, na.rm = TRUE)

###Estimated STDE is the sqrt of the variance of the error bar lengths at each gauge. 
STDE_estimate_all = apply(Total_estimated_error_all,1, FUN = sd.p)

###Estimated rmse is the rmse of the error bar lengths at each gauge. 
rmse_estimate_all = Total_estimated_error_all^2
rmse_estimate_all = sqrt(apply(rmse_estimate_all, 1, FUN = mean, na.rm = TRUE))
###########################################################################################################################
##Convert to Tall dataframe and plot. 
###########################################################################################################################
training_select = gage_stats%>%dplyr::select(Bias, STDE, RMSE)
training_select$Bias = abs(training_select$Bias)
training_tall = gather(training_select, key = "Variable", value = "value")
tall_cal = training_tall
tall_cal$calibration = "Calibrated"
tall_uncal = training_tall
tall_uncal$calibration = "Uncalibrated"
tall = bind_rows(tall_uncal, tall_cal)
tall$calibration = as.factor(tall$calibration)
yVals = c(abs(Bias_estimate_all), STDE_estimate_all, rmse_estimate_all,abs(Bias_estimate_all), STDE_estimate_all, rmse_estimate_all)
tall$y = yVals
my.formula = y~x+0
labs = c("Uncalibrated", "Calibrated", "test")

##Plot
p1 = ggplot(data = tall[tall$calibration =="Uncalibrated",], aes(x = value, y = y, colour = factor(Variable, levels = c("Bias", "STDE", "RMSE"))))+
  geom_point(alpha =0.25)+facet_grid(factor(calibration, levels = c("Uncalibrated", "Calibrated"))~factor(Variable, levels = c("Bias", "STDE", "RMSE")))+coord_fixed(xlim = c(mn, mx), ylim = c(mn, mx))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")+geom_smooth(method = "lm", se = FALSE, formula = my.formula)+
  stat_poly_eq(data = tall[tall$calibration =="Uncalibrated",], formula = my.formula, 
               aes(label = paste(..eq.label..,"X", ..rr.label.., sep = "~~~")), 
               parse = TRUE)+ylab("Estimated Error")+xlab("Actual Error")+
  coord_trans(x="log10", y="log10", xlim = c(mn, mx), ylim = c(mn,mx))+
  guides(color = guide_legend(override.aes = list(size = 3) ) )+
  scale_x_continuous("Actual Error", breaks = c(.1,1,10,100,1000,10000), labels = trans_format('log10',math_format(10^.x)))+
  scale_y_continuous("Estimated Error", breaks = c(.1,1,10,100,1000,10000), labels = trans_format('log10',math_format(10^.x)))+
  theme_bw()+ggtitle("Uncalibrated")+
  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank(),
        strip.text.y = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.title.x = element_text(colour = "white"), legend.justification = "top",
        axis.title.y = element_blank(), legend.title = element_blank())+
  scale_color_manual(values =error_colors[1:3])
p1
training_select = testing%>%dplyr::select(Bias, STDE, RMSE)
training_select$Bias = abs(training_select$Bias)
training_tall = gather(training_select, key = "Variable", value = "value")
tall_cal = training_tall
tall_cal$calibration = "Calibrated"
tall_uncal = training_tall
tall_uncal$calibration = "Uncalibrated"
tall = bind_rows(tall_uncal, tall_cal)
tall$calibration = as.factor(tall$calibration)
yVals = c(abs(Bias_estimate), STDE_estimate, rmse_estimate,corr_bias, corr_stde, corr_rmse)
tall$y = yVals
my.formula = y~x+0
labs = c("Uncalibrated", "Calibrated", "test")
p2 = ggplot(data = tall[tall$calibration =="Calibrated",], aes(x = value, y = y, colour = factor(Variable, levels = c("Bias", "STDE", "RMSE"))))+
  geom_point(alpha = 0.25)+facet_grid(factor(calibration, levels = c("Uncalibrated", "Calibrated"))~factor(Variable, levels = c("Bias", "STDE", "RMSE")))+coord_fixed(xlim = c(mn, mx), ylim = c(mn, mx))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")+geom_smooth(method = "lm", se = FALSE, formula = my.formula)+
  stat_poly_eq(data = tall[tall$calibration =="Calibrated",], formula = my.formula, 
               aes(label = paste(..eq.label..,"X", ..rr.label.., sep = "~~~")), 
               parse = TRUE)+ylab("Estimated Error")+xlab("Actual Error")+
  guides(color = guide_legend(override.aes = list(size = 3) ) )+
  coord_trans(x="log10", y="log10", xlim = c(mn, mx), ylim = c(mn,mx))+
  scale_x_continuous("Actual Error", breaks = c(.1,1,10,100,1000,10000), labels = trans_format('log10',math_format(10^.x)))+
  scale_y_continuous("Estimated Error", breaks = c(.1,1,10,100,1000,10000), labels = trans_format('log10',math_format(10^.x)))+
  theme_bw()+ggtitle("Calibrated")+
  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank(),
        strip.text.y = element_blank(), strip.text.x = element_blank(), legend.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())+
  scale_color_manual(values =error_colors[1:3])
p2
y.grob = textGrob("Uncertainty", rot = 90, gp = gpar(cex = 1.25))
x.grob = textGrob("Error", gp = gpar(cex = 1.25))
fig = ggarrange(p1, p2, common.legend = T, legend = "right", align = "v", ncol= 1)
annotate_figure(fig, left = y.grob, bottom = x.grob)












