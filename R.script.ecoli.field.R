
# Packages
require(ggplot2)
require(sf)
require(mosaic)
require(proxy)
require(cowplot)
library(car)
require(tmap)
require(corrplot)
require(spdep)
require(lmtest)
require(osmdata)
require(tmaptools)
require(OpenStreetMap)
require(gridExtra)
require(MASS)
require(dplyr)
require(mgcv)
require(PMCMRplus)
require(effects)
require(vegan)
require(gridExtra)
require(tmap)
require(corrplot)
require(spdep)
require(lmtest)
require(osmdata)
require(tmaptools)
require(OpenStreetMap)
require(gridExtra)
require(FSA)
require(reshape2)

results_fen = read.csv("results_f_e_n.csv")
results_fen$Site = as.factor(results_fen$Site)
is.factor(results_fen$Site)

#Ecoli stats
hist(results_fen$Average_E_Coli_per_100ml)#big right skew
favstats(results_fen$Average_E_Coli_per_100ml) #Large range. Median far below mean
results_fen$ecoli_logged = log(results_fen$Average_E_Coli_per_100ml)
hist(results_fen$ecoli_logged)#not normally distributed
shapiro.test(results_fen$ecoli_logged) #0.7W and p value < 0.05. So still not normal
#Nutrients stats
favstats(results_fen$TON_mg_l) 
favstats(results_fen$NH4_ug_l) 
favstats(results_fen$PO4_ug_l)
favstats(results_fen$TN_mg_l)
favstats(results_fen$TP_ug_l)
favstats(results_fen$TDN_mg_l)
favstats(results_fen$PON_mg_l)
#TON and TN have much smaller ranges than NH4. PON slightly higher range than TDN. PO4 and TP both have large  ranges. 
hist(results_fen$TON_mg_l)
hist(results_fen$NH4_ug_l)
hist(results_fen$PO4_ug_l)
hist(results_fen$TN_mg_l)
hist(results_fen$TP_ug_l)
hist(results_fen$TDN_mg_l)
hist(results_fen$PON_mg_l)
shapiro.test(results_fen$TON_mg_l)
shapiro.test(results_fen$NH4_ug_l)
shapiro.test(results_fen$PO4_ug_l)
shapiro.test(results_fen$TN_mg_l)
shapiro.test(results_fen$TP_ug_l)
shapiro.test(results_fen$TDN_mg_l)
shapiro.test(results_fen$PON_mg_l)
#Only PON and TN are normally distirbuted
hist(results_fen$DO..mg.l.)
hist(results_fen$EC..μS.)
hist(results_fen$pH)
hist(results_fen$Temperature...C.)
shapiro.test(results_fen$DO..mg.l.)
shapiro.test(results_fen$EC..μS.)
shapiro.test(results_fen$pH)
shapiro.test(results_fen$Temperature...C.)

#Investigating how E Coli changes between the sites
kruskal.test(Average_E_Coli_per_100ml ~ Site, data = results_fen) #siginicant differences
#Poisson Generalized Linear Model 
poisson_glm_ecoli = glm(Average_E_Coli_per_100ml ~ Site, family = poisson(link = "log"), data = results_fen)
summary(poisson_glm_ecoli)
null_dev_pe = poisson_glm_ecoli$null.deviance
resid_dev_pe = poisson_glm_ecoli$deviance
100 * (null_dev_pe - resid_dev_pe) / null_dev_pe #the explained deviance (r2) is 98% 
resid_df_pe = df.residual(poisson_glm_ecoli) 
resid_dev_pe / resid_df_pe # very high dispersion (3629) = unreliable model
# Negative Binomial
nb_model_ecoli = glm.nb(Average_E_Coli_per_100ml ~ Site, data = results_fen)
resid_df_nbec = df.residual(nb_model_ecoli)
resid_dev_nbec = nb_model_ecoli$deviance
resid_df_nbec / resid_dev_nbec
null_dev_nbec = nb_model_ecoli$null.deviance
100 * (null_dev_nbec - resid_dev_nbec) / null_dev_nbec
# Compare AIC
AIC(nb_model_ecoli, poisson_glm_ecoli) #chosing nb as it is significantly smaller
plot(fitted(nb_model_ecoli), residuals(nb_model_ecoli, type = "pearson"),
     xlab = "Fitted values", ylab = "Pearson Residuals")+ 
  abline(h = 0, col = "red") # centered around 0 and not really a pattern which is good
qqnorm(residuals(nb_model_ecoli, type = "pearson"))
qqline(residuals(nb_model_ecoli, type = "pearson"), col = "red") #normally distributed qq plot
hist(residuals(nb_model_ecoli), xlab = "Residuals") #mostly normally distributed
cooksd <- cooks.distance(nb_model_ecoli)
# Plot Cook's Distance
plot(cooksd, main = "Cook's Distance for Observations",
     xlab = "Observation", ylab = "Cook's Distance", type = "h")
abline(h = 1, col = "red")  # Add a reference line at 1. No outliers 
which(cooksd > 1)  # 0 influential points dispraportionatley affecting the model
#Can now proceed with confidence after checking model
summary(nb_model_ecoli)

#Plotting the model
ecoli_means = read.csv("ecoli_means.csv")
ecoli_means$logged_means = log(ecoli_means$E.coli)
results_fen$ecoli_logged = log(results_fen$Average_E_Coli_per_100ml)
ecoli_boxplot = ggplot(results_fen, aes(x = factor(Site), y = ecoli_logged)) +
  geom_boxplot() +
  labs(x = "Site", y = expression("Logged "*italic("E. coli")*" (CFU / 100 mL)")) +
  theme_minimal()

log_ecoli_bar_chart = ggplot(ecoli_means, aes(x = factor(Site), y = logged_means)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Site",
       y = expression("Logged "*italic("E. coli")*" (CFU / 100 mL)")) +
  theme_minimal()
grid.arrange(log_ecoli_bar_chart, ecoli_boxplot, ncol = 2)


#Comparing sites before and after Cumberland basin.
results_pre_basin = results_fen[1:9, ] 
results_post_basin = results_fen[10:27, ]
results_4 = results_fen[10:12, ]
favstats(results_pre_basin$Average_E_Coli_per_100ml)#mean is 262,569/ml median is 216,800/ml IQR 216,800
favstats(results_post_basin$Average_E_Coli_per_100ml)#mean is 19,439 median is 195 IQR 140
wilcox.test(results_pre_basin$Average_E_Coli_per_100ml, results_post_basin$Average_E_Coli_per_100ml) # Mann Whit U test Confirms statistical difference as P value v small.
wilcox.test(results_pre_basin$Average_E_Coli_per_100ml, results_fen$Average_E_Coli_per_100ml[10:12])#Site 4 is jsut about different to Sites 1,2,3 but only just so bare in mind its v similar.
wilcox.test(results_post_basin$Average_E_Coli_per_100ml, results_fen$Average_E_Coli_per_100ml[10:12]) #Site 4 statistically different to 5, 6, 7, 8
# Combining to get limits to have consistent scale
all_logged = c(results_4$ecoli_logged, 
                results_pre_basin$ecoli_logged, 
                results_post_basin_no_4$ecoli_logged)
x_limits = range(all_logged*1.1, na.rm = TRUE)

boxplot_pre_basin = bwplot(results_pre_basin$ecoli_logged, xlab = expression("Logged "*italic("E. coli")*" (CFU / 100 mL)"), main = "Sites 1, 2, and 3", xlim = x_limits)
boxplot_4 = bwplot(results_4$ecoli_logged, xlab = expression("Logged "*italic("E. coli")*" (CFU / 100 mL)"), main = "Site 4", xlim = x_limits)
boxplot_post_basin = bwplot(results_post_basin$ecoli_logged, xlab = "Logged E.Coli/100 mL", main = "Sites 4, 5, 6,7, 8 and 9")
boxplot_post_basin_no_4 = bwplot(results_post_basin_no_4$ecoli_logged, xlab = expression("Logged "*italic("E. coli")*" (CFU / 100 mL)"), main = "Sites 5, 6,7, 8 and 9", xlim = x_limits)
boxplots_basins = plot_grid(boxplot_pre_basin, boxplot_4, boxplot_post_basin_no_4, nrow = 3, align = "v", rel_heights = c(1, 1, 1))
boxplots_basins #Shows the range has decreased signficantly without Site 4 

#field measurements tests
wilcox.test(results_post_basin$pH, results_pre_basin$pH)
wilcox.test(results_post_basin$DO..mg.l., results_pre_basin$DO..mg.l.) 
wilcox.test(results_post_basin$EC..μS., results_pre_basin$EC..μS.)
wilcox.test(results_post_basin$Temperature...C., results_pre_basin$Temperature...C.)
#Shows statistical with only DO  

field_measurement_model_ecoli = glm.nb(Average_E_Coli_per_100ml ~ Temperature...C. + DO..mg.l. + EC..μS. + pH, data = results_fen)
resid_df_field = df.residual(field_measurement_model_ecoli)
resid_dev_field = field_measurement_model_ecoli$deviance
resid_df_field / resid_dev_field
null_dev_field = nb_model_ecoli$null.deviance
100 * (null_dev_field - resid_dev_field) / null_dev_field
plot(fitted(field_measurement_model_ecoli), residuals(field_measurement_model_ecoli, type = "pearson"),
     xlab = "Fitted values", ylab = "Pearson Residuals")+
  abline(h = 0, col = "red") # slight clustering at right end
qqnorm(residuals(field_measurement_model_ecoli, type = "pearson"))
qqline(residuals(field_measurement_model_ecoli, type = "pearson"), col = "red") # not too bad
hist(residuals(field_measurement_model_ecoli), xlab = "Residuals")# not that normal
summary(field_measurement_model_ecoli) #only significant effect on e.coli is DO
#Plotting the DO
ggplot(results_fen, aes(x = ecoli_logged, y = `DO..mg.l.`)) +
         geom_point() +
         geom_smooth(method = "lm", se = FALSE, color = "steelblue") +
         labs(
           x = expression("Logged "*italic("E. coli")*" (CFU / 100 mL)"),
           y = expression("DO (mg "*L^{-1}*")"))+
         theme_minimal()
  

cor(results_fen$ecoli_logged, results_fen$TN_mg_l, method = "spearman") #very weak relationship between e.coli and TN which is expected

cor(results_fen$ecoli_logged, results_fen$TON_mg_l, method = "spearman")
cor(results_fen$ecoli_logged, results_fen$PO4_mg_l, method = "spearman")
cor(results_fen$ecoli_logged, results_fen$NH4_mg_l, method = "spearman")
#TON and NH4 were -0.69 and 0.69 respectively. po4 not significant  (0.3). 
# Generalized Additive Model (GAM) with a smooth term for po4 to check potential non-linear fit
gam_model_po4_ecoli = gam(ecoli_logged ~ s(PO4_ug_l), data = results_fen)
plot(fitted(gam_model_po4_ecoli), residuals(gam_model_po4_ecoli, type = "pearson"),
     xlab = "Fitted values", ylab = "Pearson Residuals")+
  abline(h = 0, col = "red") # slight clustering at right end
qqnorm(residuals(gam_model_po4_ecoli, type = "pearson"))
qqline(residuals(gam_model_po4_ecoli, type = "pearson"), col = "red") # acceptable
hist(residuals(gam_model_po4_ecoli), xlab = "Residuals")# acceptable
plot(fitted(gam_model_po4_ecoli), resid_gam_po4, 
     xlab = "Fitted values", ylab = "Residuals", 
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lwd = 2)
summary(gam_model_po4_ecoli)
#Insignificant (p = 0.07) with an adjusted R² of 0.332. No relationship between PO4 and E coli
#Model for TON and e coli
ton_ecoli = glm.nb(log(Average_E_Coli_per_100ml) ~ TON_mg_l, data = results_fen)
resid_df_tonec = df.residual(ton_ecoli)
resid_dev_tonec = ton_ecoli$deviance
resid_df_tonec / resid_dev_tonec
null_dev_tonec = ton_ecoli$null.deviance
100 * (null_dev_tonec - resid_dev_tonec) / null_dev_tonec 
plot(fitted(ton_ecoli), residuals(ton_ecoli, type = "pearson"),
     xlab = "Fitted values", ylab = "Pearson Residuals")+
  abline(h = 0, col = "red") # not normal really
qqnorm(residuals(ton_ecoli, type = "pearson"))
qqline(residuals(ton_ecoli, type = "pearson"), col = "red") # acceptable
hist(residuals(ton_ecoli), xlab = "Residuals")# acceptable
shapiro.test(residuals(ton_ecoli)) 
summary(ton_ecoli)# statistically significant decreases
ggplot(results_fen, aes(x = ecoli_logged, y = TON_mg_l)) +
  geom_point() +                                      
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  theme_minimal()+
  labs(x = "Average Logged E. coli CFUs per 100 mL",
     y = "TON (mg/L)")

# Model for nh4 and ecoli
nh4_ecoli = glm.nb(log(Average_E_Coli_per_100ml) ~ log(NH4_ug_l), data = results_fen)
resid_df_nhec = df.residual(nh4_ecoli)
resid_dev_nhec = nh4_ecoli$deviance
resid_df_nhec / resid_dev_nhec
null_dev_nhec = nh4_ecoli$null.deviance
100 * (null_dev_nhec - resid_dev_nhec) / null_dev_nhec 
plot(fitted(nh4_ecoli), residuals(nh4_ecoli, type = "pearson"),
     xlab = "Fitted values", ylab = "Pearson Residuals")+
  abline(h = 0, col = "red") # not normal really
qqnorm(residuals(nh4_ecoli, type = "pearson"))
qqline(residuals(nh4_ecoli, type = "pearson"), col = "red") # average
hist(residuals(nh4_ecoli), xlab = "Residuals")# right skew in residuals
summary(nh4_ecoli)# statistically significant increases
ggplot(results_fen, aes(x = ecoli_logged, y = NH4_mg_l)) +
  geom_point() +                                      
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  theme_minimal()+
  labs(x = "Average Logged E. coli CFUs per 100 mL",
       y = "NH4 (mg/L)")

#INTERACTION PLOTS
# Categorizing e coli into Low, Medium, High
results_fen$E_Coli_Category <- cut(results_fen$Average_E_Coli_per_100ml,
                                   breaks = quantile(results_fen$Average_E_Coli_per_100ml, probs = c(0, 0.33, 0.66, 1), na.rm = TRUE),
                                   labels = c("Low", "Medium", "High"),
                                   include.lowest = TRUE)
#NH4 interaction plot
nh4_interaction = ggplot(results_fen, aes(x = factor(Site), y = NH4_mg_l, color = E_Coli_Category, group = E_Coli_Category)) +
  stat_summary(fun = mean, geom = "line", size = 1) +  
  stat_summary(fun = mean, geom = "point", size = 2) + 
  labs(
    x = "Site",
    y = expression(NH[4]^"+"* " (mg "*L^{-1}*")"),
    color = expression(italic("E. coli") * " category")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    legend.position = "top"
  ) +
  scale_color_manual(values = c("blue", "red", "green"))

#TON interaction plot
ton_interaction = ggplot(results_fen, aes(x = factor(Site), y = TON_mg_l, color = E_Coli_Category, group = E_Coli_Category)) +
    stat_summary(fun = mean, geom = "line", size = 1) +  
    stat_summary(fun = mean, geom = "point", size = 2) + 
    labs(
      x = "Site",
      y = expression("TON (mg "*L^{-1}*")"),
      color = expression(italic("E. coli") * " category")
    ) +  # <— You were missing this closing parenthesis
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title = element_text(size = 12),
      legend.position = "top"
    ) +
    scale_color_manual(values = c("blue", "red", "green"))
  
   nh4_interaction
 ton_interaction
grid.arrange(ton_interaction, nh4_interaction, ncol = 2)

#Looking at EC and E.Coli relationship. Need to log EC because of the large range
results_fen$ec_logged = log(results_fen$EC..μS.)
cor(results_fen$ecoli_logged, results_fen$ec_logged, method = "spearman") #0.1
#0.1 coefficient =  no strong monotonic relationship  
#Use Generalized Additive Model (GAM) to investigate potential non-linear trend.
gam_model_ec_ecoli <- gam(ecoli_logged ~ s(ec_logged), data = results_fen)
plot(fitted(gam_model_ec_ecoli), residuals(gam_model_ec_ecoli, type = "pearson"),
     xlab = "Fitted values", ylab = "Pearson Residuals")+
  abline(h = 0, col = "red") # not normal really
qqnorm(residuals(gam_model_ec_ecoli, type = "pearson"))
qqline(residuals(gam_model_ec_ecoli, type = "pearson"), col = "red") # average
hist(residuals(gam_model_ec_ecoli), xlab = "Residuals")
shapiro.test(residuals(gam_model_ec_ecoli)) #normal
summary(gam_model_ec_ecoli)
#98% of the variability in the logged E. coli data is explained by the model.Unimodal relationship.Significant to 99%
# Plot GAM fit
ggplot(results_fen, aes(x = ec_logged, y = ecoli_logged)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x), col = "blue") +
  labs(x = expression("Logged "*italic("E. coli")*" (CFU / 100 mL)"),
       y = expression("EC ("*mu*"S cm"^-1*")")) +
  theme_minimal()


hist(results_fen$Average_Salmonella_per_100ml)
cor.test(x= results_fen$Average_Salmonella_per_100ml, y= results_fen$Average_E_Coli_per_100ml, method = 'spearman')
#0.8 and p value is very small. Very strong positive correlation

#Nutrient Callibrations
NH4_data = read.csv("NH4_data.csv")
# Subset the data for calibration
NH4_callibration = NH4_data[1:18, ]
NH4_callibration$TARGET_CONC = as.numeric(NH4_callibration$TARGET_CONC)
NH4_callibration_model <- lm(RESULT ~ TARGET_CONC, data = NH4_callibration)
R2_nh4_callibration = summary(NH4_callibration_model)$r.squared
ggplot(NH4_callibration, aes(x = TARGET_CONC, y = RESULT)) +
  geom_point() +                                      
  geom_smooth(method = "lm", se = FALSE, color = "blue") + 
  labs(x = "Target Concentration",
       y = "Reponse") +
  annotate("text", x = max(NH4_callibration$TARGET_CONC), 
           y = max(NH4_callibration$RESULT), 
           label = paste0("R² = ", round(R2_nh4_callibration, 3)), 
           hjust = 2, 
           vjust = 1, 
           size = 5, 
           color = "red")
PO4_data = read.csv("PO4_data.csv")
# Subset the data for calibration
PO4_callibration = PO4_data[1:22, ]
PO4_callibration$Target_conc = as.numeric(PO4_callibration$Target_conc)
PO4_callibration_model <- lm(Response ~ Target_conc, data = PO4_callibration)
R2_po4_callibration = summary(PO4_callibration_model)$r.squared
ggplot(PO4_callibration, aes(x = Target_conc, y = Response)) +
  geom_point() +                                      
  geom_smooth(method = "lm", se = FALSE, color = "blue") + 
  labs(x = "Target Concentration",
       y = "Reponse") +
  annotate("text", x = max(PO4_callibration$Target_conc), 
           y = max(PO4_callibration$Response), 
           label = paste0("R² = ", round(R2_po4_callibration, 3)), 
           hjust = 2, 
           vjust = 1, 
           size = 5, 
           color = "red")
TON_data = read.csv("TON_data.csv")
# Subset the data for calibration
TON_callibration = TON_data[1:10, ]
TON_callibration$Target_conc = as.numeric(TON_callibration$Target_conc)
TON_callibration_model <- lm(Response ~ Target_conc, data = TON_callibration)
R2_ton_callibration = summary(TON_callibration_model)$r.squared
ggplot(TON_callibration, aes(x = Target_conc, y = Response)) +
  geom_point() +                                      
  geom_smooth(method = "lm", se = FALSE, color = "blue") + 
  labs(x = "Target Concentration",
       y = "Reponse") +
  annotate("text", x = max(TON_callibration$Target_conc), 
           y = max(TON_callibration$Response), 
           label = paste0("R² = ", round(R2_ton_callibration, 3)), 
           hjust = 2, 
           vjust = 1, 
           size = 5, 
           color = "red")
#TN TP Data
tn_data = read.csv("TN_data.csv")
tn_data_bare = read.csv("TN_bare.csv")
tp_data = read.csv("TP_data.csv")
tp_data_bare = read.csv("TP_bare.csv")
tn_data_bare$Site = as.factor(tn_data_bare$Site)
tp_data_bare$Site = as.factor(tp_data_bare$Site)

#TNTP Callibration
tn_callibration = read.csv("TN_callibration.csv")
tn_callibration = tn_callibration[1:5, ]
TN_callibration_model <- lm(Response ~ Target_conc_ppm, data = tn_callibration)
R2_tn_callibration = summary(TN_callibration_model)$r.squared
ggplot(tn_callibration, aes(x = Target_conc_ppm, y = Response)) +
  geom_point() +                                      
  geom_smooth(method = "lm", se = FALSE, color = "blue") + 
  labs(x = "Target Concentration",
       y = "Reponse") +
  annotate("text", x = max(tn_callibration$Target_conc_ppm), 
           y = max(tn_callibration$Response), 
           label = paste0("R² = ", round(R2_tn_callibration, 3)), 
           hjust = 2, 
           vjust = 1, 
           size = 5, 
           color = "red")
tp_callibration = read.csv("TP_callibration.csv")
tp_callibration = tp_callibration[1:13, ]
TP_callibration_model <- lm(Response ~ Target_conc_ppb, data = tp_callibration)
R2_tp_callibration = summary(TP_callibration_model)$r.squared
ggplot(tp_callibration, aes(x = Target_conc_ppb, y = Response)) +
  geom_point() +                                      
  geom_smooth(method = "lm", se = FALSE, color = "blue") + 
  labs(x = "Target Concentration",
       y = "Reponse") +
  annotate("text", x = max(tp_callibration$Target_conc_ppb), 
           y = max(tp_callibration$Response), 
           label = paste0("R² = ", round(R2_tp_callibration, 3)), 
           hjust = 2, 
           vjust = 1, 
           size = 5, 
           color = "red")

ton_data_bare = read.csv("TON_data_bare.csv")
po4_data_bare = read.csv("po4_data_bare.csv")
nh4_data_bare = read.csv("nh4_data_bare.csv")
ton_data_bare$Site = as.factor(ton_data_bare$Site)
po4_data_bare$Site = as.factor(po4_data_bare$Site)
nh4_data_bare$Site = as.factor(nh4_data_bare$Site)



kruskal.test(Result ~ Site, data = ton_data_bare) #Differences between sites with TON
kruskal.test(Result ~ Site, data = po4_data_bare) #No differences between sites with po4
kruskal.test(Result ~ Site, data = nh4_data_bare) #Differences between sites with NH4

#Modelling above nutrients
#TON
favstats(ton_data_bare$Result)
lm_ton = lm(Result ~ Site, data = ton_data_bare) #Nutrients increase as sites go up
summary(lm_ton) #0.93 r2
plot(fitted(lm_ton), residuals(lm_ton, type = "pearson"),
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Pearson Residuals")
abline(h = 0, col = "red") # No pattern. A few higher than 0
qqnorm(residuals(lm_ton, type = "pearson"))
qqline(residuals(lm_ton, type = "pearson"), col = "red")#Reasonably normal
#NH4
favstats(nh4_data_bare$Result)
lm_nh4 = lm(Result ~ Site, data = nh4_data_bare) 
summary(lm_nh4) #0.86 r2
plot(fitted(lm_nh4), residuals(lm_nh4, type = "pearson"),
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Pearson Residuals")
abline(h = 0, col = "red") # No pattern. A few higher than 0
qqnorm(residuals(lm_nh4, type = "pearson"))
qqline(residuals(lm_nh4, type = "pearson"), col = "red")#Reasonably normal

lm_po4 = lm(Result ~ Site, data = po4_data_bare)  
summary(lm_po4) # Not significant p > 0.05

ggplot(nh4_data_bare, aes(x = Site, y = Result)) +
  geom_boxplot() +
  labs(x = "Site",
       y = "NH4 (mg/L)") +
  theme_minimal()
ggplot(ton_data_bare, aes(x = Site, y = Result)) +
  geom_boxplot() +
  labs(x = "Site",
       y = "TON (mg/L)") +
  theme_minimal()
ggplot(po4_data_bare, aes(x = Site, y = Result)) +
  geom_boxplot() +
  labs(x = "Site",
       y = "PO4 (mg/L)") +
  theme_minimal()


cor.test(x= results_fen$DO..mg.l., y= results_fen$TON_mg_l, method = 'spearman') #not significant
cor.test(x= results_fen$DO..mg.l., y= results_fen$NH4_mg_l, method = 'spearman') #nearly significant (0.06) 0.3 rho
cor.test(x= results_fen$DO..mg.l., y= results_fen$PO4_mg_l, method = 'spearman') #significant weak relationship (0.5 rho)
cor.test(x= results_fen$DO..mg.l., y= results_fen$TDN_mg_l, method = 'spearman') #not significant
cor.test(x= results_fen$DO..mg.l., y= results_fen$PON_mg_l, method = 'spearman') #not significant


cor.test(x= results_fen$EC..μS., y= results_fen$TON_mg_l, method = 'spearman') 
cor.test(x= results_fen$EC..μS., y= results_fen$NH4_mg_l, method = 'spearman')
cor.test(x= results_fen$EC..μS., y= results_fen$PO4_mg_l, method = 'spearman')
cor.test(x= results_fen$EC..μS., y= results_fen$TDN_mg_l, method = 'spearman')
cor.test(x= results_fen$EC..μS., y= results_fen$PON_mg_l, method = 'spearman')
#None show statistical significance 

cor.test(x= results_fen$Temperature...C., y= results_fen$TON_mg_l, method = 'spearman') 
cor.test(x= results_fen$Temperature...C., y= results_fen$NH4_mg_l, method = 'spearman')
cor.test(x= results_fen$Temperature...C., y= results_fen$PO4_mg_l, method = 'spearman')
cor.test(x= results_fen$Temperature...C., y= results_fen$TDN_mg_l, method = 'spearman')
cor.test(x= results_fen$Temperature...C., y= results_fen$PON_mg_l, method = 'spearman')
#None show statistical significance

cor.test(x= results_fen$pH, y= results_fen$TON_mg_l, method = 'spearman') 
cor.test(x= results_fen$pH, y= results_fen$NH4_mg_l, method = 'spearman')
cor.test(x= results_fen$pH, y= results_fen$PO4_mg_l, method = 'spearman')
cor.test(x= results_fen$pH, y= results_fen$TDN_mg_l, method = 'spearman')
cor.test(x= results_fen$pH, y= results_fen$PON_mg_l, method = 'spearman')
#None show statistical significance

pon_proportions = read.csv("pon_proportions.csv")
pon_proportions$site = as.factor(pon_proportions$site)
kruskal.test(proportion ~ site, data = pon_proportions) #not statistical differences. p is 0.4
is.factor(group1_proportions$site)
group1_proportions = pon_proportions[1:4, ]
group2_proportions = pon_proportions[5:9, ]
wilcox.test(group1_proportions$proportion, group2_proportions$proportion) # 0.01 p so they are statistically different
ggplot(pon_proportions, aes(x = as.factor(site), y = proportion, fill = as.factor(site))) +
  geom_col(fill ="steelblue") +  
  theme_minimal() +
  labs(x = "Site",
       y = "PON as a % TN") +
  theme(legend.position = "none")

#Pie chart for PON proportions
data_nitrogen = read.csv("nitrogen_proportions.csv")
data_long = melt(data_nitrogen, id.vars = "Site", measure.vars = c("TDN", "PON"), 
                  variable.name = "Nitrogen_Type", value.name = "Value")
data_long = data_long %>%
  group_by(Site) %>%
  mutate(Percentage = Value / sum(Value) * 100) #makes each pie full
ggplot(data_long, aes(x = "", y = Percentage, fill =  Nitrogen_Type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  facet_wrap(~ Site) +
  labs(title = expression("Nitrogen Proportions Across Sites (mg "*L^{-1}*")"), x = NULL, y = NULL, fill = "") +  
  theme_minimal() +
  theme(axis.text.x = element_blank(), 
        axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#1f78b4", "#33a02c"))

#Morans I
results_coords = read.csv("results_f_e_n.csv")
results_sf = st_as_sf(results_coords, coords = c("Longitude", "Latitude"), crs = 4326) #converted data to spatial data using long and lat
coords = st_coordinates(results_sf)
nb = knn2nb(knearneigh(coords, k = 4))
lw = nb2listw(nb, style = "W")
moran_test_ecoli <- moran.test(results_sf$Average_E_Coli_per_100ml, lw)
print(moran_test_ecoli) #alot of clustering

#Comparing E. coli with other studies plot
other_means = read.csv("means.csv")
other_means = other_means %>%
  mutate(Study = gsub("_", " ", Study))
ggplot(other_means, aes(x = Study, y = Mean)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 900, color = "red", linetype = "dashed", size = 1) +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(
    y = expression("Logged "*italic("E. coli")*" (CFU / 100 mL)"),
    x = "Study"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    panel.grid.minor = element_blank()
  )
other_means_2 = read.csv("other_means_2.csv")
other_means_2 = other_means_2 %>%
  mutate(Study = gsub("_", " ", Study))
ggplot(other_means_2, aes(x = Study, y = Mean)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 900, color = "red", linetype = "dashed", size = 1) +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(
    y = expression("Logged "*italic("E. coli")*" (CFU / 100 mL)"),
    x = "Study"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    panel.grid.minor = element_blank()
  )


cor.test(x= results_fen$NH4_mg_l, y= results_fen$TON_mg_l, method = 'spearman')
cor.test(x= results_fen$NH4_mg_l, y= results_fen$PON_mg_l, method = 'spearman')
cor.test(x= results_fen$NH4_mg_l, y= results_fen$TDN_mg_l, method = 'spearman')
cor.test(x= results_fen$PON_mg_l, y= results_fen$TON_mg_l, method = 'spearman')
cor.test(x= results_fen$PO4_mg_l, y= results_fen$TDN_mg_l, method = 'spearman')
cor.test(x= results_fen$PO4_mg_l, y= results_fen$PON_mg_l, method = 'spearman')
cor.test(x= results_fen$PO4_mg_l, y= results_fen$NH4_mg_l, method = 'spearman')
cor.test(x= results_fen$PO4_mg_l, y= results_fen$TON_mg_l, method = 'spearman')
#All siginificant, strongest between NH4 and TON, then PON and TON, then slight positive correlation between PON and NH4

#PO4
cor.test(x= results_fen$NH4_ug_l, y= results_fen$PO4_ug_l, method = 'spearman') #No correlation
wilcox.test(results_pre_basin$PO4_ug_l, results_post_basin$PO4_ug_l) #No differneces
po4_means = read.csv("po4_means.csv")
#No changes in PO4
ggplot(po4_means, aes(x = factor(Site), y = PO4_mg_l)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 0.1, color = "red", linetype = "dashed", size = 1) +
  labs(x = "Site",
       y = expression(PO[4]^"3-"* " (mg "*L^{-1}*")")) +
  theme_minimal()

ton_means = read.csv("ton_means.csv")
ggplot(ton_means, aes(x = factor(Site), y = TON_mg_l)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 5.6, color = "red", linetype = "dashed", size = 1) +
  labs(x = "Site",
       y = expression("TON (mg "*L^{-1}*")")) +
  theme_minimal()

nh4_means = read.csv("nh4_means.csv")
ggplot(nh4_means, aes(x = factor(Site), y = NH4_mg_l)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 0.2, color = "red", linetype = "dashed", size = 1) +
  labs(x = "Site",
       y = expression(NH[4] ^"+"* " (mg "*L^{-1}*")")) +
  theme_minimal()

do_means = read.csv("do_means.csv")
ggplot(do_means, aes(x = factor(Site), y = DO)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    x = "Site",
    y = expression("DO (mg "*L^{-1}*")")
  ) +
  theme_minimal()

ec_means = read.csv("ec_means.csv") 
ggplot(ec_means, aes(x = factor(Site), y = EC)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Site", y = expression("Electrical Conductivity ("*mu*"S cm"^-1*")")+
         theme_minimal()

       ggplot(ec_means, aes(x = factor(Site), y = EC)) +
         geom_bar(stat = "identity", fill = "steelblue") +
         labs(x = "Site", 
              y = expression("EC ("*mu*"S cm"^-1*")")) +
         theme_minimal()
  
