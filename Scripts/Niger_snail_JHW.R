#' ### Prepare the environment

# Remove all items in environment, for a clean start
rm(list = ls(all=TRUE))  

# Lod the necessary packages
library(tidyverse)
library(lubridate)
library(glmmTMB)
library(ggthemes)
library(DHARMa) # for model diagnostics
library(emmeans)  # for post hoc test
library(DataExplorer)   # for data exploration
library(gridExtra)

# Change font size and ggthemes globally
theme_set(ggthemes::theme_few(base_size = 8))

# Change globally how numbers are displayed
options("scipen"=100, "digits"=4)

# Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  

# Import some helpful functions
source("HighstatLibV10.R")
setwd("..")
test <- fulldf %>% group_by(site_type) %>% 
      summarize(mT = mean(Temp_Water, na.rm = T ), minT = min(Temp_Water, na.rm = T ), maxT = max(Temp_Water, na.rm = T ))

#' ### Import and prepare the data set
fulldf <- read.csv("Data/Niger_snail_survey_cref_FTA_counts_2019.csv") %>% 
  dplyr::select(filter.min, coll_date, month, locality, site_irn, BP_tot, BP_pos_tot, BF_tot, BF_pos_tot, 
                BT_tot, BT_pos_tot, BG_tot, BS_tot, BT_prev, BF_prev, Bulinus_tot, Bulinus_pos_tot,
                bp_pres, bt_pres, bf_pres, visit_no, site_type, L_tot, 
                Temp_Air,Temp_Water, water_speed_ms, water_depth, pH, Cond, PPM, Latitude.v, Longitude.v, Stagnante,
                wmo_min_temp, wmo_max_temp, wmo_av_temp, wmo_prec, seas_wmo, duration, Heure,
                water_level.v, water_level3) %>% 
  # Remove NA values in predictor variables
  dplyr::filter(Bulinus_tot <= 1000,
                !is.na(site_irn)) %>% 
  mutate(coll_date = lubridate::dmy(coll_date), tz = "Africa/Niger",
         duration = as.numeric(as.character(duration)),
         month = as.factor(month),
         year = lubridate::year(coll_date),
         site_irn = as.factor(site_irn), 
         visit_no = as.factor(visit_no),
         bp_pres = as.factor(as.character(bp_pres)),
         bf_pres = as.factor(as.character(bf_pres)),
         Heure = as.factor(as.character(Heure))) %>% 
  # Re-arrange the columns
  dplyr::select(locality, site_irn, visit_no, Bulinus_tot, Bulinus_pos_tot, coll_date, everything()) %>% 
  arrange(locality, site_irn) %>% 
  # This site has no measurements in the wet season, and a lot in the dry season. Later in the models, this leads to 
  # very high variance, when we're looking at interactions containing the season variable.
  # For now, remove this site from the analysis.
  filter(filter.min != "y") %>% 
      # Rescale variables with large values, change water depth unit to meters
      mutate(Cond = scale(Cond, center = 0),
             PPM = scale(PPM, center = 0),
             wmo_prec = scale(wmo_prec, center = 0),
             Temp_Water = scale(Temp_Water, center = 0),
             water_speed_ms = scale(water_speed_ms, center = 0),
             water_depth = water_depth/100,
             water_depth = scale(water_depth, center = 0),
             pH = scale(pH, center = 0)) %>% 
      group_by(site_irn) %>% 
      mutate(Latitude = mean(Latitude.v, na.rm = TRUE),
             Longitude = mean(Longitude.v, na.rm = TRUE))
fulldf$Heure[fulldf$Heure == ""] <- NA

#' Overview over missing values
plot_missing(fulldf) + theme_few()
# For 33.5% of the data, information on the duration of sampling is missing. How is sampling duration distributed?
# To decide how to deal with the missing values for duration, we need to know if the values are missing at complete random 
# (in that case, removing them would not be an issue)

# Is the duration more likely to be missing, given certain other values?
fulldf$dur_missing <- ifelse(is.na(fulldf$duration), 1, 0)
miss_model <- glm(dur_missing ~ Temp_Air + BP_tot + BT_tot + BF_tot +
                        site_type, 
                  family = binomial,
                  data = fulldf)
summary(miss_model) # Not really any strong effects



################################################################################################################################
# Remove NA values in duration, as well as predictors
fulldf <- fulldf %>% 
      filter(Bulinus_tot < 1000, !is.na(duration))
################################################################################################################################

#' Overview over the data structure
plot_str(fulldf) + theme_few()
introduce(fulldf)

#' Overview over predictor variables
plot_histogram(fulldf)  # Outliers in Cond, PPM, water_depth, wmo_prec, water_speed. These outliers will drive results. Look at them critically

#' Look closer at the outliers

# 1: Cond
plot(fulldf$Cond, fulldf$Bulinus_tot)

# Where do these high values occur?
fulldf$site_irn[fulldf$Cond > 6]   
View(fulldf[fulldf$site_irn == 382877,] %>% dplyr::select(Cond, site_irn, water_depth, pH, PPM, everything()))
View(fulldf[fulldf$site_irn == 487402,] %>% dplyr::select(Cond, site_irn, water_depth, pH, PPM, everything()))

# High Cond seem to coincide with high PPM, so it's not very likely that these are measurement outliers

# 2: water_speed_ms
plot(fulldf$water_speed_ms, fulldf$Bulinus_tot)  # One strong outlier

# Where do these high values occur?
fulldf$site_irn[fulldf$water_speed_ms > 3]   
View(fulldf[fulldf$site_irn == 382867,] %>% dplyr::select(water_speed_ms, Cond, site_irn, water_depth, pH, PPM, everything()))

# Is it generally an outlier, looking at all site types?
plot(fulldf$site_type, fulldf$water_speed_ms)   # General outlier

# The measurement of 5 is possible, but a general outlier. It was the single measurement higher than 2, throughout
# the whole study. Even though it is most likely a real measurement, the complete lack of other measurements
# in that range means that it is hard to say that this is an actual representation of counts at high water speeds,
# yet the measurement looks like it's a highly influential outliers, due to its position and high leverage.
# It might be best to exclude this observation from the data.

# 2: wmo_prec
plot(fulldf$wmo_prec, fulldf$Bulinus_tot)  # Many observations at most extreme precipitation

# Where do these high values occur?
View(fulldf[fulldf$wmo_prec > 60,] %>% dplyr::select(wmo_prec, seas_wmo, site_irn, water_depth, pH, PPM, everything()))  

# All outliers here come from the same, very rainy day at Lata Kabia (2014-08-02)
ggplot(fulldf[fulldf$locality == "Lata Kabia",]) +
      geom_point(aes(x = as.factor(month), y = wmo_prec))

#' Make a data frame with outliers removed
fulldf_out <- fulldf %>% 
      dplyr::filter(water_speed_ms < 2,
                    Cond < 8)

# Make data frame with all missing values remove in water variables
chemdf <- fulldf %>% 
      dplyr::filter(!is.na(pH), !is.na(Cond), !is.na(Temp_Water), !is.na(PPM), !is.na(pH),
                    !is.na(water_depth), !is.na(water_speed_ms), water_level.v != "")
chemdf_out <- chemdf %>% 
      dplyr::filter(water_speed_ms < 2,
                    Cond < 8)
#' Are there any variance inflation factors (multicollinearity)? Check using a function from Zuur et al. 2010

pairs(fulldf[,c("Temp_Air", "Temp_Water", "water_speed_ms", "water_depth", "pH", "Cond", "PPM",
                 "wmo_av_temp", "Bulinus_tot", "wmo_prec")],
      lower.panel = panel.cor)

corvif(data.table::as.data.table(fulldf)[, c("Temp_Air", "Temp_Water", "water_speed_ms", "water_depth", "pH", "Cond", "PPM",
                                  "wmo_av_temp", "wmo_prec"), with=FALSE])

# Cond and PPM have high GVIF values (10). For values of higher than 4, only one of the two variables should be used in models, 
# to avoid multicollinearity. 


#' Look at general correlation matrix
plot_correlation(na.omit(fulldfc), maxcat = 5L)   


# Study design:
# I have count data of snails per date, counted over many dates at sites, nested in localities.
# So, in each locality the snail counts come from several different sites, repeatedly sampled on different dates.

# Goal:
# Test if snail counts differ between localities, and test influence of environmental factors (e.g. water pH)

# Things to account for:
# A: Sites in localities might show variation in intercepts due to higher initial snail abundance
# B: Sampling duration differed (5-33 minutes), which will most likely influence counts

# How to account for it:
# A: Include site as a random intercept, to account for variation in counts between the sites
# B: Include sampling duration as an offset, to account for differences in sampling effort


#'<br><br>
#'
#'--------------------------------------------------------------------------------------------
#'--------------------------------------------------------------------------------------------
#'
#'
#'--------------------------------------------------------------------------------------------
#'--------------------------------------------------------------------------------------------
#'
#'
#########################################################################################
############################      Bulinus total      ####################################  To investigate water attribute effects specifically
#########################################################################################

# There are generally many zeroes. This doesn't mean that the data is zero-inflated, but it might be worth checking for zero-inflation anyways

#' ### Make a GLMM

#' Make a maximum model. glmmTMB is a new package by Ben Bolker, that fits models faster, and allows to include arguments
# to account for zero-inflation, if needed.
# Include the sampling duration as an offset. It needs to be specified as log(), since we are using a family distribution
# with log link (nbinom2)

# Decide for a family. It's count count data, looking very overdispersed, so negative binomial is most likely appropriate
poiss <- glmmTMB(Bulinus_tot ~ (1|locality/site_irn/visit_no)  + 
                       locality + pH + water_speed_ms + water_depth + water_level.v + Cond + wmo_prec +
                       Temp_Water + site_type + seas_wmo + 
                       locality*seas_wmo + site_type*seas_wmo + site_type*wmo_prec + 
                       site_type*Cond + site_type*Temp_Water + site_type*pH  +
                   offset(log(duration)),
                 data=chemdf,
                 family=poisson)

model <- glmmTMB(Bulinus_tot ~ (1|locality/site_irn/visit_no) + 
                       locality + pH + water_speed_ms + water_depth + water_level.v + Cond + wmo_prec +
                       Temp_Water + site_type + seas_wmo + 
                       locality*seas_wmo + site_type*seas_wmo + site_type*wmo_prec + 
                       site_type*Cond + site_type*Temp_Water + site_type*pH +
                      offset(log(duration)),
                 data=chemdf,
                 family=nbinom1)

model1 <- glmmTMB(Bulinus_tot ~ (1|locality/site_irn) + (1|coll_date) +
                        Temp_Water + pH + water_speed_ms + water_depth + Cond + wmo_prec +
                        site_type + Bulinus_pos_tot + 
                        site_type*Temp_Water + site_type*pH + site_type*Cond + site_type*wmo_prec + 
                        offset(log(duration)),
                  data=chemdf,
                  family=nbinom2)
Anova.glmmTMB(model1)
#' Model diagnostics
# Visually check the model fit using the DHARMa package, a package for model diagnostics
# https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html#formal-goodness-of-fit-tests-on-the-scaled-residuals
sim_residuals1 <- DHARMa::simulateResiduals(model, 1000)  # Ignore warnings
# Plot the residuals to visually test for over-/underdispersion
plot(sim_residuals1)  # Good fit

sim_residuals2 <- DHARMa::simulateResiduals(model1, 1000) 
plot(sim_residuals2)  # Virtually the same. Difference between nbinom1 and nbinom2 seems neglible

testZeroInflation(sim_residuals2)  # There is no evidence for zero-inflation

# Look at the model

summary(model1)  # Next to no variance is coming from locality, so technically we could exclude this from the model.
                 # Aesthetically, we can keep it in, to highlight the nested design of the study
                 
# Standard errors for Kohan Garantche:seas_wmowet are huge. Try to find out why
View(chemdf[chemdf$locality == "Kohan Garantche",] %>% dplyr::select(seas_wmo, locality, Bulinus_tot, everything())) 
# There were many snails in the dry season, but zero in the wet season. This inflates the standard error hugely.
# One way to potentially deal with this is to artificially add a count of 1 for one day of the wet season
# Put count of 1 for visit_no 40
bulinus_df <- chemdf
bulinus_df$Bulinus_tot[bulinus_df$locality == "Kohan Garantche" & bulinus_df$visit_no == 40] <- 1
model2 <- glmmTMB(Bulinus_tot ~ (1 |locality/site_irn/visit_no) + 
                        locality + pH + water_speed_ms + water_depth + Cond + wmo_prec +
                        Temp_Water + site_type + month + Bulinus_pos_tot + 
                        site_type*wmo_prec + 
                        site_type*Cond + site_type*Temp_Water + site_type*pH +
                        offset(log(duration)),
                  data=bulinus_df,
                  family=nbinom2)
summary(model2)  # Standard errors are now normal


#'<br><br>
#'
#'--------------------------------------------------------------------------------------------
#'--------------------------------------------------------------------------------------------
#'
#'
#'--------------------------------------------------------------------------------------------
#'--------------------------------------------------------------------------------------------
#'
#'
#########################################################################################
########################      Bulinus truncatus total      ##############################  To look at species site preferences and seasonal trends
#########################################################################################
mod_1 <- glmmTMB(BT_tot ~ (1|locality/site_irn) + (1|coll_date) + wmo_prec +
                        site_type + month + BT_pos_tot + 
                        BP_tot + BF_tot + L_tot +
                        offset(log(duration)),
                  data=fulldf,
                  family=nbinom2)

Anova.glmmTMB(mod_1)
sim_res_mod1 <- DHARMa::simulateResiduals(mod_1, 1000)
plot(sim_res_mod1)  # Not overdispersed.
testZeroInflation(sim_res_mod1)  # No zero-inflation. 
summary(mod_1)  
View(chemdf %>% dplyr::select(BT_tot, site_type, wmo_prec, Temp_Water, everything()))

#'<br><br>
#'
#'--------------------------------------------------------------------------------------------
#'--------------------------------------------------------------------------------------------
#'
#'
#'--------------------------------------------------------------------------------------------
#'--------------------------------------------------------------------------------------------
#'
#'
#########################################################################################
########################      Bulinus forskalii total      ##############################
#########################################################################################
mod_A <- glmmTMB(BF_tot ~ (1|locality/site_irn) + (1|coll_date) + wmo_prec +
              site_type + month + BF_pos_tot + 
              BP_tot + BT_tot + L_tot +
              offset(log(duration)),
        data=fulldf,
        family=nbinom2)

sim_res_modA <- DHARMa::simulateResiduals(mod_A, 1000)
plot(sim_res_modA)  # Not overdispersed.
testZeroInflation(sim_res_modA)  # No zero-inflation. 

summary(mod_A)
Anova.glmmTMB(mod_A)

#'<br><br>
#'
#'--------------------------------------------------------------------------------------------
#'--------------------------------------------------------------------------------------------
#'
#'
#'--------------------------------------------------------------------------------------------
#'--------------------------------------------------------------------------------------------
#'
#'
#########################################################################################
##############################      Limnea total      ###################################
#########################################################################################
mod_a <- glmmTMB(L_tot ~ (1|locality/site_irn) + (1|coll_date) + wmo_prec +
                       site_type + month + 
                       BP_tot + BT_tot + BF_tot +
                       offset(log(duration)),
                 data= fulldf, 
                 family=nbinom2)
summary(mod_a)
glmmTMB::Anova.glmmTMB(mod_a)
sim_res_moda <- DHARMa::simulateResiduals(mod_a, 1000)
plot(sim_res_moda)  # Not overdispersed.
testZeroInflation(sim_res_moda)  # No zero-inflation. 

View(chemdf[chemdf$site_type == "spillway",] %>% dplyr::select(L_tot, Temp_Water, site_type, everything()))
#########################################################################################
##############################      Anova tables      ###################################
#########################################################################################
glmmTMB::Anova.glmmTMB(model1)
summary(model1)
glmmTMB::Anova.glmmTMB(mod_1)   # BT
summary(mod_1)
glmmTMB::Anova.glmmTMB(mod_A)   # BF
summary(mod_A)
glmmTMB::Anova.glmmTMB(mod_a)   # L


#########################################################################################
############################      Plot of emmeans      ##################################
#########################################################################################

#' Temp_Water:site_type

Temp_Water_site_type_df <-  data.frame(emtrends(mod_A,  ~ site_type, var = "Temp_Water",
                                                lmer.df = "tukey", type = "response"))
ggplot(Temp_Water_site_type_df) +
      geom_errorbar(aes(x = reorder(site_type, Temp_Water.trend), ymin = Temp_Water.trend-SE, ymax = Temp_Water.trend+SE),
                    width = .3) +
      geom_point(aes(x = reorder(site_type, Temp_Water.trend), y = Temp_Water.trend), pch = 23, cex = 4,
                 fill = "white") +
      ylab("Estimated slope with Temp_Water") +
      xlab("")


# Month
month_BT  <- data.frame(emmeans(mod_1,  ~ month, type = "response")) %>% mutate(species = "Bulinus truncatus")
month_BF  <- data.frame(emmeans(mod_A,  ~ month, type = "response")) %>% mutate(species = "Bulinus forskalii")
month_L  <- data.frame(emmeans(mod_a,  ~ month, type = "response")) %>% mutate(species = "Limnea")
month_all <- rbind(month_BT, month_BF, month_L)

ggplot(month_all) +
      geom_bar(aes(x = month, y = rate), col = "black", 
               position = position_dodge(), stat = "identity") +
      geom_errorbar(aes(x = month, ymin = rate-SE, ymax = rate+SE, group = species), position = position_dodge()) +
      facet_grid(rows = vars(species), scales = "free")


pairs(emmeans(mod_1,  ~ month, type = "response"))


