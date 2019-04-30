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
library(lme4)

# Change font size and ggthemes globally
theme_set(ggthemes::theme_few(base_size = 8))

# Change globally how numbers are displayed
options("scipen"=100, "digits"=4)

# Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  

# Import some helpful functions
source("HighstatLibV10.R")
setwd("..")

#' ### Import and prepare the data set
fulldf <- read.csv("Data/Niger_snail_survey_Bulinus_comb_2019-04-28.csv") %>% 
  dplyr::select(filter.min, coll_date, month, locality, site_irn, BP_tot, BP_pos_tot, BF_tot, BF_pos_tot, 
                BT_tot, BT_pos_tot, BT_prev, BF_prev, Bulinus_tot, Bulinus_pos_tot,
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
subdf <- fulldf %>% 
      filter(!is.na(duration))
################################################################################################################################

#' Overview over the data structure
plot_str(subdf) + theme_few()
introduce(subdf)

#' Overview over predictor variables
plot_histogram(subdf)  # Outliers in Cond, PPM, water_depth, wmo_prec, water_speed. These outliers will drive results. Look at them critically

#' Look closer at the outliers

# 1: Cond
plot(subdf$Cond, subdf$Bulinus_tot)

# Where do these high values occur?
subdf$site_irn[subdf$Cond > 6]   
View(subdf[subdf$site_irn == 382877,] %>% dplyr::select(Cond, site_irn, water_depth, pH, PPM, everything()))
View(subdf[subdf$site_irn == 487402,] %>% dplyr::select(Cond, site_irn, water_depth, pH, PPM, everything()))

# High Cond seem to coincide with high PPM, so it's not very likely that these are measurement outliers

# 2: water_speed_ms
plot(subdf$water_speed_ms, subdf$Bulinus_tot)  # One strong outlier

# Where do these high values occur?
subdf$site_irn[subdf$water_speed_ms > 3]   
View(subdf[subdf$site_irn == 382867,] %>% dplyr::select(water_speed_ms, Cond, site_irn, water_depth, pH, PPM, everything()))

# Is it generally an outlier, looking at all site types?
plot(subdf$site_type, subdf$water_speed_ms)   # General outlier

# The measurement of 5 is possible, but a general outlier. It was the single measurement higher than 2, throughout
# the whole study. Even though it is most likely a real measurement, the complete lack of other measurements
# in that range means that it is hard to say that this is an actual representation of counts at high water speeds,
# yet the measurement looks like it's a highly influential outliers, due to its position and high leverage.
# It might be best to exclude this observation from the data.

# 2: wmo_prec
plot(subdf$wmo_prec, subdf$Bulinus_tot)  # Many observations at most extreme precipitation

# Where do these high values occur?
View(subdf[subdf$wmo_prec > 60,] %>% dplyr::select(wmo_prec, seas_wmo, site_irn, water_depth, pH, PPM, everything()))  

# All outliers here come from the same, very rainy day at Lata Kabia (2014-08-02)
ggplot(subdf[subdf$locality == "Lata Kabia",]) +
      geom_point(aes(x = as.factor(month), y = wmo_prec))

#' Make a data frame with outliers removed
subdf_out <- subdf %>% 
      dplyr::filter(water_speed_ms < 2,
                    Cond < 8)

#' Are there any variance inflation factors (multicollinearity)? Check using a function from Zuur et al. 2010

pairs(subdf[,c("Temp_Air", "Temp_Water", "water_speed_ms", "water_depth", "pH", "Cond", "PPM",
                 "wmo_av_temp", "Bulinus_tot", "wmo_prec")],
      lower.panel = panel.cor)

corvif(data.table::as.data.table(subdf)[, c("Temp_Air", "Temp_Water", "water_speed_ms", "water_depth", "pH", "Cond", "PPM",
                                  "wmo_av_temp", "wmo_prec"), with=FALSE])

# Cond and PPM have high GVIF values (10). For values of higher than 4, only one of the two variables should be used in models, 
# to avoid multicollinearity. 


#' Look at general correlation matrix
plot_correlation(na.omit(subdfc), maxcat = 5L)   


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
Bulinus_poiss <- glmmTMB(Bulinus_tot ~ (1|locality/site_irn/visit_no)  + 
                       Temp_Water + pH + water_speed_ms + water_depth + Cond + wmo_prec +
                       locality + site_type + Bulinus_pos_tot + 
                       site_type*Temp_Water + site_type*pH + site_type*Cond + site_type*wmo_prec + 
                       offset(log(duration)),
                 data=subdf,
                 family=poisson)

Bulinus_m <- glmmTMB(Bulinus_tot ~ (1|locality/site_irn/visit_no) + 
                       Temp_Water + pH + water_speed_ms + water_depth + Cond + wmo_prec +
                       locality + site_type + Bulinus_pos_tot + 
                       site_type*Temp_Water + site_type*pH + site_type*Cond + site_type*wmo_prec + 
                       offset(log(duration)),
                 data=subdf,
                 family=nbinom1)

Bulinus_m1 <- glmmTMB(Bulinus_tot ~ (1|locality/site_irn/visit_no) + 
                        Temp_Water + pH + water_speed_ms + water_depth + Cond + wmo_prec +
                        locality + site_type + Bulinus_pos_tot + 
                        site_type*Temp_Water + site_type*pH + site_type*Cond + site_type*wmo_prec + 
                        offset(log(duration)),
                  data=subdf,
                  family=nbinom2)

#' Model diagnostics
# Visually check the model fit using the DHARMa package, a package for model diagnostics
# https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html#formal-goodness-of-fit-tests-on-the-scaled-residuals
sim_residualsPoiss <- DHARMa::simulateResiduals(Bulinus_poiss, 1000)  # Ignore warnings
# Plot the residuals to visually test for over-/underdispersion
plot(sim_residualsPoiss)  # significant deviation

sim_residuals1 <- DHARMa::simulateResiduals(Bulinus_m, 1000)  # Ignore warnings
# Plot the residuals to visually test for over-/underdispersion
plot(sim_residuals1)  # Good fit

sim_residuals2 <- DHARMa::simulateResiduals(Bulinus_m1, 1000) 
plot(sim_residuals2)  # Virtually the same. Difference between nbinom1 and nbinom2 seems neglible

testZeroInflation(sim_residuals2)  # There is no evidence for zero-inflation

# Model results
summary(Bulinus_m1)  # Next to no variance is coming from locality, so technically we could exclude this from the model.
                 # Aesthetically, we can keep it in, to highlight the nested design of the study
Anova.glmmTMB(Bulinus_m1)




#########################################################################################
########################      Bulinus truncatus total      ##############################  
#########################################################################################
BT_m1 <- glmmTMB(BT_tot ~ (1|locality/site_irn/visit_no)  + wmo_prec +
                        site_type + month + BT_pos_tot + 
                        BP_tot + BF_tot + L_tot +
                        offset(log(duration)),
                  data=subdf,
                  family=nbinom2)

sim_res_BT_m1 <- DHARMa::simulateResiduals(BT_m1, 1000)
plot(sim_res_BT_m1)  # Not overdispersed.
testZeroInflation(sim_res_BT_m1)  # No zero-inflation. 

# Model results
summary(BT_m1)  
Anova.glmmTMB(BT_m1)




#########################################################################################
########################      Bulinus forskalii total      ##############################
#########################################################################################
BF_m1 <- glmmTMB(BF_tot ~ (1|locality/site_irn) + (1|coll_date) + wmo_prec +
              site_type + month + BF_pos_tot + 
              BP_tot + BT_tot + L_tot +
              offset(log(duration)),
        data=subdf,
        family=nbinom2)

sim_res_BF_m1 <- DHARMa::simulateResiduals(BF_m1, 1000)
plot(sim_res_BF_m1)  # Not overdispersed.
testZeroInflation(sim_res_BF_m1)  # No zero-inflation. 

# Model results
summary(BF_m1)
Anova.glmmTMB(BF_m1)




#########################################################################################
##############################      Limnea total      ###################################
#########################################################################################
L_m1 <- glmmTMB(L_tot ~ (1|locality/site_irn) + (1|coll_date) + wmo_prec +
                       site_type + month + 
                       BP_tot + BT_tot + BF_tot +
                       offset(log(duration)),
                 data= subdf, 
                 family=nbinom2)

sim_res_L_m1 <- DHARMa::simulateResiduals(L_m1, 1000)
plot(sim_res_L_m1)  # Not overdispersed.
testZeroInflation(sim_res_L_m1)  # No zero-inflation. 

# Model results
summary(L_m1)
Anova.glmmTMB(L_m1)


#########################################################################################
###############################      Prevalence      ####################################
#########################################################################################
df_monthly <- fulldf %>% group_by(month) %>% 
      mutate(av_prec = mean(wmo_prec)) %>% ungroup() %>% 
      group_by(month, locality, site_type, seas_wmo) %>% 
      summarize(BT_pos_tot = as.numeric(sum(BT_pos_tot)),
                BT_tot = as.numeric(sum(BT_tot)),
                BP_pos_tot = as.numeric(sum(BP_pos_tot)),
                BP_tot = as.numeric(sum(BP_tot)),
                BF_pos_tot = as.numeric(sum(BF_pos_tot)),
                BF_tot = as.numeric(sum(BF_tot)),
                av_prec = mean(av_prec))


# test <- df_monthly %>% 
#       gather("species", "value", -c(month, locality, site_type, seas_wmo, av_prec)) %>% 
#       separate(col = "species", into = c("species", "var"), sep = "_") %>% 
#       spread(key = var, value = value)

BT_prev_m1 <- glm(BT_pos_tot/BT_tot ~ site_type + locality + month,
                      weights = BT_tot,
                      data= df_monthly[df_monthly$locality != "Gantchi Bassarou" & df_monthly$locality != "Tagabati" & 
                                             df_monthly$locality != "Tiaguirire" & df_monthly$locality != "Yoreize Koira" & 
                                             df_monthly$site_type != "stream",],
                      family= binomial)

sim_residuals_BT_prev <- DHARMa::simulateResiduals(BT_prev_m1, 1000)  
plot(sim_residuals_BT_prev) 
DHARMa::testDispersion(sim_residuals_BT_prev)

car::Anova(BT_prev_m1)
summary(BT_prev_m1)

BP_prev_m1 <- glm(BP_pos_tot/BP_tot ~ site_type + locality + month,
                          weights = BP_tot,
                          data=subset(df_monthly, locality %in% c("Namari Goungou", "Diambala")),
                          family= binomial)

sim_residuals_BP_prev <- DHARMa::simulateResiduals(BP_prev_m1, 1000)  
plot(sim_residuals_BP_prev) 
DHARMa::testDispersion(sim_residuals_BP_prev)

car::Anova(BP_prev_m1)
summary(BP_prev_m1)


#########################################################################################
##############################      Anova tables      ###################################
#########################################################################################
glmmTMB::Anova.glmmTMB(Bulinus_m1)
summary(Bulinus_m1)
glmmTMB::Anova.glmmTMB(BT_m1)   
summary(BT_m1)
glmmTMB::Anova.glmmTMB(BF_m1)  
summary(BF_m1)
glmmTMB::Anova.glmmTMB(L_m1)   
summary(L_m1)
car::Anova(BT_prev_m1)
summary(BT_prev_m1)
car::Anova(BP_prev_m1)
summary(BP_prev_m1)

#########################################################################################
############################      Plot of emmeans      ##################################
#########################################################################################

#' Temp_Water:site_type
Temp_Water_site_type_df <-  data.frame(emtrends(Bulinus_m1,  ~ site_type, var = "Temp_Water", type = "response"))
ggplot(Temp_Water_site_type_df) +
      geom_errorbar(aes(x = reorder(site_type, Temp_Water.trend), ymin = Temp_Water.trend-SE, ymax = Temp_Water.trend+SE),
                    width = .3) +
      geom_point(aes(x = reorder(site_type, Temp_Water.trend), y = Temp_Water.trend), pch = 23, cex = 4,
                 fill = "white") +
      ylab("Estimated slope with Temp_Water") +
      xlab("")

ggplot(Temp_Water_site_type_df) +
      geom_bar(aes(x = reorder(site_type, Temp_Water.trend), y = Temp_Water.trend), col = "black",
               fill = "white", position = position_dodge(), stat = "identity") +
      geom_errorbar(aes(x = reorder(site_type, Temp_Water.trend), ymin = Temp_Water.trend-SE, ymax = Temp_Water.trend+SE),
                    width = .3) +
      ylab("Estimated slope with Temp_Water") +
      xlab("")

#' Cond:site_type
Cond_site_type_df <-  data.frame(emtrends(Bulinus_m1,  ~ site_type, var = "Cond", type = "response"))
ggplot(Cond_site_type_df) +
      geom_errorbar(aes(x = reorder(site_type, Cond.trend), ymin = Cond.trend-SE, ymax = Cond.trend+SE),
                    width = .3) +
      geom_point(aes(x = reorder(site_type, Cond.trend), y = Cond.trend), pch = 23, cex = 4,
                 fill = "white") +
      ylab("Estimated slope with Cond") +
      xlab("")

#' Bulinus_pos_tot
emmip(Bulinus_m1, ~ Bulinus_pos_tot, cov.reduce = range)

#' locality
locality_df <-  data.frame(emmeans(Bulinus_m1,  ~ locality, type = "response"))
ggplot(locality_df) +
      geom_bar(aes(x = reorder(locality, rate), y = rate), col = "black",
               fill = "white", position = position_dodge(), stat = "identity") +
      geom_errorbar(aes(x = reorder(locality, rate), ymin = rate-SE, ymax = rate+SE),
                    width = .3) +
      ylab("Estimated marginal mean") +
      xlab("")

#' water_speed_ms
emmip(Bulinus_m1, ~ water_speed_ms, cov.reduce = range)


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

test  <- data.frame(emmeans(BT_prev_m1,  ~ month, type = "response"))
ggplot(test) +
      geom_bar(aes(x = month, y = prob), col = "black", 
               position = position_dodge(), stat = "identity") +
      geom_errorbar(aes(x = month, ymin = prob-SE, ymax = prob+SE), position = position_dodge()) 

#' Prevalence site_type
prev_site_type_BT <-  data.frame(emmeans(BT_prev_m1,  ~ site_type, type = "response")) %>% mutate(species = "Bulinus truncatus")
prev_site_type_BP <-  data.frame(emmeans(BP_prev_m1,  ~ site_type, type = "response")) %>% mutate(species = "Biomphilaria pfeifferi")
prev_site_type_all <- rbind(prev_site_type_BT, prev_site_type_BP)

ggplot(prev_site_type_all) +
      geom_bar(aes(x = site_type, y = prob, fill = species), col = "black", 
               position = position_dodge(), stat = "identity") +
      geom_errorbar(aes(x = site_type, ymin = prob-SE, ymax = prob+SE, group = species), position = position_dodge()) 

#' Prevalence month
prev_month_BT <-  data.frame(emmeans(BT_prev_m1,  ~ month, type = "response")) %>% mutate(species = "Bulinus truncatus")
prev_month_BP <-  data.frame(emmeans(BP_prev_m1,  ~ month, type = "response")) %>% mutate(species = "Biomphilaria pfeifferi")
prev_month_all <- rbind(prev_month_BT, prev_month_BP)

ggplot(prev_month_all) +
      geom_bar(aes(x = month, y = prob, fill = species), col = "black", 
               position = position_dodge(), stat = "identity") +
      geom_errorbar(aes(x = month, ymin = prob-SE, ymax = prob+SE, group = species), position = position_dodge()) 

#' Prevalence locality
prev_locality_BT <-  data.frame(emmeans(BT_prev_m1,  ~ locality, type = "response")) %>% mutate(species = "Bulinus truncatus")
prev_locality_BP <-  data.frame(emmeans(BP_prev_m1,  ~ locality, type = "response")) %>% mutate(species = "Biomphilaria pfeifferi")
prev_locality_all <- rbind(prev_locality_BT, prev_locality_BP)

tiff(file="Figures/prev_locality.tiff", width = 170, height = 150, units = "mm", res = 200)
ggplot(prev_locality_BT) +
      geom_bar(aes(x = reorder(locality, prob), y = prob), col = "black", fill = "lightgrey",
               position = position_dodge(), stat = "identity") +
      geom_errorbar(aes(x = locality, ymin = prob-SE, ymax = prob+SE), position = position_dodge()) +
      ggtitle("Bulinus truncatus")
dev.off()
