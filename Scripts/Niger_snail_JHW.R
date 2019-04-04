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
fulldf <- read.csv("Data/Niger_snail_survey_ALL_2019.csv") %>% 
  dplyr::select(filter.min, coll_date, month, locality, site_irn, BP_tot, BP_pos_tot, BF_tot, BF_pos_tot, 
                BT_tot, BT_pos_tot, BG_tot, BS_tot, BT_prev, BF_prev, Bulinus_tot, Bulinus_pos_tot,
                bp_pres, bt_pres, bf_pres, visit_no, site_type,
                Temp_Air,Temp_Water, water_speed_ms, water_depth, pH, Cond, PPM, Latitude, Longitude, Stagnante,
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
             water_depth = water_depth/100) %>% 
      group_by(site_irn) %>% 
      mutate(Latitude = mean(Latitude, na.rm = TRUE),
             Longitude = mean(Longitude, na.rm = TRUE))
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
############################      Bulinus total      ####################################
#########################################################################################

# There are generally many zeroes. This doesn't mean that the data is zero-inflated, but it might be worth checking for zero-inflation anyways

#' ### Make a GLMM

#' Make a maximum model. glmmTMB is a new package by Ben Bolker, that fits models faster, and allows to include arguments
# to account for zero-inflation, if needed.
# Include the sampling duration as an offset. It needs to be specified as log(), since we are using a family distribution
# with log link (nbinom2)

# Decide for a family. It's count count data, looking very overdispersed, so negative binomial is most likely appropriate
poiss <- glmmTMB(Bulinus_tot ~ (1|locality/site_irn/visit_no)  + 
                       locality + pH + water_speed_ms + water_depth + water_level + Cond + wmo_prec +
                       Temp_Water + site_type + seas_wmo + 
                       locality*seas_wmo + site_type*seas_wmo + site_type*wmo_prec + 
                       site_type*Cond + site_type*Temp_Water + site_type*pH  +
                   offset(log(duration)),
                 data=chemdf,
                 family=poisson)

model <- glmmTMB(Bulinus_tot ~ (1|locality/site_irn/visit_no) + 
                       locality + pH + water_speed_ms + water_depth + water_level + Cond + wmo_prec +
                       Temp_Water + site_type + seas_wmo + 
                       locality*seas_wmo + site_type*seas_wmo + site_type*wmo_prec + 
                       site_type*Cond + site_type*Temp_Water + site_type*pH +
                      offset(log(duration)),
                 data=chemdf,
                 family=nbinom1)

model1 <- glmmTMB(Bulinus_tot ~ (1 |locality/site_irn/visit_no) + 
                        locality + pH + water_speed_ms + water_depth + water_level + Cond + wmo_prec +
                        Temp_Water + site_type + seas_wmo + 
                        locality*seas_wmo + site_type*seas_wmo + site_type*wmo_prec + 
                        site_type*Cond + site_type*Temp_Water + site_type*pH +
                             offset(log(duration)),
                           data=chemdf,
                           family=nbinom2)

# Check which model fits best, based on AIC scores
anova(poiss, model, model1) # model and model1 have the same AICs, better than poiss

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
View(fulldf[fulldf$locality == "Kohan Garantche",] %>% dplyr::select(seas_wmo, locality, Bulinus_tot, everything()))
# There were many snails in the dry season, but zero in the wet season. This inflates the standard error hugely.
# Model fits can still be used, but standard errors will be inflated.

#' Test significance of predictors
# Use chisquare test on maximum model, to decide which term to remove first

drop1(model1, test = "Chisq")        # Temp_Water:site_type  7 10408  3.398 0.8458894  

model2 <- glmmTMB(Bulinus_tot ~ (1 |locality/site_irn/visit_no) + 
                        locality + pH + water_speed_ms + water_depth + water_level + Cond + wmo_prec +
                        Temp_Water + site_type + seas_wmo + 
                        locality*seas_wmo + site_type*seas_wmo + site_type*wmo_prec + 
                        site_type*Cond + site_type*pH +
                        offset(log(duration)),
                  data=chemdf,
                  family=nbinom2)
drop1(model2, test = "Chisq")        # pH:site_type        7 10399  5.205 0.6349038 

model3 <- glmmTMB(Bulinus_tot ~ (1 |locality/site_irn/visit_no) + 
                        locality + pH + water_speed_ms + water_depth + water_level + Cond + wmo_prec +
                        Temp_Water + site_type + seas_wmo + 
                        locality*seas_wmo + site_type*seas_wmo + site_type*wmo_prec + 
                        site_type*Cond +
                        offset(log(duration)),
                  data=chemdf,
                  family=nbinom2)
drop1(model3, test = "Chisq")        # pH                  1 10397  0.090 0.7635879   

model4 <- glmmTMB(Bulinus_tot ~ (1 |locality/site_irn/visit_no) + 
                        locality + water_speed_ms + water_depth + water_level + Cond + wmo_prec +
                        Temp_Water + site_type + seas_wmo + 
                        locality*seas_wmo + site_type*seas_wmo + site_type*wmo_prec + 
                        site_type*Cond +
                        offset(log(duration)),
                  data=chemdf,
                  family=nbinom2)
drop1(model4, test = "Chisq")        # water_depth         1 10396  1.202 0.2729720  

model5 <- glmmTMB(Bulinus_tot ~ (1 |locality/site_irn/visit_no) + 
                        locality + water_speed_ms + water_level + Cond + wmo_prec +
                        Temp_Water + site_type + seas_wmo +
                        locality*seas_wmo + site_type*seas_wmo + site_type*wmo_prec + 
                        site_type*Cond +
                        offset(log(duration)),
                  data=chemdf,
                  family=nbinom2)
drop1(model5, test = "Chisq")      # Temp_Water          1 10396  1.799  0.179807  
MASS::dropterm(model5, test = "Chisq")

model6 <- glmmTMB(Bulinus_tot ~ (1 |locality/site_irn/visit_no) + 
                        locality + water_speed_ms + water_level + Cond + wmo_prec +
                        site_type + seas_wmo +
                        locality*seas_wmo + site_type*seas_wmo + site_type*wmo_prec + 
                        site_type*Cond +
                        offset(log(duration)),
                  data=chemdf,
                  family=nbinom2)
drop1(model6, test = "Chisq")      # water_speed_ms      1 10397  3.103 0.0781317 .  

model7 <- glmmTMB(Bulinus_tot ~ (1 |locality/site_irn/visit_no) +
                        locality + water_level + Cond + wmo_prec +
                        site_type + seas_wmo +
                        locality*seas_wmo + site_type*seas_wmo + site_type*wmo_prec + 
                        site_type*Cond +
                        offset(log(duration)),
                  data=chemdf[chemdf$locality != "Kohan Garantche",],
                  family=nbinom2)

#                          Df   AIC  LRT    Pr(>Chi)    
#       <none>                10232                     
#       water_level         2 10267 39.2 0.000000003 ***
#       locality:seas_wmo  18 10234 37.6      0.0044 ** 
#       site_type:seas_wmo  7 10241 23.3      0.0015 ** 
#       wmo_prec:site_type  7 10233 15.3      0.0329 *  
#       Cond:site_type      7 10238 20.2      0.0052 ** 
  
summary(model7)

#' Interaction analysis using emmeans
glmmTMB::Anova.glmmTMB(model7, test.statistic = "Chisq")

sim_residuals7 <- DHARMa::simulateResiduals(model7, 1000) 
plot(sim_residuals7)

emmip(model7, site_type ~ seas_wmo)
test <- data.frame(emmeans(model7, pairwise ~ seas_wmo | site_type))
test2 <- allEffects(model7)
plot(test2$water_level)
plot(test2$`locality:seas_wmo`)
df <- data.frame(test2$`site_type:seas_wmo`)


confint(test, adjust = "tukey")
df <- data.frame(Effect.glmmTMB(model7, focal.predictors = c("site_type", "seas_wmo")))

ggplot(df) +
      geom_bar(aes(x = reorder(site_type, fit), y = fit, fill = seas_wmo), col = "black",
               position = position_dodge(), stat = "identity") +
      geom_errorbar(aes(x = site_type, ymin = fit-se, ymax = fit+se,
                        group = seas_wmo), position = position_dodge())


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
########################      Bulinus positive total      ###############################
#########################################################################################
# There is insufficient data, if we want to include the same environmental factors as in the 
# previous model. Use a data frame that excludes fewer rows, by using fewer preictors

# Hurdle models would not be appropriate (see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3238139/pdf/nihms-342982.pdf).
# They assume that that 'all zero data are from one “structural” source.' E.g. Study with number of cigarettes smoked,
# where only non-smokers will smoke zero cigarettes, and smokers will score some positive (non-zero) number of cigarettes.
# Hence the zero observations can come from only one “structural” source, the non-smokers.
# If a subject is considered a smoker, they do not have the ‘ability’ to score zero cigarettes and will always score 
# a positive number of cigarettes.
# 
# Zero-inflated models assume mix of structural and sampling zero data (some zeros by chance, some zeros due to specific
# structure in the data). 
# This is appropriate. Some zero positive counts might be because some snails just didn't get infected, but some zero counts 
# might be because there were no schistosomes in the water, so no chance of infection to begin with.

mod_a <- glmmTMB(Bulinus_pos_tot ~ (1|locality/site_irn/visit_no) + locality + site_type + seas_wmo + 
                       locality*seas_wmo + site_type*seas_wmo +
                   offset(log(duration)),
                 data=fulldf,
                 family=poisson)

# Model convergence problem; non-positive-definite Hessian matrix.
# This is problematic, and could mean that a random-effect variance is estimated to be zero, or
# the model is overparameterized (data does not contain enough information to estimate the parameters reliably).

#' Diagnose the convergence problem
# Try different family
mod_b <- glmmTMB(Bulinus_pos_tot ~ (1|locality/site_irn/visit_no) + locality + site_type + seas_wmo + 
                       locality*seas_wmo + site_type*seas_wmo +
                       offset(log(duration)),
                 data=fulldf,
                 family=nbinom1) # Many convergence problems

mod_c <- glmmTMB(Bulinus_pos_tot ~ (1|locality/site_irn/visit_no) + locality + site_type + seas_wmo + 
                       locality*seas_wmo + site_type*seas_wmo +
                       offset(log(duration)),
                 data=fulldf,
                 family=nbinom2) # Some convergence problems

# Check overdispersion
sim_res_moda <- DHARMa::simulateResiduals(mod_a, 1000)
sim_res_modb <- DHARMa::simulateResiduals(mod_b, 1000)
sim_res_modc <- DHARMa::simulateResiduals(mod_c, 1000)


plot(sim_res_moda) # Deviation significant
plot(sim_res_modb) # Deviation not significant
plot(sim_res_modc) # Deviation not significant

testDispersion(sim_res_moda) # Overdispersed
testDispersion(sim_res_modb) # Overdispersed
testDispersion(sim_res_modc) # Not overdispersed

# Check zero inflation
testZeroInflation(sim_res_moda)  # Zero-inflated.
testZeroInflation(sim_res_modb)  # Not zero-inflated. 
testZeroInflation(sim_res_modc)  # Not zero-inflated.

# Use mod_c, remove terms untl no more convergence problems
mod_d <- glmmTMB(Bulinus_pos_tot ~ (1|locality/site_irn/visit_no) + locality +
                       offset(log(duration)),
                 data=fulldf,
                 family=nbinom2) # No convergence problems
drop1(mod_d, test = "Chisq")  # locality 18 861.01 37.543  0.00445 **
summary(mod_d)
mod_e <- glmmTMB(Bulinus_pos_tot ~ (1|locality/site_irn/visit_no) + seas_wmo + site_type*seas_wmo +
                       offset(log(duration)),
                 data=fulldf,
                 family=nbinom2) # No convergence problems
drop1(mod_e, test = "Chisq")  # seas_wmo:site_type  7 867 9.35     0.23

mod_f <- glmmTMB(Bulinus_pos_tot ~ (1|locality/site_irn/visit_no) + seas_wmo +
                       offset(log(duration)),
                 data=fulldf,
                 family=nbinom2) # No convergence problems
drop1(mod_f, test = "Chisq")  # seas_wmo  1 862 0.148      0.7

mod_g <- glmmTMB(Bulinus_pos_tot ~ (1|locality/site_irn/visit_no) + site_type + 
                       offset(log(duration)),
                 data=snaildf2,
                 family=nbinom2) # No convergence problems
drop1(mod_g, test = "Chisq")  # site_type  7 861.01 10.324   0.1709

mod_h <- glmmTMB(Bulinus_pos_tot ~ (1|locality/site_irn/visit_no) + Cond + 
                       offset(log(duration)),
                 data=snaildf,
                 family=nbinom2) # No convergence problems
drop1(mod_h, test = "Chisq")  # Cond    1 787.42 0.99221   0.3192

mod_i <- glmmTMB(Bulinus_pos_tot ~ (1|locality/site_irn/visit_no) + pH +
                       offset(log(duration)),
                 data=chemdf,
                 family=nbinom2) # Convergence problems
plot(snaildf$pH, snaildf$Bulinus_pos_tot)    # Maybe nonlinear relationship?

mod_j <- glmmTMB(Bulinus_pos_tot ~ (1|locality/site_irn/visit_no) + PPM + 
                       offset(log(duration)),
                 data=snaildf,
                 family=nbinom2) # No convergence problems
drop1(mod_j, test = "Chisq")  # PPM     1 787.42 0.64707   0.4212

mod_k <- glmmTMB(Bulinus_pos_tot ~ (1|locality/site_irn/visit_no) + wmo_prec + 
                       offset(log(duration)),
                 data=snaildf,
                 family=nbinom2) # No convergence problems
drop1(mod_j, test = "Chisq")  # wmo_prec  1 787.42 1.198   0.2737
plot(snaildf$Temp_Water, snaildf$Bulinus_pos_tot)   # Maybe nonlinear relationship?

mod_k <- glmmTMB(Bulinus_pos_tot ~ (1|locality/site_irn/visit_no) + Stagnante + 
                       offset(log(duration)),
                 data=snaildf,
                 family=nbinom2) # No convergence problems
drop1(mod_k, test = "Chisq")  # Stagnante  1 787.42 0.48501   0.4862

mod_l <- glmmTMB(Bulinus_pos_tot ~ (1|locality/site_irn/visit_no) + water_depth + 
                       offset(log(duration)),
                 data=snaildf,
                 family=nbinom2) # No convergence problems
drop1(mod_k, test = "Chisq")  # water_depth  1 787.42 0.098649   0.7535

mod_l <- glmmTMB(Bulinus_pos_tot ~ (1|locality/site_irn/visit_no) + water_level + 
                       offset(log(duration)),
                 data=snaildf,
                 family=nbinom2) # No convergence problems
drop1(mod_l, test = "Chisq")  # water_depth  1 787.42 0.098649   0.7535

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
########################      Bulinus truncatus total      ##############################
#########################################################################################
# Include other species present yes/no
mod_1 <- glmmTMB(BT_tot ~ (1 |locality/site_irn/visit_no) + 
locality + pH + water_speed_ms + water_depth + water_level + Cond + wmo_prec +
      Temp_Water + site_type + seas_wmo + BP_tot + BF_tot + BP_tot*BF_tot +
      locality*seas_wmo + site_type*seas_wmo + site_type*Temp_Water + site_type*pH +
      offset(log(duration)),
data=chemdf,
family=nbinom2)

summary(mod_1)
sim_res_mod1 <- DHARMa::simulateResiduals(mod_1, 1000)
plot(sim_res_mod1)  # Not overdispersed.
testZeroInflation(sim_res_mod1)  # No zero-inflation. 

drop1(mod_1, test = "Chisq")   # pH:site_type          7 8215  5.1    0.65174

mod_2 <- glmmTMB(BT_tot ~ (1 |locality/site_irn/visit_no) + 
                       locality + pH + water_speed_ms + water_depth + water_level + Cond + wmo_prec +
                       Temp_Water + site_type + seas_wmo + BP_tot + BF_tot + BP_tot*BF_tot +
                       locality*seas_wmo + site_type*seas_wmo + site_type*Temp_Water +
                       offset(log(duration)),
                 data=chemdf,
                 family=nbinom2)
drop1(mod_2, test = "Chisq")   # pH                    1 8216  2.3   0.12689    

mod_3 <- glmmTMB(BT_tot ~ (1 |locality/site_irn/visit_no) + 
                       locality + water_speed_ms + water_depth + water_level + Cond + wmo_prec +
                       Temp_Water + site_type + seas_wmo + BP_tot + BF_tot + BP_tot*BF_tot +
                       locality*seas_wmo + site_type*seas_wmo + site_type*Temp_Water +
                       offset(log(duration)),
                 data=chemdf,
                 family=nbinom2)
drop1(mod_3, test = "Chisq")   # Temp_Water:site_type  7 8214 12.5  0.08452 .       

mod_4 <- glmmTMB(BT_tot ~ (1 |locality/site_irn/visit_no) + 
                       locality + water_speed_ms + water_depth + water_level + Cond + wmo_prec +
                       Temp_Water + site_type + seas_wmo + BP_tot + BF_tot + BP_tot*BF_tot +
                       locality*seas_wmo + site_type*seas_wmo +
                       offset(log(duration)),
                 data=chemdf,
                 family=nbinom2)
drop1(mod_4, test = "Chisq")   # Temp_Water          1 8212  0.1    0.7109  

mod_5 <- glmmTMB(BT_tot ~ (1 |locality/site_irn/visit_no) + 
                       locality + water_speed_ms + water_depth + water_level + Cond + wmo_prec +
                       site_type + seas_wmo + BP_tot + BF_tot + BP_tot*BF_tot +
                       locality*seas_wmo + site_type*seas_wmo +
                       offset(log(duration)),
                 data=chemdf,
                 family=nbinom2)
drop1(mod_5, test = "Chisq")   # nothing to remove
mod_5 <- glmmTMB(BT_tot ~ (1 |locality/site_irn/visit_no) + 
                       locality + water_speed_ms + water_depth + water_level + Cond + wmo_prec +
                       site_type + seas_wmo + BP_tot + BF_tot + BP_tot*BF_tot +
                       locality*seas_wmo + site_type*seas_wmo +
                       offset(log(duration)),
                 data=chemdf[chemdf$locality != "Kohan Garantche" & chemdf$locality != "Gantchi Bassarou", ],
                 family=nbinom2)

df2 <- data.frame(Effect.glmmTMB(mod_5, focal.predictors = c("site_type", "seas_wmo")))

ggplot(df2) +
      geom_bar(aes(x = site_type, y = fit, fill = seas_wmo), col = "black",
               position = position_dodge(), stat = "identity") +
      geom_errorbar(aes(x = site_type, ymin = fit-se, ymax = fit+se,
                        group = seas_wmo), position = position_dodge())
summary(mod_5)

plot(chemdf$site_type)

test <- allEffects(mod_5)
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
mod_A <- glmmTMB(BF_tot ~ (1 |locality/site_irn/visit_no) + 
                       locality + pH + water_speed_ms + water_depth + water_level + Cond + wmo_prec +
                       Temp_Water + site_type + seas_wmo + BP_tot + BT_tot + BP_tot*BT_tot +
                       locality*seas_wmo + site_type*seas_wmo + site_type*Temp_Water + site_type*pH +
                       offset(log(duration)),
                 data=chemdf[chemdf$locality != "Kohan Garantche", ],
                 family=nbinom2)
summary(mod_A)
sim_res_modA <- DHARMa::simulateResiduals(mod_A, 1000)
plot(sim_res_modA)  # Not overdispersed.
testZeroInflation(sim_res_modA)  # No zero-inflation. 

drop1(mod_A, test = "Chisq")   # Cond                  1 5948  0.0   0.9540  

mod_B <- glmmTMB(BF_tot ~ (1 |locality/site_irn/visit_no) + 
                       locality + pH + water_speed_ms + water_depth + water_level + wmo_prec +
                       Temp_Water + site_type + seas_wmo + BP_tot + BT_tot + BP_tot*BT_tot +
                       locality*seas_wmo + site_type*seas_wmo + site_type*Temp_Water + site_type*pH +
                       offset(log(duration)),
                 data=chemdf[chemdf$locality != "Kohan Garantche", ],
                 family=nbinom2)
drop1(mod_B, test = "Chisq")   # Temp_Water:site_type  7 5939  4.7   0.6934  

mod_C <- glmmTMB(BF_tot ~ (1 |locality/site_irn/visit_no) + 
                       locality + pH + water_speed_ms + water_depth + water_level + wmo_prec +
                       Temp_Water + site_type + seas_wmo + BP_tot + BT_tot + BP_tot*BT_tot +
                       locality*seas_wmo + site_type*seas_wmo + site_type*pH +
                       offset(log(duration)),
                 data=chemdf[chemdf$locality != "Kohan Garantche", ],
                 family=nbinom2)
drop1(mod_C, test = "Chisq")   # wmo_prec            1 5937  0.1  0.80410 

mod_D <- glmmTMB(BF_tot ~ (1 |locality/site_irn/visit_no) + 
                       locality + pH + water_speed_ms + water_depth + water_level + 
                       Temp_Water + site_type + seas_wmo + BP_tot + BT_tot + BP_tot*BT_tot +
                       locality*seas_wmo + site_type*seas_wmo + site_type*pH +
                       offset(log(duration)),
                 data=chemdf[chemdf$locality != "Kohan Garantche", ],
                 family=nbinom2)
drop1(mod_D, test = "Chisq")   # water_depth         1 5935  0.3  0.61265   

mod_E <- glmmTMB(BF_tot ~ (1 |locality/site_irn/visit_no) + 
                       locality + pH + water_speed_ms + water_level + 
                       Temp_Water + site_type + seas_wmo + BP_tot + BT_tot + BP_tot*BT_tot +
                       locality*seas_wmo + site_type*seas_wmo + site_type*pH +
                       offset(log(duration)),
                 data=chemdf[chemdf$locality != "Kohan Garantche", ],
                 family=nbinom2)
drop1(mod_E, test = "Chisq")   # BP_tot:BT_tot       1 5934  0.4   0.54452  

mod_F <- glmmTMB(BF_tot ~ (1 |locality/site_irn/visit_no) + 
                       locality + pH + water_speed_ms + water_level + 
                       Temp_Water + site_type + seas_wmo + BP_tot + BT_tot +
                       locality*seas_wmo + site_type*seas_wmo + site_type*pH +
                       offset(log(duration)),
                 data=chemdf[chemdf$locality != "Kohan Garantche", ],
                 family=nbinom2)
drop1(mod_F, test = "Chisq")   # pH:site_type        7 5930 10.2     0.17847  

mod_G <- glmmTMB(BF_tot ~ (1 |locality/site_irn/visit_no) + 
                       locality + pH + water_speed_ms + water_level + 
                       Temp_Water + site_type + seas_wmo + BP_tot + BT_tot +
                       locality*seas_wmo + site_type*seas_wmo +
                       offset(log(duration)),
                 data=chemdf[chemdf$locality != "Kohan Garantche", ],
                 family=nbinom2)
drop1(mod_G, test = "Chisq")   # pH                  1 5929  1.3     0.26017    

mod_H <- glmmTMB(BF_tot ~ (1 |locality/site_irn/visit_no) + 
                       locality + water_speed_ms + water_level + 
                       Temp_Water + site_type + seas_wmo + BP_tot + BT_tot +
                       locality*seas_wmo + site_type*seas_wmo +
                       offset(log(duration)),
                 data=chemdf[chemdf$locality != "Kohan Garantche", ],
                 family=nbinom2)
summary(mod_H)

df3 <- data.frame(Effect.glmmTMB(mod_H, focal.predictors = c("site_type", "seas_wmo")))

ggplot(df3) +
      geom_bar(aes(x = site_type, y = fit, fill = seas_wmo), col = "black",
               position = position_dodge(), stat = "identity") +
      geom_errorbar(aes(x = site_type, ymin = fit-se, ymax = fit+se,
                        group = seas_wmo), position = position_dodge())
summary(mod_H)



#########################################################################################
########################      Summary of final models      ##############################
#########################################################################################
model7 <- glmmTMB(Bulinus_tot ~ (1 |locality/site_irn/visit_no) +
                        locality + water_level.v + Cond + wmo_prec +
                        site_type + seas_wmo +
                        locality*seas_wmo + site_type*seas_wmo + site_type*wmo_prec + 
                        site_type*Cond +
                        offset(log(duration)),
                  data=chemdf[chemdf$locality != "Kohan Garantche", ],
                  family=nbinom2)

mod_d <- glmmTMB(Bulinus_pos_tot ~ (1|locality/site_irn/visit_no) + locality +
                       offset(log(duration)),
                 data=chemdf,
                 family=nbinom2) 

mod_5 <- glmmTMB(BT_tot ~ (1 |locality/site_irn/visit_no) + 
                       locality + water_speed_ms + water_depth + water_level + Cond + wmo_prec +
                       site_type + seas_wmo + BP_tot + BF_tot + BP_tot*BF_tot +
                       locality*seas_wmo + site_type*seas_wmo +
                       offset(log(duration)),
                 data=chemdf[chemdf$locality != "Kohan Garantche" & chemdf$locality != "Gantchi Bassarou", ],
                 family=nbinom2)

mod_H <- glmmTMB(BF_tot ~ (1 |locality/site_irn/visit_no) + 
                       locality + water_speed_ms + water_level + 
                       Temp_Water + site_type + seas_wmo + BP_tot + BT_tot +
                       locality*seas_wmo + site_type*seas_wmo +
                       offset(log(duration)),
                 data=chemdf[chemdf$locality != "Kohan Garantche", ],
                 family=nbinom2)

#########################################################################################
############################      Plot of emmeans      ##################################
#########################################################################################

#' Bulinus_tot
# locality:seas_wmo
#  emmeans
emmip(model7, locality ~ seas_wmo)

locality_seas_wmo_df <- data.frame(emmeans(model7,  ~ locality*seas_wmo, lmer.df = "tukey"))
pairs(emmeans(model7,  ~ locality*seas_wmo), simple = "seas_wmo")
ggplot(locality_seas_wmo_df) +
      geom_bar(aes(x = reorder(locality, emmean), y = emmean, fill = seas_wmo), col = "black",
               position = position_dodge(), stat = "identity") +
      geom_errorbar(aes(x = locality, ymin = emmean-SE, ymax = emmean+SE,
                        group = seas_wmo), position = position_dodge()) +
      scale_fill_manual(values = c("lightgrey", "darkgrey")) 

summary(model7)

# effects
df <- data.frame(glmmTMB::Effect.glmmTMB(model7, focal.predictors = c("locality", "seas_wmo")))
ggplot(df) +
      geom_bar(aes(x = reorder(locality, fit), y = fit, fill = seas_wmo), col = "black",
               position = position_dodge(), stat = "identity") +
      geom_errorbar(aes(x = locality, ymin = fit-se, ymax = fit+se,
                        group = seas_wmo), position = position_dodge()) +
      scale_fill_manual(values = c("lightgrey", "darkgrey")) 

# site_type:seas_wmo
#  emmeans
emmip(model7, site_type ~ seas_wmo)

site_type_seas_wmo_df <- data.frame(emmeans(model7,  ~ site_type*seas_wmo, lmer.df = "tukey"))
pairs(emmeans(model7,  ~ site_type*seas_wmo), simple = "seas_wmo")
ggplot(site_type_seas_wmo_df) +
      geom_bar(aes(x = reorder(site_type, emmean), y = emmean, fill = seas_wmo), col = "black",
               position = position_dodge(), stat = "identity") +
      geom_errorbar(aes(x = site_type, ymin = emmean-SE, ymax = emmean+SE,
                        group = seas_wmo), position = position_dodge()) +
      scale_fill_manual(values = c("lightgrey", "darkgrey")) 


# effects
df <- data.frame(glmmTMB::Effect.glmmTMB(model7, focal.predictors = c("site_type", "seas_wmo")))
ggplot(df) +
      geom_bar(aes(x = reorder(site_type, fit), y = fit, fill = seas_wmo), col = "black",
               position = position_dodge(), stat = "identity") +
      geom_errorbar(aes(x = site_type, ymin = fit-se, ymax = fit+se,
                        group = seas_wmo), position = position_dodge()) +
      scale_fill_manual(values = c("lightgrey", "darkgrey")) 


