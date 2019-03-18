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



# Change font size and ggthemes globally
theme_set(ggthemes::theme_few(base_size = 8))

# Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  

# Import some helpful functions
source("HighstatLibV10.R")
setwd("..")

#' ### Import and prepare the data set
snaildf <- read.csv("Data/survey_para_merge_ALL_2018-10-13.csv") %>% 
  dplyr::select(filter.min, coll_date, month, passage_no, locality, site_irn, site_type_clean, BP_tot, BP_pos_tot, BF_tot, BF_pos_tot, 
                BT_tot, BT_pos_tot, BG_tot, BS_tot, BT_prev, BF_prev, Bulinus_tot, Bulinus_pos_tot,
                bp_pres, bt_pres, bf_pres,
                Temp_Air, Temp_Eau,w_level_ed, water_speed_ms, water_depth, pH, Cond, PPM, Latitude, Longitude, Stagnante,
                wmo_min_temp, wmo_max_temp, wmo_av_temp, wmo_prec, seas_wmo, duration, Heure) %>% 
  rename(visit_no = passage_no, Temp_Water = Temp_Eau, water_level = w_level_ed, site_type = site_type_clean) %>% 
  # Remove NA values in predictor variables
  dplyr::filter(Bulinus_tot <= 1000
                ,!is.na(BT_tot),
                !is.na(Temp_Air), !is.na(water_speed_ms), !is.na(water_depth), !is.na(pH), !is.na(PPM),
                !is.na(wmo_av_temp), !is.na(wmo_prec), !is.na(Stagnante), !is.na(site_type), !is.na(site_irn),
                !is.na(Temp_Water)
                , duration != "") %>% 
  mutate(coll_date = lubridate::dmy(coll_date), tz = "Africa/Niger",
         duration = as.numeric(as.character(duration)),
         site_irn = as.factor(site_irn), 
         visit_no = as.factor(visit_no),
         bp_pres = as.factor(as.character(bp_pres)),
         bf_pres = as.factor(as.character(bf_pres))) %>% 
  filter(!is.na(duration)) %>% 
  # Re-arrange the columns
  dplyr::select(locality, site_irn, visit_no, Bulinus_tot, Bulinus_pos_tot, coll_date, everything()) %>% 
  arrange(locality, site_irn) %>% 
  # This site has no measurements in the wet season, and a lot in the dry season. Later in the models, this leads to 
  # very high variance, when we're looking at interactions containing the season variable.
  # For now, remove this site from the analysis.
  filter(filter.min != "y") %>% 
      # Rescale variables with large values
      mutate(Cond = scale(Cond, center = 0),
             water_depth = scale(water_depth, center = 0))

#' ### Explore the data
#' How uniform was the sampling effort (duration of sampling)?
tiff(file="Figures/sampling_effort.tiff", width = 166, height = 166, units = "mm", res = 150)
ggplot(snaildf) +
  geom_line(aes(x = as.factor(visit_no), y = duration, group = locality), alpha = .2) +
  geom_point(aes(x = as.factor(visit_no), y = duration, group = locality), pch = 23, fill = "white", cex = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~locality) +
  ggtitle("Sampling effort")
dev.off()

# There was a lot of variation. Since sampling time most likely influences counts, this needs to be integrated in the analylsis

#' Are there any variance inflation factors (multicollinearity)? Check using a function from Zuur et al. 2010

pairs(snaildf[,c("Temp_Air", "Temp_Water", "water_speed_ms", "water_depth", "pH", "Cond", "PPM",
                 "wmo_av_temp", "Bulinus_tot", "wmo_prec", "site")],
      lower.panel = panel.cor)

corvif(as.data.table(snaildf)[, c("Temp_Air", "Temp_Water", "water_speed_ms", "water_depth", "pH", "Cond", "PPM",
                                  "wmo_av_temp", "wmo_prec"), with=FALSE])
# Cond and PPM have high GVIF values (10). For values of higher than 4, only one of the two variables should be used in models, 
# to avoid multicollinearity. 

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

#' How is the data distributed?
ggplot(snaildf) +
  geom_histogram(aes(Bulinus_tot)) +
  theme_few() +
  xlab("Snail count") +
  facet_wrap(~locality, scales = "free", nrow = 2)

# There are generally many zeroes. This doesn't mean that the data is zero-inflated, but it might be worth checking for zero-inflation anyways

#' ### Make a GLMM

#' Make a maximum model. glmmTMB is a new package by Ben Bolker, that fits models faster, and allows to include arguments
# to account for zero-inflation, if needed

# Decide for a family. It's count count data, looking very overdispersed, so negative binomial is most likely appropriate
poiss <- glmmTMB(Bulinus_tot ~ (1|locality/site_irn/visit_no) + locality + pH + water_speed_ms + water_depth + Cond + wmo_prec +
                   Temp_Air + Temp_Water + site_type + seas_wmo + month +
                   locality*seas_wmo + wmo_prec*seas_wmo + locality*month + site_type*seas_wmo +
                   offset(log(duration)),
                 data=snaildf,
                 family=poisson)

model <- glmmTMB(Bulinus_tot ~ (1|locality/site_irn/visit_no) + locality + pH + water_speed_ms + water_depth + Cond + wmo_prec +
                   Temp_Air + Temp_Water + site_type + seas_wmo + month +
                   locality*seas_wmo + wmo_prec*seas_wmo + locality*month + site_type*seas_wmo +
                      offset(log(duration)),
                 data=snaildf,
                 family=nbinom1)

model1 <- glmmTMB(Bulinus_tot ~ (1|locality/site_irn/visit_no) + locality + pH + water_speed_ms + water_depth + Cond + wmo_prec +
                             Temp_Air + Temp_Water + site_type + seas_wmo + month +
                             locality*seas_wmo + wmo_prec*seas_wmo + locality*month + site_type*seas_wmo +
                             offset(log(duration)),
                           data=snaildf,
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

testZeroInflation(sim_residuals1)  # There is no evidence for zero-inflation

# Look at the model

summary(model1)  # Next to no variance is coming from locality, so technically we could exclude this from the model.
                 # Aesthetically, we can keep it in, to highlight the nested design of the study

#' Test significance of predictors
# Use chisquare test on maximum model, to decide which term to remove first

drop1(model1, test = "Chisq")        # wmo_prec:seas_wmo   1 10668  0.064 0.8009513  

model2 <- glmmTMB(Bulinus_tot ~ (1|locality/site_irn/visit_no) + locality + pH + water_speed_ms + water_depth + Cond + wmo_prec +
                    Temp_Air + Temp_Water + site_type + seas_wmo + month +
                    locality*seas_wmo + locality*month + site_type*seas_wmo +
                    offset(log(duration)),
                  data=snaildf,
                  family=nbinom2)
drop1(model2, test = "Chisq")        # pH                  1 10666  0.293 0.5880312

model3 <- glmmTMB(Bulinus_tot ~ (1|locality/site_irn/visit_no) + locality + water_speed_ms + water_depth + Cond + wmo_prec +
                    Temp_Air + Temp_Water + site_type + seas_wmo + month +
                    locality*seas_wmo + locality*month + site_type*seas_wmo +
                    offset(log(duration)),
                  data=snaildf,
                  family=nbinom2)
drop1(model3, test = "Chisq")        # wmo_prec            1 10666  1.534 0.2155185

model4 <- glmmTMB(Bulinus_tot ~ (1|locality/site_irn/visit_no) + locality + water_speed_ms + water_depth + Cond +
                    Temp_Air + Temp_Water + site_type + seas_wmo + month +
                    locality*seas_wmo + locality*month + site_type*seas_wmo +
                    offset(log(duration)),
                  data=snaildf,
                  family=nbinom2)
drop1(model4, test = "Chisq")        # Cond                1 10666  1.705 0.1916683 

model5 <- glmmTMB(Bulinus_tot ~ (1|locality/site_irn/visit_no) + locality + water_speed_ms + water_depth +
                    Temp_Air + Temp_Water + site_type + seas_wmo + month +
                    locality*seas_wmo + locality*month + site_type*seas_wmo +
                    offset(log(duration)),
                  data=snaildf,
                  family=nbinom2)
drop1(model5, test = "Chisq")        # water_depth         1 10666  2.248 0.1338204  

model6 <- glmmTMB(Bulinus_tot ~ (1|locality/site_irn/visit_no) + locality + water_speed_ms +
                    Temp_Air + Temp_Water + site_type + seas_wmo + month +
                    locality*seas_wmo + locality*month + site_type*seas_wmo +
                    offset(log(duration)),
                  data=snaildf,
                  family=nbinom2)
drop1(model6, test = "Chisq")        # water_depth         1 10666  2.248 0.1338204  

  # water_speed_ms      1 10670  6.018 0.0141583 *  
  # Temp_Air            1 10670  6.326 0.0118976 *  
  # Temp_Water          1 10677 12.727 0.0003603 ***
  # locality:seas_wmo  18 10668 37.863 0.0040386 ** 
  # locality:month     18 10664 33.961 0.0127358 *  
  # site_type:seas_wmo  7 10671 18.633 0.0094164 ** 
  
summary(model6)
glmmTMB::Anova.glmmTMB(model6)
sim_residuals6 <- DHARMa::simulateResiduals(model6, 1000) 
plot(sim_residuals6)
# Could the effect of water temperature have something to do with water depth? Consider interaction.

model7 <- glmmTMB(Bulinus_tot ~ (1|locality/site_irn/visit_no) + locality + water_speed_ms + water_depth +
                    Temp_Air + Temp_Water + site_type + seas_wmo + month +
                    locality*seas_wmo + locality*month + site_type*seas_wmo + water_speed_ms*Temp_Water +
                    offset(log(duration)),
                  data=snaildf,
                  family=nbinom2)
drop1(model7, test = "Chisq")        # 

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

snaildf2 <- read.csv("Data/survey_para_merge_ALL_2018-10-13.csv") %>% 
      dplyr::select(filter.min, coll_date, month, passage_no, locality, site_irn, site_type_clean, BP_tot, BP_pos_tot, BF_tot, BF_pos_tot, 
                    BT_tot, BT_pos_tot, BG_tot, BS_tot, BT_prev, BF_prev, Bulinus_tot, Bulinus_pos_tot,
                    Temp_Air, Temp_Eau,w_level_ed, water_speed_ms, water_depth, pH, Cond, PPM, Latitude, Longitude, Stagnante,
                    wmo_min_temp, wmo_max_temp, wmo_av_temp, wmo_prec, seas_wmo, duration, Heure) %>% 
      rename(visit_no = passage_no, Temp_Water = Temp_Eau, water_level = w_level_ed, site_type = site_type_clean) %>% 
      # Remove NA values in predictor variables
      dplyr::filter(Bulinus_tot <= 1000,!is.na(site_type), !is.na(site_irn), !is.na(seas_wmo)
                    , duration != "") %>% 
      mutate(coll_date = lubridate::dmy(coll_date), tz = "Africa/Niger",
             duration = as.numeric(as.character(duration)),
             site_irn = as.factor(site_irn), 
             visit_no = as.factor(visit_no)) %>% 
      filter(!is.na(duration)) %>% 
      # Re-arrange the columns
      dplyr::select(Bulinus_pos_tot, locality, site_irn, visit_no, coll_date, seas_wmo, everything()) %>% 
      arrange(locality, site_irn) %>% 
      # This site has no measurements in the wet season, and a lot in the dry season. Later in the models, this leads to 
      # very high variance, when we're looking at interactions containing the season variable.
      # For now, remove this site from the analysis.
      filter(locality != "Kohan Garantche", filter.min != "y")


ggplot(snaildf2) +
  geom_histogram(aes(Bulinus_pos_tot), binwidth = .5, position = position_dodge()) +
  theme_few() +
  xlab("Snail count") +
  facet_wrap(~locality, scales = "free", nrow = 2)



mod_a <- glmmTMB(Bulinus_pos_tot ~ (1|locality/site_irn/visit_no) + locality + site_type + seas_wmo + 
                       locality*seas_wmo + site_type*seas_wmo +
                   offset(log(duration)),
                 data=snaildf2,
                 family=poisson)

# Model convergence problem; non-positive-definite Hessian matrix.
# This is problematic, and could mean that a random-effect variance is estimated to be zero, or
# the model is overparameterized (data does not contain enough information to estimate the parameters reliably).

#' Diagnose the convergence problem
# Try different family
mod_b <- glmmTMB(Bulinus_pos_tot ~ (1|locality/site_irn/visit_no) + locality + site_type + seas_wmo + 
                       locality*seas_wmo + site_type*seas_wmo +
                       offset(log(duration)),
                 data=snaildf2,
                 family=nbinom1) # Many convergence problems

mod_c <- glmmTMB(Bulinus_pos_tot ~ (1|locality/site_irn/visit_no) + locality + site_type + seas_wmo + 
                       locality*seas_wmo + site_type*seas_wmo +
                       offset(log(duration)),
                 data=snaildf2,
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
                 data=snaildf2,
                 family=nbinom2) # No convergence problems
drop1(mod_d, test = "Chisq")  # locality 18 861.01 37.543  0.00445 **
summary(mod_d)
mod_e <- glmmTMB(Bulinus_pos_tot ~ (1|locality/site_irn/visit_no) + seas_wmo + 
                       offset(log(duration)),
                 data=snaildf2,
                 family=nbinom2) # No convergence problems
drop1(mod_e, test = "Chisq")  # seas_wmo  1 861.01 0.15887   0.6902

mod_f <- glmmTMB(Bulinus_pos_tot ~ (1|locality/site_irn/visit_no) + site_type + 
                       offset(log(duration)),
                 data=snaildf2,
                 family=nbinom2) # No convergence problems
drop1(mod_f, test = "Chisq")  # site_type  7 861.01 10.324   0.1709

mod_g <- glmmTMB(Bulinus_pos_tot ~ (1|locality/site_irn/visit_no) + Cond + 
                       offset(log(duration)),
                 data=snaildf,
                 family=nbinom2) # No convergence problems
drop1(mod_g, test = "Chisq")  # Cond    1 787.42 0.99221   0.3192

mod_h <- glmmTMB(Bulinus_pos_tot ~ (1|locality/site_irn/visit_no) + pH +
                       offset(log(duration)),
                 data=snaildf,
                 family=nbinom2) # Convergence problems
plot(snaildf$pH, snaildf$Bulinus_pos_tot)    # Maybe nonlinear relationship?

mod_i <- glmmTMB(Bulinus_pos_tot ~ (1|locality/site_irn/visit_no) + PPM + 
                       offset(log(duration)),
                 data=snaildf,
                 family=nbinom2) # No convergence problems
drop1(mod_i, test = "Chisq")  # PPM     1 787.42 0.64707   0.4212

mod_j <- glmmTMB(Bulinus_pos_tot ~ (1|locality/site_irn/visit_no) + wmo_prec + 
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
mod_1 <- glmmTMB(BT_tot ~ (1|locality/site_irn/visit_no) + locality + pH + water_speed_ms + water_depth + Cond + wmo_prec +
                   Temp_Air + Temp_Water + site_type + seas_wmo + month + bp_pres + bf_pres + bp_pres*bf_pres +
                   locality*seas_wmo + wmo_prec*seas_wmo + locality*month + site_type*seas_wmo + water_speed_ms*Temp_Water +
                   water_depth*water_speed_ms +
                   offset(log(duration)),
                 data=snaildf,
                 family=nbinom2) # Model convergence problem; extreme or very small eigen values detected

# Removing 'month' from the model resolves this problem
mod_2 <- glmmTMB(BT_tot ~ (1|locality/site_irn/visit_no) + locality + pH + water_speed_ms + water_depth + Cond + wmo_prec +
                       Temp_Air + Temp_Water + site_type + seas_wmo + bp_pres + bf_pres + bp_pres*bf_pres +
                       locality*seas_wmo + wmo_prec*seas_wmo + site_type*seas_wmo + water_speed_ms*Temp_Water +
                       water_depth*water_speed_ms +
                       offset(log(duration)),
                 data=snaildf,
                 family=nbinom2) # No more model convergence problems


sim_res_mod1 <- DHARMa::simulateResiduals(mod_2, 1000)
plot(sim_res_mod1)  # Not overdispersed.
testZeroInflation(sim_res_mod1)  # No zero-inflation. 

drop1(mod_2, test = "Chisq")   # Temp_Air                    1 8229.9  0.005 0.943678  

mod_3 <- glmmTMB(BT_tot ~ (1|locality/site_irn/visit_no) + locality + pH + water_speed_ms + water_depth + Cond + wmo_prec +
                       Temp_Water + site_type + seas_wmo + bp_pres + bf_pres + bp_pres*bf_pres +
                       locality*seas_wmo + wmo_prec*seas_wmo + site_type*seas_wmo + water_speed_ms*Temp_Water +
                       water_depth*water_speed_ms +
                       offset(log(duration)),
                 data=snaildf,
                 family=nbinom2)
drop1(mod_3, test = "Chisq")   # water_speed_ms:Temp_Water   1 8229.7  1.804 0.179283   

mod_4 <- glmmTMB(BT_tot ~ (1|locality/site_irn/visit_no) + locality + pH + water_speed_ms + water_depth + Cond + wmo_prec +
                       Temp_Water + site_type + seas_wmo + bp_pres + bf_pres + bp_pres*bf_pres +
                       locality*seas_wmo + wmo_prec*seas_wmo + site_type*seas_wmo + 
                       water_depth*water_speed_ms +
                       offset(log(duration)),
                 data=snaildf,
                 family=nbinom2)
drop1(mod_4, test = "Chisq")  # Temp_Water                  1 8228.4  0.713 0.398490   

mod_5 <- glmmTMB(BT_tot ~ (1|locality/site_irn/visit_no) + locality + pH + water_speed_ms + water_depth + Cond + wmo_prec +
                       site_type + seas_wmo + bp_pres + bf_pres + bp_pres*bf_pres +
                       locality*seas_wmo + wmo_prec*seas_wmo + site_type*seas_wmo + 
                       water_depth*water_speed_ms +
                       offset(log(duration)),
                 data=snaildf,
                 family=nbinom2)
drop1(mod_5, test = "Chisq")  # wmo_prec:seas_wmo           1 8228.4  1.969 0.160567  

mod_6 <- glmmTMB(BT_tot ~ (1|locality/site_irn/visit_no) + locality + pH + water_speed_ms + water_depth + Cond + wmo_prec +
                       site_type + seas_wmo + bp_pres + bf_pres + bp_pres*bf_pres +
                       locality*seas_wmo + site_type*seas_wmo + 
                       water_depth*water_speed_ms +
                       offset(log(duration)),
                 data=snaildf,
                 family=nbinom2)
drop1(mod_6, test = "Chisq")  # pH                          1 8229.4  2.988 0.083878 . 

mod_7 <- glmmTMB(BT_tot ~ (1|locality/site_irn/visit_no) + locality + water_speed_ms + water_depth + Cond + wmo_prec +
                       site_type + seas_wmo + bp_pres + bf_pres + bp_pres*bf_pres +
                       locality*seas_wmo + site_type*seas_wmo + 
                       water_depth*water_speed_ms +
                       offset(log(duration)),
                 data=snaildf,
                 family=nbinom2) 
drop1(mod_7, test = "Chisq")  # water_speed_ms:water_depth  1 8230.7  3.330 0.068028 . 

mod_8 <- glmmTMB(BT_tot ~ (1|locality/site_irn/visit_no) + locality + water_speed_ms + water_depth + Cond + wmo_prec +
                       site_type + seas_wmo + bp_pres + bf_pres + bp_pres*bf_pres +
                       locality*seas_wmo + site_type*seas_wmo +
                       offset(log(duration)),
                 data=snaildf,
                 family=nbinom2) 
drop1(mod_8, test = "Chisq")  # water_depth         1 8229.4  0.674 0.411710   

mod_9 <- glmer.nb(BT_tot ~ (1|locality/site_irn/visit_no) + locality + water_speed_ms + Cond + wmo_prec +
                       site_type + seas_wmo + bp_pres + bf_pres + bp_pres*bf_pres +
                       locality*seas_wmo + site_type*seas_wmo,
                       offset = log(duration),
                 data=snaildf)    # Model convergence problem; extreme or very small eigen values detected.
drop1(mod_9, test = "Chisq")  # Nothing to remove

# Since model convergence failed at mod_9, it might be best to use mod_7 for the summary
library(lme4)
summary(mod_7)
Anova.glmmTMB(mod_7)
TukeyHSD(mod_7)
emmeans(mod_7, ~ (locality | seas_wmo))
