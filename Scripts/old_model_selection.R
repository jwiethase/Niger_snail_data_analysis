
#' Test significance of predictors
# Use chisquare test on maximum model, to decide which term to remove first

drop1(model1, test = "Chisq")        # Temp_Water:site_type  7 10408  3.398 0.8458894  

model2 <- glmmTMB(Bulinus_tot ~ (1 |locality/site_irn/visit_no) + 
                        locality + pH + water_speed_ms + water_depth + water_level.v + Cond + wmo_prec +
                        Temp_Water + site_type + seas_wmo + 
                        locality*seas_wmo + site_type*seas_wmo + site_type*wmo_prec + 
                        site_type*Cond + site_type*pH +
                        offset(log(duration)),
                  data=bulinus_df,
                  family=nbinom2)
drop1(model2, test = "Chisq")        # pH:site_type        7 10399  5.205 0.6349038 

model3 <- glmmTMB(Bulinus_tot ~ (1 |locality/site_irn/visit_no) + 
                        locality + pH + water_speed_ms + water_depth + water_level.v + Cond + wmo_prec +
                        Temp_Water + site_type + seas_wmo + 
                        locality*seas_wmo + site_type*seas_wmo + site_type*wmo_prec + 
                        site_type*Cond +
                        offset(log(duration)),
                  data=bulinus_df,
                  family=nbinom2)
drop1(model3, test = "Chisq")        # pH                  1 10397  0.090 0.7635879   

model4 <- glmmTMB(Bulinus_tot ~ (1 |locality/site_irn/visit_no) + 
                        locality + water_speed_ms + water_depth + water_level.v + Cond + wmo_prec +
                        Temp_Water + site_type + seas_wmo + 
                        locality*seas_wmo + site_type*seas_wmo + site_type*wmo_prec + 
                        site_type*Cond +
                        offset(log(duration)),
                  data=bulinus_df,
                  family=nbinom2)
drop1(model4, test = "Chisq")        # water_depth         1 10396  1.202 0.2729720  

model5 <- glmmTMB(Bulinus_tot ~ (1 |locality/site_irn/visit_no) + 
                        locality + water_speed_ms + water_level.v + Cond + wmo_prec +
                        Temp_Water + site_type + seas_wmo +
                        locality*seas_wmo + site_type*seas_wmo + site_type*wmo_prec + 
                        site_type*Cond +
                        offset(log(duration)),
                  data=bulinus_df,
                  family=nbinom2)
drop1(model5, test = "Chisq")      # Temp_Water          1 10396  1.799  0.179807  
MASS::dropterm(model5, test = "Chisq")

model6 <- glmmTMB(Bulinus_tot ~ (1 |locality/site_irn/visit_no) + 
                        locality + water_speed_ms + water_level.v + Cond + wmo_prec +
                        site_type + seas_wmo +
                        locality*seas_wmo + site_type*seas_wmo + site_type*wmo_prec + 
                        site_type*Cond +
                        offset(log(duration)),
                  data=bulinus_df,
                  family=nbinom2)
drop1(model6, test = "Chisq")      # water_speed_ms      1 10634  2.7     0.10083    

model7 <- glmmTMB(Bulinus_tot ~ (1 |locality/site_irn/visit_no) +
                        locality + water_level.v + Cond + wmo_prec +
                        site_type + seas_wmo +
                        locality*seas_wmo + site_type*seas_wmo + site_type*wmo_prec + 
                        site_type*Cond +
                        offset(log(duration)),
                  data=bulinus_df,
                  family=nbinom2)
drop1(model7, test = "Chisq")

# Add month
model8 <- glmmTMB(Bulinus_tot ~ (1 |locality/site_irn/visit_no) +
                        locality + water_level.v + Cond + wmo_prec +
                        site_type + seas_wmo + month +
                        locality*seas_wmo + site_type*seas_wmo + site_type*wmo_prec + 
                        site_type*Cond +
                        offset(log(duration)),
                  data=bulinus_df,
                  family=nbinom2)
anova(model7, model8)  # model8 89 10533 11010  -5177    10355   129     11 <0.0000000000000002 ***

# Add Bulinus_pos_tot
model9 <- glmmTMB(Bulinus_tot ~ (1 |locality/site_irn/visit_no) +
                        locality + water_level.v + Cond + wmo_prec +
                        site_type + seas_wmo + month + Bulinus_pos_tot +
                        locality*seas_wmo + site_type*seas_wmo + site_type*wmo_prec + 
                        site_type*Cond +
                        offset(log(duration)),
                  data=bulinus_df,
                  family=nbinom2)
anova(model8, model9)  # model9 90 10494 10977  -5157    10314  40.2      1 0.00000000022 ***

# Add interaction of conductivity and bulinus positive
model10 <- glmmTMB(Bulinus_tot ~ (1 |locality/site_irn/visit_no) + 
                         locality + water_level.v + Cond + wmo_prec + 
                         site_type + seas_wmo + month + Bulinus_pos_tot + Bulinus_pos_tot*Cond +
                         locality*seas_wmo + site_type*seas_wmo + site_type*wmo_prec + 
                         site_type*Cond +
                         offset(log(duration)),
                   data=bulinus_df,
                   family=nbinom2)
anova(model9, model10)
drop1(model10, test = "Chisq")
summary(model10)

model11 <- glmmTMB(Bulinus_tot ~ (1 |locality/site_irn/visit_no) + 
                         locality + water_level.v + Cond + wmo_prec +
                         site_type + seas_wmo + month + Bulinus_pos_tot + Bulinus_pos_tot*Cond +
                         locality*seas_wmo + site_type*seas_wmo + site_type*wmo_prec + 
                         site_type*Cond + seas_wmo*water_level.v +
                         offset(log(duration)),
                   data=bulinus_df,
                   family=nbinom2)

#                          Df   AIC  LRT    Pr(>Chi)    
#       <none>                10232                     
#       water_level.v         2 10267 39.2 0.000000003 ***
#       locality:seas_wmo  18 10234 37.6      0.0044 ** 
#       site_type:seas_wmo  7 10241 23.3      0.0015 ** 
#       wmo_prec:site_type  7 10233 15.3      0.0329 *  
#       Cond:site_type      7 10238 20.2      0.0052 ** 


drop1(mod_1, test = "Chisq")   # pH:site_type          7 8449  1.2              0.98990    

mod_2 <- glmmTMB(BT_tot ~ (1 |locality/site_irn/visit_no) + locality + pH + water_speed_ms + 
                       water_depth + water_level.v + Cond + wmo_prec +
                       Temp_Water + site_type + seas_wmo + BP_tot + BF_tot + BT_pos_tot +
                       locality*seas_wmo + site_type*seas_wmo + site_type*Temp_Water + 
                       offset(log(duration)),
                 data=BT_df,
                 family=nbinom2)
drop1(mod_2, test = "Chisq")   # pH                    1 8446  0.0              0.86229  

mod_3 <- glmmTMB(BT_tot ~ (1 |locality/site_irn/visit_no) + locality + water_speed_ms + 
                       water_depth + water_level.v + Cond + wmo_prec +
                       Temp_Water + site_type + seas_wmo + BP_tot + BF_tot + BT_pos_tot +
                       locality*seas_wmo + site_type*seas_wmo + site_type*Temp_Water + 
                       offset(log(duration)),
                 data=BT_df,
                 family=nbinom2)
drop1(mod_3, test = "Chisq")  # BP_tot                1 8444  0.6               0.4361  

mod_4 <- glmmTMB(BT_tot ~ (1 |locality/site_irn/visit_no) + locality + water_speed_ms + 
                       water_depth + water_level.v + Cond + wmo_prec +
                       Temp_Water + site_type + seas_wmo + BF_tot + BT_pos_tot +
                       locality*seas_wmo + site_type*seas_wmo + site_type*Temp_Water + 
                       offset(log(duration)),
                 data=BT_df,
                 family=nbinom2)
drop1(mod_4, test = "Chisq") # nothing to remove

summary(mod_4)


drop1(mod_A, test = "Chisq")   # Cond                  1 6203  0.0  0.96744 

mod_B <- glmmTMB(BF_tot ~ (1 |locality/site_irn/visit_no) + 
                       locality + pH + water_speed_ms + water_depth + water_level.v + wmo_prec +
                       Temp_Water + site_type + seas_wmo + BP_tot + BT_tot + BP_tot*BT_tot + BF_pos_tot + 
                       locality*seas_wmo + site_type*seas_wmo + site_type*Temp_Water + site_type*pH +
                       offset(log(duration)),
                 data=chemdf,
                 family=nbinom2)
drop1(mod_B, test = "Chisq")   # wmo_prec              1 6196  0.2   0.6434 

mod_C <- glmmTMB(BF_tot ~ (1 |locality/site_irn/visit_no) + 
                       locality + pH + water_speed_ms + water_depth + water_level.v + 
                       Temp_Water + site_type + seas_wmo + BP_tot + BT_tot + BP_tot*BT_tot + BF_pos_tot + 
                       locality*seas_wmo + site_type*seas_wmo + site_type*Temp_Water + site_type*pH +
                       offset(log(duration)),
                 data=chemdf,
                 family=nbinom2)
drop1(mod_C, test = "Chisq")  # Temp_Water:site_type  7 6187  5.5  0.59771 

mod_D <-glmmTMB(BF_tot ~ (1 |locality/site_irn/visit_no) + 
                      locality + pH + water_speed_ms + water_depth + water_level.v + 
                      Temp_Water + site_type + seas_wmo + BP_tot + BT_tot + BP_tot*BT_tot + BF_pos_tot + 
                      locality*seas_wmo + site_type*seas_wmo + site_type*pH +
                      offset(log(duration)),
                data=chemdf,
                family=nbinom2)
drop1(mod_C, test = "Chisq") 

####FINISH