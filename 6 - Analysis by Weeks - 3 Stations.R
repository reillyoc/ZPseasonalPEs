sim_start <- read.csv("simcoe_manual_subset_weeks_3.csv", header = T)

library(tidyverse)
library(reshape2)

set.seed(018)

simcoe <- sim_start %>% select(Year, Week, Station_ID, Species, Density)

##### Shannon Diversity & Evenness #####

library(vegan)

simcoe_diversity <- simcoe %>% group_by(Year, Species, Station_ID) %>%
  summarize(Density = sum(Density))

simcoe_diversity <- dcast(simcoe_diversity, Year + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")


simcoe_diversity$Shannon_Div_Local <- diversity(simcoe_diversity[,c(3:50)])
simcoe_diversity$Evenness_Local <- diversity(simcoe_diversity[,c(3:50)])/log(specnumber(simcoe_diversity[,c(3:50)]))

simcoe_diversity <- simcoe_diversity %>% select(Year, Shannon_Div_Local, Evenness_Local) %>%
  group_by(Year) %>%
  summarize(Shannon_Div_Local = mean(Shannon_Div_Local),
            Evenness_Local = mean(Evenness_Local))


#Calculating regional diversity (# of metapopulations per year)
simcoe_richness <- simcoe %>% group_by(Year, Species) %>%
  summarize(Density = sum(Density))

simcoe_richness <- dcast(simcoe_richness, Year  ~ Species, fun.aggregate = sum, value.var = "Density")

simcoe_richness$Species_Richness <- specnumber(simcoe_richness[,c(2:49)])

simcoe_richness <- simcoe_richness %>% select(Year, Species_Richness)

detach("package:vegan", unload=TRUE)

simcoe_diversity <- merge(simcoe_diversity, simcoe_richness, by = c("Year"))

write.csv(simcoe_diversity, "simcoe_final_3_weeks_diversity")

##### Variability, Stabilization, & Asynchrony Calculations as per Hammond et al. 2020 #####


##### 1986 #####
simcoe_1986 <- subset(simcoe, Year == '1986')


#Transpose dataframe to add 0s within years
simcoe_1986_t1 <- dcast(simcoe_1986, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_1986_m <- melt(simcoe_1986_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_1986 <- subset(simcoe_1986_m, Station_ID == 'C9')
simcoe_c9_1986$Species <- paste("C9", simcoe_c9_1986$Species, sep=" ")
simcoe_k42_1986 <- subset(simcoe_1986_m, Station_ID == 'K42')
simcoe_k42_1986$Species <- paste("K42", simcoe_k42_1986$Species, sep=" ")
simcoe_k45_1986 <- subset(simcoe_1986_m, Station_ID == 'K45')
simcoe_k45_1986$Species <- paste("K45", simcoe_k45_1986$Species, sep=" ")

#recombine dataframes
simcoe_spec_1986 <- rbind(simcoe_c9_1986, simcoe_k42_1986, simcoe_k45_1986)
simcoe_spec_1986 <- simcoe_spec_1986 %>% select(- Station_ID)

#transpose dataframe
simcoe_1986_t2 <- dcast(simcoe_spec_1986, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_1986_t2 <- log(simcoe_1986_t2 + 1)
simcoe_1986_t2 <- simcoe_1986_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_1986_cv <- simcoe_1986_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_1986 = mean(Density),
            pop_sd_1986 = sd(Density)) %>%
  mutate(uw_pop_cv_1986 = pop_sd_1986/pop_mean_1986)

simcoe_1986_cv <- simcoe_1986_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_1986),
         mean_pop_cv = mean(uw_pop_cv_1986, na.rm = T),
         mean_pop_density = sum(pop_mean_1986),
         mean_pop_variance = sum(pop_sd_1986))

simcoe_1986_cv <- simcoe_1986_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_1986/mc_sum_mean_density)%>%
  mutate(w_pop_cv_1986 = (pop_mean_1986/mc_sum_mean_density)*uw_pop_cv_1986) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_1986_wa_pop_var_1986 <- simcoe_1986_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_1986 = sum(w_pop_cv_1986, na.rm = T),
            lcv = sum(w_pop_cv_1986, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_1986_cor_list_ts <- as.dist(round(cor(simcoe_1986_t2[]),2))
simcoe_1986_cor_ts <- stack(simcoe_1986_cor_list_ts, dim.names = TRUE)
simcoe_1986_cor_ts <- simcoe_1986_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_1986_cor - simcoe_1986_cv
simcoe_c9_1986_w_cv_ts <- subset(simcoe_1986_cv, Station_ID == 'C9')
simcoe_c9_1986_w_cv_ts$Species <- paste("C9", simcoe_c9_1986_w_cv_ts$Species, sep=" ")
simcoe_c9_1986_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_1986_w_cv_ts)
simcoe_k42_1986_w_cv_ts <- subset(simcoe_1986_cv, Station_ID == 'K42')
simcoe_k42_1986_w_cv_ts$Species <- paste("K42", simcoe_k42_1986_w_cv_ts$Species, sep=" ")
simcoe_k42_1986_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_1986_w_cv_ts)
simcoe_k45_1986_w_cv_ts <- subset(simcoe_1986_cv, Station_ID == 'K45')
simcoe_k45_1986_w_cv_ts$Species <- paste("K45", simcoe_k45_1986_w_cv_ts$Species, sep=" ")
simcoe_k45_1986_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_1986_w_cv_ts)

simcoe_spec_1986_w_sp_ts <- rbind(simcoe_c9_1986_w_cv_ts, simcoe_k42_1986_w_cv_ts, simcoe_k45_1986_w_cv_ts)
simcoe_spec_1986_w_sp_ts <- simcoe_spec_1986_w_sp_ts %>% select(Species, w_pop_cv_1986, Species_ID, Station_ID)
simcoe_spec_1986_w_sp1_ts <- simcoe_spec_1986_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_1986_sp1 = w_pop_cv_1986) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_1986_sp1, Species_ID, station_1)

simcoe_1986_cor_cv_ts <- merge(simcoe_1986_cor_ts, simcoe_spec_1986_w_sp1_ts, by = "species_1")

simcoe_spec_1986_w_sp2_ts <- simcoe_spec_1986_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_1986_sp2 = w_pop_cv_1986) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_1986_sp2, Species_ID, station_2)

simcoe_1986_cor_cv_ts <- merge(simcoe_1986_cor_cv_ts, simcoe_spec_1986_w_sp2_ts, by = "species_2")

simcoe_1986_cor_cv_ts_omit <- na.omit(simcoe_1986_cor_cv_ts)

simcoe_1986_ind_W_ts <- simcoe_1986_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_1986 = w_pop_cv_1986_sp1*w_pop_cv_1986_sp2,
         ind_W_1986 = (1 - corr)*(w_pop_cv_1986_sp1*w_pop_cv_1986_sp2),
         Year = 1986)

simcoe_1986_ind_W_ts$number_species_pairs <- nrow(simcoe_1986_ind_W_ts)

simcoe_mean_W_sp_pairs_1986 <- simcoe_1986_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_1986, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_1986_W_ts <- simcoe_1986_ind_W_ts %>% group_by(Year) %>%
  summarize(W_1986 = sum(ind_W_1986, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_1986_ind_W_ts_c9 <- filter(simcoe_1986_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_1986_ind_W_ts_k42 <- filter(simcoe_1986_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_1986_ind_W_ts_k45 <- filter(simcoe_1986_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_1986_c9 <- na.omit(simcoe_1986_ind_W_ts_c9)

simcoe_species_pairs_count_1986_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_1986_c9)

simcoe_mean_W_sp_pairs_1986_c9 <- simcoe_species_pairs_count_1986_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_1986, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1986_W_ts_c9 <- simcoe_1986_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_1986 = sum(ind_W_1986, na.rm = TRUE))

#K42
simcoe_species_pairs_count_1986_k42 <- na.omit(simcoe_1986_ind_W_ts_k42)

simcoe_species_pairs_count_1986_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_1986_k42)

simcoe_mean_W_sp_pairs_1986_k42 <- simcoe_species_pairs_count_1986_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_1986, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1986_W_ts_k42 <- simcoe_1986_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_1986 = sum(ind_W_1986, na.rm = TRUE))

#K45
simcoe_species_pairs_count_1986_k45 <- na.omit(simcoe_1986_ind_W_ts_k45)

simcoe_species_pairs_count_1986_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_1986_k45)

simcoe_mean_W_sp_pairs_1986_k45 <- simcoe_species_pairs_count_1986_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_1986, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1986_W_ts_k45 <- simcoe_1986_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_1986 = sum(ind_W_1986, na.rm = TRUE))

#Step 4
simcoe_1986_ind_W_as1 <- rbind(simcoe_1986_W_ts_c9, simcoe_1986_W_ts_k42, simcoe_1986_W_ts_k45)

simcoe_1986_W_as1 <- simcoe_1986_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_1986 = sum(W_1986, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_1986 <- rbind(simcoe_mean_W_sp_pairs_1986_c9, simcoe_mean_W_sp_pairs_1986_k42, simcoe_mean_W_sp_pairs_1986_k45)

simcoe_lc_sp_pairs_1986 <- simcoe_mean_lc_sp_pairs_1986 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_1986_ind_W_as4 <- filter(simcoe_1986_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_1986_mp  <- na.omit(simcoe_1986_ind_W_as4)

simcoe_species_pairs_count_1986_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_1986_mp)

simcoe_mean_mp_sp_pairs_1986  <- simcoe_species_pairs_count_1986_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_1986, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_1986_sum_W_as4 <- simcoe_1986_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_1986 = sum(ind_W_1986, na.rm = TRUE))


#Step 3
simcoe_1986_sum_W_as4 <- simcoe_1986_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_1986 = sum(as4_1986, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_1986_ind_W_as5 <- filter(simcoe_1986_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_1986_cc <- na.omit(simcoe_1986_ind_W_as5)

simcoe_species_pairs_count_1986_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_1986_cc)

simcoe_mean_cc_sp_pairs_1986  <- simcoe_species_pairs_count_1986_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_1986, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_1986_sum_W_as5 <- simcoe_1986_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_1986 = sum(ind_W_1986, na.rm = TRUE))


#Step 3
simcoe_1986_sum_W_as5 <- simcoe_1986_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_1986 = sum(as5_1986, na.rm = TRUE))

simcoe_1986_W <- merge(simcoe_1986_sum_W_as4, simcoe_1986_W_ts, by = c('Year'))
simcoe_1986_W <- merge(simcoe_1986_W, simcoe_1986_W_as1, by = c('Year'))
simcoe_1986_W <- merge(simcoe_1986_W, simcoe_1986_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_1986 <- merge(simcoe_mean_W_sp_pairs_1986, simcoe_lc_sp_pairs_1986, by=c("Year"))
simcoe_sp_pair_ave_stab_1986 <- merge(simcoe_sp_pair_ave_stab_1986, simcoe_mean_mp_sp_pairs_1986, by=c("Year"))
simcoe_sp_pair_ave_stab_1986 <- merge(simcoe_sp_pair_ave_stab_1986, simcoe_mean_cc_sp_pairs_1986, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_1986_final <- merge(simcoe_1986_W, simcoe_1986_wa_pop_var_1986, by = c('Year'))
simcoe_1986_final <- simcoe_1986_final %>% group_by(Year) %>%
  mutate(gcv_1986 = lcv - W_1986) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_1986, 
            lcv = lcv,
            W = W_1986,
            lc_stab = as2_1986,
            mp_stab = as4_1986,
            cc_stab = as5_1986,
            mc_cv = sqrt(gcv_1986),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_1986/lcv + as4_1986/lcv + as5_1986/lcv,
            lc_asynchrony = as2_1986/lcv,
            mp_asynchrony = as4_1986/lcv,
            cc_asynchrony= as5_1986/lcv)

simcoe_1986_cv <- simcoe_1986_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_1986_final$mc_sum_mean_density <- simcoe_1986_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 1986 #####
simcoe_mean_corr_1986_c9 <- simcoe_species_pairs_count_1986_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1986, ind_W_1986)
simcoe_mean_corr_1986_k42 <- simcoe_species_pairs_count_1986_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1986, ind_W_1986)
simcoe_mean_corr_1986_k45 <- simcoe_species_pairs_count_1986_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1986, ind_W_1986)
simcoe_mean_corr_1986_lc <- rbind(simcoe_mean_corr_1986_c9, simcoe_mean_corr_1986_k42, simcoe_mean_corr_1986_k45)
simcoe_mean_corr_1986_mp <- simcoe_species_pairs_count_1986_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_1986, ind_W_1986)
simcoe_mean_corr_1986_cc <- simcoe_species_pairs_count_1986_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_1986, ind_W_1986)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_1986_lc[sample(1:nrow(simcoe_mean_corr_1986_lc),  nrow(simcoe_mean_corr_1986_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_1986),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_1986),
              mean_lc_stab = mean(ind_W_1986),
              sd_lc_stab = sd(ind_W_1986))
  
  lc_corr_1986 <- as.data.frame(mean_lc)
}
lc_corr_1986

mp_corr_1986 <- simcoe_mean_corr_1986_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_1986),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_1986),
                                                       mean_mp_stab = mean(ind_W_1986),
                                                       sd_mp_stab = sd(ind_W_1986))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_1986_cc[sample(1:nrow(simcoe_mean_corr_1986_cc), nrow(simcoe_mean_corr_1986_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_1986),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_1986),
              mean_cc_stab = mean(ind_W_1986),
              sd_cc_stab = sd(ind_W_1986))
  
  cc_corr_1986 <- as.data.frame(mean_cc)
}
cc_corr_1986


corr_1986 <- cbind(lc_corr_1986, mp_corr_1986)
corr_1986 <- cbind(corr_1986, cc_corr_1986)
corr_1986 <- corr_1986 %>% mutate(Year = 1986)





##### 1987 #####
simcoe_1987 <- subset(simcoe, Year == '1987')


#Transpose dataframe to add 0s within years
simcoe_1987_t1 <- dcast(simcoe_1987, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_1987_m <- melt(simcoe_1987_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_1987 <- subset(simcoe_1987_m, Station_ID == 'C9')
simcoe_c9_1987$Species <- paste("C9", simcoe_c9_1987$Species, sep=" ")
simcoe_k42_1987 <- subset(simcoe_1987_m, Station_ID == 'K42')
simcoe_k42_1987$Species <- paste("K42", simcoe_k42_1987$Species, sep=" ")
simcoe_k45_1987 <- subset(simcoe_1987_m, Station_ID == 'K45')
simcoe_k45_1987$Species <- paste("K45", simcoe_k45_1987$Species, sep=" ")

#recombine dataframes
simcoe_spec_1987 <- rbind(simcoe_c9_1987, simcoe_k42_1987, simcoe_k45_1987)
simcoe_spec_1987 <- simcoe_spec_1987 %>% select(- Station_ID)

#transpose dataframe
simcoe_1987_t2 <- dcast(simcoe_spec_1987, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_1987_t2 <- log(simcoe_1987_t2 + 1)
simcoe_1987_t2 <- simcoe_1987_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_1987_cv <- simcoe_1987_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_1987 = mean(Density),
            pop_sd_1987 = sd(Density)) %>%
  mutate(uw_pop_cv_1987 = pop_sd_1987/pop_mean_1987)

simcoe_1987_cv <- simcoe_1987_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_1987),
         mean_pop_cv = mean(uw_pop_cv_1987, na.rm = T),
         mean_pop_density = sum(pop_mean_1987),
         mean_pop_variance = sum(pop_sd_1987))

simcoe_1987_cv <- simcoe_1987_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_1987/mc_sum_mean_density)%>%
  mutate(w_pop_cv_1987 = (pop_mean_1987/mc_sum_mean_density)*uw_pop_cv_1987) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_1987_wa_pop_var_1987 <- simcoe_1987_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_1987 = sum(w_pop_cv_1987, na.rm = T),
            lcv = sum(w_pop_cv_1987, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_1987_cor_list_ts <- as.dist(round(cor(simcoe_1987_t2[]),2))
simcoe_1987_cor_ts <- stack(simcoe_1987_cor_list_ts, dim.names = TRUE)
simcoe_1987_cor_ts <- simcoe_1987_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_1987_cor - simcoe_1987_cv
simcoe_c9_1987_w_cv_ts <- subset(simcoe_1987_cv, Station_ID == 'C9')
simcoe_c9_1987_w_cv_ts$Species <- paste("C9", simcoe_c9_1987_w_cv_ts$Species, sep=" ")
simcoe_c9_1987_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_1987_w_cv_ts)
simcoe_k42_1987_w_cv_ts <- subset(simcoe_1987_cv, Station_ID == 'K42')
simcoe_k42_1987_w_cv_ts$Species <- paste("K42", simcoe_k42_1987_w_cv_ts$Species, sep=" ")
simcoe_k42_1987_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_1987_w_cv_ts)
simcoe_k45_1987_w_cv_ts <- subset(simcoe_1987_cv, Station_ID == 'K45')
simcoe_k45_1987_w_cv_ts$Species <- paste("K45", simcoe_k45_1987_w_cv_ts$Species, sep=" ")
simcoe_k45_1987_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_1987_w_cv_ts)

simcoe_spec_1987_w_sp_ts <- rbind(simcoe_c9_1987_w_cv_ts, simcoe_k42_1987_w_cv_ts, simcoe_k45_1987_w_cv_ts)
simcoe_spec_1987_w_sp_ts <- simcoe_spec_1987_w_sp_ts %>% select(Species, w_pop_cv_1987, Species_ID, Station_ID)
simcoe_spec_1987_w_sp1_ts <- simcoe_spec_1987_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_1987_sp1 = w_pop_cv_1987) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_1987_sp1, Species_ID, station_1)

simcoe_1987_cor_cv_ts <- merge(simcoe_1987_cor_ts, simcoe_spec_1987_w_sp1_ts, by = "species_1")

simcoe_spec_1987_w_sp2_ts <- simcoe_spec_1987_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_1987_sp2 = w_pop_cv_1987) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_1987_sp2, Species_ID, station_2)

simcoe_1987_cor_cv_ts <- merge(simcoe_1987_cor_cv_ts, simcoe_spec_1987_w_sp2_ts, by = "species_2")

simcoe_1987_cor_cv_ts_omit <- na.omit(simcoe_1987_cor_cv_ts)

simcoe_1987_ind_W_ts <- simcoe_1987_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_1987 = w_pop_cv_1987_sp1*w_pop_cv_1987_sp2,
         ind_W_1987 = (1 - corr)*(w_pop_cv_1987_sp1*w_pop_cv_1987_sp2),
         Year = 1987)

simcoe_1987_ind_W_ts$number_species_pairs <- nrow(simcoe_1987_ind_W_ts)

simcoe_mean_W_sp_pairs_1987 <- simcoe_1987_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_1987, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_1987_W_ts <- simcoe_1987_ind_W_ts %>% group_by(Year) %>%
  summarize(W_1987 = sum(ind_W_1987, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_1987_ind_W_ts_c9 <- filter(simcoe_1987_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_1987_ind_W_ts_k42 <- filter(simcoe_1987_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_1987_ind_W_ts_k45 <- filter(simcoe_1987_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_1987_c9 <- na.omit(simcoe_1987_ind_W_ts_c9)

simcoe_species_pairs_count_1987_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_1987_c9)

simcoe_mean_W_sp_pairs_1987_c9 <- simcoe_species_pairs_count_1987_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_1987, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1987_W_ts_c9 <- simcoe_1987_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_1987 = sum(ind_W_1987, na.rm = TRUE))

#K42
simcoe_species_pairs_count_1987_k42 <- na.omit(simcoe_1987_ind_W_ts_k42)

simcoe_species_pairs_count_1987_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_1987_k42)

simcoe_mean_W_sp_pairs_1987_k42 <- simcoe_species_pairs_count_1987_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_1987, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1987_W_ts_k42 <- simcoe_1987_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_1987 = sum(ind_W_1987, na.rm = TRUE))

#K45
simcoe_species_pairs_count_1987_k45 <- na.omit(simcoe_1987_ind_W_ts_k45)

simcoe_species_pairs_count_1987_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_1987_k45)

simcoe_mean_W_sp_pairs_1987_k45 <- simcoe_species_pairs_count_1987_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_1987, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1987_W_ts_k45 <- simcoe_1987_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_1987 = sum(ind_W_1987, na.rm = TRUE))

#Step 4
simcoe_1987_ind_W_as1 <- rbind(simcoe_1987_W_ts_c9, simcoe_1987_W_ts_k42, simcoe_1987_W_ts_k45)

simcoe_1987_W_as1 <- simcoe_1987_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_1987 = sum(W_1987, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_1987 <- rbind(simcoe_mean_W_sp_pairs_1987_c9, simcoe_mean_W_sp_pairs_1987_k42, simcoe_mean_W_sp_pairs_1987_k45)

simcoe_lc_sp_pairs_1987 <- simcoe_mean_lc_sp_pairs_1987 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_1987_ind_W_as4 <- filter(simcoe_1987_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_1987_mp  <- na.omit(simcoe_1987_ind_W_as4)

simcoe_species_pairs_count_1987_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_1987_mp)

simcoe_mean_mp_sp_pairs_1987  <- simcoe_species_pairs_count_1987_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_1987, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_1987_sum_W_as4 <- simcoe_1987_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_1987 = sum(ind_W_1987, na.rm = TRUE))


#Step 3
simcoe_1987_sum_W_as4 <- simcoe_1987_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_1987 = sum(as4_1987, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_1987_ind_W_as5 <- filter(simcoe_1987_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_1987_cc <- na.omit(simcoe_1987_ind_W_as5)

simcoe_species_pairs_count_1987_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_1987_cc)

simcoe_mean_cc_sp_pairs_1987  <- simcoe_species_pairs_count_1987_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_1987, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_1987_sum_W_as5 <- simcoe_1987_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_1987 = sum(ind_W_1987, na.rm = TRUE))


#Step 3
simcoe_1987_sum_W_as5 <- simcoe_1987_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_1987 = sum(as5_1987, na.rm = TRUE))

simcoe_1987_W <- merge(simcoe_1987_sum_W_as4, simcoe_1987_W_ts, by = c('Year'))
simcoe_1987_W <- merge(simcoe_1987_W, simcoe_1987_W_as1, by = c('Year'))
simcoe_1987_W <- merge(simcoe_1987_W, simcoe_1987_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_1987 <- merge(simcoe_mean_W_sp_pairs_1987, simcoe_lc_sp_pairs_1987, by=c("Year"))
simcoe_sp_pair_ave_stab_1987 <- merge(simcoe_sp_pair_ave_stab_1987, simcoe_mean_mp_sp_pairs_1987, by=c("Year"))
simcoe_sp_pair_ave_stab_1987 <- merge(simcoe_sp_pair_ave_stab_1987, simcoe_mean_cc_sp_pairs_1987, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_1987_final <- merge(simcoe_1987_W, simcoe_1987_wa_pop_var_1987, by = c('Year'))
simcoe_1987_final <- simcoe_1987_final %>% group_by(Year) %>%
  mutate(gcv_1987 = lcv - W_1987) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_1987, 
            lcv = lcv,
            W = W_1987,
            lc_stab = as2_1987,
            mp_stab = as4_1987,
            cc_stab = as5_1987,
            mc_cv = sqrt(gcv_1987),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_1987/lcv + as4_1987/lcv + as5_1987/lcv,
            lc_asynchrony = as2_1987/lcv,
            mp_asynchrony = as4_1987/lcv,
            cc_asynchrony= as5_1987/lcv)

simcoe_1987_cv <- simcoe_1987_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_1987_final$mc_sum_mean_density <- simcoe_1987_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 1987 #####
simcoe_mean_corr_1987_c9 <- simcoe_species_pairs_count_1987_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1987, ind_W_1987)
simcoe_mean_corr_1987_k42 <- simcoe_species_pairs_count_1987_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1987, ind_W_1987)
simcoe_mean_corr_1987_k45 <- simcoe_species_pairs_count_1987_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1987, ind_W_1987)
simcoe_mean_corr_1987_lc <- rbind(simcoe_mean_corr_1987_c9, simcoe_mean_corr_1987_k42, simcoe_mean_corr_1987_k45)
simcoe_mean_corr_1987_mp <- simcoe_species_pairs_count_1987_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_1987, ind_W_1987)
simcoe_mean_corr_1987_cc <- simcoe_species_pairs_count_1987_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_1987, ind_W_1987)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_1987_lc[sample(1:nrow(simcoe_mean_corr_1987_lc),  nrow(simcoe_mean_corr_1987_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_1987),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_1987),
              mean_lc_stab = mean(ind_W_1987),
              sd_lc_stab = sd(ind_W_1987))
  
  lc_corr_1987 <- as.data.frame(mean_lc)
}
lc_corr_1987

mp_corr_1987 <- simcoe_mean_corr_1987_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_1987),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_1987),
                                                       mean_mp_stab = mean(ind_W_1987),
                                                       sd_mp_stab = sd(ind_W_1987))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_1987_cc[sample(1:nrow(simcoe_mean_corr_1987_cc), nrow(simcoe_mean_corr_1987_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_1987),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_1987),
              mean_cc_stab = mean(ind_W_1987),
              sd_cc_stab = sd(ind_W_1987))
  
  cc_corr_1987 <- as.data.frame(mean_cc)
}
cc_corr_1987


corr_1987 <- cbind(lc_corr_1987, mp_corr_1987)
corr_1987 <- cbind(corr_1987, cc_corr_1987)
corr_1987 <- corr_1987 %>% mutate(Year = 1987)






##### 1988 #####
simcoe_1988 <- subset(simcoe, Year == '1988')


#Transpose dataframe to add 0s within years
simcoe_1988_t1 <- dcast(simcoe_1988, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_1988_m <- melt(simcoe_1988_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_1988 <- subset(simcoe_1988_m, Station_ID == 'C9')
simcoe_c9_1988$Species <- paste("C9", simcoe_c9_1988$Species, sep=" ")
simcoe_k42_1988 <- subset(simcoe_1988_m, Station_ID == 'K42')
simcoe_k42_1988$Species <- paste("K42", simcoe_k42_1988$Species, sep=" ")
simcoe_k45_1988 <- subset(simcoe_1988_m, Station_ID == 'K45')
simcoe_k45_1988$Species <- paste("K45", simcoe_k45_1988$Species, sep=" ")

#recombine dataframes
simcoe_spec_1988 <- rbind(simcoe_c9_1988, simcoe_k42_1988, simcoe_k45_1988)
simcoe_spec_1988 <- simcoe_spec_1988 %>% select(- Station_ID)

#transpose dataframe
simcoe_1988_t2 <- dcast(simcoe_spec_1988, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_1988_t2 <- log(simcoe_1988_t2 + 1)
simcoe_1988_t2 <- simcoe_1988_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_1988_cv <- simcoe_1988_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_1988 = mean(Density),
            pop_sd_1988 = sd(Density)) %>%
  mutate(uw_pop_cv_1988 = pop_sd_1988/pop_mean_1988)

simcoe_1988_cv <- simcoe_1988_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_1988),
         mean_pop_cv = mean(uw_pop_cv_1988, na.rm = T),
         mean_pop_density = sum(pop_mean_1988),
         mean_pop_variance = sum(pop_sd_1988))

simcoe_1988_cv <- simcoe_1988_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_1988/mc_sum_mean_density)%>%
  mutate(w_pop_cv_1988 = (pop_mean_1988/mc_sum_mean_density)*uw_pop_cv_1988) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_1988_wa_pop_var_1988 <- simcoe_1988_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_1988 = sum(w_pop_cv_1988, na.rm = T),
            lcv = sum(w_pop_cv_1988, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_1988_cor_list_ts <- as.dist(round(cor(simcoe_1988_t2[]),2))
simcoe_1988_cor_ts <- stack(simcoe_1988_cor_list_ts, dim.names = TRUE)
simcoe_1988_cor_ts <- simcoe_1988_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_1988_cor - simcoe_1988_cv
simcoe_c9_1988_w_cv_ts <- subset(simcoe_1988_cv, Station_ID == 'C9')
simcoe_c9_1988_w_cv_ts$Species <- paste("C9", simcoe_c9_1988_w_cv_ts$Species, sep=" ")
simcoe_c9_1988_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_1988_w_cv_ts)
simcoe_k42_1988_w_cv_ts <- subset(simcoe_1988_cv, Station_ID == 'K42')
simcoe_k42_1988_w_cv_ts$Species <- paste("K42", simcoe_k42_1988_w_cv_ts$Species, sep=" ")
simcoe_k42_1988_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_1988_w_cv_ts)
simcoe_k45_1988_w_cv_ts <- subset(simcoe_1988_cv, Station_ID == 'K45')
simcoe_k45_1988_w_cv_ts$Species <- paste("K45", simcoe_k45_1988_w_cv_ts$Species, sep=" ")
simcoe_k45_1988_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_1988_w_cv_ts)

simcoe_spec_1988_w_sp_ts <- rbind(simcoe_c9_1988_w_cv_ts, simcoe_k42_1988_w_cv_ts, simcoe_k45_1988_w_cv_ts)
simcoe_spec_1988_w_sp_ts <- simcoe_spec_1988_w_sp_ts %>% select(Species, w_pop_cv_1988, Species_ID, Station_ID)
simcoe_spec_1988_w_sp1_ts <- simcoe_spec_1988_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_1988_sp1 = w_pop_cv_1988) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_1988_sp1, Species_ID, station_1)

simcoe_1988_cor_cv_ts <- merge(simcoe_1988_cor_ts, simcoe_spec_1988_w_sp1_ts, by = "species_1")

simcoe_spec_1988_w_sp2_ts <- simcoe_spec_1988_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_1988_sp2 = w_pop_cv_1988) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_1988_sp2, Species_ID, station_2)

simcoe_1988_cor_cv_ts <- merge(simcoe_1988_cor_cv_ts, simcoe_spec_1988_w_sp2_ts, by = "species_2")

simcoe_1988_cor_cv_ts_omit <- na.omit(simcoe_1988_cor_cv_ts)

simcoe_1988_ind_W_ts <- simcoe_1988_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_1988 = w_pop_cv_1988_sp1*w_pop_cv_1988_sp2,
         ind_W_1988 = (1 - corr)*(w_pop_cv_1988_sp1*w_pop_cv_1988_sp2),
         Year = 1988)

simcoe_1988_ind_W_ts$number_species_pairs <- nrow(simcoe_1988_ind_W_ts)

simcoe_mean_W_sp_pairs_1988 <- simcoe_1988_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_1988, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_1988_W_ts <- simcoe_1988_ind_W_ts %>% group_by(Year) %>%
  summarize(W_1988 = sum(ind_W_1988, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_1988_ind_W_ts_c9 <- filter(simcoe_1988_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_1988_ind_W_ts_k42 <- filter(simcoe_1988_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_1988_ind_W_ts_k45 <- filter(simcoe_1988_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_1988_c9 <- na.omit(simcoe_1988_ind_W_ts_c9)

simcoe_species_pairs_count_1988_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_1988_c9)

simcoe_mean_W_sp_pairs_1988_c9 <- simcoe_species_pairs_count_1988_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_1988, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1988_W_ts_c9 <- simcoe_1988_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_1988 = sum(ind_W_1988, na.rm = TRUE))

#K42
simcoe_species_pairs_count_1988_k42 <- na.omit(simcoe_1988_ind_W_ts_k42)

simcoe_species_pairs_count_1988_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_1988_k42)

simcoe_mean_W_sp_pairs_1988_k42 <- simcoe_species_pairs_count_1988_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_1988, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1988_W_ts_k42 <- simcoe_1988_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_1988 = sum(ind_W_1988, na.rm = TRUE))

#K45
simcoe_species_pairs_count_1988_k45 <- na.omit(simcoe_1988_ind_W_ts_k45)

simcoe_species_pairs_count_1988_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_1988_k45)

simcoe_mean_W_sp_pairs_1988_k45 <- simcoe_species_pairs_count_1988_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_1988, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1988_W_ts_k45 <- simcoe_1988_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_1988 = sum(ind_W_1988, na.rm = TRUE))

#Step 4
simcoe_1988_ind_W_as1 <- rbind(simcoe_1988_W_ts_c9, simcoe_1988_W_ts_k42, simcoe_1988_W_ts_k45)

simcoe_1988_W_as1 <- simcoe_1988_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_1988 = sum(W_1988, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_1988 <- rbind(simcoe_mean_W_sp_pairs_1988_c9, simcoe_mean_W_sp_pairs_1988_k42, simcoe_mean_W_sp_pairs_1988_k45)

simcoe_lc_sp_pairs_1988 <- simcoe_mean_lc_sp_pairs_1988 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_1988_ind_W_as4 <- filter(simcoe_1988_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_1988_mp  <- na.omit(simcoe_1988_ind_W_as4)

simcoe_species_pairs_count_1988_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_1988_mp)

simcoe_mean_mp_sp_pairs_1988  <- simcoe_species_pairs_count_1988_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_1988, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_1988_sum_W_as4 <- simcoe_1988_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_1988 = sum(ind_W_1988, na.rm = TRUE))


#Step 3
simcoe_1988_sum_W_as4 <- simcoe_1988_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_1988 = sum(as4_1988, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_1988_ind_W_as5 <- filter(simcoe_1988_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_1988_cc <- na.omit(simcoe_1988_ind_W_as5)

simcoe_species_pairs_count_1988_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_1988_cc)

simcoe_mean_cc_sp_pairs_1988  <- simcoe_species_pairs_count_1988_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_1988, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_1988_sum_W_as5 <- simcoe_1988_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_1988 = sum(ind_W_1988, na.rm = TRUE))


#Step 3
simcoe_1988_sum_W_as5 <- simcoe_1988_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_1988 = sum(as5_1988, na.rm = TRUE))

simcoe_1988_W <- merge(simcoe_1988_sum_W_as4, simcoe_1988_W_ts, by = c('Year'))
simcoe_1988_W <- merge(simcoe_1988_W, simcoe_1988_W_as1, by = c('Year'))
simcoe_1988_W <- merge(simcoe_1988_W, simcoe_1988_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_1988 <- merge(simcoe_mean_W_sp_pairs_1988, simcoe_lc_sp_pairs_1988, by=c("Year"))
simcoe_sp_pair_ave_stab_1988 <- merge(simcoe_sp_pair_ave_stab_1988, simcoe_mean_mp_sp_pairs_1988, by=c("Year"))
simcoe_sp_pair_ave_stab_1988 <- merge(simcoe_sp_pair_ave_stab_1988, simcoe_mean_cc_sp_pairs_1988, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_1988_final <- merge(simcoe_1988_W, simcoe_1988_wa_pop_var_1988, by = c('Year'))
simcoe_1988_final <- simcoe_1988_final %>% group_by(Year) %>%
  mutate(gcv_1988 = lcv - W_1988) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_1988, 
            lcv = lcv,
            W = W_1988,
            lc_stab = as2_1988,
            mp_stab = as4_1988,
            cc_stab = as5_1988,
            mc_cv = sqrt(gcv_1988),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_1988/lcv + as4_1988/lcv + as5_1988/lcv,
            lc_asynchrony = as2_1988/lcv,
            mp_asynchrony = as4_1988/lcv,
            cc_asynchrony= as5_1988/lcv)

simcoe_1988_cv <- simcoe_1988_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_1988_final$mc_sum_mean_density <- simcoe_1988_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 1988 #####
simcoe_mean_corr_1988_c9 <- simcoe_species_pairs_count_1988_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1988, ind_W_1988)
simcoe_mean_corr_1988_k42 <- simcoe_species_pairs_count_1988_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1988, ind_W_1988)
simcoe_mean_corr_1988_k45 <- simcoe_species_pairs_count_1988_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1988, ind_W_1988)
simcoe_mean_corr_1988_lc <- rbind(simcoe_mean_corr_1988_c9, simcoe_mean_corr_1988_k42, simcoe_mean_corr_1988_k45)
simcoe_mean_corr_1988_mp <- simcoe_species_pairs_count_1988_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_1988, ind_W_1988)
simcoe_mean_corr_1988_cc <- simcoe_species_pairs_count_1988_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_1988, ind_W_1988)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_1988_lc[sample(1:nrow(simcoe_mean_corr_1988_lc),  nrow(simcoe_mean_corr_1988_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_1988),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_1988),
              mean_lc_stab = mean(ind_W_1988),
              sd_lc_stab = sd(ind_W_1988))
  
  lc_corr_1988 <- as.data.frame(mean_lc)
}
lc_corr_1988

mp_corr_1988 <- simcoe_mean_corr_1988_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_1988),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_1988),
                                                       mean_mp_stab = mean(ind_W_1988),
                                                       sd_mp_stab = sd(ind_W_1988))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_1988_cc[sample(1:nrow(simcoe_mean_corr_1988_cc), nrow(simcoe_mean_corr_1988_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_1988),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_1988),
              mean_cc_stab = mean(ind_W_1988),
              sd_cc_stab = sd(ind_W_1988))
  
  cc_corr_1988 <- as.data.frame(mean_cc)
}
cc_corr_1988


corr_1988 <- cbind(lc_corr_1988, mp_corr_1988)
corr_1988 <- cbind(corr_1988, cc_corr_1988)
corr_1988 <- corr_1988 %>% mutate(Year = 1988)




##### 1989 #####
simcoe_1989 <- subset(simcoe, Year == '1989')


#Transpose dataframe to add 0s within years
simcoe_1989_t1 <- dcast(simcoe_1989, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_1989_m <- melt(simcoe_1989_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_1989 <- subset(simcoe_1989_m, Station_ID == 'C9')
simcoe_c9_1989$Species <- paste("C9", simcoe_c9_1989$Species, sep=" ")
simcoe_k42_1989 <- subset(simcoe_1989_m, Station_ID == 'K42')
simcoe_k42_1989$Species <- paste("K42", simcoe_k42_1989$Species, sep=" ")
simcoe_k45_1989 <- subset(simcoe_1989_m, Station_ID == 'K45')
simcoe_k45_1989$Species <- paste("K45", simcoe_k45_1989$Species, sep=" ")

#recombine dataframes
simcoe_spec_1989 <- rbind(simcoe_c9_1989, simcoe_k42_1989, simcoe_k45_1989)
simcoe_spec_1989 <- simcoe_spec_1989 %>% select(- Station_ID)

#transpose dataframe
simcoe_1989_t2 <- dcast(simcoe_spec_1989, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_1989_t2 <- log(simcoe_1989_t2 + 1)
simcoe_1989_t2 <- simcoe_1989_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_1989_cv <- simcoe_1989_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_1989 = mean(Density),
            pop_sd_1989 = sd(Density)) %>%
  mutate(uw_pop_cv_1989 = pop_sd_1989/pop_mean_1989)

simcoe_1989_cv <- simcoe_1989_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_1989),
         mean_pop_cv = mean(uw_pop_cv_1989, na.rm = T),
         mean_pop_density = sum(pop_mean_1989),
         mean_pop_variance = sum(pop_sd_1989))

simcoe_1989_cv <- simcoe_1989_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_1989/mc_sum_mean_density)%>%
  mutate(w_pop_cv_1989 = (pop_mean_1989/mc_sum_mean_density)*uw_pop_cv_1989) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_1989_wa_pop_var_1989 <- simcoe_1989_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_1989 = sum(w_pop_cv_1989, na.rm = T),
            lcv = sum(w_pop_cv_1989, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_1989_cor_list_ts <- as.dist(round(cor(simcoe_1989_t2[]),2))
simcoe_1989_cor_ts <- stack(simcoe_1989_cor_list_ts, dim.names = TRUE)
simcoe_1989_cor_ts <- simcoe_1989_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_1989_cor - simcoe_1989_cv
simcoe_c9_1989_w_cv_ts <- subset(simcoe_1989_cv, Station_ID == 'C9')
simcoe_c9_1989_w_cv_ts$Species <- paste("C9", simcoe_c9_1989_w_cv_ts$Species, sep=" ")
simcoe_c9_1989_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_1989_w_cv_ts)
simcoe_k42_1989_w_cv_ts <- subset(simcoe_1989_cv, Station_ID == 'K42')
simcoe_k42_1989_w_cv_ts$Species <- paste("K42", simcoe_k42_1989_w_cv_ts$Species, sep=" ")
simcoe_k42_1989_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_1989_w_cv_ts)
simcoe_k45_1989_w_cv_ts <- subset(simcoe_1989_cv, Station_ID == 'K45')
simcoe_k45_1989_w_cv_ts$Species <- paste("K45", simcoe_k45_1989_w_cv_ts$Species, sep=" ")
simcoe_k45_1989_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_1989_w_cv_ts)

simcoe_spec_1989_w_sp_ts <- rbind(simcoe_c9_1989_w_cv_ts, simcoe_k42_1989_w_cv_ts, simcoe_k45_1989_w_cv_ts)
simcoe_spec_1989_w_sp_ts <- simcoe_spec_1989_w_sp_ts %>% select(Species, w_pop_cv_1989, Species_ID, Station_ID)
simcoe_spec_1989_w_sp1_ts <- simcoe_spec_1989_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_1989_sp1 = w_pop_cv_1989) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_1989_sp1, Species_ID, station_1)

simcoe_1989_cor_cv_ts <- merge(simcoe_1989_cor_ts, simcoe_spec_1989_w_sp1_ts, by = "species_1")

simcoe_spec_1989_w_sp2_ts <- simcoe_spec_1989_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_1989_sp2 = w_pop_cv_1989) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_1989_sp2, Species_ID, station_2)

simcoe_1989_cor_cv_ts <- merge(simcoe_1989_cor_cv_ts, simcoe_spec_1989_w_sp2_ts, by = "species_2")

simcoe_1989_cor_cv_ts_omit <- na.omit(simcoe_1989_cor_cv_ts)

simcoe_1989_ind_W_ts <- simcoe_1989_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_1989 = w_pop_cv_1989_sp1*w_pop_cv_1989_sp2,
         ind_W_1989 = (1 - corr)*(w_pop_cv_1989_sp1*w_pop_cv_1989_sp2),
         Year = 1989)

simcoe_1989_ind_W_ts$number_species_pairs <- nrow(simcoe_1989_ind_W_ts)

simcoe_mean_W_sp_pairs_1989 <- simcoe_1989_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_1989, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_1989_W_ts <- simcoe_1989_ind_W_ts %>% group_by(Year) %>%
  summarize(W_1989 = sum(ind_W_1989, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_1989_ind_W_ts_c9 <- filter(simcoe_1989_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_1989_ind_W_ts_k42 <- filter(simcoe_1989_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_1989_ind_W_ts_k45 <- filter(simcoe_1989_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_1989_c9 <- na.omit(simcoe_1989_ind_W_ts_c9)

simcoe_species_pairs_count_1989_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_1989_c9)

simcoe_mean_W_sp_pairs_1989_c9 <- simcoe_species_pairs_count_1989_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_1989, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1989_W_ts_c9 <- simcoe_1989_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_1989 = sum(ind_W_1989, na.rm = TRUE))

#K42
simcoe_species_pairs_count_1989_k42 <- na.omit(simcoe_1989_ind_W_ts_k42)

simcoe_species_pairs_count_1989_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_1989_k42)

simcoe_mean_W_sp_pairs_1989_k42 <- simcoe_species_pairs_count_1989_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_1989, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1989_W_ts_k42 <- simcoe_1989_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_1989 = sum(ind_W_1989, na.rm = TRUE))

#K45
simcoe_species_pairs_count_1989_k45 <- na.omit(simcoe_1989_ind_W_ts_k45)

simcoe_species_pairs_count_1989_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_1989_k45)

simcoe_mean_W_sp_pairs_1989_k45 <- simcoe_species_pairs_count_1989_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_1989, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1989_W_ts_k45 <- simcoe_1989_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_1989 = sum(ind_W_1989, na.rm = TRUE))

#Step 4
simcoe_1989_ind_W_as1 <- rbind(simcoe_1989_W_ts_c9, simcoe_1989_W_ts_k42, simcoe_1989_W_ts_k45)

simcoe_1989_W_as1 <- simcoe_1989_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_1989 = sum(W_1989, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_1989 <- rbind(simcoe_mean_W_sp_pairs_1989_c9, simcoe_mean_W_sp_pairs_1989_k42, simcoe_mean_W_sp_pairs_1989_k45)

simcoe_lc_sp_pairs_1989 <- simcoe_mean_lc_sp_pairs_1989 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_1989_ind_W_as4 <- filter(simcoe_1989_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_1989_mp  <- na.omit(simcoe_1989_ind_W_as4)

simcoe_species_pairs_count_1989_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_1989_mp)

simcoe_mean_mp_sp_pairs_1989  <- simcoe_species_pairs_count_1989_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_1989, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_1989_sum_W_as4 <- simcoe_1989_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_1989 = sum(ind_W_1989, na.rm = TRUE))


#Step 3
simcoe_1989_sum_W_as4 <- simcoe_1989_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_1989 = sum(as4_1989, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_1989_ind_W_as5 <- filter(simcoe_1989_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_1989_cc <- na.omit(simcoe_1989_ind_W_as5)

simcoe_species_pairs_count_1989_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_1989_cc)

simcoe_mean_cc_sp_pairs_1989  <- simcoe_species_pairs_count_1989_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_1989, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_1989_sum_W_as5 <- simcoe_1989_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_1989 = sum(ind_W_1989, na.rm = TRUE))


#Step 3
simcoe_1989_sum_W_as5 <- simcoe_1989_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_1989 = sum(as5_1989, na.rm = TRUE))

simcoe_1989_W <- merge(simcoe_1989_sum_W_as4, simcoe_1989_W_ts, by = c('Year'))
simcoe_1989_W <- merge(simcoe_1989_W, simcoe_1989_W_as1, by = c('Year'))
simcoe_1989_W <- merge(simcoe_1989_W, simcoe_1989_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_1989 <- merge(simcoe_mean_W_sp_pairs_1989, simcoe_lc_sp_pairs_1989, by=c("Year"))
simcoe_sp_pair_ave_stab_1989 <- merge(simcoe_sp_pair_ave_stab_1989, simcoe_mean_mp_sp_pairs_1989, by=c("Year"))
simcoe_sp_pair_ave_stab_1989 <- merge(simcoe_sp_pair_ave_stab_1989, simcoe_mean_cc_sp_pairs_1989, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_1989_final <- merge(simcoe_1989_W, simcoe_1989_wa_pop_var_1989, by = c('Year'))
simcoe_1989_final <- simcoe_1989_final %>% group_by(Year) %>%
  mutate(gcv_1989 = lcv - W_1989) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_1989, 
            lcv = lcv,
            W = W_1989,
            lc_stab = as2_1989,
            mp_stab = as4_1989,
            cc_stab = as5_1989,
            mc_cv = sqrt(gcv_1989),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_1989/lcv + as4_1989/lcv + as5_1989/lcv,
            lc_asynchrony = as2_1989/lcv,
            mp_asynchrony = as4_1989/lcv,
            cc_asynchrony= as5_1989/lcv)

simcoe_1989_cv <- simcoe_1989_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_1989_final$mc_sum_mean_density <- simcoe_1989_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 1989 #####
simcoe_mean_corr_1989_c9 <- simcoe_species_pairs_count_1989_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1989, ind_W_1989)
simcoe_mean_corr_1989_k42 <- simcoe_species_pairs_count_1989_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1989, ind_W_1989)
simcoe_mean_corr_1989_k45 <- simcoe_species_pairs_count_1989_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1989, ind_W_1989)
simcoe_mean_corr_1989_lc <- rbind(simcoe_mean_corr_1989_c9, simcoe_mean_corr_1989_k42, simcoe_mean_corr_1989_k45)
simcoe_mean_corr_1989_mp <- simcoe_species_pairs_count_1989_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_1989, ind_W_1989)
simcoe_mean_corr_1989_cc <- simcoe_species_pairs_count_1989_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_1989, ind_W_1989)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_1989_lc[sample(1:nrow(simcoe_mean_corr_1989_lc),  nrow(simcoe_mean_corr_1989_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_1989),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_1989),
              mean_lc_stab = mean(ind_W_1989),
              sd_lc_stab = sd(ind_W_1989))
  
  lc_corr_1989 <- as.data.frame(mean_lc)
}
lc_corr_1989

mp_corr_1989 <- simcoe_mean_corr_1989_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_1989),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_1989),
                                                       mean_mp_stab = mean(ind_W_1989),
                                                       sd_mp_stab = sd(ind_W_1989))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_1989_cc[sample(1:nrow(simcoe_mean_corr_1989_cc), nrow(simcoe_mean_corr_1989_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_1989),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_1989),
              mean_cc_stab = mean(ind_W_1989),
              sd_cc_stab = sd(ind_W_1989))
  
  cc_corr_1989 <- as.data.frame(mean_cc)
}
cc_corr_1989


corr_1989 <- cbind(lc_corr_1989, mp_corr_1989)
corr_1989 <- cbind(corr_1989, cc_corr_1989)
corr_1989 <- corr_1989 %>% mutate(Year = 1989)






##### 1990 #####
simcoe_1990 <- subset(simcoe, Year == '1990')


#Transpose dataframe to add 0s within years
simcoe_1990_t1 <- dcast(simcoe_1990, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_1990_m <- melt(simcoe_1990_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_1990 <- subset(simcoe_1990_m, Station_ID == 'C9')
simcoe_c9_1990$Species <- paste("C9", simcoe_c9_1990$Species, sep=" ")
simcoe_k42_1990 <- subset(simcoe_1990_m, Station_ID == 'K42')
simcoe_k42_1990$Species <- paste("K42", simcoe_k42_1990$Species, sep=" ")
simcoe_k45_1990 <- subset(simcoe_1990_m, Station_ID == 'K45')
simcoe_k45_1990$Species <- paste("K45", simcoe_k45_1990$Species, sep=" ")

#recombine dataframes
simcoe_spec_1990 <- rbind(simcoe_c9_1990, simcoe_k42_1990, simcoe_k45_1990)
simcoe_spec_1990 <- simcoe_spec_1990 %>% select(- Station_ID)

#transpose dataframe
simcoe_1990_t2 <- dcast(simcoe_spec_1990, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_1990_t2 <- log(simcoe_1990_t2 + 1)
simcoe_1990_t2 <- simcoe_1990_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_1990_cv <- simcoe_1990_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_1990 = mean(Density),
            pop_sd_1990 = sd(Density)) %>%
  mutate(uw_pop_cv_1990 = pop_sd_1990/pop_mean_1990)

simcoe_1990_cv <- simcoe_1990_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_1990),
         mean_pop_cv = mean(uw_pop_cv_1990, na.rm = T),
         mean_pop_density = sum(pop_mean_1990),
         mean_pop_variance = sum(pop_sd_1990))

simcoe_1990_cv <- simcoe_1990_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_1990/mc_sum_mean_density)%>%
  mutate(w_pop_cv_1990 = (pop_mean_1990/mc_sum_mean_density)*uw_pop_cv_1990) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_1990_wa_pop_var_1990 <- simcoe_1990_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_1990 = sum(w_pop_cv_1990, na.rm = T),
            lcv = sum(w_pop_cv_1990, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_1990_cor_list_ts <- as.dist(round(cor(simcoe_1990_t2[]),2))
simcoe_1990_cor_ts <- stack(simcoe_1990_cor_list_ts, dim.names = TRUE)
simcoe_1990_cor_ts <- simcoe_1990_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_1990_cor - simcoe_1990_cv
simcoe_c9_1990_w_cv_ts <- subset(simcoe_1990_cv, Station_ID == 'C9')
simcoe_c9_1990_w_cv_ts$Species <- paste("C9", simcoe_c9_1990_w_cv_ts$Species, sep=" ")
simcoe_c9_1990_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_1990_w_cv_ts)
simcoe_k42_1990_w_cv_ts <- subset(simcoe_1990_cv, Station_ID == 'K42')
simcoe_k42_1990_w_cv_ts$Species <- paste("K42", simcoe_k42_1990_w_cv_ts$Species, sep=" ")
simcoe_k42_1990_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_1990_w_cv_ts)
simcoe_k45_1990_w_cv_ts <- subset(simcoe_1990_cv, Station_ID == 'K45')
simcoe_k45_1990_w_cv_ts$Species <- paste("K45", simcoe_k45_1990_w_cv_ts$Species, sep=" ")
simcoe_k45_1990_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_1990_w_cv_ts)

simcoe_spec_1990_w_sp_ts <- rbind(simcoe_c9_1990_w_cv_ts, simcoe_k42_1990_w_cv_ts, simcoe_k45_1990_w_cv_ts)
simcoe_spec_1990_w_sp_ts <- simcoe_spec_1990_w_sp_ts %>% select(Species, w_pop_cv_1990, Species_ID, Station_ID)
simcoe_spec_1990_w_sp1_ts <- simcoe_spec_1990_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_1990_sp1 = w_pop_cv_1990) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_1990_sp1, Species_ID, station_1)

simcoe_1990_cor_cv_ts <- merge(simcoe_1990_cor_ts, simcoe_spec_1990_w_sp1_ts, by = "species_1")

simcoe_spec_1990_w_sp2_ts <- simcoe_spec_1990_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_1990_sp2 = w_pop_cv_1990) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_1990_sp2, Species_ID, station_2)

simcoe_1990_cor_cv_ts <- merge(simcoe_1990_cor_cv_ts, simcoe_spec_1990_w_sp2_ts, by = "species_2")

simcoe_1990_cor_cv_ts_omit <- na.omit(simcoe_1990_cor_cv_ts)

simcoe_1990_ind_W_ts <- simcoe_1990_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_1990 = w_pop_cv_1990_sp1*w_pop_cv_1990_sp2,
         ind_W_1990 = (1 - corr)*(w_pop_cv_1990_sp1*w_pop_cv_1990_sp2),
         Year = 1990)

simcoe_1990_ind_W_ts$number_species_pairs <- nrow(simcoe_1990_ind_W_ts)

simcoe_mean_W_sp_pairs_1990 <- simcoe_1990_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_1990, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_1990_W_ts <- simcoe_1990_ind_W_ts %>% group_by(Year) %>%
  summarize(W_1990 = sum(ind_W_1990, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_1990_ind_W_ts_c9 <- filter(simcoe_1990_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_1990_ind_W_ts_k42 <- filter(simcoe_1990_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_1990_ind_W_ts_k45 <- filter(simcoe_1990_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_1990_c9 <- na.omit(simcoe_1990_ind_W_ts_c9)

simcoe_species_pairs_count_1990_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_1990_c9)

simcoe_mean_W_sp_pairs_1990_c9 <- simcoe_species_pairs_count_1990_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_1990, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1990_W_ts_c9 <- simcoe_1990_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_1990 = sum(ind_W_1990, na.rm = TRUE))

#K42
simcoe_species_pairs_count_1990_k42 <- na.omit(simcoe_1990_ind_W_ts_k42)

simcoe_species_pairs_count_1990_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_1990_k42)

simcoe_mean_W_sp_pairs_1990_k42 <- simcoe_species_pairs_count_1990_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_1990, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1990_W_ts_k42 <- simcoe_1990_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_1990 = sum(ind_W_1990, na.rm = TRUE))

#K45
simcoe_species_pairs_count_1990_k45 <- na.omit(simcoe_1990_ind_W_ts_k45)

simcoe_species_pairs_count_1990_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_1990_k45)

simcoe_mean_W_sp_pairs_1990_k45 <- simcoe_species_pairs_count_1990_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_1990, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1990_W_ts_k45 <- simcoe_1990_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_1990 = sum(ind_W_1990, na.rm = TRUE))

#Step 4
simcoe_1990_ind_W_as1 <- rbind(simcoe_1990_W_ts_c9, simcoe_1990_W_ts_k42, simcoe_1990_W_ts_k45)

simcoe_1990_W_as1 <- simcoe_1990_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_1990 = sum(W_1990, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_1990 <- rbind(simcoe_mean_W_sp_pairs_1990_c9, simcoe_mean_W_sp_pairs_1990_k42, simcoe_mean_W_sp_pairs_1990_k45)

simcoe_lc_sp_pairs_1990 <- simcoe_mean_lc_sp_pairs_1990 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_1990_ind_W_as4 <- filter(simcoe_1990_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_1990_mp  <- na.omit(simcoe_1990_ind_W_as4)

simcoe_species_pairs_count_1990_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_1990_mp)

simcoe_mean_mp_sp_pairs_1990  <- simcoe_species_pairs_count_1990_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_1990, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_1990_sum_W_as4 <- simcoe_1990_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_1990 = sum(ind_W_1990, na.rm = TRUE))


#Step 3
simcoe_1990_sum_W_as4 <- simcoe_1990_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_1990 = sum(as4_1990, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_1990_ind_W_as5 <- filter(simcoe_1990_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_1990_cc <- na.omit(simcoe_1990_ind_W_as5)

simcoe_species_pairs_count_1990_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_1990_cc)

simcoe_mean_cc_sp_pairs_1990  <- simcoe_species_pairs_count_1990_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_1990, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_1990_sum_W_as5 <- simcoe_1990_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_1990 = sum(ind_W_1990, na.rm = TRUE))


#Step 3
simcoe_1990_sum_W_as5 <- simcoe_1990_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_1990 = sum(as5_1990, na.rm = TRUE))

simcoe_1990_W <- merge(simcoe_1990_sum_W_as4, simcoe_1990_W_ts, by = c('Year'))
simcoe_1990_W <- merge(simcoe_1990_W, simcoe_1990_W_as1, by = c('Year'))
simcoe_1990_W <- merge(simcoe_1990_W, simcoe_1990_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_1990 <- merge(simcoe_mean_W_sp_pairs_1990, simcoe_lc_sp_pairs_1990, by=c("Year"))
simcoe_sp_pair_ave_stab_1990 <- merge(simcoe_sp_pair_ave_stab_1990, simcoe_mean_mp_sp_pairs_1990, by=c("Year"))
simcoe_sp_pair_ave_stab_1990 <- merge(simcoe_sp_pair_ave_stab_1990, simcoe_mean_cc_sp_pairs_1990, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_1990_final <- merge(simcoe_1990_W, simcoe_1990_wa_pop_var_1990, by = c('Year'))
simcoe_1990_final <- simcoe_1990_final %>% group_by(Year) %>%
  mutate(gcv_1990 = lcv - W_1990) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_1990, 
            lcv = lcv,
            W = W_1990,
            lc_stab = as2_1990,
            mp_stab = as4_1990,
            cc_stab = as5_1990,
            mc_cv = sqrt(gcv_1990),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_1990/lcv + as4_1990/lcv + as5_1990/lcv,
            lc_asynchrony = as2_1990/lcv,
            mp_asynchrony = as4_1990/lcv,
            cc_asynchrony= as5_1990/lcv)

simcoe_1990_cv <- simcoe_1990_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_1990_final$mc_sum_mean_density <- simcoe_1990_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 1990 #####
simcoe_mean_corr_1990_c9 <- simcoe_species_pairs_count_1990_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1990, ind_W_1990)
simcoe_mean_corr_1990_k42 <- simcoe_species_pairs_count_1990_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1990, ind_W_1990)
simcoe_mean_corr_1990_k45 <- simcoe_species_pairs_count_1990_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1990, ind_W_1990)
simcoe_mean_corr_1990_lc <- rbind(simcoe_mean_corr_1990_c9, simcoe_mean_corr_1990_k42, simcoe_mean_corr_1990_k45)
simcoe_mean_corr_1990_mp <- simcoe_species_pairs_count_1990_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_1990, ind_W_1990)
simcoe_mean_corr_1990_cc <- simcoe_species_pairs_count_1990_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_1990, ind_W_1990)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_1990_lc[sample(1:nrow(simcoe_mean_corr_1990_lc),  nrow(simcoe_mean_corr_1990_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_1990),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_1990),
              mean_lc_stab = mean(ind_W_1990),
              sd_lc_stab = sd(ind_W_1990))
  
  lc_corr_1990 <- as.data.frame(mean_lc)
}
lc_corr_1990

mp_corr_1990 <- simcoe_mean_corr_1990_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_1990),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_1990),
                                                       mean_mp_stab = mean(ind_W_1990),
                                                       sd_mp_stab = sd(ind_W_1990))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_1990_cc[sample(1:nrow(simcoe_mean_corr_1990_cc), nrow(simcoe_mean_corr_1990_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_1990),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_1990),
              mean_cc_stab = mean(ind_W_1990),
              sd_cc_stab = sd(ind_W_1990))
  
  cc_corr_1990 <- as.data.frame(mean_cc)
}
cc_corr_1990


corr_1990 <- cbind(lc_corr_1990, mp_corr_1990)
corr_1990 <- cbind(corr_1990, cc_corr_1990)
corr_1990 <- corr_1990 %>% mutate(Year = 1990)





##### 1991 #####
simcoe_1991 <- subset(simcoe, Year == '1991')


#Transpose dataframe to add 0s within years
simcoe_1991_t1 <- dcast(simcoe_1991, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_1991_m <- melt(simcoe_1991_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_1991 <- subset(simcoe_1991_m, Station_ID == 'C9')
simcoe_c9_1991$Species <- paste("C9", simcoe_c9_1991$Species, sep=" ")
simcoe_k42_1991 <- subset(simcoe_1991_m, Station_ID == 'K42')
simcoe_k42_1991$Species <- paste("K42", simcoe_k42_1991$Species, sep=" ")
simcoe_k45_1991 <- subset(simcoe_1991_m, Station_ID == 'K45')
simcoe_k45_1991$Species <- paste("K45", simcoe_k45_1991$Species, sep=" ")

#recombine dataframes
simcoe_spec_1991 <- rbind(simcoe_c9_1991, simcoe_k42_1991, simcoe_k45_1991)
simcoe_spec_1991 <- simcoe_spec_1991 %>% select(- Station_ID)

#transpose dataframe
simcoe_1991_t2 <- dcast(simcoe_spec_1991, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_1991_t2 <- log(simcoe_1991_t2 + 1)
simcoe_1991_t2 <- simcoe_1991_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_1991_cv <- simcoe_1991_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_1991 = mean(Density),
            pop_sd_1991 = sd(Density)) %>%
  mutate(uw_pop_cv_1991 = pop_sd_1991/pop_mean_1991)

simcoe_1991_cv <- simcoe_1991_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_1991),
         mean_pop_cv = mean(uw_pop_cv_1991, na.rm = T),
         mean_pop_density = sum(pop_mean_1991),
         mean_pop_variance = sum(pop_sd_1991))

simcoe_1991_cv <- simcoe_1991_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_1991/mc_sum_mean_density)%>%
  mutate(w_pop_cv_1991 = (pop_mean_1991/mc_sum_mean_density)*uw_pop_cv_1991) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_1991_wa_pop_var_1991 <- simcoe_1991_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_1991 = sum(w_pop_cv_1991, na.rm = T),
            lcv = sum(w_pop_cv_1991, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_1991_cor_list_ts <- as.dist(round(cor(simcoe_1991_t2[]),2))
simcoe_1991_cor_ts <- stack(simcoe_1991_cor_list_ts, dim.names = TRUE)
simcoe_1991_cor_ts <- simcoe_1991_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_1991_cor - simcoe_1991_cv
simcoe_c9_1991_w_cv_ts <- subset(simcoe_1991_cv, Station_ID == 'C9')
simcoe_c9_1991_w_cv_ts$Species <- paste("C9", simcoe_c9_1991_w_cv_ts$Species, sep=" ")
simcoe_c9_1991_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_1991_w_cv_ts)
simcoe_k42_1991_w_cv_ts <- subset(simcoe_1991_cv, Station_ID == 'K42')
simcoe_k42_1991_w_cv_ts$Species <- paste("K42", simcoe_k42_1991_w_cv_ts$Species, sep=" ")
simcoe_k42_1991_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_1991_w_cv_ts)
simcoe_k45_1991_w_cv_ts <- subset(simcoe_1991_cv, Station_ID == 'K45')
simcoe_k45_1991_w_cv_ts$Species <- paste("K45", simcoe_k45_1991_w_cv_ts$Species, sep=" ")
simcoe_k45_1991_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_1991_w_cv_ts)

simcoe_spec_1991_w_sp_ts <- rbind(simcoe_c9_1991_w_cv_ts, simcoe_k42_1991_w_cv_ts, simcoe_k45_1991_w_cv_ts)
simcoe_spec_1991_w_sp_ts <- simcoe_spec_1991_w_sp_ts %>% select(Species, w_pop_cv_1991, Species_ID, Station_ID)
simcoe_spec_1991_w_sp1_ts <- simcoe_spec_1991_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_1991_sp1 = w_pop_cv_1991) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_1991_sp1, Species_ID, station_1)

simcoe_1991_cor_cv_ts <- merge(simcoe_1991_cor_ts, simcoe_spec_1991_w_sp1_ts, by = "species_1")

simcoe_spec_1991_w_sp2_ts <- simcoe_spec_1991_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_1991_sp2 = w_pop_cv_1991) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_1991_sp2, Species_ID, station_2)

simcoe_1991_cor_cv_ts <- merge(simcoe_1991_cor_cv_ts, simcoe_spec_1991_w_sp2_ts, by = "species_2")

simcoe_1991_cor_cv_ts_omit <- na.omit(simcoe_1991_cor_cv_ts)

simcoe_1991_ind_W_ts <- simcoe_1991_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_1991 = w_pop_cv_1991_sp1*w_pop_cv_1991_sp2,
         ind_W_1991 = (1 - corr)*(w_pop_cv_1991_sp1*w_pop_cv_1991_sp2),
         Year = 1991)

simcoe_1991_ind_W_ts$number_species_pairs <- nrow(simcoe_1991_ind_W_ts)

simcoe_mean_W_sp_pairs_1991 <- simcoe_1991_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_1991, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_1991_W_ts <- simcoe_1991_ind_W_ts %>% group_by(Year) %>%
  summarize(W_1991 = sum(ind_W_1991, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_1991_ind_W_ts_c9 <- filter(simcoe_1991_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_1991_ind_W_ts_k42 <- filter(simcoe_1991_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_1991_ind_W_ts_k45 <- filter(simcoe_1991_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_1991_c9 <- na.omit(simcoe_1991_ind_W_ts_c9)

simcoe_species_pairs_count_1991_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_1991_c9)

simcoe_mean_W_sp_pairs_1991_c9 <- simcoe_species_pairs_count_1991_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_1991, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1991_W_ts_c9 <- simcoe_1991_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_1991 = sum(ind_W_1991, na.rm = TRUE))

#K42
simcoe_species_pairs_count_1991_k42 <- na.omit(simcoe_1991_ind_W_ts_k42)

simcoe_species_pairs_count_1991_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_1991_k42)

simcoe_mean_W_sp_pairs_1991_k42 <- simcoe_species_pairs_count_1991_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_1991, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1991_W_ts_k42 <- simcoe_1991_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_1991 = sum(ind_W_1991, na.rm = TRUE))

#K45
simcoe_species_pairs_count_1991_k45 <- na.omit(simcoe_1991_ind_W_ts_k45)

simcoe_species_pairs_count_1991_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_1991_k45)

simcoe_mean_W_sp_pairs_1991_k45 <- simcoe_species_pairs_count_1991_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_1991, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1991_W_ts_k45 <- simcoe_1991_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_1991 = sum(ind_W_1991, na.rm = TRUE))

#Step 4
simcoe_1991_ind_W_as1 <- rbind(simcoe_1991_W_ts_c9, simcoe_1991_W_ts_k42, simcoe_1991_W_ts_k45)

simcoe_1991_W_as1 <- simcoe_1991_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_1991 = sum(W_1991, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_1991 <- rbind(simcoe_mean_W_sp_pairs_1991_c9, simcoe_mean_W_sp_pairs_1991_k42, simcoe_mean_W_sp_pairs_1991_k45)

simcoe_lc_sp_pairs_1991 <- simcoe_mean_lc_sp_pairs_1991 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_1991_ind_W_as4 <- filter(simcoe_1991_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_1991_mp  <- na.omit(simcoe_1991_ind_W_as4)

simcoe_species_pairs_count_1991_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_1991_mp)

simcoe_mean_mp_sp_pairs_1991  <- simcoe_species_pairs_count_1991_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_1991, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_1991_sum_W_as4 <- simcoe_1991_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_1991 = sum(ind_W_1991, na.rm = TRUE))


#Step 3
simcoe_1991_sum_W_as4 <- simcoe_1991_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_1991 = sum(as4_1991, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_1991_ind_W_as5 <- filter(simcoe_1991_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_1991_cc <- na.omit(simcoe_1991_ind_W_as5)

simcoe_species_pairs_count_1991_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_1991_cc)

simcoe_mean_cc_sp_pairs_1991  <- simcoe_species_pairs_count_1991_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_1991, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_1991_sum_W_as5 <- simcoe_1991_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_1991 = sum(ind_W_1991, na.rm = TRUE))


#Step 3
simcoe_1991_sum_W_as5 <- simcoe_1991_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_1991 = sum(as5_1991, na.rm = TRUE))

simcoe_1991_W <- merge(simcoe_1991_sum_W_as4, simcoe_1991_W_ts, by = c('Year'))
simcoe_1991_W <- merge(simcoe_1991_W, simcoe_1991_W_as1, by = c('Year'))
simcoe_1991_W <- merge(simcoe_1991_W, simcoe_1991_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_1991 <- merge(simcoe_mean_W_sp_pairs_1991, simcoe_lc_sp_pairs_1991, by=c("Year"))
simcoe_sp_pair_ave_stab_1991 <- merge(simcoe_sp_pair_ave_stab_1991, simcoe_mean_mp_sp_pairs_1991, by=c("Year"))
simcoe_sp_pair_ave_stab_1991 <- merge(simcoe_sp_pair_ave_stab_1991, simcoe_mean_cc_sp_pairs_1991, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_1991_final <- merge(simcoe_1991_W, simcoe_1991_wa_pop_var_1991, by = c('Year'))
simcoe_1991_final <- simcoe_1991_final %>% group_by(Year) %>%
  mutate(gcv_1991 = lcv - W_1991) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_1991, 
            lcv = lcv,
            W = W_1991,
            lc_stab = as2_1991,
            mp_stab = as4_1991,
            cc_stab = as5_1991,
            mc_cv = sqrt(gcv_1991),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_1991/lcv + as4_1991/lcv + as5_1991/lcv,
            lc_asynchrony = as2_1991/lcv,
            mp_asynchrony = as4_1991/lcv,
            cc_asynchrony= as5_1991/lcv)

simcoe_1991_cv <- simcoe_1991_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_1991_final$mc_sum_mean_density <- simcoe_1991_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 1991 #####
simcoe_mean_corr_1991_c9 <- simcoe_species_pairs_count_1991_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1991, ind_W_1991)
simcoe_mean_corr_1991_k42 <- simcoe_species_pairs_count_1991_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1991, ind_W_1991)
simcoe_mean_corr_1991_k45 <- simcoe_species_pairs_count_1991_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1991, ind_W_1991)
simcoe_mean_corr_1991_lc <- rbind(simcoe_mean_corr_1991_c9, simcoe_mean_corr_1991_k42, simcoe_mean_corr_1991_k45)
simcoe_mean_corr_1991_mp <- simcoe_species_pairs_count_1991_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_1991, ind_W_1991)
simcoe_mean_corr_1991_cc <- simcoe_species_pairs_count_1991_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_1991, ind_W_1991)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_1991_lc[sample(1:nrow(simcoe_mean_corr_1991_lc),  nrow(simcoe_mean_corr_1991_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_1991),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_1991),
              mean_lc_stab = mean(ind_W_1991),
              sd_lc_stab = sd(ind_W_1991))
  
  lc_corr_1991 <- as.data.frame(mean_lc)
}
lc_corr_1991

mp_corr_1991 <- simcoe_mean_corr_1991_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_1991),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_1991),
                                                       mean_mp_stab = mean(ind_W_1991),
                                                       sd_mp_stab = sd(ind_W_1991))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_1991_cc[sample(1:nrow(simcoe_mean_corr_1991_cc), nrow(simcoe_mean_corr_1991_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_1991),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_1991),
              mean_cc_stab = mean(ind_W_1991),
              sd_cc_stab = sd(ind_W_1991))
  
  cc_corr_1991 <- as.data.frame(mean_cc)
}
cc_corr_1991


corr_1991 <- cbind(lc_corr_1991, mp_corr_1991)
corr_1991 <- cbind(corr_1991, cc_corr_1991)
corr_1991 <- corr_1991 %>% mutate(Year = 1991)





##### 1992 #####
simcoe_1992 <- subset(simcoe, Year == '1992')


#Transpose dataframe to add 0s within years
simcoe_1992_t1 <- dcast(simcoe_1992, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_1992_m <- melt(simcoe_1992_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_1992 <- subset(simcoe_1992_m, Station_ID == 'C9')
simcoe_c9_1992$Species <- paste("C9", simcoe_c9_1992$Species, sep=" ")
simcoe_k42_1992 <- subset(simcoe_1992_m, Station_ID == 'K42')
simcoe_k42_1992$Species <- paste("K42", simcoe_k42_1992$Species, sep=" ")
simcoe_k45_1992 <- subset(simcoe_1992_m, Station_ID == 'K45')
simcoe_k45_1992$Species <- paste("K45", simcoe_k45_1992$Species, sep=" ")

#recombine dataframes
simcoe_spec_1992 <- rbind(simcoe_c9_1992, simcoe_k42_1992, simcoe_k45_1992)
simcoe_spec_1992 <- simcoe_spec_1992 %>% select(- Station_ID)

#transpose dataframe
simcoe_1992_t2 <- dcast(simcoe_spec_1992, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_1992_t2 <- log(simcoe_1992_t2 + 1)
simcoe_1992_t2 <- simcoe_1992_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_1992_cv <- simcoe_1992_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_1992 = mean(Density),
            pop_sd_1992 = sd(Density)) %>%
  mutate(uw_pop_cv_1992 = pop_sd_1992/pop_mean_1992)

simcoe_1992_cv <- simcoe_1992_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_1992),
         mean_pop_cv = mean(uw_pop_cv_1992, na.rm = T),
         mean_pop_density = sum(pop_mean_1992),
         mean_pop_variance = sum(pop_sd_1992))

simcoe_1992_cv <- simcoe_1992_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_1992/mc_sum_mean_density)%>%
  mutate(w_pop_cv_1992 = (pop_mean_1992/mc_sum_mean_density)*uw_pop_cv_1992) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_1992_wa_pop_var_1992 <- simcoe_1992_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_1992 = sum(w_pop_cv_1992, na.rm = T),
            lcv = sum(w_pop_cv_1992, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_1992_cor_list_ts <- as.dist(round(cor(simcoe_1992_t2[]),2))
simcoe_1992_cor_ts <- stack(simcoe_1992_cor_list_ts, dim.names = TRUE)
simcoe_1992_cor_ts <- simcoe_1992_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_1992_cor - simcoe_1992_cv
simcoe_c9_1992_w_cv_ts <- subset(simcoe_1992_cv, Station_ID == 'C9')
simcoe_c9_1992_w_cv_ts$Species <- paste("C9", simcoe_c9_1992_w_cv_ts$Species, sep=" ")
simcoe_c9_1992_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_1992_w_cv_ts)
simcoe_k42_1992_w_cv_ts <- subset(simcoe_1992_cv, Station_ID == 'K42')
simcoe_k42_1992_w_cv_ts$Species <- paste("K42", simcoe_k42_1992_w_cv_ts$Species, sep=" ")
simcoe_k42_1992_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_1992_w_cv_ts)
simcoe_k45_1992_w_cv_ts <- subset(simcoe_1992_cv, Station_ID == 'K45')
simcoe_k45_1992_w_cv_ts$Species <- paste("K45", simcoe_k45_1992_w_cv_ts$Species, sep=" ")
simcoe_k45_1992_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_1992_w_cv_ts)

simcoe_spec_1992_w_sp_ts <- rbind(simcoe_c9_1992_w_cv_ts, simcoe_k42_1992_w_cv_ts, simcoe_k45_1992_w_cv_ts)
simcoe_spec_1992_w_sp_ts <- simcoe_spec_1992_w_sp_ts %>% select(Species, w_pop_cv_1992, Species_ID, Station_ID)
simcoe_spec_1992_w_sp1_ts <- simcoe_spec_1992_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_1992_sp1 = w_pop_cv_1992) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_1992_sp1, Species_ID, station_1)

simcoe_1992_cor_cv_ts <- merge(simcoe_1992_cor_ts, simcoe_spec_1992_w_sp1_ts, by = "species_1")

simcoe_spec_1992_w_sp2_ts <- simcoe_spec_1992_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_1992_sp2 = w_pop_cv_1992) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_1992_sp2, Species_ID, station_2)

simcoe_1992_cor_cv_ts <- merge(simcoe_1992_cor_cv_ts, simcoe_spec_1992_w_sp2_ts, by = "species_2")

simcoe_1992_cor_cv_ts_omit <- na.omit(simcoe_1992_cor_cv_ts)

simcoe_1992_ind_W_ts <- simcoe_1992_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_1992 = w_pop_cv_1992_sp1*w_pop_cv_1992_sp2,
         ind_W_1992 = (1 - corr)*(w_pop_cv_1992_sp1*w_pop_cv_1992_sp2),
         Year = 1992)

simcoe_1992_ind_W_ts$number_species_pairs <- nrow(simcoe_1992_ind_W_ts)

simcoe_mean_W_sp_pairs_1992 <- simcoe_1992_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_1992, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_1992_W_ts <- simcoe_1992_ind_W_ts %>% group_by(Year) %>%
  summarize(W_1992 = sum(ind_W_1992, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_1992_ind_W_ts_c9 <- filter(simcoe_1992_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_1992_ind_W_ts_k42 <- filter(simcoe_1992_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_1992_ind_W_ts_k45 <- filter(simcoe_1992_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_1992_c9 <- na.omit(simcoe_1992_ind_W_ts_c9)

simcoe_species_pairs_count_1992_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_1992_c9)

simcoe_mean_W_sp_pairs_1992_c9 <- simcoe_species_pairs_count_1992_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_1992, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1992_W_ts_c9 <- simcoe_1992_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_1992 = sum(ind_W_1992, na.rm = TRUE))

#K42
simcoe_species_pairs_count_1992_k42 <- na.omit(simcoe_1992_ind_W_ts_k42)

simcoe_species_pairs_count_1992_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_1992_k42)

simcoe_mean_W_sp_pairs_1992_k42 <- simcoe_species_pairs_count_1992_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_1992, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1992_W_ts_k42 <- simcoe_1992_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_1992 = sum(ind_W_1992, na.rm = TRUE))

#K45
simcoe_species_pairs_count_1992_k45 <- na.omit(simcoe_1992_ind_W_ts_k45)

simcoe_species_pairs_count_1992_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_1992_k45)

simcoe_mean_W_sp_pairs_1992_k45 <- simcoe_species_pairs_count_1992_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_1992, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1992_W_ts_k45 <- simcoe_1992_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_1992 = sum(ind_W_1992, na.rm = TRUE))

#Step 4
simcoe_1992_ind_W_as1 <- rbind(simcoe_1992_W_ts_c9, simcoe_1992_W_ts_k42, simcoe_1992_W_ts_k45)

simcoe_1992_W_as1 <- simcoe_1992_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_1992 = sum(W_1992, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_1992 <- rbind(simcoe_mean_W_sp_pairs_1992_c9, simcoe_mean_W_sp_pairs_1992_k42, simcoe_mean_W_sp_pairs_1992_k45)

simcoe_lc_sp_pairs_1992 <- simcoe_mean_lc_sp_pairs_1992 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_1992_ind_W_as4 <- filter(simcoe_1992_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_1992_mp  <- na.omit(simcoe_1992_ind_W_as4)

simcoe_species_pairs_count_1992_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_1992_mp)

simcoe_mean_mp_sp_pairs_1992  <- simcoe_species_pairs_count_1992_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_1992, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_1992_sum_W_as4 <- simcoe_1992_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_1992 = sum(ind_W_1992, na.rm = TRUE))


#Step 3
simcoe_1992_sum_W_as4 <- simcoe_1992_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_1992 = sum(as4_1992, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_1992_ind_W_as5 <- filter(simcoe_1992_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_1992_cc <- na.omit(simcoe_1992_ind_W_as5)

simcoe_species_pairs_count_1992_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_1992_cc)

simcoe_mean_cc_sp_pairs_1992  <- simcoe_species_pairs_count_1992_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_1992, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_1992_sum_W_as5 <- simcoe_1992_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_1992 = sum(ind_W_1992, na.rm = TRUE))


#Step 3
simcoe_1992_sum_W_as5 <- simcoe_1992_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_1992 = sum(as5_1992, na.rm = TRUE))

simcoe_1992_W <- merge(simcoe_1992_sum_W_as4, simcoe_1992_W_ts, by = c('Year'))
simcoe_1992_W <- merge(simcoe_1992_W, simcoe_1992_W_as1, by = c('Year'))
simcoe_1992_W <- merge(simcoe_1992_W, simcoe_1992_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_1992 <- merge(simcoe_mean_W_sp_pairs_1992, simcoe_lc_sp_pairs_1992, by=c("Year"))
simcoe_sp_pair_ave_stab_1992 <- merge(simcoe_sp_pair_ave_stab_1992, simcoe_mean_mp_sp_pairs_1992, by=c("Year"))
simcoe_sp_pair_ave_stab_1992 <- merge(simcoe_sp_pair_ave_stab_1992, simcoe_mean_cc_sp_pairs_1992, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_1992_final <- merge(simcoe_1992_W, simcoe_1992_wa_pop_var_1992, by = c('Year'))
simcoe_1992_final <- simcoe_1992_final %>% group_by(Year) %>%
  mutate(gcv_1992 = lcv - W_1992) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_1992, 
            lcv = lcv,
            W = W_1992,
            lc_stab = as2_1992,
            mp_stab = as4_1992,
            cc_stab = as5_1992,
            mc_cv = sqrt(gcv_1992),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_1992/lcv + as4_1992/lcv + as5_1992/lcv,
            lc_asynchrony = as2_1992/lcv,
            mp_asynchrony = as4_1992/lcv,
            cc_asynchrony= as5_1992/lcv)

simcoe_1992_cv <- simcoe_1992_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_1992_final$mc_sum_mean_density <- simcoe_1992_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 1992 #####
simcoe_mean_corr_1992_c9 <- simcoe_species_pairs_count_1992_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1992, ind_W_1992)
simcoe_mean_corr_1992_k42 <- simcoe_species_pairs_count_1992_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1992, ind_W_1992)
simcoe_mean_corr_1992_k45 <- simcoe_species_pairs_count_1992_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1992, ind_W_1992)
simcoe_mean_corr_1992_lc <- rbind(simcoe_mean_corr_1992_c9, simcoe_mean_corr_1992_k42, simcoe_mean_corr_1992_k45)
simcoe_mean_corr_1992_mp <- simcoe_species_pairs_count_1992_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_1992, ind_W_1992)
simcoe_mean_corr_1992_cc <- simcoe_species_pairs_count_1992_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_1992, ind_W_1992)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_1992_lc[sample(1:nrow(simcoe_mean_corr_1992_lc),  nrow(simcoe_mean_corr_1992_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_1992),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_1992),
              mean_lc_stab = mean(ind_W_1992),
              sd_lc_stab = sd(ind_W_1992))
  
  lc_corr_1992 <- as.data.frame(mean_lc)
}
lc_corr_1992

mp_corr_1992 <- simcoe_mean_corr_1992_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_1992),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_1992),
                                                       mean_mp_stab = mean(ind_W_1992),
                                                       sd_mp_stab = sd(ind_W_1992))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_1992_cc[sample(1:nrow(simcoe_mean_corr_1992_cc), nrow(simcoe_mean_corr_1992_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_1992),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_1992),
              mean_cc_stab = mean(ind_W_1992),
              sd_cc_stab = sd(ind_W_1992))
  
  cc_corr_1992 <- as.data.frame(mean_cc)
}
cc_corr_1992


corr_1992 <- cbind(lc_corr_1992, mp_corr_1992)
corr_1992 <- cbind(corr_1992, cc_corr_1992)
corr_1992 <- corr_1992 %>% mutate(Year = 1992)





##### 1993 #####
simcoe_1993 <- subset(simcoe, Year == '1993')


#Transpose dataframe to add 0s within years
simcoe_1993_t1 <- dcast(simcoe_1993, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_1993_m <- melt(simcoe_1993_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_1993 <- subset(simcoe_1993_m, Station_ID == 'C9')
simcoe_c9_1993$Species <- paste("C9", simcoe_c9_1993$Species, sep=" ")
simcoe_k42_1993 <- subset(simcoe_1993_m, Station_ID == 'K42')
simcoe_k42_1993$Species <- paste("K42", simcoe_k42_1993$Species, sep=" ")
simcoe_k45_1993 <- subset(simcoe_1993_m, Station_ID == 'K45')
simcoe_k45_1993$Species <- paste("K45", simcoe_k45_1993$Species, sep=" ")

#recombine dataframes
simcoe_spec_1993 <- rbind(simcoe_c9_1993, simcoe_k42_1993, simcoe_k45_1993)
simcoe_spec_1993 <- simcoe_spec_1993 %>% select(- Station_ID)

#transpose dataframe
simcoe_1993_t2 <- dcast(simcoe_spec_1993, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_1993_t2 <- log(simcoe_1993_t2 + 1)
simcoe_1993_t2 <- simcoe_1993_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_1993_cv <- simcoe_1993_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_1993 = mean(Density),
            pop_sd_1993 = sd(Density)) %>%
  mutate(uw_pop_cv_1993 = pop_sd_1993/pop_mean_1993)

simcoe_1993_cv <- simcoe_1993_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_1993),
         mean_pop_cv = mean(uw_pop_cv_1993, na.rm = T),
         mean_pop_density = sum(pop_mean_1993),
         mean_pop_variance = sum(pop_sd_1993))

simcoe_1993_cv <- simcoe_1993_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_1993/mc_sum_mean_density)%>%
  mutate(w_pop_cv_1993 = (pop_mean_1993/mc_sum_mean_density)*uw_pop_cv_1993) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_1993_wa_pop_var_1993 <- simcoe_1993_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_1993 = sum(w_pop_cv_1993, na.rm = T),
            lcv = sum(w_pop_cv_1993, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_1993_cor_list_ts <- as.dist(round(cor(simcoe_1993_t2[]),2))
simcoe_1993_cor_ts <- stack(simcoe_1993_cor_list_ts, dim.names = TRUE)
simcoe_1993_cor_ts <- simcoe_1993_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_1993_cor - simcoe_1993_cv
simcoe_c9_1993_w_cv_ts <- subset(simcoe_1993_cv, Station_ID == 'C9')
simcoe_c9_1993_w_cv_ts$Species <- paste("C9", simcoe_c9_1993_w_cv_ts$Species, sep=" ")
simcoe_c9_1993_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_1993_w_cv_ts)
simcoe_k42_1993_w_cv_ts <- subset(simcoe_1993_cv, Station_ID == 'K42')
simcoe_k42_1993_w_cv_ts$Species <- paste("K42", simcoe_k42_1993_w_cv_ts$Species, sep=" ")
simcoe_k42_1993_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_1993_w_cv_ts)
simcoe_k45_1993_w_cv_ts <- subset(simcoe_1993_cv, Station_ID == 'K45')
simcoe_k45_1993_w_cv_ts$Species <- paste("K45", simcoe_k45_1993_w_cv_ts$Species, sep=" ")
simcoe_k45_1993_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_1993_w_cv_ts)

simcoe_spec_1993_w_sp_ts <- rbind(simcoe_c9_1993_w_cv_ts, simcoe_k42_1993_w_cv_ts, simcoe_k45_1993_w_cv_ts)
simcoe_spec_1993_w_sp_ts <- simcoe_spec_1993_w_sp_ts %>% select(Species, w_pop_cv_1993, Species_ID, Station_ID)
simcoe_spec_1993_w_sp1_ts <- simcoe_spec_1993_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_1993_sp1 = w_pop_cv_1993) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_1993_sp1, Species_ID, station_1)

simcoe_1993_cor_cv_ts <- merge(simcoe_1993_cor_ts, simcoe_spec_1993_w_sp1_ts, by = "species_1")

simcoe_spec_1993_w_sp2_ts <- simcoe_spec_1993_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_1993_sp2 = w_pop_cv_1993) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_1993_sp2, Species_ID, station_2)

simcoe_1993_cor_cv_ts <- merge(simcoe_1993_cor_cv_ts, simcoe_spec_1993_w_sp2_ts, by = "species_2")

simcoe_1993_cor_cv_ts_omit <- na.omit(simcoe_1993_cor_cv_ts)

simcoe_1993_ind_W_ts <- simcoe_1993_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_1993 = w_pop_cv_1993_sp1*w_pop_cv_1993_sp2,
         ind_W_1993 = (1 - corr)*(w_pop_cv_1993_sp1*w_pop_cv_1993_sp2),
         Year = 1993)

simcoe_1993_ind_W_ts$number_species_pairs <- nrow(simcoe_1993_ind_W_ts)

simcoe_mean_W_sp_pairs_1993 <- simcoe_1993_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_1993, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_1993_W_ts <- simcoe_1993_ind_W_ts %>% group_by(Year) %>%
  summarize(W_1993 = sum(ind_W_1993, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_1993_ind_W_ts_c9 <- filter(simcoe_1993_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_1993_ind_W_ts_k42 <- filter(simcoe_1993_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_1993_ind_W_ts_k45 <- filter(simcoe_1993_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_1993_c9 <- na.omit(simcoe_1993_ind_W_ts_c9)

simcoe_species_pairs_count_1993_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_1993_c9)

simcoe_mean_W_sp_pairs_1993_c9 <- simcoe_species_pairs_count_1993_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_1993, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1993_W_ts_c9 <- simcoe_1993_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_1993 = sum(ind_W_1993, na.rm = TRUE))

#K42
simcoe_species_pairs_count_1993_k42 <- na.omit(simcoe_1993_ind_W_ts_k42)

simcoe_species_pairs_count_1993_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_1993_k42)

simcoe_mean_W_sp_pairs_1993_k42 <- simcoe_species_pairs_count_1993_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_1993, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1993_W_ts_k42 <- simcoe_1993_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_1993 = sum(ind_W_1993, na.rm = TRUE))

#K45
simcoe_species_pairs_count_1993_k45 <- na.omit(simcoe_1993_ind_W_ts_k45)

simcoe_species_pairs_count_1993_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_1993_k45)

simcoe_mean_W_sp_pairs_1993_k45 <- simcoe_species_pairs_count_1993_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_1993, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1993_W_ts_k45 <- simcoe_1993_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_1993 = sum(ind_W_1993, na.rm = TRUE))

#Step 4
simcoe_1993_ind_W_as1 <- rbind(simcoe_1993_W_ts_c9, simcoe_1993_W_ts_k42, simcoe_1993_W_ts_k45)

simcoe_1993_W_as1 <- simcoe_1993_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_1993 = sum(W_1993, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_1993 <- rbind(simcoe_mean_W_sp_pairs_1993_c9, simcoe_mean_W_sp_pairs_1993_k42, simcoe_mean_W_sp_pairs_1993_k45)

simcoe_lc_sp_pairs_1993 <- simcoe_mean_lc_sp_pairs_1993 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_1993_ind_W_as4 <- filter(simcoe_1993_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_1993_mp  <- na.omit(simcoe_1993_ind_W_as4)

simcoe_species_pairs_count_1993_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_1993_mp)

simcoe_mean_mp_sp_pairs_1993  <- simcoe_species_pairs_count_1993_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_1993, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_1993_sum_W_as4 <- simcoe_1993_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_1993 = sum(ind_W_1993, na.rm = TRUE))


#Step 3
simcoe_1993_sum_W_as4 <- simcoe_1993_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_1993 = sum(as4_1993, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_1993_ind_W_as5 <- filter(simcoe_1993_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_1993_cc <- na.omit(simcoe_1993_ind_W_as5)

simcoe_species_pairs_count_1993_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_1993_cc)

simcoe_mean_cc_sp_pairs_1993  <- simcoe_species_pairs_count_1993_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_1993, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_1993_sum_W_as5 <- simcoe_1993_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_1993 = sum(ind_W_1993, na.rm = TRUE))


#Step 3
simcoe_1993_sum_W_as5 <- simcoe_1993_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_1993 = sum(as5_1993, na.rm = TRUE))

simcoe_1993_W <- merge(simcoe_1993_sum_W_as4, simcoe_1993_W_ts, by = c('Year'))
simcoe_1993_W <- merge(simcoe_1993_W, simcoe_1993_W_as1, by = c('Year'))
simcoe_1993_W <- merge(simcoe_1993_W, simcoe_1993_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_1993 <- merge(simcoe_mean_W_sp_pairs_1993, simcoe_lc_sp_pairs_1993, by=c("Year"))
simcoe_sp_pair_ave_stab_1993 <- merge(simcoe_sp_pair_ave_stab_1993, simcoe_mean_mp_sp_pairs_1993, by=c("Year"))
simcoe_sp_pair_ave_stab_1993 <- merge(simcoe_sp_pair_ave_stab_1993, simcoe_mean_cc_sp_pairs_1993, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_1993_final <- merge(simcoe_1993_W, simcoe_1993_wa_pop_var_1993, by = c('Year'))
simcoe_1993_final <- simcoe_1993_final %>% group_by(Year) %>%
  mutate(gcv_1993 = lcv - W_1993) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_1993, 
            lcv = lcv,
            W = W_1993,
            lc_stab = as2_1993,
            mp_stab = as4_1993,
            cc_stab = as5_1993,
            mc_cv = sqrt(gcv_1993),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_1993/lcv + as4_1993/lcv + as5_1993/lcv,
            lc_asynchrony = as2_1993/lcv,
            mp_asynchrony = as4_1993/lcv,
            cc_asynchrony= as5_1993/lcv)

simcoe_1993_cv <- simcoe_1993_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_1993_final$mc_sum_mean_density <- simcoe_1993_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 1993 #####
simcoe_mean_corr_1993_c9 <- simcoe_species_pairs_count_1993_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1993, ind_W_1993)
simcoe_mean_corr_1993_k42 <- simcoe_species_pairs_count_1993_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1993, ind_W_1993)
simcoe_mean_corr_1993_k45 <- simcoe_species_pairs_count_1993_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1993, ind_W_1993)
simcoe_mean_corr_1993_lc <- rbind(simcoe_mean_corr_1993_c9, simcoe_mean_corr_1993_k42, simcoe_mean_corr_1993_k45)
simcoe_mean_corr_1993_mp <- simcoe_species_pairs_count_1993_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_1993, ind_W_1993)
simcoe_mean_corr_1993_cc <- simcoe_species_pairs_count_1993_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_1993, ind_W_1993)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_1993_lc[sample(1:nrow(simcoe_mean_corr_1993_lc),  nrow(simcoe_mean_corr_1993_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_1993),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_1993),
              mean_lc_stab = mean(ind_W_1993),
              sd_lc_stab = sd(ind_W_1993))
  
  lc_corr_1993 <- as.data.frame(mean_lc)
}
lc_corr_1993

mp_corr_1993 <- simcoe_mean_corr_1993_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_1993),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_1993),
                                                       mean_mp_stab = mean(ind_W_1993),
                                                       sd_mp_stab = sd(ind_W_1993))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_1993_cc[sample(1:nrow(simcoe_mean_corr_1993_cc), nrow(simcoe_mean_corr_1993_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_1993),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_1993),
              mean_cc_stab = mean(ind_W_1993),
              sd_cc_stab = sd(ind_W_1993))
  
  cc_corr_1993 <- as.data.frame(mean_cc)
}
cc_corr_1993


corr_1993 <- cbind(lc_corr_1993, mp_corr_1993)
corr_1993 <- cbind(corr_1993, cc_corr_1993)
corr_1993 <- corr_1993 %>% mutate(Year = 1993)





##### 1994 #####
simcoe_1994 <- subset(simcoe, Year == '1994')


#Transpose dataframe to add 0s within years
simcoe_1994_t1 <- dcast(simcoe_1994, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_1994_m <- melt(simcoe_1994_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_1994 <- subset(simcoe_1994_m, Station_ID == 'C9')
simcoe_c9_1994$Species <- paste("C9", simcoe_c9_1994$Species, sep=" ")
simcoe_k42_1994 <- subset(simcoe_1994_m, Station_ID == 'K42')
simcoe_k42_1994$Species <- paste("K42", simcoe_k42_1994$Species, sep=" ")
simcoe_k45_1994 <- subset(simcoe_1994_m, Station_ID == 'K45')
simcoe_k45_1994$Species <- paste("K45", simcoe_k45_1994$Species, sep=" ")

#recombine dataframes
simcoe_spec_1994 <- rbind(simcoe_c9_1994, simcoe_k42_1994, simcoe_k45_1994)
simcoe_spec_1994 <- simcoe_spec_1994 %>% select(- Station_ID)

#transpose dataframe
simcoe_1994_t2 <- dcast(simcoe_spec_1994, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_1994_t2 <- log(simcoe_1994_t2 + 1)
simcoe_1994_t2 <- simcoe_1994_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_1994_cv <- simcoe_1994_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_1994 = mean(Density),
            pop_sd_1994 = sd(Density)) %>%
  mutate(uw_pop_cv_1994 = pop_sd_1994/pop_mean_1994)

simcoe_1994_cv <- simcoe_1994_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_1994),
         mean_pop_cv = mean(uw_pop_cv_1994, na.rm = T),
         mean_pop_density = sum(pop_mean_1994),
         mean_pop_variance = sum(pop_sd_1994))

simcoe_1994_cv <- simcoe_1994_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_1994/mc_sum_mean_density)%>%
  mutate(w_pop_cv_1994 = (pop_mean_1994/mc_sum_mean_density)*uw_pop_cv_1994) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_1994_wa_pop_var_1994 <- simcoe_1994_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_1994 = sum(w_pop_cv_1994, na.rm = T),
            lcv = sum(w_pop_cv_1994, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_1994_cor_list_ts <- as.dist(round(cor(simcoe_1994_t2[]),2))
simcoe_1994_cor_ts <- stack(simcoe_1994_cor_list_ts, dim.names = TRUE)
simcoe_1994_cor_ts <- simcoe_1994_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_1994_cor - simcoe_1994_cv
simcoe_c9_1994_w_cv_ts <- subset(simcoe_1994_cv, Station_ID == 'C9')
simcoe_c9_1994_w_cv_ts$Species <- paste("C9", simcoe_c9_1994_w_cv_ts$Species, sep=" ")
simcoe_c9_1994_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_1994_w_cv_ts)
simcoe_k42_1994_w_cv_ts <- subset(simcoe_1994_cv, Station_ID == 'K42')
simcoe_k42_1994_w_cv_ts$Species <- paste("K42", simcoe_k42_1994_w_cv_ts$Species, sep=" ")
simcoe_k42_1994_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_1994_w_cv_ts)
simcoe_k45_1994_w_cv_ts <- subset(simcoe_1994_cv, Station_ID == 'K45')
simcoe_k45_1994_w_cv_ts$Species <- paste("K45", simcoe_k45_1994_w_cv_ts$Species, sep=" ")
simcoe_k45_1994_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_1994_w_cv_ts)

simcoe_spec_1994_w_sp_ts <- rbind(simcoe_c9_1994_w_cv_ts, simcoe_k42_1994_w_cv_ts, simcoe_k45_1994_w_cv_ts)
simcoe_spec_1994_w_sp_ts <- simcoe_spec_1994_w_sp_ts %>% select(Species, w_pop_cv_1994, Species_ID, Station_ID)
simcoe_spec_1994_w_sp1_ts <- simcoe_spec_1994_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_1994_sp1 = w_pop_cv_1994) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_1994_sp1, Species_ID, station_1)

simcoe_1994_cor_cv_ts <- merge(simcoe_1994_cor_ts, simcoe_spec_1994_w_sp1_ts, by = "species_1")

simcoe_spec_1994_w_sp2_ts <- simcoe_spec_1994_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_1994_sp2 = w_pop_cv_1994) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_1994_sp2, Species_ID, station_2)

simcoe_1994_cor_cv_ts <- merge(simcoe_1994_cor_cv_ts, simcoe_spec_1994_w_sp2_ts, by = "species_2")

simcoe_1994_cor_cv_ts_omit <- na.omit(simcoe_1994_cor_cv_ts)

simcoe_1994_ind_W_ts <- simcoe_1994_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_1994 = w_pop_cv_1994_sp1*w_pop_cv_1994_sp2,
         ind_W_1994 = (1 - corr)*(w_pop_cv_1994_sp1*w_pop_cv_1994_sp2),
         Year = 1994)

simcoe_1994_ind_W_ts$number_species_pairs <- nrow(simcoe_1994_ind_W_ts)

simcoe_mean_W_sp_pairs_1994 <- simcoe_1994_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_1994, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_1994_W_ts <- simcoe_1994_ind_W_ts %>% group_by(Year) %>%
  summarize(W_1994 = sum(ind_W_1994, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_1994_ind_W_ts_c9 <- filter(simcoe_1994_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_1994_ind_W_ts_k42 <- filter(simcoe_1994_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_1994_ind_W_ts_k45 <- filter(simcoe_1994_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_1994_c9 <- na.omit(simcoe_1994_ind_W_ts_c9)

simcoe_species_pairs_count_1994_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_1994_c9)

simcoe_mean_W_sp_pairs_1994_c9 <- simcoe_species_pairs_count_1994_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_1994, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1994_W_ts_c9 <- simcoe_1994_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_1994 = sum(ind_W_1994, na.rm = TRUE))

#K42
simcoe_species_pairs_count_1994_k42 <- na.omit(simcoe_1994_ind_W_ts_k42)

simcoe_species_pairs_count_1994_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_1994_k42)

simcoe_mean_W_sp_pairs_1994_k42 <- simcoe_species_pairs_count_1994_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_1994, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1994_W_ts_k42 <- simcoe_1994_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_1994 = sum(ind_W_1994, na.rm = TRUE))

#K45
simcoe_species_pairs_count_1994_k45 <- na.omit(simcoe_1994_ind_W_ts_k45)

simcoe_species_pairs_count_1994_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_1994_k45)

simcoe_mean_W_sp_pairs_1994_k45 <- simcoe_species_pairs_count_1994_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_1994, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1994_W_ts_k45 <- simcoe_1994_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_1994 = sum(ind_W_1994, na.rm = TRUE))

#Step 4
simcoe_1994_ind_W_as1 <- rbind(simcoe_1994_W_ts_c9, simcoe_1994_W_ts_k42, simcoe_1994_W_ts_k45)

simcoe_1994_W_as1 <- simcoe_1994_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_1994 = sum(W_1994, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_1994 <- rbind(simcoe_mean_W_sp_pairs_1994_c9, simcoe_mean_W_sp_pairs_1994_k42, simcoe_mean_W_sp_pairs_1994_k45)

simcoe_lc_sp_pairs_1994 <- simcoe_mean_lc_sp_pairs_1994 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_1994_ind_W_as4 <- filter(simcoe_1994_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_1994_mp  <- na.omit(simcoe_1994_ind_W_as4)

simcoe_species_pairs_count_1994_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_1994_mp)

simcoe_mean_mp_sp_pairs_1994  <- simcoe_species_pairs_count_1994_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_1994, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_1994_sum_W_as4 <- simcoe_1994_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_1994 = sum(ind_W_1994, na.rm = TRUE))


#Step 3
simcoe_1994_sum_W_as4 <- simcoe_1994_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_1994 = sum(as4_1994, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_1994_ind_W_as5 <- filter(simcoe_1994_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_1994_cc <- na.omit(simcoe_1994_ind_W_as5)

simcoe_species_pairs_count_1994_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_1994_cc)

simcoe_mean_cc_sp_pairs_1994  <- simcoe_species_pairs_count_1994_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_1994, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_1994_sum_W_as5 <- simcoe_1994_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_1994 = sum(ind_W_1994, na.rm = TRUE))


#Step 3
simcoe_1994_sum_W_as5 <- simcoe_1994_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_1994 = sum(as5_1994, na.rm = TRUE))

simcoe_1994_W <- merge(simcoe_1994_sum_W_as4, simcoe_1994_W_ts, by = c('Year'))
simcoe_1994_W <- merge(simcoe_1994_W, simcoe_1994_W_as1, by = c('Year'))
simcoe_1994_W <- merge(simcoe_1994_W, simcoe_1994_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_1994 <- merge(simcoe_mean_W_sp_pairs_1994, simcoe_lc_sp_pairs_1994, by=c("Year"))
simcoe_sp_pair_ave_stab_1994 <- merge(simcoe_sp_pair_ave_stab_1994, simcoe_mean_mp_sp_pairs_1994, by=c("Year"))
simcoe_sp_pair_ave_stab_1994 <- merge(simcoe_sp_pair_ave_stab_1994, simcoe_mean_cc_sp_pairs_1994, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_1994_final <- merge(simcoe_1994_W, simcoe_1994_wa_pop_var_1994, by = c('Year'))
simcoe_1994_final <- simcoe_1994_final %>% group_by(Year) %>%
  mutate(gcv_1994 = lcv - W_1994) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_1994, 
            lcv = lcv,
            W = W_1994,
            lc_stab = as2_1994,
            mp_stab = as4_1994,
            cc_stab = as5_1994,
            mc_cv = sqrt(gcv_1994),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_1994/lcv + as4_1994/lcv + as5_1994/lcv,
            lc_asynchrony = as2_1994/lcv,
            mp_asynchrony = as4_1994/lcv,
            cc_asynchrony= as5_1994/lcv)

simcoe_1994_cv <- simcoe_1994_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_1994_final$mc_sum_mean_density <- simcoe_1994_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 1994 #####
simcoe_mean_corr_1994_c9 <- simcoe_species_pairs_count_1994_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1994, ind_W_1994)
simcoe_mean_corr_1994_k42 <- simcoe_species_pairs_count_1994_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1994, ind_W_1994)
simcoe_mean_corr_1994_k45 <- simcoe_species_pairs_count_1994_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1994, ind_W_1994)
simcoe_mean_corr_1994_lc <- rbind(simcoe_mean_corr_1994_c9, simcoe_mean_corr_1994_k42, simcoe_mean_corr_1994_k45)
simcoe_mean_corr_1994_mp <- simcoe_species_pairs_count_1994_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_1994, ind_W_1994)
simcoe_mean_corr_1994_cc <- simcoe_species_pairs_count_1994_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_1994, ind_W_1994)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_1994_lc[sample(1:nrow(simcoe_mean_corr_1994_lc),  nrow(simcoe_mean_corr_1994_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_1994),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_1994),
              mean_lc_stab = mean(ind_W_1994),
              sd_lc_stab = sd(ind_W_1994))
  
  lc_corr_1994 <- as.data.frame(mean_lc)
}
lc_corr_1994

mp_corr_1994 <- simcoe_mean_corr_1994_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_1994),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_1994),
                                                       mean_mp_stab = mean(ind_W_1994),
                                                       sd_mp_stab = sd(ind_W_1994))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_1994_cc[sample(1:nrow(simcoe_mean_corr_1994_cc), nrow(simcoe_mean_corr_1994_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_1994),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_1994),
              mean_cc_stab = mean(ind_W_1994),
              sd_cc_stab = sd(ind_W_1994))
  
  cc_corr_1994 <- as.data.frame(mean_cc)
}
cc_corr_1994


corr_1994 <- cbind(lc_corr_1994, mp_corr_1994)
corr_1994 <- cbind(corr_1994, cc_corr_1994)
corr_1994 <- corr_1994 %>% mutate(Year = 1994)





##### 1995 #####
simcoe_1995 <- subset(simcoe, Year == '1995')


#Transpose dataframe to add 0s within years
simcoe_1995_t1 <- dcast(simcoe_1995, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_1995_m <- melt(simcoe_1995_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_1995 <- subset(simcoe_1995_m, Station_ID == 'C9')
simcoe_c9_1995$Species <- paste("C9", simcoe_c9_1995$Species, sep=" ")
simcoe_k42_1995 <- subset(simcoe_1995_m, Station_ID == 'K42')
simcoe_k42_1995$Species <- paste("K42", simcoe_k42_1995$Species, sep=" ")
simcoe_k45_1995 <- subset(simcoe_1995_m, Station_ID == 'K45')
simcoe_k45_1995$Species <- paste("K45", simcoe_k45_1995$Species, sep=" ")

#recombine dataframes
simcoe_spec_1995 <- rbind(simcoe_c9_1995, simcoe_k42_1995, simcoe_k45_1995)
simcoe_spec_1995 <- simcoe_spec_1995 %>% select(- Station_ID)

#transpose dataframe
simcoe_1995_t2 <- dcast(simcoe_spec_1995, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_1995_t2 <- log(simcoe_1995_t2 + 1)
simcoe_1995_t2 <- simcoe_1995_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_1995_cv <- simcoe_1995_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_1995 = mean(Density),
            pop_sd_1995 = sd(Density)) %>%
  mutate(uw_pop_cv_1995 = pop_sd_1995/pop_mean_1995)

simcoe_1995_cv <- simcoe_1995_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_1995),
         mean_pop_cv = mean(uw_pop_cv_1995, na.rm = T),
         mean_pop_density = sum(pop_mean_1995),
         mean_pop_variance = sum(pop_sd_1995))

simcoe_1995_cv <- simcoe_1995_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_1995/mc_sum_mean_density)%>%
  mutate(w_pop_cv_1995 = (pop_mean_1995/mc_sum_mean_density)*uw_pop_cv_1995) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_1995_wa_pop_var_1995 <- simcoe_1995_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_1995 = sum(w_pop_cv_1995, na.rm = T),
            lcv = sum(w_pop_cv_1995, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_1995_cor_list_ts <- as.dist(round(cor(simcoe_1995_t2[]),2))
simcoe_1995_cor_ts <- stack(simcoe_1995_cor_list_ts, dim.names = TRUE)
simcoe_1995_cor_ts <- simcoe_1995_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_1995_cor - simcoe_1995_cv
simcoe_c9_1995_w_cv_ts <- subset(simcoe_1995_cv, Station_ID == 'C9')
simcoe_c9_1995_w_cv_ts$Species <- paste("C9", simcoe_c9_1995_w_cv_ts$Species, sep=" ")
simcoe_c9_1995_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_1995_w_cv_ts)
simcoe_k42_1995_w_cv_ts <- subset(simcoe_1995_cv, Station_ID == 'K42')
simcoe_k42_1995_w_cv_ts$Species <- paste("K42", simcoe_k42_1995_w_cv_ts$Species, sep=" ")
simcoe_k42_1995_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_1995_w_cv_ts)
simcoe_k45_1995_w_cv_ts <- subset(simcoe_1995_cv, Station_ID == 'K45')
simcoe_k45_1995_w_cv_ts$Species <- paste("K45", simcoe_k45_1995_w_cv_ts$Species, sep=" ")
simcoe_k45_1995_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_1995_w_cv_ts)

simcoe_spec_1995_w_sp_ts <- rbind(simcoe_c9_1995_w_cv_ts, simcoe_k42_1995_w_cv_ts, simcoe_k45_1995_w_cv_ts)
simcoe_spec_1995_w_sp_ts <- simcoe_spec_1995_w_sp_ts %>% select(Species, w_pop_cv_1995, Species_ID, Station_ID)
simcoe_spec_1995_w_sp1_ts <- simcoe_spec_1995_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_1995_sp1 = w_pop_cv_1995) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_1995_sp1, Species_ID, station_1)

simcoe_1995_cor_cv_ts <- merge(simcoe_1995_cor_ts, simcoe_spec_1995_w_sp1_ts, by = "species_1")

simcoe_spec_1995_w_sp2_ts <- simcoe_spec_1995_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_1995_sp2 = w_pop_cv_1995) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_1995_sp2, Species_ID, station_2)

simcoe_1995_cor_cv_ts <- merge(simcoe_1995_cor_cv_ts, simcoe_spec_1995_w_sp2_ts, by = "species_2")

simcoe_1995_cor_cv_ts_omit <- na.omit(simcoe_1995_cor_cv_ts)

simcoe_1995_ind_W_ts <- simcoe_1995_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_1995 = w_pop_cv_1995_sp1*w_pop_cv_1995_sp2,
         ind_W_1995 = (1 - corr)*(w_pop_cv_1995_sp1*w_pop_cv_1995_sp2),
         Year = 1995)

simcoe_1995_ind_W_ts$number_species_pairs <- nrow(simcoe_1995_ind_W_ts)

simcoe_mean_W_sp_pairs_1995 <- simcoe_1995_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_1995, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_1995_W_ts <- simcoe_1995_ind_W_ts %>% group_by(Year) %>%
  summarize(W_1995 = sum(ind_W_1995, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_1995_ind_W_ts_c9 <- filter(simcoe_1995_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_1995_ind_W_ts_k42 <- filter(simcoe_1995_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_1995_ind_W_ts_k45 <- filter(simcoe_1995_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_1995_c9 <- na.omit(simcoe_1995_ind_W_ts_c9)

simcoe_species_pairs_count_1995_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_1995_c9)

simcoe_mean_W_sp_pairs_1995_c9 <- simcoe_species_pairs_count_1995_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_1995, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1995_W_ts_c9 <- simcoe_1995_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_1995 = sum(ind_W_1995, na.rm = TRUE))

#K42
simcoe_species_pairs_count_1995_k42 <- na.omit(simcoe_1995_ind_W_ts_k42)

simcoe_species_pairs_count_1995_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_1995_k42)

simcoe_mean_W_sp_pairs_1995_k42 <- simcoe_species_pairs_count_1995_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_1995, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1995_W_ts_k42 <- simcoe_1995_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_1995 = sum(ind_W_1995, na.rm = TRUE))

#K45
simcoe_species_pairs_count_1995_k45 <- na.omit(simcoe_1995_ind_W_ts_k45)

simcoe_species_pairs_count_1995_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_1995_k45)

simcoe_mean_W_sp_pairs_1995_k45 <- simcoe_species_pairs_count_1995_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_1995, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1995_W_ts_k45 <- simcoe_1995_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_1995 = sum(ind_W_1995, na.rm = TRUE))

#Step 4
simcoe_1995_ind_W_as1 <- rbind(simcoe_1995_W_ts_c9, simcoe_1995_W_ts_k42, simcoe_1995_W_ts_k45)

simcoe_1995_W_as1 <- simcoe_1995_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_1995 = sum(W_1995, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_1995 <- rbind(simcoe_mean_W_sp_pairs_1995_c9, simcoe_mean_W_sp_pairs_1995_k42, simcoe_mean_W_sp_pairs_1995_k45)

simcoe_lc_sp_pairs_1995 <- simcoe_mean_lc_sp_pairs_1995 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_1995_ind_W_as4 <- filter(simcoe_1995_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_1995_mp  <- na.omit(simcoe_1995_ind_W_as4)

simcoe_species_pairs_count_1995_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_1995_mp)

simcoe_mean_mp_sp_pairs_1995  <- simcoe_species_pairs_count_1995_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_1995, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_1995_sum_W_as4 <- simcoe_1995_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_1995 = sum(ind_W_1995, na.rm = TRUE))


#Step 3
simcoe_1995_sum_W_as4 <- simcoe_1995_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_1995 = sum(as4_1995, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_1995_ind_W_as5 <- filter(simcoe_1995_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_1995_cc <- na.omit(simcoe_1995_ind_W_as5)

simcoe_species_pairs_count_1995_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_1995_cc)

simcoe_mean_cc_sp_pairs_1995  <- simcoe_species_pairs_count_1995_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_1995, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_1995_sum_W_as5 <- simcoe_1995_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_1995 = sum(ind_W_1995, na.rm = TRUE))


#Step 3
simcoe_1995_sum_W_as5 <- simcoe_1995_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_1995 = sum(as5_1995, na.rm = TRUE))

simcoe_1995_W <- merge(simcoe_1995_sum_W_as4, simcoe_1995_W_ts, by = c('Year'))
simcoe_1995_W <- merge(simcoe_1995_W, simcoe_1995_W_as1, by = c('Year'))
simcoe_1995_W <- merge(simcoe_1995_W, simcoe_1995_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_1995 <- merge(simcoe_mean_W_sp_pairs_1995, simcoe_lc_sp_pairs_1995, by=c("Year"))
simcoe_sp_pair_ave_stab_1995 <- merge(simcoe_sp_pair_ave_stab_1995, simcoe_mean_mp_sp_pairs_1995, by=c("Year"))
simcoe_sp_pair_ave_stab_1995 <- merge(simcoe_sp_pair_ave_stab_1995, simcoe_mean_cc_sp_pairs_1995, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_1995_final <- merge(simcoe_1995_W, simcoe_1995_wa_pop_var_1995, by = c('Year'))
simcoe_1995_final <- simcoe_1995_final %>% group_by(Year) %>%
  mutate(gcv_1995 = lcv - W_1995) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_1995, 
            lcv = lcv,
            W = W_1995,
            lc_stab = as2_1995,
            mp_stab = as4_1995,
            cc_stab = as5_1995,
            mc_cv = sqrt(gcv_1995),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_1995/lcv + as4_1995/lcv + as5_1995/lcv,
            lc_asynchrony = as2_1995/lcv,
            mp_asynchrony = as4_1995/lcv,
            cc_asynchrony= as5_1995/lcv)

simcoe_1995_cv <- simcoe_1995_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_1995_final$mc_sum_mean_density <- simcoe_1995_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 1995 #####
simcoe_mean_corr_1995_c9 <- simcoe_species_pairs_count_1995_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1995, ind_W_1995)
simcoe_mean_corr_1995_k42 <- simcoe_species_pairs_count_1995_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1995, ind_W_1995)
simcoe_mean_corr_1995_k45 <- simcoe_species_pairs_count_1995_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1995, ind_W_1995)
simcoe_mean_corr_1995_lc <- rbind(simcoe_mean_corr_1995_c9, simcoe_mean_corr_1995_k42, simcoe_mean_corr_1995_k45)
simcoe_mean_corr_1995_mp <- simcoe_species_pairs_count_1995_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_1995, ind_W_1995)
simcoe_mean_corr_1995_cc <- simcoe_species_pairs_count_1995_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_1995, ind_W_1995)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_1995_lc[sample(1:nrow(simcoe_mean_corr_1995_lc),  nrow(simcoe_mean_corr_1995_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_1995),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_1995),
              mean_lc_stab = mean(ind_W_1995),
              sd_lc_stab = sd(ind_W_1995))
  
  lc_corr_1995 <- as.data.frame(mean_lc)
}
lc_corr_1995

mp_corr_1995 <- simcoe_mean_corr_1995_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_1995),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_1995),
                                                       mean_mp_stab = mean(ind_W_1995),
                                                       sd_mp_stab = sd(ind_W_1995))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_1995_cc[sample(1:nrow(simcoe_mean_corr_1995_cc), nrow(simcoe_mean_corr_1995_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_1995),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_1995),
              mean_cc_stab = mean(ind_W_1995),
              sd_cc_stab = sd(ind_W_1995))
  
  cc_corr_1995 <- as.data.frame(mean_cc)
}
cc_corr_1995


corr_1995 <- cbind(lc_corr_1995, mp_corr_1995)
corr_1995 <- cbind(corr_1995, cc_corr_1995)
corr_1995 <- corr_1995 %>% mutate(Year = 1995)




##### 1996 #####
simcoe_1996 <- subset(simcoe, Year == '1996')


#Transpose dataframe to add 0s within years
simcoe_1996_t1 <- dcast(simcoe_1996, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_1996_m <- melt(simcoe_1996_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_1996 <- subset(simcoe_1996_m, Station_ID == 'C9')
simcoe_c9_1996$Species <- paste("C9", simcoe_c9_1996$Species, sep=" ")
simcoe_k42_1996 <- subset(simcoe_1996_m, Station_ID == 'K42')
simcoe_k42_1996$Species <- paste("K42", simcoe_k42_1996$Species, sep=" ")
simcoe_k45_1996 <- subset(simcoe_1996_m, Station_ID == 'K45')
simcoe_k45_1996$Species <- paste("K45", simcoe_k45_1996$Species, sep=" ")

#recombine dataframes
simcoe_spec_1996 <- rbind(simcoe_c9_1996, simcoe_k42_1996, simcoe_k45_1996)
simcoe_spec_1996 <- simcoe_spec_1996 %>% select(- Station_ID)

#transpose dataframe
simcoe_1996_t2 <- dcast(simcoe_spec_1996, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_1996_t2 <- log(simcoe_1996_t2 + 1)
simcoe_1996_t2 <- simcoe_1996_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_1996_cv <- simcoe_1996_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_1996 = mean(Density),
            pop_sd_1996 = sd(Density)) %>%
  mutate(uw_pop_cv_1996 = pop_sd_1996/pop_mean_1996)

simcoe_1996_cv <- simcoe_1996_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_1996),
         mean_pop_cv = mean(uw_pop_cv_1996, na.rm = T),
         mean_pop_density = sum(pop_mean_1996),
         mean_pop_variance = sum(pop_sd_1996))

simcoe_1996_cv <- simcoe_1996_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_1996/mc_sum_mean_density)%>%
  mutate(w_pop_cv_1996 = (pop_mean_1996/mc_sum_mean_density)*uw_pop_cv_1996) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_1996_wa_pop_var_1996 <- simcoe_1996_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_1996 = sum(w_pop_cv_1996, na.rm = T),
            lcv = sum(w_pop_cv_1996, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_1996_cor_list_ts <- as.dist(round(cor(simcoe_1996_t2[]),2))
simcoe_1996_cor_ts <- stack(simcoe_1996_cor_list_ts, dim.names = TRUE)
simcoe_1996_cor_ts <- simcoe_1996_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_1996_cor - simcoe_1996_cv
simcoe_c9_1996_w_cv_ts <- subset(simcoe_1996_cv, Station_ID == 'C9')
simcoe_c9_1996_w_cv_ts$Species <- paste("C9", simcoe_c9_1996_w_cv_ts$Species, sep=" ")
simcoe_c9_1996_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_1996_w_cv_ts)
simcoe_k42_1996_w_cv_ts <- subset(simcoe_1996_cv, Station_ID == 'K42')
simcoe_k42_1996_w_cv_ts$Species <- paste("K42", simcoe_k42_1996_w_cv_ts$Species, sep=" ")
simcoe_k42_1996_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_1996_w_cv_ts)
simcoe_k45_1996_w_cv_ts <- subset(simcoe_1996_cv, Station_ID == 'K45')
simcoe_k45_1996_w_cv_ts$Species <- paste("K45", simcoe_k45_1996_w_cv_ts$Species, sep=" ")
simcoe_k45_1996_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_1996_w_cv_ts)

simcoe_spec_1996_w_sp_ts <- rbind(simcoe_c9_1996_w_cv_ts, simcoe_k42_1996_w_cv_ts, simcoe_k45_1996_w_cv_ts)
simcoe_spec_1996_w_sp_ts <- simcoe_spec_1996_w_sp_ts %>% select(Species, w_pop_cv_1996, Species_ID, Station_ID)
simcoe_spec_1996_w_sp1_ts <- simcoe_spec_1996_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_1996_sp1 = w_pop_cv_1996) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_1996_sp1, Species_ID, station_1)

simcoe_1996_cor_cv_ts <- merge(simcoe_1996_cor_ts, simcoe_spec_1996_w_sp1_ts, by = "species_1")

simcoe_spec_1996_w_sp2_ts <- simcoe_spec_1996_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_1996_sp2 = w_pop_cv_1996) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_1996_sp2, Species_ID, station_2)

simcoe_1996_cor_cv_ts <- merge(simcoe_1996_cor_cv_ts, simcoe_spec_1996_w_sp2_ts, by = "species_2")

simcoe_1996_cor_cv_ts_omit <- na.omit(simcoe_1996_cor_cv_ts)

simcoe_1996_ind_W_ts <- simcoe_1996_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_1996 = w_pop_cv_1996_sp1*w_pop_cv_1996_sp2,
         ind_W_1996 = (1 - corr)*(w_pop_cv_1996_sp1*w_pop_cv_1996_sp2),
         Year = 1996)

simcoe_1996_ind_W_ts$number_species_pairs <- nrow(simcoe_1996_ind_W_ts)

simcoe_mean_W_sp_pairs_1996 <- simcoe_1996_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_1996, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_1996_W_ts <- simcoe_1996_ind_W_ts %>% group_by(Year) %>%
  summarize(W_1996 = sum(ind_W_1996, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_1996_ind_W_ts_c9 <- filter(simcoe_1996_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_1996_ind_W_ts_k42 <- filter(simcoe_1996_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_1996_ind_W_ts_k45 <- filter(simcoe_1996_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_1996_c9 <- na.omit(simcoe_1996_ind_W_ts_c9)

simcoe_species_pairs_count_1996_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_1996_c9)

simcoe_mean_W_sp_pairs_1996_c9 <- simcoe_species_pairs_count_1996_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_1996, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1996_W_ts_c9 <- simcoe_1996_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_1996 = sum(ind_W_1996, na.rm = TRUE))

#K42
simcoe_species_pairs_count_1996_k42 <- na.omit(simcoe_1996_ind_W_ts_k42)

simcoe_species_pairs_count_1996_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_1996_k42)

simcoe_mean_W_sp_pairs_1996_k42 <- simcoe_species_pairs_count_1996_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_1996, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1996_W_ts_k42 <- simcoe_1996_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_1996 = sum(ind_W_1996, na.rm = TRUE))

#K45
simcoe_species_pairs_count_1996_k45 <- na.omit(simcoe_1996_ind_W_ts_k45)

simcoe_species_pairs_count_1996_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_1996_k45)

simcoe_mean_W_sp_pairs_1996_k45 <- simcoe_species_pairs_count_1996_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_1996, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1996_W_ts_k45 <- simcoe_1996_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_1996 = sum(ind_W_1996, na.rm = TRUE))

#Step 4
simcoe_1996_ind_W_as1 <- rbind(simcoe_1996_W_ts_c9, simcoe_1996_W_ts_k42, simcoe_1996_W_ts_k45)

simcoe_1996_W_as1 <- simcoe_1996_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_1996 = sum(W_1996, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_1996 <- rbind(simcoe_mean_W_sp_pairs_1996_c9, simcoe_mean_W_sp_pairs_1996_k42, simcoe_mean_W_sp_pairs_1996_k45)

simcoe_lc_sp_pairs_1996 <- simcoe_mean_lc_sp_pairs_1996 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_1996_ind_W_as4 <- filter(simcoe_1996_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_1996_mp  <- na.omit(simcoe_1996_ind_W_as4)

simcoe_species_pairs_count_1996_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_1996_mp)

simcoe_mean_mp_sp_pairs_1996  <- simcoe_species_pairs_count_1996_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_1996, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_1996_sum_W_as4 <- simcoe_1996_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_1996 = sum(ind_W_1996, na.rm = TRUE))


#Step 3
simcoe_1996_sum_W_as4 <- simcoe_1996_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_1996 = sum(as4_1996, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_1996_ind_W_as5 <- filter(simcoe_1996_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_1996_cc <- na.omit(simcoe_1996_ind_W_as5)

simcoe_species_pairs_count_1996_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_1996_cc)

simcoe_mean_cc_sp_pairs_1996  <- simcoe_species_pairs_count_1996_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_1996, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_1996_sum_W_as5 <- simcoe_1996_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_1996 = sum(ind_W_1996, na.rm = TRUE))


#Step 3
simcoe_1996_sum_W_as5 <- simcoe_1996_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_1996 = sum(as5_1996, na.rm = TRUE))

simcoe_1996_W <- merge(simcoe_1996_sum_W_as4, simcoe_1996_W_ts, by = c('Year'))
simcoe_1996_W <- merge(simcoe_1996_W, simcoe_1996_W_as1, by = c('Year'))
simcoe_1996_W <- merge(simcoe_1996_W, simcoe_1996_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_1996 <- merge(simcoe_mean_W_sp_pairs_1996, simcoe_lc_sp_pairs_1996, by=c("Year"))
simcoe_sp_pair_ave_stab_1996 <- merge(simcoe_sp_pair_ave_stab_1996, simcoe_mean_mp_sp_pairs_1996, by=c("Year"))
simcoe_sp_pair_ave_stab_1996 <- merge(simcoe_sp_pair_ave_stab_1996, simcoe_mean_cc_sp_pairs_1996, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_1996_final <- merge(simcoe_1996_W, simcoe_1996_wa_pop_var_1996, by = c('Year'))
simcoe_1996_final <- simcoe_1996_final %>% group_by(Year) %>%
  mutate(gcv_1996 = lcv - W_1996) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_1996, 
            lcv = lcv,
            W = W_1996,
            lc_stab = as2_1996,
            mp_stab = as4_1996,
            cc_stab = as5_1996,
            mc_cv = sqrt(gcv_1996),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_1996/lcv + as4_1996/lcv + as5_1996/lcv,
            lc_asynchrony = as2_1996/lcv,
            mp_asynchrony = as4_1996/lcv,
            cc_asynchrony= as5_1996/lcv)

simcoe_1996_cv <- simcoe_1996_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_1996_final$mc_sum_mean_density <- simcoe_1996_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 1996 #####
simcoe_mean_corr_1996_c9 <- simcoe_species_pairs_count_1996_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1996, ind_W_1996)
simcoe_mean_corr_1996_k42 <- simcoe_species_pairs_count_1996_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1996, ind_W_1996)
simcoe_mean_corr_1996_k45 <- simcoe_species_pairs_count_1996_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1996, ind_W_1996)
simcoe_mean_corr_1996_lc <- rbind(simcoe_mean_corr_1996_c9, simcoe_mean_corr_1996_k42, simcoe_mean_corr_1996_k45)
simcoe_mean_corr_1996_mp <- simcoe_species_pairs_count_1996_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_1996, ind_W_1996)
simcoe_mean_corr_1996_cc <- simcoe_species_pairs_count_1996_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_1996, ind_W_1996)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_1996_lc[sample(1:nrow(simcoe_mean_corr_1996_lc),  nrow(simcoe_mean_corr_1996_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_1996),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_1996),
              mean_lc_stab = mean(ind_W_1996),
              sd_lc_stab = sd(ind_W_1996))
  
  lc_corr_1996 <- as.data.frame(mean_lc)
}
lc_corr_1996

mp_corr_1996 <- simcoe_mean_corr_1996_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_1996),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_1996),
                                                       mean_mp_stab = mean(ind_W_1996),
                                                       sd_mp_stab = sd(ind_W_1996))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_1996_cc[sample(1:nrow(simcoe_mean_corr_1996_cc), nrow(simcoe_mean_corr_1996_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_1996),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_1996),
              mean_cc_stab = mean(ind_W_1996),
              sd_cc_stab = sd(ind_W_1996))
  
  cc_corr_1996 <- as.data.frame(mean_cc)
}
cc_corr_1996


corr_1996 <- cbind(lc_corr_1996, mp_corr_1996)
corr_1996 <- cbind(corr_1996, cc_corr_1996)
corr_1996 <- corr_1996 %>% mutate(Year = 1996)





##### 1997 #####
simcoe_1997 <- subset(simcoe, Year == '1997')


#Transpose dataframe to add 0s within years
simcoe_1997_t1 <- dcast(simcoe_1997, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_1997_m <- melt(simcoe_1997_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_1997 <- subset(simcoe_1997_m, Station_ID == 'C9')
simcoe_c9_1997$Species <- paste("C9", simcoe_c9_1997$Species, sep=" ")
simcoe_k42_1997 <- subset(simcoe_1997_m, Station_ID == 'K42')
simcoe_k42_1997$Species <- paste("K42", simcoe_k42_1997$Species, sep=" ")
simcoe_k45_1997 <- subset(simcoe_1997_m, Station_ID == 'K45')
simcoe_k45_1997$Species <- paste("K45", simcoe_k45_1997$Species, sep=" ")

#recombine dataframes
simcoe_spec_1997 <- rbind(simcoe_c9_1997, simcoe_k42_1997, simcoe_k45_1997)
simcoe_spec_1997 <- simcoe_spec_1997 %>% select(- Station_ID)

#transpose dataframe
simcoe_1997_t2 <- dcast(simcoe_spec_1997, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_1997_t2 <- log(simcoe_1997_t2 + 1)
simcoe_1997_t2 <- simcoe_1997_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_1997_cv <- simcoe_1997_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_1997 = mean(Density),
            pop_sd_1997 = sd(Density)) %>%
  mutate(uw_pop_cv_1997 = pop_sd_1997/pop_mean_1997)

simcoe_1997_cv <- simcoe_1997_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_1997),
         mean_pop_cv = mean(uw_pop_cv_1997, na.rm = T),
         mean_pop_density = sum(pop_mean_1997),
         mean_pop_variance = sum(pop_sd_1997))

simcoe_1997_cv <- simcoe_1997_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_1997/mc_sum_mean_density)%>%
  mutate(w_pop_cv_1997 = (pop_mean_1997/mc_sum_mean_density)*uw_pop_cv_1997) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_1997_wa_pop_var_1997 <- simcoe_1997_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_1997 = sum(w_pop_cv_1997, na.rm = T),
            lcv = sum(w_pop_cv_1997, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_1997_cor_list_ts <- as.dist(round(cor(simcoe_1997_t2[]),2))
simcoe_1997_cor_ts <- stack(simcoe_1997_cor_list_ts, dim.names = TRUE)
simcoe_1997_cor_ts <- simcoe_1997_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_1997_cor - simcoe_1997_cv
simcoe_c9_1997_w_cv_ts <- subset(simcoe_1997_cv, Station_ID == 'C9')
simcoe_c9_1997_w_cv_ts$Species <- paste("C9", simcoe_c9_1997_w_cv_ts$Species, sep=" ")
simcoe_c9_1997_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_1997_w_cv_ts)
simcoe_k42_1997_w_cv_ts <- subset(simcoe_1997_cv, Station_ID == 'K42')
simcoe_k42_1997_w_cv_ts$Species <- paste("K42", simcoe_k42_1997_w_cv_ts$Species, sep=" ")
simcoe_k42_1997_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_1997_w_cv_ts)
simcoe_k45_1997_w_cv_ts <- subset(simcoe_1997_cv, Station_ID == 'K45')
simcoe_k45_1997_w_cv_ts$Species <- paste("K45", simcoe_k45_1997_w_cv_ts$Species, sep=" ")
simcoe_k45_1997_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_1997_w_cv_ts)

simcoe_spec_1997_w_sp_ts <- rbind(simcoe_c9_1997_w_cv_ts, simcoe_k42_1997_w_cv_ts, simcoe_k45_1997_w_cv_ts)
simcoe_spec_1997_w_sp_ts <- simcoe_spec_1997_w_sp_ts %>% select(Species, w_pop_cv_1997, Species_ID, Station_ID)
simcoe_spec_1997_w_sp1_ts <- simcoe_spec_1997_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_1997_sp1 = w_pop_cv_1997) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_1997_sp1, Species_ID, station_1)

simcoe_1997_cor_cv_ts <- merge(simcoe_1997_cor_ts, simcoe_spec_1997_w_sp1_ts, by = "species_1")

simcoe_spec_1997_w_sp2_ts <- simcoe_spec_1997_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_1997_sp2 = w_pop_cv_1997) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_1997_sp2, Species_ID, station_2)

simcoe_1997_cor_cv_ts <- merge(simcoe_1997_cor_cv_ts, simcoe_spec_1997_w_sp2_ts, by = "species_2")

simcoe_1997_cor_cv_ts_omit <- na.omit(simcoe_1997_cor_cv_ts)

simcoe_1997_ind_W_ts <- simcoe_1997_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_1997 = w_pop_cv_1997_sp1*w_pop_cv_1997_sp2,
         ind_W_1997 = (1 - corr)*(w_pop_cv_1997_sp1*w_pop_cv_1997_sp2),
         Year = 1997)

simcoe_1997_ind_W_ts$number_species_pairs <- nrow(simcoe_1997_ind_W_ts)

simcoe_mean_W_sp_pairs_1997 <- simcoe_1997_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_1997, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_1997_W_ts <- simcoe_1997_ind_W_ts %>% group_by(Year) %>%
  summarize(W_1997 = sum(ind_W_1997, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_1997_ind_W_ts_c9 <- filter(simcoe_1997_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_1997_ind_W_ts_k42 <- filter(simcoe_1997_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_1997_ind_W_ts_k45 <- filter(simcoe_1997_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_1997_c9 <- na.omit(simcoe_1997_ind_W_ts_c9)

simcoe_species_pairs_count_1997_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_1997_c9)

simcoe_mean_W_sp_pairs_1997_c9 <- simcoe_species_pairs_count_1997_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_1997, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1997_W_ts_c9 <- simcoe_1997_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_1997 = sum(ind_W_1997, na.rm = TRUE))

#K42
simcoe_species_pairs_count_1997_k42 <- na.omit(simcoe_1997_ind_W_ts_k42)

simcoe_species_pairs_count_1997_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_1997_k42)

simcoe_mean_W_sp_pairs_1997_k42 <- simcoe_species_pairs_count_1997_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_1997, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1997_W_ts_k42 <- simcoe_1997_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_1997 = sum(ind_W_1997, na.rm = TRUE))

#K45
simcoe_species_pairs_count_1997_k45 <- na.omit(simcoe_1997_ind_W_ts_k45)

simcoe_species_pairs_count_1997_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_1997_k45)

simcoe_mean_W_sp_pairs_1997_k45 <- simcoe_species_pairs_count_1997_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_1997, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1997_W_ts_k45 <- simcoe_1997_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_1997 = sum(ind_W_1997, na.rm = TRUE))

#Step 4
simcoe_1997_ind_W_as1 <- rbind(simcoe_1997_W_ts_c9, simcoe_1997_W_ts_k42, simcoe_1997_W_ts_k45)

simcoe_1997_W_as1 <- simcoe_1997_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_1997 = sum(W_1997, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_1997 <- rbind(simcoe_mean_W_sp_pairs_1997_c9, simcoe_mean_W_sp_pairs_1997_k42, simcoe_mean_W_sp_pairs_1997_k45)

simcoe_lc_sp_pairs_1997 <- simcoe_mean_lc_sp_pairs_1997 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_1997_ind_W_as4 <- filter(simcoe_1997_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_1997_mp  <- na.omit(simcoe_1997_ind_W_as4)

simcoe_species_pairs_count_1997_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_1997_mp)

simcoe_mean_mp_sp_pairs_1997  <- simcoe_species_pairs_count_1997_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_1997, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_1997_sum_W_as4 <- simcoe_1997_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_1997 = sum(ind_W_1997, na.rm = TRUE))


#Step 3
simcoe_1997_sum_W_as4 <- simcoe_1997_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_1997 = sum(as4_1997, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_1997_ind_W_as5 <- filter(simcoe_1997_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_1997_cc <- na.omit(simcoe_1997_ind_W_as5)

simcoe_species_pairs_count_1997_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_1997_cc)

simcoe_mean_cc_sp_pairs_1997  <- simcoe_species_pairs_count_1997_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_1997, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_1997_sum_W_as5 <- simcoe_1997_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_1997 = sum(ind_W_1997, na.rm = TRUE))


#Step 3
simcoe_1997_sum_W_as5 <- simcoe_1997_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_1997 = sum(as5_1997, na.rm = TRUE))

simcoe_1997_W <- merge(simcoe_1997_sum_W_as4, simcoe_1997_W_ts, by = c('Year'))
simcoe_1997_W <- merge(simcoe_1997_W, simcoe_1997_W_as1, by = c('Year'))
simcoe_1997_W <- merge(simcoe_1997_W, simcoe_1997_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_1997 <- merge(simcoe_mean_W_sp_pairs_1997, simcoe_lc_sp_pairs_1997, by=c("Year"))
simcoe_sp_pair_ave_stab_1997 <- merge(simcoe_sp_pair_ave_stab_1997, simcoe_mean_mp_sp_pairs_1997, by=c("Year"))
simcoe_sp_pair_ave_stab_1997 <- merge(simcoe_sp_pair_ave_stab_1997, simcoe_mean_cc_sp_pairs_1997, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_1997_final <- merge(simcoe_1997_W, simcoe_1997_wa_pop_var_1997, by = c('Year'))
simcoe_1997_final <- simcoe_1997_final %>% group_by(Year) %>%
  mutate(gcv_1997 = lcv - W_1997) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_1997, 
            lcv = lcv,
            W = W_1997,
            lc_stab = as2_1997,
            mp_stab = as4_1997,
            cc_stab = as5_1997,
            mc_cv = sqrt(gcv_1997),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_1997/lcv + as4_1997/lcv + as5_1997/lcv,
            lc_asynchrony = as2_1997/lcv,
            mp_asynchrony = as4_1997/lcv,
            cc_asynchrony= as5_1997/lcv)

simcoe_1997_cv <- simcoe_1997_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_1997_final$mc_sum_mean_density <- simcoe_1997_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 1997 #####
simcoe_mean_corr_1997_c9 <- simcoe_species_pairs_count_1997_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1997, ind_W_1997)
simcoe_mean_corr_1997_k42 <- simcoe_species_pairs_count_1997_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1997, ind_W_1997)
simcoe_mean_corr_1997_k45 <- simcoe_species_pairs_count_1997_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1997, ind_W_1997)
simcoe_mean_corr_1997_lc <- rbind(simcoe_mean_corr_1997_c9, simcoe_mean_corr_1997_k42, simcoe_mean_corr_1997_k45)
simcoe_mean_corr_1997_mp <- simcoe_species_pairs_count_1997_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_1997, ind_W_1997)
simcoe_mean_corr_1997_cc <- simcoe_species_pairs_count_1997_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_1997, ind_W_1997)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_1997_lc[sample(1:nrow(simcoe_mean_corr_1997_lc),  nrow(simcoe_mean_corr_1997_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_1997),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_1997),
              mean_lc_stab = mean(ind_W_1997),
              sd_lc_stab = sd(ind_W_1997))
  
  lc_corr_1997 <- as.data.frame(mean_lc)
}
lc_corr_1997

mp_corr_1997 <- simcoe_mean_corr_1997_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_1997),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_1997),
                                                       mean_mp_stab = mean(ind_W_1997),
                                                       sd_mp_stab = sd(ind_W_1997))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_1997_cc[sample(1:nrow(simcoe_mean_corr_1997_cc), nrow(simcoe_mean_corr_1997_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_1997),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_1997),
              mean_cc_stab = mean(ind_W_1997),
              sd_cc_stab = sd(ind_W_1997))
  
  cc_corr_1997 <- as.data.frame(mean_cc)
}
cc_corr_1997


corr_1997 <- cbind(lc_corr_1997, mp_corr_1997)
corr_1997 <- cbind(corr_1997, cc_corr_1997)
corr_1997 <- corr_1997 %>% mutate(Year = 1997)






##### 1998 #####
simcoe_1998 <- subset(simcoe, Year == '1998')


#Transpose dataframe to add 0s within years
simcoe_1998_t1 <- dcast(simcoe_1998, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_1998_m <- melt(simcoe_1998_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_1998 <- subset(simcoe_1998_m, Station_ID == 'C9')
simcoe_c9_1998$Species <- paste("C9", simcoe_c9_1998$Species, sep=" ")
simcoe_k42_1998 <- subset(simcoe_1998_m, Station_ID == 'K42')
simcoe_k42_1998$Species <- paste("K42", simcoe_k42_1998$Species, sep=" ")
simcoe_k45_1998 <- subset(simcoe_1998_m, Station_ID == 'K45')
simcoe_k45_1998$Species <- paste("K45", simcoe_k45_1998$Species, sep=" ")

#recombine dataframes
simcoe_spec_1998 <- rbind(simcoe_c9_1998, simcoe_k42_1998, simcoe_k45_1998)
simcoe_spec_1998 <- simcoe_spec_1998 %>% select(- Station_ID)

#transpose dataframe
simcoe_1998_t2 <- dcast(simcoe_spec_1998, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_1998_t2 <- log(simcoe_1998_t2 + 1)
simcoe_1998_t2 <- simcoe_1998_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_1998_cv <- simcoe_1998_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_1998 = mean(Density),
            pop_sd_1998 = sd(Density)) %>%
  mutate(uw_pop_cv_1998 = pop_sd_1998/pop_mean_1998)

simcoe_1998_cv <- simcoe_1998_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_1998),
         mean_pop_cv = mean(uw_pop_cv_1998, na.rm = T),
         mean_pop_density = sum(pop_mean_1998),
         mean_pop_variance = sum(pop_sd_1998))

simcoe_1998_cv <- simcoe_1998_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_1998/mc_sum_mean_density)%>%
  mutate(w_pop_cv_1998 = (pop_mean_1998/mc_sum_mean_density)*uw_pop_cv_1998) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_1998_wa_pop_var_1998 <- simcoe_1998_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_1998 = sum(w_pop_cv_1998, na.rm = T),
            lcv = sum(w_pop_cv_1998, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_1998_cor_list_ts <- as.dist(round(cor(simcoe_1998_t2[]),2))
simcoe_1998_cor_ts <- stack(simcoe_1998_cor_list_ts, dim.names = TRUE)
simcoe_1998_cor_ts <- simcoe_1998_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_1998_cor - simcoe_1998_cv
simcoe_c9_1998_w_cv_ts <- subset(simcoe_1998_cv, Station_ID == 'C9')
simcoe_c9_1998_w_cv_ts$Species <- paste("C9", simcoe_c9_1998_w_cv_ts$Species, sep=" ")
simcoe_c9_1998_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_1998_w_cv_ts)
simcoe_k42_1998_w_cv_ts <- subset(simcoe_1998_cv, Station_ID == 'K42')
simcoe_k42_1998_w_cv_ts$Species <- paste("K42", simcoe_k42_1998_w_cv_ts$Species, sep=" ")
simcoe_k42_1998_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_1998_w_cv_ts)
simcoe_k45_1998_w_cv_ts <- subset(simcoe_1998_cv, Station_ID == 'K45')
simcoe_k45_1998_w_cv_ts$Species <- paste("K45", simcoe_k45_1998_w_cv_ts$Species, sep=" ")
simcoe_k45_1998_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_1998_w_cv_ts)

simcoe_spec_1998_w_sp_ts <- rbind(simcoe_c9_1998_w_cv_ts, simcoe_k42_1998_w_cv_ts, simcoe_k45_1998_w_cv_ts)
simcoe_spec_1998_w_sp_ts <- simcoe_spec_1998_w_sp_ts %>% select(Species, w_pop_cv_1998, Species_ID, Station_ID)
simcoe_spec_1998_w_sp1_ts <- simcoe_spec_1998_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_1998_sp1 = w_pop_cv_1998) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_1998_sp1, Species_ID, station_1)

simcoe_1998_cor_cv_ts <- merge(simcoe_1998_cor_ts, simcoe_spec_1998_w_sp1_ts, by = "species_1")

simcoe_spec_1998_w_sp2_ts <- simcoe_spec_1998_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_1998_sp2 = w_pop_cv_1998) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_1998_sp2, Species_ID, station_2)

simcoe_1998_cor_cv_ts <- merge(simcoe_1998_cor_cv_ts, simcoe_spec_1998_w_sp2_ts, by = "species_2")

simcoe_1998_cor_cv_ts_omit <- na.omit(simcoe_1998_cor_cv_ts)

simcoe_1998_ind_W_ts <- simcoe_1998_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_1998 = w_pop_cv_1998_sp1*w_pop_cv_1998_sp2,
         ind_W_1998 = (1 - corr)*(w_pop_cv_1998_sp1*w_pop_cv_1998_sp2),
         Year = 1998)

simcoe_1998_ind_W_ts$number_species_pairs <- nrow(simcoe_1998_ind_W_ts)

simcoe_mean_W_sp_pairs_1998 <- simcoe_1998_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_1998, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_1998_W_ts <- simcoe_1998_ind_W_ts %>% group_by(Year) %>%
  summarize(W_1998 = sum(ind_W_1998, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_1998_ind_W_ts_c9 <- filter(simcoe_1998_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_1998_ind_W_ts_k42 <- filter(simcoe_1998_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_1998_ind_W_ts_k45 <- filter(simcoe_1998_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_1998_c9 <- na.omit(simcoe_1998_ind_W_ts_c9)

simcoe_species_pairs_count_1998_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_1998_c9)

simcoe_mean_W_sp_pairs_1998_c9 <- simcoe_species_pairs_count_1998_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_1998, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1998_W_ts_c9 <- simcoe_1998_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_1998 = sum(ind_W_1998, na.rm = TRUE))

#K42
simcoe_species_pairs_count_1998_k42 <- na.omit(simcoe_1998_ind_W_ts_k42)

simcoe_species_pairs_count_1998_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_1998_k42)

simcoe_mean_W_sp_pairs_1998_k42 <- simcoe_species_pairs_count_1998_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_1998, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1998_W_ts_k42 <- simcoe_1998_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_1998 = sum(ind_W_1998, na.rm = TRUE))

#K45
simcoe_species_pairs_count_1998_k45 <- na.omit(simcoe_1998_ind_W_ts_k45)

simcoe_species_pairs_count_1998_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_1998_k45)

simcoe_mean_W_sp_pairs_1998_k45 <- simcoe_species_pairs_count_1998_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_1998, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1998_W_ts_k45 <- simcoe_1998_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_1998 = sum(ind_W_1998, na.rm = TRUE))

#Step 4
simcoe_1998_ind_W_as1 <- rbind(simcoe_1998_W_ts_c9, simcoe_1998_W_ts_k42, simcoe_1998_W_ts_k45)

simcoe_1998_W_as1 <- simcoe_1998_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_1998 = sum(W_1998, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_1998 <- rbind(simcoe_mean_W_sp_pairs_1998_c9, simcoe_mean_W_sp_pairs_1998_k42, simcoe_mean_W_sp_pairs_1998_k45)

simcoe_lc_sp_pairs_1998 <- simcoe_mean_lc_sp_pairs_1998 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_1998_ind_W_as4 <- filter(simcoe_1998_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_1998_mp  <- na.omit(simcoe_1998_ind_W_as4)

simcoe_species_pairs_count_1998_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_1998_mp)

simcoe_mean_mp_sp_pairs_1998  <- simcoe_species_pairs_count_1998_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_1998, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_1998_sum_W_as4 <- simcoe_1998_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_1998 = sum(ind_W_1998, na.rm = TRUE))


#Step 3
simcoe_1998_sum_W_as4 <- simcoe_1998_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_1998 = sum(as4_1998, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_1998_ind_W_as5 <- filter(simcoe_1998_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_1998_cc <- na.omit(simcoe_1998_ind_W_as5)

simcoe_species_pairs_count_1998_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_1998_cc)

simcoe_mean_cc_sp_pairs_1998  <- simcoe_species_pairs_count_1998_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_1998, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_1998_sum_W_as5 <- simcoe_1998_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_1998 = sum(ind_W_1998, na.rm = TRUE))


#Step 3
simcoe_1998_sum_W_as5 <- simcoe_1998_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_1998 = sum(as5_1998, na.rm = TRUE))

simcoe_1998_W <- merge(simcoe_1998_sum_W_as4, simcoe_1998_W_ts, by = c('Year'))
simcoe_1998_W <- merge(simcoe_1998_W, simcoe_1998_W_as1, by = c('Year'))
simcoe_1998_W <- merge(simcoe_1998_W, simcoe_1998_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_1998 <- merge(simcoe_mean_W_sp_pairs_1998, simcoe_lc_sp_pairs_1998, by=c("Year"))
simcoe_sp_pair_ave_stab_1998 <- merge(simcoe_sp_pair_ave_stab_1998, simcoe_mean_mp_sp_pairs_1998, by=c("Year"))
simcoe_sp_pair_ave_stab_1998 <- merge(simcoe_sp_pair_ave_stab_1998, simcoe_mean_cc_sp_pairs_1998, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_1998_final <- merge(simcoe_1998_W, simcoe_1998_wa_pop_var_1998, by = c('Year'))
simcoe_1998_final <- simcoe_1998_final %>% group_by(Year) %>%
  mutate(gcv_1998 = lcv - W_1998) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_1998, 
            lcv = lcv,
            W = W_1998,
            lc_stab = as2_1998,
            mp_stab = as4_1998,
            cc_stab = as5_1998,
            mc_cv = sqrt(gcv_1998),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_1998/lcv + as4_1998/lcv + as5_1998/lcv,
            lc_asynchrony = as2_1998/lcv,
            mp_asynchrony = as4_1998/lcv,
            cc_asynchrony= as5_1998/lcv)

simcoe_1998_cv <- simcoe_1998_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_1998_final$mc_sum_mean_density <- simcoe_1998_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 1998 #####
simcoe_mean_corr_1998_c9 <- simcoe_species_pairs_count_1998_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1998, ind_W_1998)
simcoe_mean_corr_1998_k42 <- simcoe_species_pairs_count_1998_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1998, ind_W_1998)
simcoe_mean_corr_1998_k45 <- simcoe_species_pairs_count_1998_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1998, ind_W_1998)
simcoe_mean_corr_1998_lc <- rbind(simcoe_mean_corr_1998_c9, simcoe_mean_corr_1998_k42, simcoe_mean_corr_1998_k45)
simcoe_mean_corr_1998_mp <- simcoe_species_pairs_count_1998_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_1998, ind_W_1998)
simcoe_mean_corr_1998_cc <- simcoe_species_pairs_count_1998_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_1998, ind_W_1998)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_1998_lc[sample(1:nrow(simcoe_mean_corr_1998_lc),  nrow(simcoe_mean_corr_1998_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_1998),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_1998),
              mean_lc_stab = mean(ind_W_1998),
              sd_lc_stab = sd(ind_W_1998))
  
  lc_corr_1998 <- as.data.frame(mean_lc)
}
lc_corr_1998

mp_corr_1998 <- simcoe_mean_corr_1998_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_1998),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_1998),
                                                       mean_mp_stab = mean(ind_W_1998),
                                                       sd_mp_stab = sd(ind_W_1998))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_1998_cc[sample(1:nrow(simcoe_mean_corr_1998_cc), nrow(simcoe_mean_corr_1998_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_1998),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_1998),
              mean_cc_stab = mean(ind_W_1998),
              sd_cc_stab = sd(ind_W_1998))
  
  cc_corr_1998 <- as.data.frame(mean_cc)
}
cc_corr_1998


corr_1998 <- cbind(lc_corr_1998, mp_corr_1998)
corr_1998 <- cbind(corr_1998, cc_corr_1998)
corr_1998 <- corr_1998 %>% mutate(Year = 1998)




##### 1999 #####
simcoe_1999 <- subset(simcoe, Year == '1999')


#Transpose dataframe to add 0s within years
simcoe_1999_t1 <- dcast(simcoe_1999, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_1999_m <- melt(simcoe_1999_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_1999 <- subset(simcoe_1999_m, Station_ID == 'C9')
simcoe_c9_1999$Species <- paste("C9", simcoe_c9_1999$Species, sep=" ")
simcoe_k42_1999 <- subset(simcoe_1999_m, Station_ID == 'K42')
simcoe_k42_1999$Species <- paste("K42", simcoe_k42_1999$Species, sep=" ")
simcoe_k45_1999 <- subset(simcoe_1999_m, Station_ID == 'K45')
simcoe_k45_1999$Species <- paste("K45", simcoe_k45_1999$Species, sep=" ")

#recombine dataframes
simcoe_spec_1999 <- rbind(simcoe_c9_1999, simcoe_k42_1999, simcoe_k45_1999)
simcoe_spec_1999 <- simcoe_spec_1999 %>% select(- Station_ID)

#transpose dataframe
simcoe_1999_t2 <- dcast(simcoe_spec_1999, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_1999_t2 <- log(simcoe_1999_t2 + 1)
simcoe_1999_t2 <- simcoe_1999_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_1999_cv <- simcoe_1999_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_1999 = mean(Density),
            pop_sd_1999 = sd(Density)) %>%
  mutate(uw_pop_cv_1999 = pop_sd_1999/pop_mean_1999)

simcoe_1999_cv <- simcoe_1999_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_1999),
         mean_pop_cv = mean(uw_pop_cv_1999, na.rm = T),
         mean_pop_density = sum(pop_mean_1999),
         mean_pop_variance = sum(pop_sd_1999))

simcoe_1999_cv <- simcoe_1999_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_1999/mc_sum_mean_density)%>%
  mutate(w_pop_cv_1999 = (pop_mean_1999/mc_sum_mean_density)*uw_pop_cv_1999) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_1999_wa_pop_var_1999 <- simcoe_1999_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_1999 = sum(w_pop_cv_1999, na.rm = T),
            lcv = sum(w_pop_cv_1999, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_1999_cor_list_ts <- as.dist(round(cor(simcoe_1999_t2[]),2))
simcoe_1999_cor_ts <- stack(simcoe_1999_cor_list_ts, dim.names = TRUE)
simcoe_1999_cor_ts <- simcoe_1999_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_1999_cor - simcoe_1999_cv
simcoe_c9_1999_w_cv_ts <- subset(simcoe_1999_cv, Station_ID == 'C9')
simcoe_c9_1999_w_cv_ts$Species <- paste("C9", simcoe_c9_1999_w_cv_ts$Species, sep=" ")
simcoe_c9_1999_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_1999_w_cv_ts)
simcoe_k42_1999_w_cv_ts <- subset(simcoe_1999_cv, Station_ID == 'K42')
simcoe_k42_1999_w_cv_ts$Species <- paste("K42", simcoe_k42_1999_w_cv_ts$Species, sep=" ")
simcoe_k42_1999_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_1999_w_cv_ts)
simcoe_k45_1999_w_cv_ts <- subset(simcoe_1999_cv, Station_ID == 'K45')
simcoe_k45_1999_w_cv_ts$Species <- paste("K45", simcoe_k45_1999_w_cv_ts$Species, sep=" ")
simcoe_k45_1999_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_1999_w_cv_ts)

simcoe_spec_1999_w_sp_ts <- rbind(simcoe_c9_1999_w_cv_ts, simcoe_k42_1999_w_cv_ts, simcoe_k45_1999_w_cv_ts)
simcoe_spec_1999_w_sp_ts <- simcoe_spec_1999_w_sp_ts %>% select(Species, w_pop_cv_1999, Species_ID, Station_ID)
simcoe_spec_1999_w_sp1_ts <- simcoe_spec_1999_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_1999_sp1 = w_pop_cv_1999) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_1999_sp1, Species_ID, station_1)

simcoe_1999_cor_cv_ts <- merge(simcoe_1999_cor_ts, simcoe_spec_1999_w_sp1_ts, by = "species_1")

simcoe_spec_1999_w_sp2_ts <- simcoe_spec_1999_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_1999_sp2 = w_pop_cv_1999) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_1999_sp2, Species_ID, station_2)

simcoe_1999_cor_cv_ts <- merge(simcoe_1999_cor_cv_ts, simcoe_spec_1999_w_sp2_ts, by = "species_2")

simcoe_1999_cor_cv_ts_omit <- na.omit(simcoe_1999_cor_cv_ts)

simcoe_1999_ind_W_ts <- simcoe_1999_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_1999 = w_pop_cv_1999_sp1*w_pop_cv_1999_sp2,
         ind_W_1999 = (1 - corr)*(w_pop_cv_1999_sp1*w_pop_cv_1999_sp2),
         Year = 1999)

simcoe_1999_ind_W_ts$number_species_pairs <- nrow(simcoe_1999_ind_W_ts)

simcoe_mean_W_sp_pairs_1999 <- simcoe_1999_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_1999, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_1999_W_ts <- simcoe_1999_ind_W_ts %>% group_by(Year) %>%
  summarize(W_1999 = sum(ind_W_1999, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_1999_ind_W_ts_c9 <- filter(simcoe_1999_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_1999_ind_W_ts_k42 <- filter(simcoe_1999_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_1999_ind_W_ts_k45 <- filter(simcoe_1999_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_1999_c9 <- na.omit(simcoe_1999_ind_W_ts_c9)

simcoe_species_pairs_count_1999_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_1999_c9)

simcoe_mean_W_sp_pairs_1999_c9 <- simcoe_species_pairs_count_1999_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_1999, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1999_W_ts_c9 <- simcoe_1999_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_1999 = sum(ind_W_1999, na.rm = TRUE))

#K42
simcoe_species_pairs_count_1999_k42 <- na.omit(simcoe_1999_ind_W_ts_k42)

simcoe_species_pairs_count_1999_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_1999_k42)

simcoe_mean_W_sp_pairs_1999_k42 <- simcoe_species_pairs_count_1999_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_1999, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1999_W_ts_k42 <- simcoe_1999_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_1999 = sum(ind_W_1999, na.rm = TRUE))

#K45
simcoe_species_pairs_count_1999_k45 <- na.omit(simcoe_1999_ind_W_ts_k45)

simcoe_species_pairs_count_1999_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_1999_k45)

simcoe_mean_W_sp_pairs_1999_k45 <- simcoe_species_pairs_count_1999_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_1999, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_1999_W_ts_k45 <- simcoe_1999_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_1999 = sum(ind_W_1999, na.rm = TRUE))

#Step 4
simcoe_1999_ind_W_as1 <- rbind(simcoe_1999_W_ts_c9, simcoe_1999_W_ts_k42, simcoe_1999_W_ts_k45)

simcoe_1999_W_as1 <- simcoe_1999_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_1999 = sum(W_1999, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_1999 <- rbind(simcoe_mean_W_sp_pairs_1999_c9, simcoe_mean_W_sp_pairs_1999_k42, simcoe_mean_W_sp_pairs_1999_k45)

simcoe_lc_sp_pairs_1999 <- simcoe_mean_lc_sp_pairs_1999 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_1999_ind_W_as4 <- filter(simcoe_1999_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_1999_mp  <- na.omit(simcoe_1999_ind_W_as4)

simcoe_species_pairs_count_1999_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_1999_mp)

simcoe_mean_mp_sp_pairs_1999  <- simcoe_species_pairs_count_1999_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_1999, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_1999_sum_W_as4 <- simcoe_1999_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_1999 = sum(ind_W_1999, na.rm = TRUE))


#Step 3
simcoe_1999_sum_W_as4 <- simcoe_1999_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_1999 = sum(as4_1999, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_1999_ind_W_as5 <- filter(simcoe_1999_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_1999_cc <- na.omit(simcoe_1999_ind_W_as5)

simcoe_species_pairs_count_1999_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_1999_cc)

simcoe_mean_cc_sp_pairs_1999  <- simcoe_species_pairs_count_1999_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_1999, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_1999_sum_W_as5 <- simcoe_1999_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_1999 = sum(ind_W_1999, na.rm = TRUE))


#Step 3
simcoe_1999_sum_W_as5 <- simcoe_1999_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_1999 = sum(as5_1999, na.rm = TRUE))

simcoe_1999_W <- merge(simcoe_1999_sum_W_as4, simcoe_1999_W_ts, by = c('Year'))
simcoe_1999_W <- merge(simcoe_1999_W, simcoe_1999_W_as1, by = c('Year'))
simcoe_1999_W <- merge(simcoe_1999_W, simcoe_1999_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_1999 <- merge(simcoe_mean_W_sp_pairs_1999, simcoe_lc_sp_pairs_1999, by=c("Year"))
simcoe_sp_pair_ave_stab_1999 <- merge(simcoe_sp_pair_ave_stab_1999, simcoe_mean_mp_sp_pairs_1999, by=c("Year"))
simcoe_sp_pair_ave_stab_1999 <- merge(simcoe_sp_pair_ave_stab_1999, simcoe_mean_cc_sp_pairs_1999, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_1999_final <- merge(simcoe_1999_W, simcoe_1999_wa_pop_var_1999, by = c('Year'))
simcoe_1999_final <- simcoe_1999_final %>% group_by(Year) %>%
  mutate(gcv_1999 = lcv - W_1999) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_1999, 
            lcv = lcv,
            W = W_1999,
            lc_stab = as2_1999,
            mp_stab = as4_1999,
            cc_stab = as5_1999,
            mc_cv = sqrt(gcv_1999),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_1999/lcv + as4_1999/lcv + as5_1999/lcv,
            lc_asynchrony = as2_1999/lcv,
            mp_asynchrony = as4_1999/lcv,
            cc_asynchrony= as5_1999/lcv)

simcoe_1999_cv <- simcoe_1999_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_1999_final$mc_sum_mean_density <- simcoe_1999_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 1999 #####
simcoe_mean_corr_1999_c9 <- simcoe_species_pairs_count_1999_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1999, ind_W_1999)
simcoe_mean_corr_1999_k42 <- simcoe_species_pairs_count_1999_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1999, ind_W_1999)
simcoe_mean_corr_1999_k45 <- simcoe_species_pairs_count_1999_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_1999, ind_W_1999)
simcoe_mean_corr_1999_lc <- rbind(simcoe_mean_corr_1999_c9, simcoe_mean_corr_1999_k42, simcoe_mean_corr_1999_k45)
simcoe_mean_corr_1999_mp <- simcoe_species_pairs_count_1999_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_1999, ind_W_1999)
simcoe_mean_corr_1999_cc <- simcoe_species_pairs_count_1999_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_1999, ind_W_1999)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_1999_lc[sample(1:nrow(simcoe_mean_corr_1999_lc),  nrow(simcoe_mean_corr_1999_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_1999),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_1999),
              mean_lc_stab = mean(ind_W_1999),
              sd_lc_stab = sd(ind_W_1999))
  
  lc_corr_1999 <- as.data.frame(mean_lc)
}
lc_corr_1999

mp_corr_1999 <- simcoe_mean_corr_1999_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_1999),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_1999),
                                                       mean_mp_stab = mean(ind_W_1999),
                                                       sd_mp_stab = sd(ind_W_1999))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_1999_cc[sample(1:nrow(simcoe_mean_corr_1999_cc), nrow(simcoe_mean_corr_1999_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_1999),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_1999),
              mean_cc_stab = mean(ind_W_1999),
              sd_cc_stab = sd(ind_W_1999))
  
  cc_corr_1999 <- as.data.frame(mean_cc)
}
cc_corr_1999


corr_1999 <- cbind(lc_corr_1999, mp_corr_1999)
corr_1999 <- cbind(corr_1999, cc_corr_1999)
corr_1999 <- corr_1999 %>% mutate(Year = 1999)





##### 2000 #####
simcoe_2000 <- subset(simcoe, Year == '2000')


#Transpose dataframe to add 0s within years
simcoe_2000_t1 <- dcast(simcoe_2000, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2000_m <- melt(simcoe_2000_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_2000 <- subset(simcoe_2000_m, Station_ID == 'C9')
simcoe_c9_2000$Species <- paste("C9", simcoe_c9_2000$Species, sep=" ")
simcoe_k42_2000 <- subset(simcoe_2000_m, Station_ID == 'K42')
simcoe_k42_2000$Species <- paste("K42", simcoe_k42_2000$Species, sep=" ")
simcoe_k45_2000 <- subset(simcoe_2000_m, Station_ID == 'K45')
simcoe_k45_2000$Species <- paste("K45", simcoe_k45_2000$Species, sep=" ")

#recombine dataframes
simcoe_spec_2000 <- rbind(simcoe_c9_2000, simcoe_k42_2000, simcoe_k45_2000)
simcoe_spec_2000 <- simcoe_spec_2000 %>% select(- Station_ID)

#transpose dataframe
simcoe_2000_t2 <- dcast(simcoe_spec_2000, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2000_t2 <- log(simcoe_2000_t2 + 1)
simcoe_2000_t2 <- simcoe_2000_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_2000_cv <- simcoe_2000_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_2000 = mean(Density),
            pop_sd_2000 = sd(Density)) %>%
  mutate(uw_pop_cv_2000 = pop_sd_2000/pop_mean_2000)

simcoe_2000_cv <- simcoe_2000_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_2000),
         mean_pop_cv = mean(uw_pop_cv_2000, na.rm = T),
         mean_pop_density = sum(pop_mean_2000),
         mean_pop_variance = sum(pop_sd_2000))

simcoe_2000_cv <- simcoe_2000_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_2000/mc_sum_mean_density)%>%
  mutate(w_pop_cv_2000 = (pop_mean_2000/mc_sum_mean_density)*uw_pop_cv_2000) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_2000_wa_pop_var_2000 <- simcoe_2000_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_2000 = sum(w_pop_cv_2000, na.rm = T),
            lcv = sum(w_pop_cv_2000, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_2000_cor_list_ts <- as.dist(round(cor(simcoe_2000_t2[]),2))
simcoe_2000_cor_ts <- stack(simcoe_2000_cor_list_ts, dim.names = TRUE)
simcoe_2000_cor_ts <- simcoe_2000_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_2000_cor - simcoe_2000_cv
simcoe_c9_2000_w_cv_ts <- subset(simcoe_2000_cv, Station_ID == 'C9')
simcoe_c9_2000_w_cv_ts$Species <- paste("C9", simcoe_c9_2000_w_cv_ts$Species, sep=" ")
simcoe_c9_2000_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_2000_w_cv_ts)
simcoe_k42_2000_w_cv_ts <- subset(simcoe_2000_cv, Station_ID == 'K42')
simcoe_k42_2000_w_cv_ts$Species <- paste("K42", simcoe_k42_2000_w_cv_ts$Species, sep=" ")
simcoe_k42_2000_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_2000_w_cv_ts)
simcoe_k45_2000_w_cv_ts <- subset(simcoe_2000_cv, Station_ID == 'K45')
simcoe_k45_2000_w_cv_ts$Species <- paste("K45", simcoe_k45_2000_w_cv_ts$Species, sep=" ")
simcoe_k45_2000_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_2000_w_cv_ts)

simcoe_spec_2000_w_sp_ts <- rbind(simcoe_c9_2000_w_cv_ts, simcoe_k42_2000_w_cv_ts, simcoe_k45_2000_w_cv_ts)
simcoe_spec_2000_w_sp_ts <- simcoe_spec_2000_w_sp_ts %>% select(Species, w_pop_cv_2000, Species_ID, Station_ID)
simcoe_spec_2000_w_sp1_ts <- simcoe_spec_2000_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_2000_sp1 = w_pop_cv_2000) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_2000_sp1, Species_ID, station_1)

simcoe_2000_cor_cv_ts <- merge(simcoe_2000_cor_ts, simcoe_spec_2000_w_sp1_ts, by = "species_1")

simcoe_spec_2000_w_sp2_ts <- simcoe_spec_2000_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_2000_sp2 = w_pop_cv_2000) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_2000_sp2, Species_ID, station_2)

simcoe_2000_cor_cv_ts <- merge(simcoe_2000_cor_cv_ts, simcoe_spec_2000_w_sp2_ts, by = "species_2")

simcoe_2000_cor_cv_ts_omit <- na.omit(simcoe_2000_cor_cv_ts)

simcoe_2000_ind_W_ts <- simcoe_2000_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_2000 = w_pop_cv_2000_sp1*w_pop_cv_2000_sp2,
         ind_W_2000 = (1 - corr)*(w_pop_cv_2000_sp1*w_pop_cv_2000_sp2),
         Year = 2000)

simcoe_2000_ind_W_ts$number_species_pairs <- nrow(simcoe_2000_ind_W_ts)

simcoe_mean_W_sp_pairs_2000 <- simcoe_2000_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_2000, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_2000_W_ts <- simcoe_2000_ind_W_ts %>% group_by(Year) %>%
  summarize(W_2000 = sum(ind_W_2000, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_2000_ind_W_ts_c9 <- filter(simcoe_2000_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_2000_ind_W_ts_k42 <- filter(simcoe_2000_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_2000_ind_W_ts_k45 <- filter(simcoe_2000_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_2000_c9 <- na.omit(simcoe_2000_ind_W_ts_c9)

simcoe_species_pairs_count_2000_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_2000_c9)

simcoe_mean_W_sp_pairs_2000_c9 <- simcoe_species_pairs_count_2000_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_2000, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2000_W_ts_c9 <- simcoe_2000_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_2000 = sum(ind_W_2000, na.rm = TRUE))

#K42
simcoe_species_pairs_count_2000_k42 <- na.omit(simcoe_2000_ind_W_ts_k42)

simcoe_species_pairs_count_2000_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_2000_k42)

simcoe_mean_W_sp_pairs_2000_k42 <- simcoe_species_pairs_count_2000_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_2000, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2000_W_ts_k42 <- simcoe_2000_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_2000 = sum(ind_W_2000, na.rm = TRUE))

#K45
simcoe_species_pairs_count_2000_k45 <- na.omit(simcoe_2000_ind_W_ts_k45)

simcoe_species_pairs_count_2000_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_2000_k45)

simcoe_mean_W_sp_pairs_2000_k45 <- simcoe_species_pairs_count_2000_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_2000, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2000_W_ts_k45 <- simcoe_2000_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_2000 = sum(ind_W_2000, na.rm = TRUE))

#Step 4
simcoe_2000_ind_W_as1 <- rbind(simcoe_2000_W_ts_c9, simcoe_2000_W_ts_k42, simcoe_2000_W_ts_k45)

simcoe_2000_W_as1 <- simcoe_2000_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_2000 = sum(W_2000, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_2000 <- rbind(simcoe_mean_W_sp_pairs_2000_c9, simcoe_mean_W_sp_pairs_2000_k42, simcoe_mean_W_sp_pairs_2000_k45)

simcoe_lc_sp_pairs_2000 <- simcoe_mean_lc_sp_pairs_2000 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_2000_ind_W_as4 <- filter(simcoe_2000_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_2000_mp  <- na.omit(simcoe_2000_ind_W_as4)

simcoe_species_pairs_count_2000_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_2000_mp)

simcoe_mean_mp_sp_pairs_2000  <- simcoe_species_pairs_count_2000_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_2000, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2000_sum_W_as4 <- simcoe_2000_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_2000 = sum(ind_W_2000, na.rm = TRUE))


#Step 3
simcoe_2000_sum_W_as4 <- simcoe_2000_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_2000 = sum(as4_2000, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_2000_ind_W_as5 <- filter(simcoe_2000_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_2000_cc <- na.omit(simcoe_2000_ind_W_as5)

simcoe_species_pairs_count_2000_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_2000_cc)

simcoe_mean_cc_sp_pairs_2000  <- simcoe_species_pairs_count_2000_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_2000, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2000_sum_W_as5 <- simcoe_2000_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_2000 = sum(ind_W_2000, na.rm = TRUE))


#Step 3
simcoe_2000_sum_W_as5 <- simcoe_2000_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_2000 = sum(as5_2000, na.rm = TRUE))

simcoe_2000_W <- merge(simcoe_2000_sum_W_as4, simcoe_2000_W_ts, by = c('Year'))
simcoe_2000_W <- merge(simcoe_2000_W, simcoe_2000_W_as1, by = c('Year'))
simcoe_2000_W <- merge(simcoe_2000_W, simcoe_2000_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_2000 <- merge(simcoe_mean_W_sp_pairs_2000, simcoe_lc_sp_pairs_2000, by=c("Year"))
simcoe_sp_pair_ave_stab_2000 <- merge(simcoe_sp_pair_ave_stab_2000, simcoe_mean_mp_sp_pairs_2000, by=c("Year"))
simcoe_sp_pair_ave_stab_2000 <- merge(simcoe_sp_pair_ave_stab_2000, simcoe_mean_cc_sp_pairs_2000, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_2000_final <- merge(simcoe_2000_W, simcoe_2000_wa_pop_var_2000, by = c('Year'))
simcoe_2000_final <- simcoe_2000_final %>% group_by(Year) %>%
  mutate(gcv_2000 = lcv - W_2000) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_2000, 
            lcv = lcv,
            W = W_2000,
            lc_stab = as2_2000,
            mp_stab = as4_2000,
            cc_stab = as5_2000,
            mc_cv = sqrt(gcv_2000),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_2000/lcv + as4_2000/lcv + as5_2000/lcv,
            lc_asynchrony = as2_2000/lcv,
            mp_asynchrony = as4_2000/lcv,
            cc_asynchrony= as5_2000/lcv)

simcoe_2000_cv <- simcoe_2000_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_2000_final$mc_sum_mean_density <- simcoe_2000_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 2000 #####
simcoe_mean_corr_2000_c9 <- simcoe_species_pairs_count_2000_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2000, ind_W_2000)
simcoe_mean_corr_2000_k42 <- simcoe_species_pairs_count_2000_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2000, ind_W_2000)
simcoe_mean_corr_2000_k45 <- simcoe_species_pairs_count_2000_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2000, ind_W_2000)
simcoe_mean_corr_2000_lc <- rbind(simcoe_mean_corr_2000_c9, simcoe_mean_corr_2000_k42, simcoe_mean_corr_2000_k45)
simcoe_mean_corr_2000_mp <- simcoe_species_pairs_count_2000_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_2000, ind_W_2000)
simcoe_mean_corr_2000_cc <- simcoe_species_pairs_count_2000_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_2000, ind_W_2000)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_2000_lc[sample(1:nrow(simcoe_mean_corr_2000_lc),  nrow(simcoe_mean_corr_2000_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_2000),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_2000),
              mean_lc_stab = mean(ind_W_2000),
              sd_lc_stab = sd(ind_W_2000))
  
  lc_corr_2000 <- as.data.frame(mean_lc)
}
lc_corr_2000

mp_corr_2000 <- simcoe_mean_corr_2000_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_2000),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_2000),
                                                       mean_mp_stab = mean(ind_W_2000),
                                                       sd_mp_stab = sd(ind_W_2000))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_2000_cc[sample(1:nrow(simcoe_mean_corr_2000_cc), nrow(simcoe_mean_corr_2000_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_2000),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_2000),
              mean_cc_stab = mean(ind_W_2000),
              sd_cc_stab = sd(ind_W_2000))
  
  cc_corr_2000 <- as.data.frame(mean_cc)
}
cc_corr_2000


corr_2000 <- cbind(lc_corr_2000, mp_corr_2000)
corr_2000 <- cbind(corr_2000, cc_corr_2000)
corr_2000 <- corr_2000 %>% mutate(Year = 2000)





##### 2001 #####
simcoe_2001 <- subset(simcoe, Year == '2001')


#Transpose dataframe to add 0s within years
simcoe_2001_t1 <- dcast(simcoe_2001, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2001_m <- melt(simcoe_2001_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_2001 <- subset(simcoe_2001_m, Station_ID == 'C9')
simcoe_c9_2001$Species <- paste("C9", simcoe_c9_2001$Species, sep=" ")
simcoe_k42_2001 <- subset(simcoe_2001_m, Station_ID == 'K42')
simcoe_k42_2001$Species <- paste("K42", simcoe_k42_2001$Species, sep=" ")
simcoe_k45_2001 <- subset(simcoe_2001_m, Station_ID == 'K45')
simcoe_k45_2001$Species <- paste("K45", simcoe_k45_2001$Species, sep=" ")

#recombine dataframes
simcoe_spec_2001 <- rbind(simcoe_c9_2001, simcoe_k42_2001, simcoe_k45_2001)
simcoe_spec_2001 <- simcoe_spec_2001 %>% select(- Station_ID)

#transpose dataframe
simcoe_2001_t2 <- dcast(simcoe_spec_2001, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2001_t2 <- log(simcoe_2001_t2 + 1)
simcoe_2001_t2 <- simcoe_2001_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_2001_cv <- simcoe_2001_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_2001 = mean(Density),
            pop_sd_2001 = sd(Density)) %>%
  mutate(uw_pop_cv_2001 = pop_sd_2001/pop_mean_2001)

simcoe_2001_cv <- simcoe_2001_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_2001),
         mean_pop_cv = mean(uw_pop_cv_2001, na.rm = T),
         mean_pop_density = sum(pop_mean_2001),
         mean_pop_variance = sum(pop_sd_2001))

simcoe_2001_cv <- simcoe_2001_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_2001/mc_sum_mean_density)%>%
  mutate(w_pop_cv_2001 = (pop_mean_2001/mc_sum_mean_density)*uw_pop_cv_2001) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_2001_wa_pop_var_2001 <- simcoe_2001_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_2001 = sum(w_pop_cv_2001, na.rm = T),
            lcv = sum(w_pop_cv_2001, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_2001_cor_list_ts <- as.dist(round(cor(simcoe_2001_t2[]),2))
simcoe_2001_cor_ts <- stack(simcoe_2001_cor_list_ts, dim.names = TRUE)
simcoe_2001_cor_ts <- simcoe_2001_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_2001_cor - simcoe_2001_cv
simcoe_c9_2001_w_cv_ts <- subset(simcoe_2001_cv, Station_ID == 'C9')
simcoe_c9_2001_w_cv_ts$Species <- paste("C9", simcoe_c9_2001_w_cv_ts$Species, sep=" ")
simcoe_c9_2001_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_2001_w_cv_ts)
simcoe_k42_2001_w_cv_ts <- subset(simcoe_2001_cv, Station_ID == 'K42')
simcoe_k42_2001_w_cv_ts$Species <- paste("K42", simcoe_k42_2001_w_cv_ts$Species, sep=" ")
simcoe_k42_2001_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_2001_w_cv_ts)
simcoe_k45_2001_w_cv_ts <- subset(simcoe_2001_cv, Station_ID == 'K45')
simcoe_k45_2001_w_cv_ts$Species <- paste("K45", simcoe_k45_2001_w_cv_ts$Species, sep=" ")
simcoe_k45_2001_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_2001_w_cv_ts)

simcoe_spec_2001_w_sp_ts <- rbind(simcoe_c9_2001_w_cv_ts, simcoe_k42_2001_w_cv_ts, simcoe_k45_2001_w_cv_ts)
simcoe_spec_2001_w_sp_ts <- simcoe_spec_2001_w_sp_ts %>% select(Species, w_pop_cv_2001, Species_ID, Station_ID)
simcoe_spec_2001_w_sp1_ts <- simcoe_spec_2001_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_2001_sp1 = w_pop_cv_2001) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_2001_sp1, Species_ID, station_1)

simcoe_2001_cor_cv_ts <- merge(simcoe_2001_cor_ts, simcoe_spec_2001_w_sp1_ts, by = "species_1")

simcoe_spec_2001_w_sp2_ts <- simcoe_spec_2001_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_2001_sp2 = w_pop_cv_2001) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_2001_sp2, Species_ID, station_2)

simcoe_2001_cor_cv_ts <- merge(simcoe_2001_cor_cv_ts, simcoe_spec_2001_w_sp2_ts, by = "species_2")

simcoe_2001_cor_cv_ts_omit <- na.omit(simcoe_2001_cor_cv_ts)

simcoe_2001_ind_W_ts <- simcoe_2001_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_2001 = w_pop_cv_2001_sp1*w_pop_cv_2001_sp2,
         ind_W_2001 = (1 - corr)*(w_pop_cv_2001_sp1*w_pop_cv_2001_sp2),
         Year = 2001)

simcoe_2001_ind_W_ts$number_species_pairs <- nrow(simcoe_2001_ind_W_ts)

simcoe_mean_W_sp_pairs_2001 <- simcoe_2001_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_2001, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_2001_W_ts <- simcoe_2001_ind_W_ts %>% group_by(Year) %>%
  summarize(W_2001 = sum(ind_W_2001, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_2001_ind_W_ts_c9 <- filter(simcoe_2001_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_2001_ind_W_ts_k42 <- filter(simcoe_2001_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_2001_ind_W_ts_k45 <- filter(simcoe_2001_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_2001_c9 <- na.omit(simcoe_2001_ind_W_ts_c9)

simcoe_species_pairs_count_2001_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_2001_c9)

simcoe_mean_W_sp_pairs_2001_c9 <- simcoe_species_pairs_count_2001_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_2001, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2001_W_ts_c9 <- simcoe_2001_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_2001 = sum(ind_W_2001, na.rm = TRUE))

#K42
simcoe_species_pairs_count_2001_k42 <- na.omit(simcoe_2001_ind_W_ts_k42)

simcoe_species_pairs_count_2001_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_2001_k42)

simcoe_mean_W_sp_pairs_2001_k42 <- simcoe_species_pairs_count_2001_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_2001, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2001_W_ts_k42 <- simcoe_2001_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_2001 = sum(ind_W_2001, na.rm = TRUE))

#K45
simcoe_species_pairs_count_2001_k45 <- na.omit(simcoe_2001_ind_W_ts_k45)

simcoe_species_pairs_count_2001_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_2001_k45)

simcoe_mean_W_sp_pairs_2001_k45 <- simcoe_species_pairs_count_2001_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_2001, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2001_W_ts_k45 <- simcoe_2001_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_2001 = sum(ind_W_2001, na.rm = TRUE))

#Step 4
simcoe_2001_ind_W_as1 <- rbind(simcoe_2001_W_ts_c9, simcoe_2001_W_ts_k42, simcoe_2001_W_ts_k45)

simcoe_2001_W_as1 <- simcoe_2001_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_2001 = sum(W_2001, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_2001 <- rbind(simcoe_mean_W_sp_pairs_2001_c9, simcoe_mean_W_sp_pairs_2001_k42, simcoe_mean_W_sp_pairs_2001_k45)

simcoe_lc_sp_pairs_2001 <- simcoe_mean_lc_sp_pairs_2001 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_2001_ind_W_as4 <- filter(simcoe_2001_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_2001_mp  <- na.omit(simcoe_2001_ind_W_as4)

simcoe_species_pairs_count_2001_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_2001_mp)

simcoe_mean_mp_sp_pairs_2001  <- simcoe_species_pairs_count_2001_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_2001, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2001_sum_W_as4 <- simcoe_2001_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_2001 = sum(ind_W_2001, na.rm = TRUE))


#Step 3
simcoe_2001_sum_W_as4 <- simcoe_2001_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_2001 = sum(as4_2001, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_2001_ind_W_as5 <- filter(simcoe_2001_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_2001_cc <- na.omit(simcoe_2001_ind_W_as5)

simcoe_species_pairs_count_2001_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_2001_cc)

simcoe_mean_cc_sp_pairs_2001  <- simcoe_species_pairs_count_2001_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_2001, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2001_sum_W_as5 <- simcoe_2001_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_2001 = sum(ind_W_2001, na.rm = TRUE))


#Step 3
simcoe_2001_sum_W_as5 <- simcoe_2001_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_2001 = sum(as5_2001, na.rm = TRUE))

simcoe_2001_W <- merge(simcoe_2001_sum_W_as4, simcoe_2001_W_ts, by = c('Year'))
simcoe_2001_W <- merge(simcoe_2001_W, simcoe_2001_W_as1, by = c('Year'))
simcoe_2001_W <- merge(simcoe_2001_W, simcoe_2001_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_2001 <- merge(simcoe_mean_W_sp_pairs_2001, simcoe_lc_sp_pairs_2001, by=c("Year"))
simcoe_sp_pair_ave_stab_2001 <- merge(simcoe_sp_pair_ave_stab_2001, simcoe_mean_mp_sp_pairs_2001, by=c("Year"))
simcoe_sp_pair_ave_stab_2001 <- merge(simcoe_sp_pair_ave_stab_2001, simcoe_mean_cc_sp_pairs_2001, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_2001_final <- merge(simcoe_2001_W, simcoe_2001_wa_pop_var_2001, by = c('Year'))
simcoe_2001_final <- simcoe_2001_final %>% group_by(Year) %>%
  mutate(gcv_2001 = lcv - W_2001) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_2001, 
            lcv = lcv,
            W = W_2001,
            lc_stab = as2_2001,
            mp_stab = as4_2001,
            cc_stab = as5_2001,
            mc_cv = sqrt(gcv_2001),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_2001/lcv + as4_2001/lcv + as5_2001/lcv,
            lc_asynchrony = as2_2001/lcv,
            mp_asynchrony = as4_2001/lcv,
            cc_asynchrony= as5_2001/lcv)

simcoe_2001_cv <- simcoe_2001_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_2001_final$mc_sum_mean_density <- simcoe_2001_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 2001 #####
simcoe_mean_corr_2001_c9 <- simcoe_species_pairs_count_2001_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2001, ind_W_2001)
simcoe_mean_corr_2001_k42 <- simcoe_species_pairs_count_2001_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2001, ind_W_2001)
simcoe_mean_corr_2001_k45 <- simcoe_species_pairs_count_2001_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2001, ind_W_2001)
simcoe_mean_corr_2001_lc <- rbind(simcoe_mean_corr_2001_c9, simcoe_mean_corr_2001_k42, simcoe_mean_corr_2001_k45)
simcoe_mean_corr_2001_mp <- simcoe_species_pairs_count_2001_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_2001, ind_W_2001)
simcoe_mean_corr_2001_cc <- simcoe_species_pairs_count_2001_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_2001, ind_W_2001)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_2001_lc[sample(1:nrow(simcoe_mean_corr_2001_lc),  nrow(simcoe_mean_corr_2001_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_2001),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_2001),
              mean_lc_stab = mean(ind_W_2001),
              sd_lc_stab = sd(ind_W_2001))
  
  lc_corr_2001 <- as.data.frame(mean_lc)
}
lc_corr_2001

mp_corr_2001 <- simcoe_mean_corr_2001_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_2001),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_2001),
                                                       mean_mp_stab = mean(ind_W_2001),
                                                       sd_mp_stab = sd(ind_W_2001))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_2001_cc[sample(1:nrow(simcoe_mean_corr_2001_cc), nrow(simcoe_mean_corr_2001_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_2001),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_2001),
              mean_cc_stab = mean(ind_W_2001),
              sd_cc_stab = sd(ind_W_2001))
  
  cc_corr_2001 <- as.data.frame(mean_cc)
}
cc_corr_2001


corr_2001 <- cbind(lc_corr_2001, mp_corr_2001)
corr_2001 <- cbind(corr_2001, cc_corr_2001)
corr_2001 <- corr_2001 %>% mutate(Year = 2001)





##### 2002 #####
simcoe_2002 <- subset(simcoe, Year == '2002')


#Transpose dataframe to add 0s within years
simcoe_2002_t1 <- dcast(simcoe_2002, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2002_m <- melt(simcoe_2002_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_2002 <- subset(simcoe_2002_m, Station_ID == 'C9')
simcoe_c9_2002$Species <- paste("C9", simcoe_c9_2002$Species, sep=" ")
simcoe_k42_2002 <- subset(simcoe_2002_m, Station_ID == 'K42')
simcoe_k42_2002$Species <- paste("K42", simcoe_k42_2002$Species, sep=" ")
simcoe_k45_2002 <- subset(simcoe_2002_m, Station_ID == 'K45')
simcoe_k45_2002$Species <- paste("K45", simcoe_k45_2002$Species, sep=" ")

#recombine dataframes
simcoe_spec_2002 <- rbind(simcoe_c9_2002, simcoe_k42_2002, simcoe_k45_2002)
simcoe_spec_2002 <- simcoe_spec_2002 %>% select(- Station_ID)

#transpose dataframe
simcoe_2002_t2 <- dcast(simcoe_spec_2002, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2002_t2 <- log(simcoe_2002_t2 + 1)
simcoe_2002_t2 <- simcoe_2002_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_2002_cv <- simcoe_2002_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_2002 = mean(Density),
            pop_sd_2002 = sd(Density)) %>%
  mutate(uw_pop_cv_2002 = pop_sd_2002/pop_mean_2002)

simcoe_2002_cv <- simcoe_2002_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_2002),
         mean_pop_cv = mean(uw_pop_cv_2002, na.rm = T),
         mean_pop_density = sum(pop_mean_2002),
         mean_pop_variance = sum(pop_sd_2002))

simcoe_2002_cv <- simcoe_2002_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_2002/mc_sum_mean_density)%>%
  mutate(w_pop_cv_2002 = (pop_mean_2002/mc_sum_mean_density)*uw_pop_cv_2002) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_2002_wa_pop_var_2002 <- simcoe_2002_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_2002 = sum(w_pop_cv_2002, na.rm = T),
            lcv = sum(w_pop_cv_2002, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_2002_cor_list_ts <- as.dist(round(cor(simcoe_2002_t2[]),2))
simcoe_2002_cor_ts <- stack(simcoe_2002_cor_list_ts, dim.names = TRUE)
simcoe_2002_cor_ts <- simcoe_2002_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_2002_cor - simcoe_2002_cv
simcoe_c9_2002_w_cv_ts <- subset(simcoe_2002_cv, Station_ID == 'C9')
simcoe_c9_2002_w_cv_ts$Species <- paste("C9", simcoe_c9_2002_w_cv_ts$Species, sep=" ")
simcoe_c9_2002_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_2002_w_cv_ts)
simcoe_k42_2002_w_cv_ts <- subset(simcoe_2002_cv, Station_ID == 'K42')
simcoe_k42_2002_w_cv_ts$Species <- paste("K42", simcoe_k42_2002_w_cv_ts$Species, sep=" ")
simcoe_k42_2002_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_2002_w_cv_ts)
simcoe_k45_2002_w_cv_ts <- subset(simcoe_2002_cv, Station_ID == 'K45')
simcoe_k45_2002_w_cv_ts$Species <- paste("K45", simcoe_k45_2002_w_cv_ts$Species, sep=" ")
simcoe_k45_2002_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_2002_w_cv_ts)

simcoe_spec_2002_w_sp_ts <- rbind(simcoe_c9_2002_w_cv_ts, simcoe_k42_2002_w_cv_ts, simcoe_k45_2002_w_cv_ts)
simcoe_spec_2002_w_sp_ts <- simcoe_spec_2002_w_sp_ts %>% select(Species, w_pop_cv_2002, Species_ID, Station_ID)
simcoe_spec_2002_w_sp1_ts <- simcoe_spec_2002_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_2002_sp1 = w_pop_cv_2002) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_2002_sp1, Species_ID, station_1)

simcoe_2002_cor_cv_ts <- merge(simcoe_2002_cor_ts, simcoe_spec_2002_w_sp1_ts, by = "species_1")

simcoe_spec_2002_w_sp2_ts <- simcoe_spec_2002_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_2002_sp2 = w_pop_cv_2002) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_2002_sp2, Species_ID, station_2)

simcoe_2002_cor_cv_ts <- merge(simcoe_2002_cor_cv_ts, simcoe_spec_2002_w_sp2_ts, by = "species_2")

simcoe_2002_cor_cv_ts_omit <- na.omit(simcoe_2002_cor_cv_ts)

simcoe_2002_ind_W_ts <- simcoe_2002_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_2002 = w_pop_cv_2002_sp1*w_pop_cv_2002_sp2,
         ind_W_2002 = (1 - corr)*(w_pop_cv_2002_sp1*w_pop_cv_2002_sp2),
         Year = 2002)

simcoe_2002_ind_W_ts$number_species_pairs <- nrow(simcoe_2002_ind_W_ts)

simcoe_mean_W_sp_pairs_2002 <- simcoe_2002_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_2002, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_2002_W_ts <- simcoe_2002_ind_W_ts %>% group_by(Year) %>%
  summarize(W_2002 = sum(ind_W_2002, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_2002_ind_W_ts_c9 <- filter(simcoe_2002_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_2002_ind_W_ts_k42 <- filter(simcoe_2002_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_2002_ind_W_ts_k45 <- filter(simcoe_2002_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_2002_c9 <- na.omit(simcoe_2002_ind_W_ts_c9)

simcoe_species_pairs_count_2002_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_2002_c9)

simcoe_mean_W_sp_pairs_2002_c9 <- simcoe_species_pairs_count_2002_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_2002, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2002_W_ts_c9 <- simcoe_2002_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_2002 = sum(ind_W_2002, na.rm = TRUE))

#K42
simcoe_species_pairs_count_2002_k42 <- na.omit(simcoe_2002_ind_W_ts_k42)

simcoe_species_pairs_count_2002_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_2002_k42)

simcoe_mean_W_sp_pairs_2002_k42 <- simcoe_species_pairs_count_2002_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_2002, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2002_W_ts_k42 <- simcoe_2002_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_2002 = sum(ind_W_2002, na.rm = TRUE))

#K45
simcoe_species_pairs_count_2002_k45 <- na.omit(simcoe_2002_ind_W_ts_k45)

simcoe_species_pairs_count_2002_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_2002_k45)

simcoe_mean_W_sp_pairs_2002_k45 <- simcoe_species_pairs_count_2002_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_2002, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2002_W_ts_k45 <- simcoe_2002_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_2002 = sum(ind_W_2002, na.rm = TRUE))

#Step 4
simcoe_2002_ind_W_as1 <- rbind(simcoe_2002_W_ts_c9, simcoe_2002_W_ts_k42, simcoe_2002_W_ts_k45)

simcoe_2002_W_as1 <- simcoe_2002_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_2002 = sum(W_2002, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_2002 <- rbind(simcoe_mean_W_sp_pairs_2002_c9, simcoe_mean_W_sp_pairs_2002_k42, simcoe_mean_W_sp_pairs_2002_k45)

simcoe_lc_sp_pairs_2002 <- simcoe_mean_lc_sp_pairs_2002 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_2002_ind_W_as4 <- filter(simcoe_2002_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_2002_mp  <- na.omit(simcoe_2002_ind_W_as4)

simcoe_species_pairs_count_2002_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_2002_mp)

simcoe_mean_mp_sp_pairs_2002  <- simcoe_species_pairs_count_2002_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_2002, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2002_sum_W_as4 <- simcoe_2002_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_2002 = sum(ind_W_2002, na.rm = TRUE))


#Step 3
simcoe_2002_sum_W_as4 <- simcoe_2002_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_2002 = sum(as4_2002, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_2002_ind_W_as5 <- filter(simcoe_2002_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_2002_cc <- na.omit(simcoe_2002_ind_W_as5)

simcoe_species_pairs_count_2002_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_2002_cc)

simcoe_mean_cc_sp_pairs_2002  <- simcoe_species_pairs_count_2002_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_2002, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2002_sum_W_as5 <- simcoe_2002_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_2002 = sum(ind_W_2002, na.rm = TRUE))


#Step 3
simcoe_2002_sum_W_as5 <- simcoe_2002_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_2002 = sum(as5_2002, na.rm = TRUE))

simcoe_2002_W <- merge(simcoe_2002_sum_W_as4, simcoe_2002_W_ts, by = c('Year'))
simcoe_2002_W <- merge(simcoe_2002_W, simcoe_2002_W_as1, by = c('Year'))
simcoe_2002_W <- merge(simcoe_2002_W, simcoe_2002_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_2002 <- merge(simcoe_mean_W_sp_pairs_2002, simcoe_lc_sp_pairs_2002, by=c("Year"))
simcoe_sp_pair_ave_stab_2002 <- merge(simcoe_sp_pair_ave_stab_2002, simcoe_mean_mp_sp_pairs_2002, by=c("Year"))
simcoe_sp_pair_ave_stab_2002 <- merge(simcoe_sp_pair_ave_stab_2002, simcoe_mean_cc_sp_pairs_2002, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_2002_final <- merge(simcoe_2002_W, simcoe_2002_wa_pop_var_2002, by = c('Year'))
simcoe_2002_final <- simcoe_2002_final %>% group_by(Year) %>%
  mutate(gcv_2002 = lcv - W_2002) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_2002, 
            lcv = lcv,
            W = W_2002,
            lc_stab = as2_2002,
            mp_stab = as4_2002,
            cc_stab = as5_2002,
            mc_cv = sqrt(gcv_2002),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_2002/lcv + as4_2002/lcv + as5_2002/lcv,
            lc_asynchrony = as2_2002/lcv,
            mp_asynchrony = as4_2002/lcv,
            cc_asynchrony= as5_2002/lcv)

simcoe_2002_cv <- simcoe_2002_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_2002_final$mc_sum_mean_density <- simcoe_2002_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 2002 #####
simcoe_mean_corr_2002_c9 <- simcoe_species_pairs_count_2002_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2002, ind_W_2002)
simcoe_mean_corr_2002_k42 <- simcoe_species_pairs_count_2002_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2002, ind_W_2002)
simcoe_mean_corr_2002_k45 <- simcoe_species_pairs_count_2002_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2002, ind_W_2002)
simcoe_mean_corr_2002_lc <- rbind(simcoe_mean_corr_2002_c9, simcoe_mean_corr_2002_k42, simcoe_mean_corr_2002_k45)
simcoe_mean_corr_2002_mp <- simcoe_species_pairs_count_2002_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_2002, ind_W_2002)
simcoe_mean_corr_2002_cc <- simcoe_species_pairs_count_2002_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_2002, ind_W_2002)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_2002_lc[sample(1:nrow(simcoe_mean_corr_2002_lc),  nrow(simcoe_mean_corr_2002_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_2002),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_2002),
              mean_lc_stab = mean(ind_W_2002),
              sd_lc_stab = sd(ind_W_2002))
  
  lc_corr_2002 <- as.data.frame(mean_lc)
}
lc_corr_2002

mp_corr_2002 <- simcoe_mean_corr_2002_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_2002),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_2002),
                                                       mean_mp_stab = mean(ind_W_2002),
                                                       sd_mp_stab = sd(ind_W_2002))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_2002_cc[sample(1:nrow(simcoe_mean_corr_2002_cc), nrow(simcoe_mean_corr_2002_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_2002),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_2002),
              mean_cc_stab = mean(ind_W_2002),
              sd_cc_stab = sd(ind_W_2002))
  
  cc_corr_2002 <- as.data.frame(mean_cc)
}
cc_corr_2002


corr_2002 <- cbind(lc_corr_2002, mp_corr_2002)
corr_2002 <- cbind(corr_2002, cc_corr_2002)
corr_2002 <- corr_2002 %>% mutate(Year = 2002)





##### 2003 #####
simcoe_2003 <- subset(simcoe, Year == '2003')


#Transpose dataframe to add 0s within years
simcoe_2003_t1 <- dcast(simcoe_2003, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2003_m <- melt(simcoe_2003_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_2003 <- subset(simcoe_2003_m, Station_ID == 'C9')
simcoe_c9_2003$Species <- paste("C9", simcoe_c9_2003$Species, sep=" ")
simcoe_k42_2003 <- subset(simcoe_2003_m, Station_ID == 'K42')
simcoe_k42_2003$Species <- paste("K42", simcoe_k42_2003$Species, sep=" ")
simcoe_k45_2003 <- subset(simcoe_2003_m, Station_ID == 'K45')
simcoe_k45_2003$Species <- paste("K45", simcoe_k45_2003$Species, sep=" ")

#recombine dataframes
simcoe_spec_2003 <- rbind(simcoe_c9_2003, simcoe_k42_2003, simcoe_k45_2003)
simcoe_spec_2003 <- simcoe_spec_2003 %>% select(- Station_ID)

#transpose dataframe
simcoe_2003_t2 <- dcast(simcoe_spec_2003, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2003_t2 <- log(simcoe_2003_t2 + 1)
simcoe_2003_t2 <- simcoe_2003_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_2003_cv <- simcoe_2003_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_2003 = mean(Density),
            pop_sd_2003 = sd(Density)) %>%
  mutate(uw_pop_cv_2003 = pop_sd_2003/pop_mean_2003)

simcoe_2003_cv <- simcoe_2003_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_2003),
         mean_pop_cv = mean(uw_pop_cv_2003, na.rm = T),
         mean_pop_density = sum(pop_mean_2003),
         mean_pop_variance = sum(pop_sd_2003))

simcoe_2003_cv <- simcoe_2003_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_2003/mc_sum_mean_density)%>%
  mutate(w_pop_cv_2003 = (pop_mean_2003/mc_sum_mean_density)*uw_pop_cv_2003) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_2003_wa_pop_var_2003 <- simcoe_2003_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_2003 = sum(w_pop_cv_2003, na.rm = T),
            lcv = sum(w_pop_cv_2003, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_2003_cor_list_ts <- as.dist(round(cor(simcoe_2003_t2[]),2))
simcoe_2003_cor_ts <- stack(simcoe_2003_cor_list_ts, dim.names = TRUE)
simcoe_2003_cor_ts <- simcoe_2003_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_2003_cor - simcoe_2003_cv
simcoe_c9_2003_w_cv_ts <- subset(simcoe_2003_cv, Station_ID == 'C9')
simcoe_c9_2003_w_cv_ts$Species <- paste("C9", simcoe_c9_2003_w_cv_ts$Species, sep=" ")
simcoe_c9_2003_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_2003_w_cv_ts)
simcoe_k42_2003_w_cv_ts <- subset(simcoe_2003_cv, Station_ID == 'K42')
simcoe_k42_2003_w_cv_ts$Species <- paste("K42", simcoe_k42_2003_w_cv_ts$Species, sep=" ")
simcoe_k42_2003_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_2003_w_cv_ts)
simcoe_k45_2003_w_cv_ts <- subset(simcoe_2003_cv, Station_ID == 'K45')
simcoe_k45_2003_w_cv_ts$Species <- paste("K45", simcoe_k45_2003_w_cv_ts$Species, sep=" ")
simcoe_k45_2003_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_2003_w_cv_ts)

simcoe_spec_2003_w_sp_ts <- rbind(simcoe_c9_2003_w_cv_ts, simcoe_k42_2003_w_cv_ts, simcoe_k45_2003_w_cv_ts)
simcoe_spec_2003_w_sp_ts <- simcoe_spec_2003_w_sp_ts %>% select(Species, w_pop_cv_2003, Species_ID, Station_ID)
simcoe_spec_2003_w_sp1_ts <- simcoe_spec_2003_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_2003_sp1 = w_pop_cv_2003) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_2003_sp1, Species_ID, station_1)

simcoe_2003_cor_cv_ts <- merge(simcoe_2003_cor_ts, simcoe_spec_2003_w_sp1_ts, by = "species_1")

simcoe_spec_2003_w_sp2_ts <- simcoe_spec_2003_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_2003_sp2 = w_pop_cv_2003) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_2003_sp2, Species_ID, station_2)

simcoe_2003_cor_cv_ts <- merge(simcoe_2003_cor_cv_ts, simcoe_spec_2003_w_sp2_ts, by = "species_2")

simcoe_2003_cor_cv_ts_omit <- na.omit(simcoe_2003_cor_cv_ts)

simcoe_2003_ind_W_ts <- simcoe_2003_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_2003 = w_pop_cv_2003_sp1*w_pop_cv_2003_sp2,
         ind_W_2003 = (1 - corr)*(w_pop_cv_2003_sp1*w_pop_cv_2003_sp2),
         Year = 2003)

simcoe_2003_ind_W_ts$number_species_pairs <- nrow(simcoe_2003_ind_W_ts)

simcoe_mean_W_sp_pairs_2003 <- simcoe_2003_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_2003, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_2003_W_ts <- simcoe_2003_ind_W_ts %>% group_by(Year) %>%
  summarize(W_2003 = sum(ind_W_2003, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_2003_ind_W_ts_c9 <- filter(simcoe_2003_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_2003_ind_W_ts_k42 <- filter(simcoe_2003_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_2003_ind_W_ts_k45 <- filter(simcoe_2003_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_2003_c9 <- na.omit(simcoe_2003_ind_W_ts_c9)

simcoe_species_pairs_count_2003_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_2003_c9)

simcoe_mean_W_sp_pairs_2003_c9 <- simcoe_species_pairs_count_2003_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_2003, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2003_W_ts_c9 <- simcoe_2003_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_2003 = sum(ind_W_2003, na.rm = TRUE))

#K42
simcoe_species_pairs_count_2003_k42 <- na.omit(simcoe_2003_ind_W_ts_k42)

simcoe_species_pairs_count_2003_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_2003_k42)

simcoe_mean_W_sp_pairs_2003_k42 <- simcoe_species_pairs_count_2003_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_2003, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2003_W_ts_k42 <- simcoe_2003_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_2003 = sum(ind_W_2003, na.rm = TRUE))

#K45
simcoe_species_pairs_count_2003_k45 <- na.omit(simcoe_2003_ind_W_ts_k45)

simcoe_species_pairs_count_2003_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_2003_k45)

simcoe_mean_W_sp_pairs_2003_k45 <- simcoe_species_pairs_count_2003_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_2003, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2003_W_ts_k45 <- simcoe_2003_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_2003 = sum(ind_W_2003, na.rm = TRUE))

#Step 4
simcoe_2003_ind_W_as1 <- rbind(simcoe_2003_W_ts_c9, simcoe_2003_W_ts_k42, simcoe_2003_W_ts_k45)

simcoe_2003_W_as1 <- simcoe_2003_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_2003 = sum(W_2003, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_2003 <- rbind(simcoe_mean_W_sp_pairs_2003_c9, simcoe_mean_W_sp_pairs_2003_k42, simcoe_mean_W_sp_pairs_2003_k45)

simcoe_lc_sp_pairs_2003 <- simcoe_mean_lc_sp_pairs_2003 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_2003_ind_W_as4 <- filter(simcoe_2003_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_2003_mp  <- na.omit(simcoe_2003_ind_W_as4)

simcoe_species_pairs_count_2003_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_2003_mp)

simcoe_mean_mp_sp_pairs_2003  <- simcoe_species_pairs_count_2003_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_2003, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2003_sum_W_as4 <- simcoe_2003_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_2003 = sum(ind_W_2003, na.rm = TRUE))


#Step 3
simcoe_2003_sum_W_as4 <- simcoe_2003_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_2003 = sum(as4_2003, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_2003_ind_W_as5 <- filter(simcoe_2003_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_2003_cc <- na.omit(simcoe_2003_ind_W_as5)

simcoe_species_pairs_count_2003_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_2003_cc)

simcoe_mean_cc_sp_pairs_2003  <- simcoe_species_pairs_count_2003_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_2003, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2003_sum_W_as5 <- simcoe_2003_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_2003 = sum(ind_W_2003, na.rm = TRUE))


#Step 3
simcoe_2003_sum_W_as5 <- simcoe_2003_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_2003 = sum(as5_2003, na.rm = TRUE))

simcoe_2003_W <- merge(simcoe_2003_sum_W_as4, simcoe_2003_W_ts, by = c('Year'))
simcoe_2003_W <- merge(simcoe_2003_W, simcoe_2003_W_as1, by = c('Year'))
simcoe_2003_W <- merge(simcoe_2003_W, simcoe_2003_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_2003 <- merge(simcoe_mean_W_sp_pairs_2003, simcoe_lc_sp_pairs_2003, by=c("Year"))
simcoe_sp_pair_ave_stab_2003 <- merge(simcoe_sp_pair_ave_stab_2003, simcoe_mean_mp_sp_pairs_2003, by=c("Year"))
simcoe_sp_pair_ave_stab_2003 <- merge(simcoe_sp_pair_ave_stab_2003, simcoe_mean_cc_sp_pairs_2003, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_2003_final <- merge(simcoe_2003_W, simcoe_2003_wa_pop_var_2003, by = c('Year'))
simcoe_2003_final <- simcoe_2003_final %>% group_by(Year) %>%
  mutate(gcv_2003 = lcv - W_2003) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_2003, 
            lcv = lcv,
            W = W_2003,
            lc_stab = as2_2003,
            mp_stab = as4_2003,
            cc_stab = as5_2003,
            mc_cv = sqrt(gcv_2003),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_2003/lcv + as4_2003/lcv + as5_2003/lcv,
            lc_asynchrony = as2_2003/lcv,
            mp_asynchrony = as4_2003/lcv,
            cc_asynchrony= as5_2003/lcv)

simcoe_2003_cv <- simcoe_2003_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_2003_final$mc_sum_mean_density <- simcoe_2003_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 2003 #####
simcoe_mean_corr_2003_c9 <- simcoe_species_pairs_count_2003_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2003, ind_W_2003)
simcoe_mean_corr_2003_k42 <- simcoe_species_pairs_count_2003_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2003, ind_W_2003)
simcoe_mean_corr_2003_k45 <- simcoe_species_pairs_count_2003_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2003, ind_W_2003)
simcoe_mean_corr_2003_lc <- rbind(simcoe_mean_corr_2003_c9, simcoe_mean_corr_2003_k42, simcoe_mean_corr_2003_k45)
simcoe_mean_corr_2003_mp <- simcoe_species_pairs_count_2003_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_2003, ind_W_2003)
simcoe_mean_corr_2003_cc <- simcoe_species_pairs_count_2003_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_2003, ind_W_2003)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_2003_lc[sample(1:nrow(simcoe_mean_corr_2003_lc),  nrow(simcoe_mean_corr_2003_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_2003),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_2003),
              mean_lc_stab = mean(ind_W_2003),
              sd_lc_stab = sd(ind_W_2003))
  
  lc_corr_2003 <- as.data.frame(mean_lc)
}
lc_corr_2003

mp_corr_2003 <- simcoe_mean_corr_2003_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_2003),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_2003),
                                                       mean_mp_stab = mean(ind_W_2003),
                                                       sd_mp_stab = sd(ind_W_2003))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_2003_cc[sample(1:nrow(simcoe_mean_corr_2003_cc), nrow(simcoe_mean_corr_2003_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_2003),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_2003),
              mean_cc_stab = mean(ind_W_2003),
              sd_cc_stab = sd(ind_W_2003))
  
  cc_corr_2003 <- as.data.frame(mean_cc)
}
cc_corr_2003


corr_2003 <- cbind(lc_corr_2003, mp_corr_2003)
corr_2003 <- cbind(corr_2003, cc_corr_2003)
corr_2003 <- corr_2003 %>% mutate(Year = 2003)





##### 2004 #####
simcoe_2004 <- subset(simcoe, Year == '2004')


#Transpose dataframe to add 0s within years
simcoe_2004_t1 <- dcast(simcoe_2004, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2004_m <- melt(simcoe_2004_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_2004 <- subset(simcoe_2004_m, Station_ID == 'C9')
simcoe_c9_2004$Species <- paste("C9", simcoe_c9_2004$Species, sep=" ")
simcoe_k42_2004 <- subset(simcoe_2004_m, Station_ID == 'K42')
simcoe_k42_2004$Species <- paste("K42", simcoe_k42_2004$Species, sep=" ")
simcoe_k45_2004 <- subset(simcoe_2004_m, Station_ID == 'K45')
simcoe_k45_2004$Species <- paste("K45", simcoe_k45_2004$Species, sep=" ")

#recombine dataframes
simcoe_spec_2004 <- rbind(simcoe_c9_2004, simcoe_k42_2004, simcoe_k45_2004)
simcoe_spec_2004 <- simcoe_spec_2004 %>% select(- Station_ID)

#transpose dataframe
simcoe_2004_t2 <- dcast(simcoe_spec_2004, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2004_t2 <- log(simcoe_2004_t2 + 1)
simcoe_2004_t2 <- simcoe_2004_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_2004_cv <- simcoe_2004_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_2004 = mean(Density),
            pop_sd_2004 = sd(Density)) %>%
  mutate(uw_pop_cv_2004 = pop_sd_2004/pop_mean_2004)

simcoe_2004_cv <- simcoe_2004_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_2004),
         mean_pop_cv = mean(uw_pop_cv_2004, na.rm = T),
         mean_pop_density = sum(pop_mean_2004),
         mean_pop_variance = sum(pop_sd_2004))

simcoe_2004_cv <- simcoe_2004_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_2004/mc_sum_mean_density)%>%
  mutate(w_pop_cv_2004 = (pop_mean_2004/mc_sum_mean_density)*uw_pop_cv_2004) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_2004_wa_pop_var_2004 <- simcoe_2004_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_2004 = sum(w_pop_cv_2004, na.rm = T),
            lcv = sum(w_pop_cv_2004, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_2004_cor_list_ts <- as.dist(round(cor(simcoe_2004_t2[]),2))
simcoe_2004_cor_ts <- stack(simcoe_2004_cor_list_ts, dim.names = TRUE)
simcoe_2004_cor_ts <- simcoe_2004_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_2004_cor - simcoe_2004_cv
simcoe_c9_2004_w_cv_ts <- subset(simcoe_2004_cv, Station_ID == 'C9')
simcoe_c9_2004_w_cv_ts$Species <- paste("C9", simcoe_c9_2004_w_cv_ts$Species, sep=" ")
simcoe_c9_2004_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_2004_w_cv_ts)
simcoe_k42_2004_w_cv_ts <- subset(simcoe_2004_cv, Station_ID == 'K42')
simcoe_k42_2004_w_cv_ts$Species <- paste("K42", simcoe_k42_2004_w_cv_ts$Species, sep=" ")
simcoe_k42_2004_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_2004_w_cv_ts)
simcoe_k45_2004_w_cv_ts <- subset(simcoe_2004_cv, Station_ID == 'K45')
simcoe_k45_2004_w_cv_ts$Species <- paste("K45", simcoe_k45_2004_w_cv_ts$Species, sep=" ")
simcoe_k45_2004_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_2004_w_cv_ts)

simcoe_spec_2004_w_sp_ts <- rbind(simcoe_c9_2004_w_cv_ts, simcoe_k42_2004_w_cv_ts, simcoe_k45_2004_w_cv_ts)
simcoe_spec_2004_w_sp_ts <- simcoe_spec_2004_w_sp_ts %>% select(Species, w_pop_cv_2004, Species_ID, Station_ID)
simcoe_spec_2004_w_sp1_ts <- simcoe_spec_2004_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_2004_sp1 = w_pop_cv_2004) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_2004_sp1, Species_ID, station_1)

simcoe_2004_cor_cv_ts <- merge(simcoe_2004_cor_ts, simcoe_spec_2004_w_sp1_ts, by = "species_1")

simcoe_spec_2004_w_sp2_ts <- simcoe_spec_2004_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_2004_sp2 = w_pop_cv_2004) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_2004_sp2, Species_ID, station_2)

simcoe_2004_cor_cv_ts <- merge(simcoe_2004_cor_cv_ts, simcoe_spec_2004_w_sp2_ts, by = "species_2")

simcoe_2004_cor_cv_ts_omit <- na.omit(simcoe_2004_cor_cv_ts)

simcoe_2004_ind_W_ts <- simcoe_2004_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_2004 = w_pop_cv_2004_sp1*w_pop_cv_2004_sp2,
         ind_W_2004 = (1 - corr)*(w_pop_cv_2004_sp1*w_pop_cv_2004_sp2),
         Year = 2004)

simcoe_2004_ind_W_ts$number_species_pairs <- nrow(simcoe_2004_ind_W_ts)

simcoe_mean_W_sp_pairs_2004 <- simcoe_2004_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_2004, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_2004_W_ts <- simcoe_2004_ind_W_ts %>% group_by(Year) %>%
  summarize(W_2004 = sum(ind_W_2004, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_2004_ind_W_ts_c9 <- filter(simcoe_2004_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_2004_ind_W_ts_k42 <- filter(simcoe_2004_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_2004_ind_W_ts_k45 <- filter(simcoe_2004_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_2004_c9 <- na.omit(simcoe_2004_ind_W_ts_c9)

simcoe_species_pairs_count_2004_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_2004_c9)

simcoe_mean_W_sp_pairs_2004_c9 <- simcoe_species_pairs_count_2004_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_2004, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2004_W_ts_c9 <- simcoe_2004_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_2004 = sum(ind_W_2004, na.rm = TRUE))

#K42
simcoe_species_pairs_count_2004_k42 <- na.omit(simcoe_2004_ind_W_ts_k42)

simcoe_species_pairs_count_2004_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_2004_k42)

simcoe_mean_W_sp_pairs_2004_k42 <- simcoe_species_pairs_count_2004_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_2004, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2004_W_ts_k42 <- simcoe_2004_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_2004 = sum(ind_W_2004, na.rm = TRUE))

#K45
simcoe_species_pairs_count_2004_k45 <- na.omit(simcoe_2004_ind_W_ts_k45)

simcoe_species_pairs_count_2004_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_2004_k45)

simcoe_mean_W_sp_pairs_2004_k45 <- simcoe_species_pairs_count_2004_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_2004, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2004_W_ts_k45 <- simcoe_2004_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_2004 = sum(ind_W_2004, na.rm = TRUE))

#Step 4
simcoe_2004_ind_W_as1 <- rbind(simcoe_2004_W_ts_c9, simcoe_2004_W_ts_k42, simcoe_2004_W_ts_k45)

simcoe_2004_W_as1 <- simcoe_2004_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_2004 = sum(W_2004, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_2004 <- rbind(simcoe_mean_W_sp_pairs_2004_c9, simcoe_mean_W_sp_pairs_2004_k42, simcoe_mean_W_sp_pairs_2004_k45)

simcoe_lc_sp_pairs_2004 <- simcoe_mean_lc_sp_pairs_2004 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_2004_ind_W_as4 <- filter(simcoe_2004_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_2004_mp  <- na.omit(simcoe_2004_ind_W_as4)

simcoe_species_pairs_count_2004_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_2004_mp)

simcoe_mean_mp_sp_pairs_2004  <- simcoe_species_pairs_count_2004_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_2004, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2004_sum_W_as4 <- simcoe_2004_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_2004 = sum(ind_W_2004, na.rm = TRUE))


#Step 3
simcoe_2004_sum_W_as4 <- simcoe_2004_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_2004 = sum(as4_2004, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_2004_ind_W_as5 <- filter(simcoe_2004_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_2004_cc <- na.omit(simcoe_2004_ind_W_as5)

simcoe_species_pairs_count_2004_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_2004_cc)

simcoe_mean_cc_sp_pairs_2004  <- simcoe_species_pairs_count_2004_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_2004, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2004_sum_W_as5 <- simcoe_2004_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_2004 = sum(ind_W_2004, na.rm = TRUE))


#Step 3
simcoe_2004_sum_W_as5 <- simcoe_2004_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_2004 = sum(as5_2004, na.rm = TRUE))

simcoe_2004_W <- merge(simcoe_2004_sum_W_as4, simcoe_2004_W_ts, by = c('Year'))
simcoe_2004_W <- merge(simcoe_2004_W, simcoe_2004_W_as1, by = c('Year'))
simcoe_2004_W <- merge(simcoe_2004_W, simcoe_2004_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_2004 <- merge(simcoe_mean_W_sp_pairs_2004, simcoe_lc_sp_pairs_2004, by=c("Year"))
simcoe_sp_pair_ave_stab_2004 <- merge(simcoe_sp_pair_ave_stab_2004, simcoe_mean_mp_sp_pairs_2004, by=c("Year"))
simcoe_sp_pair_ave_stab_2004 <- merge(simcoe_sp_pair_ave_stab_2004, simcoe_mean_cc_sp_pairs_2004, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_2004_final <- merge(simcoe_2004_W, simcoe_2004_wa_pop_var_2004, by = c('Year'))
simcoe_2004_final <- simcoe_2004_final %>% group_by(Year) %>%
  mutate(gcv_2004 = lcv - W_2004) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_2004, 
            lcv = lcv,
            W = W_2004,
            lc_stab = as2_2004,
            mp_stab = as4_2004,
            cc_stab = as5_2004,
            mc_cv = sqrt(gcv_2004),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_2004/lcv + as4_2004/lcv + as5_2004/lcv,
            lc_asynchrony = as2_2004/lcv,
            mp_asynchrony = as4_2004/lcv,
            cc_asynchrony= as5_2004/lcv)

simcoe_2004_cv <- simcoe_2004_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_2004_final$mc_sum_mean_density <- simcoe_2004_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 2004 #####
simcoe_mean_corr_2004_c9 <- simcoe_species_pairs_count_2004_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2004, ind_W_2004)
simcoe_mean_corr_2004_k42 <- simcoe_species_pairs_count_2004_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2004, ind_W_2004)
simcoe_mean_corr_2004_k45 <- simcoe_species_pairs_count_2004_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2004, ind_W_2004)
simcoe_mean_corr_2004_lc <- rbind(simcoe_mean_corr_2004_c9, simcoe_mean_corr_2004_k42, simcoe_mean_corr_2004_k45)
simcoe_mean_corr_2004_mp <- simcoe_species_pairs_count_2004_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_2004, ind_W_2004)
simcoe_mean_corr_2004_cc <- simcoe_species_pairs_count_2004_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_2004, ind_W_2004)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_2004_lc[sample(1:nrow(simcoe_mean_corr_2004_lc),  nrow(simcoe_mean_corr_2004_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_2004),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_2004),
              mean_lc_stab = mean(ind_W_2004),
              sd_lc_stab = sd(ind_W_2004))
  
  lc_corr_2004 <- as.data.frame(mean_lc)
}
lc_corr_2004

mp_corr_2004 <- simcoe_mean_corr_2004_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_2004),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_2004),
                                                       mean_mp_stab = mean(ind_W_2004),
                                                       sd_mp_stab = sd(ind_W_2004))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_2004_cc[sample(1:nrow(simcoe_mean_corr_2004_cc), nrow(simcoe_mean_corr_2004_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_2004),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_2004),
              mean_cc_stab = mean(ind_W_2004),
              sd_cc_stab = sd(ind_W_2004))
  
  cc_corr_2004 <- as.data.frame(mean_cc)
}
cc_corr_2004


corr_2004 <- cbind(lc_corr_2004, mp_corr_2004)
corr_2004 <- cbind(corr_2004, cc_corr_2004)
corr_2004 <- corr_2004 %>% mutate(Year = 2004)





##### 2005 #####
simcoe_2005 <- subset(simcoe, Year == '2005')


#Transpose dataframe to add 0s within years
simcoe_2005_t1 <- dcast(simcoe_2005, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2005_m <- melt(simcoe_2005_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_2005 <- subset(simcoe_2005_m, Station_ID == 'C9')
simcoe_c9_2005$Species <- paste("C9", simcoe_c9_2005$Species, sep=" ")
simcoe_k42_2005 <- subset(simcoe_2005_m, Station_ID == 'K42')
simcoe_k42_2005$Species <- paste("K42", simcoe_k42_2005$Species, sep=" ")
simcoe_k45_2005 <- subset(simcoe_2005_m, Station_ID == 'K45')
simcoe_k45_2005$Species <- paste("K45", simcoe_k45_2005$Species, sep=" ")

#recombine dataframes
simcoe_spec_2005 <- rbind(simcoe_c9_2005, simcoe_k42_2005, simcoe_k45_2005)
simcoe_spec_2005 <- simcoe_spec_2005 %>% select(- Station_ID)

#transpose dataframe
simcoe_2005_t2 <- dcast(simcoe_spec_2005, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2005_t2 <- log(simcoe_2005_t2 + 1)
simcoe_2005_t2 <- simcoe_2005_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_2005_cv <- simcoe_2005_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_2005 = mean(Density),
            pop_sd_2005 = sd(Density)) %>%
  mutate(uw_pop_cv_2005 = pop_sd_2005/pop_mean_2005)

simcoe_2005_cv <- simcoe_2005_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_2005),
         mean_pop_cv = mean(uw_pop_cv_2005, na.rm = T),
         mean_pop_density = sum(pop_mean_2005),
         mean_pop_variance = sum(pop_sd_2005))

simcoe_2005_cv <- simcoe_2005_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_2005/mc_sum_mean_density)%>%
  mutate(w_pop_cv_2005 = (pop_mean_2005/mc_sum_mean_density)*uw_pop_cv_2005) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_2005_wa_pop_var_2005 <- simcoe_2005_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_2005 = sum(w_pop_cv_2005, na.rm = T),
            lcv = sum(w_pop_cv_2005, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_2005_cor_list_ts <- as.dist(round(cor(simcoe_2005_t2[]),2))
simcoe_2005_cor_ts <- stack(simcoe_2005_cor_list_ts, dim.names = TRUE)
simcoe_2005_cor_ts <- simcoe_2005_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_2005_cor - simcoe_2005_cv
simcoe_c9_2005_w_cv_ts <- subset(simcoe_2005_cv, Station_ID == 'C9')
simcoe_c9_2005_w_cv_ts$Species <- paste("C9", simcoe_c9_2005_w_cv_ts$Species, sep=" ")
simcoe_c9_2005_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_2005_w_cv_ts)
simcoe_k42_2005_w_cv_ts <- subset(simcoe_2005_cv, Station_ID == 'K42')
simcoe_k42_2005_w_cv_ts$Species <- paste("K42", simcoe_k42_2005_w_cv_ts$Species, sep=" ")
simcoe_k42_2005_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_2005_w_cv_ts)
simcoe_k45_2005_w_cv_ts <- subset(simcoe_2005_cv, Station_ID == 'K45')
simcoe_k45_2005_w_cv_ts$Species <- paste("K45", simcoe_k45_2005_w_cv_ts$Species, sep=" ")
simcoe_k45_2005_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_2005_w_cv_ts)

simcoe_spec_2005_w_sp_ts <- rbind(simcoe_c9_2005_w_cv_ts, simcoe_k42_2005_w_cv_ts, simcoe_k45_2005_w_cv_ts)
simcoe_spec_2005_w_sp_ts <- simcoe_spec_2005_w_sp_ts %>% select(Species, w_pop_cv_2005, Species_ID, Station_ID)
simcoe_spec_2005_w_sp1_ts <- simcoe_spec_2005_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_2005_sp1 = w_pop_cv_2005) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_2005_sp1, Species_ID, station_1)

simcoe_2005_cor_cv_ts <- merge(simcoe_2005_cor_ts, simcoe_spec_2005_w_sp1_ts, by = "species_1")

simcoe_spec_2005_w_sp2_ts <- simcoe_spec_2005_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_2005_sp2 = w_pop_cv_2005) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_2005_sp2, Species_ID, station_2)

simcoe_2005_cor_cv_ts <- merge(simcoe_2005_cor_cv_ts, simcoe_spec_2005_w_sp2_ts, by = "species_2")

simcoe_2005_cor_cv_ts_omit <- na.omit(simcoe_2005_cor_cv_ts)

simcoe_2005_ind_W_ts <- simcoe_2005_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_2005 = w_pop_cv_2005_sp1*w_pop_cv_2005_sp2,
         ind_W_2005 = (1 - corr)*(w_pop_cv_2005_sp1*w_pop_cv_2005_sp2),
         Year = 2005)

simcoe_2005_ind_W_ts$number_species_pairs <- nrow(simcoe_2005_ind_W_ts)

simcoe_mean_W_sp_pairs_2005 <- simcoe_2005_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_2005, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_2005_W_ts <- simcoe_2005_ind_W_ts %>% group_by(Year) %>%
  summarize(W_2005 = sum(ind_W_2005, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_2005_ind_W_ts_c9 <- filter(simcoe_2005_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_2005_ind_W_ts_k42 <- filter(simcoe_2005_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_2005_ind_W_ts_k45 <- filter(simcoe_2005_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_2005_c9 <- na.omit(simcoe_2005_ind_W_ts_c9)

simcoe_species_pairs_count_2005_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_2005_c9)

simcoe_mean_W_sp_pairs_2005_c9 <- simcoe_species_pairs_count_2005_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_2005, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2005_W_ts_c9 <- simcoe_2005_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_2005 = sum(ind_W_2005, na.rm = TRUE))

#K42
simcoe_species_pairs_count_2005_k42 <- na.omit(simcoe_2005_ind_W_ts_k42)

simcoe_species_pairs_count_2005_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_2005_k42)

simcoe_mean_W_sp_pairs_2005_k42 <- simcoe_species_pairs_count_2005_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_2005, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2005_W_ts_k42 <- simcoe_2005_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_2005 = sum(ind_W_2005, na.rm = TRUE))

#K45
simcoe_species_pairs_count_2005_k45 <- na.omit(simcoe_2005_ind_W_ts_k45)

simcoe_species_pairs_count_2005_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_2005_k45)

simcoe_mean_W_sp_pairs_2005_k45 <- simcoe_species_pairs_count_2005_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_2005, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2005_W_ts_k45 <- simcoe_2005_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_2005 = sum(ind_W_2005, na.rm = TRUE))

#Step 4
simcoe_2005_ind_W_as1 <- rbind(simcoe_2005_W_ts_c9, simcoe_2005_W_ts_k42, simcoe_2005_W_ts_k45)

simcoe_2005_W_as1 <- simcoe_2005_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_2005 = sum(W_2005, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_2005 <- rbind(simcoe_mean_W_sp_pairs_2005_c9, simcoe_mean_W_sp_pairs_2005_k42, simcoe_mean_W_sp_pairs_2005_k45)

simcoe_lc_sp_pairs_2005 <- simcoe_mean_lc_sp_pairs_2005 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_2005_ind_W_as4 <- filter(simcoe_2005_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_2005_mp  <- na.omit(simcoe_2005_ind_W_as4)

simcoe_species_pairs_count_2005_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_2005_mp)

simcoe_mean_mp_sp_pairs_2005  <- simcoe_species_pairs_count_2005_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_2005, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2005_sum_W_as4 <- simcoe_2005_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_2005 = sum(ind_W_2005, na.rm = TRUE))


#Step 3
simcoe_2005_sum_W_as4 <- simcoe_2005_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_2005 = sum(as4_2005, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_2005_ind_W_as5 <- filter(simcoe_2005_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_2005_cc <- na.omit(simcoe_2005_ind_W_as5)

simcoe_species_pairs_count_2005_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_2005_cc)

simcoe_mean_cc_sp_pairs_2005  <- simcoe_species_pairs_count_2005_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_2005, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2005_sum_W_as5 <- simcoe_2005_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_2005 = sum(ind_W_2005, na.rm = TRUE))


#Step 3
simcoe_2005_sum_W_as5 <- simcoe_2005_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_2005 = sum(as5_2005, na.rm = TRUE))

simcoe_2005_W <- merge(simcoe_2005_sum_W_as4, simcoe_2005_W_ts, by = c('Year'))
simcoe_2005_W <- merge(simcoe_2005_W, simcoe_2005_W_as1, by = c('Year'))
simcoe_2005_W <- merge(simcoe_2005_W, simcoe_2005_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_2005 <- merge(simcoe_mean_W_sp_pairs_2005, simcoe_lc_sp_pairs_2005, by=c("Year"))
simcoe_sp_pair_ave_stab_2005 <- merge(simcoe_sp_pair_ave_stab_2005, simcoe_mean_mp_sp_pairs_2005, by=c("Year"))
simcoe_sp_pair_ave_stab_2005 <- merge(simcoe_sp_pair_ave_stab_2005, simcoe_mean_cc_sp_pairs_2005, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_2005_final <- merge(simcoe_2005_W, simcoe_2005_wa_pop_var_2005, by = c('Year'))
simcoe_2005_final <- simcoe_2005_final %>% group_by(Year) %>%
  mutate(gcv_2005 = lcv - W_2005) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_2005, 
            lcv = lcv,
            W = W_2005,
            lc_stab = as2_2005,
            mp_stab = as4_2005,
            cc_stab = as5_2005,
            mc_cv = sqrt(gcv_2005),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_2005/lcv + as4_2005/lcv + as5_2005/lcv,
            lc_asynchrony = as2_2005/lcv,
            mp_asynchrony = as4_2005/lcv,
            cc_asynchrony= as5_2005/lcv)

simcoe_2005_cv <- simcoe_2005_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_2005_final$mc_sum_mean_density <- simcoe_2005_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 2005 #####
simcoe_mean_corr_2005_c9 <- simcoe_species_pairs_count_2005_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2005, ind_W_2005)
simcoe_mean_corr_2005_k42 <- simcoe_species_pairs_count_2005_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2005, ind_W_2005)
simcoe_mean_corr_2005_k45 <- simcoe_species_pairs_count_2005_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2005, ind_W_2005)
simcoe_mean_corr_2005_lc <- rbind(simcoe_mean_corr_2005_c9, simcoe_mean_corr_2005_k42, simcoe_mean_corr_2005_k45)
simcoe_mean_corr_2005_mp <- simcoe_species_pairs_count_2005_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_2005, ind_W_2005)
simcoe_mean_corr_2005_cc <- simcoe_species_pairs_count_2005_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_2005, ind_W_2005)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_2005_lc[sample(1:nrow(simcoe_mean_corr_2005_lc),  nrow(simcoe_mean_corr_2005_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_2005),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_2005),
              mean_lc_stab = mean(ind_W_2005),
              sd_lc_stab = sd(ind_W_2005))
  
  lc_corr_2005 <- as.data.frame(mean_lc)
}
lc_corr_2005

mp_corr_2005 <- simcoe_mean_corr_2005_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_2005),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_2005),
                                                       mean_mp_stab = mean(ind_W_2005),
                                                       sd_mp_stab = sd(ind_W_2005))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_2005_cc[sample(1:nrow(simcoe_mean_corr_2005_cc), nrow(simcoe_mean_corr_2005_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_2005),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_2005),
              mean_cc_stab = mean(ind_W_2005),
              sd_cc_stab = sd(ind_W_2005))
  
  cc_corr_2005 <- as.data.frame(mean_cc)
}
cc_corr_2005


corr_2005 <- cbind(lc_corr_2005, mp_corr_2005)
corr_2005 <- cbind(corr_2005, cc_corr_2005)
corr_2005 <- corr_2005 %>% mutate(Year = 2005)




##### 2006 #####
simcoe_2006 <- subset(simcoe, Year == '2006')


#Transpose dataframe to add 0s within years
simcoe_2006_t1 <- dcast(simcoe_2006, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2006_m <- melt(simcoe_2006_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_2006 <- subset(simcoe_2006_m, Station_ID == 'C9')
simcoe_c9_2006$Species <- paste("C9", simcoe_c9_2006$Species, sep=" ")
simcoe_k42_2006 <- subset(simcoe_2006_m, Station_ID == 'K42')
simcoe_k42_2006$Species <- paste("K42", simcoe_k42_2006$Species, sep=" ")
simcoe_k45_2006 <- subset(simcoe_2006_m, Station_ID == 'K45')
simcoe_k45_2006$Species <- paste("K45", simcoe_k45_2006$Species, sep=" ")

#recombine dataframes
simcoe_spec_2006 <- rbind(simcoe_c9_2006, simcoe_k42_2006, simcoe_k45_2006)
simcoe_spec_2006 <- simcoe_spec_2006 %>% select(- Station_ID)

#transpose dataframe
simcoe_2006_t2 <- dcast(simcoe_spec_2006, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2006_t2 <- log(simcoe_2006_t2 + 1)
simcoe_2006_t2 <- simcoe_2006_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_2006_cv <- simcoe_2006_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_2006 = mean(Density),
            pop_sd_2006 = sd(Density)) %>%
  mutate(uw_pop_cv_2006 = pop_sd_2006/pop_mean_2006)

simcoe_2006_cv <- simcoe_2006_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_2006),
         mean_pop_cv = mean(uw_pop_cv_2006, na.rm = T),
         mean_pop_density = sum(pop_mean_2006),
         mean_pop_variance = sum(pop_sd_2006))

simcoe_2006_cv <- simcoe_2006_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_2006/mc_sum_mean_density)%>%
  mutate(w_pop_cv_2006 = (pop_mean_2006/mc_sum_mean_density)*uw_pop_cv_2006) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_2006_wa_pop_var_2006 <- simcoe_2006_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_2006 = sum(w_pop_cv_2006, na.rm = T),
            lcv = sum(w_pop_cv_2006, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_2006_cor_list_ts <- as.dist(round(cor(simcoe_2006_t2[]),2))
simcoe_2006_cor_ts <- stack(simcoe_2006_cor_list_ts, dim.names = TRUE)
simcoe_2006_cor_ts <- simcoe_2006_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_2006_cor - simcoe_2006_cv
simcoe_c9_2006_w_cv_ts <- subset(simcoe_2006_cv, Station_ID == 'C9')
simcoe_c9_2006_w_cv_ts$Species <- paste("C9", simcoe_c9_2006_w_cv_ts$Species, sep=" ")
simcoe_c9_2006_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_2006_w_cv_ts)
simcoe_k42_2006_w_cv_ts <- subset(simcoe_2006_cv, Station_ID == 'K42')
simcoe_k42_2006_w_cv_ts$Species <- paste("K42", simcoe_k42_2006_w_cv_ts$Species, sep=" ")
simcoe_k42_2006_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_2006_w_cv_ts)
simcoe_k45_2006_w_cv_ts <- subset(simcoe_2006_cv, Station_ID == 'K45')
simcoe_k45_2006_w_cv_ts$Species <- paste("K45", simcoe_k45_2006_w_cv_ts$Species, sep=" ")
simcoe_k45_2006_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_2006_w_cv_ts)

simcoe_spec_2006_w_sp_ts <- rbind(simcoe_c9_2006_w_cv_ts, simcoe_k42_2006_w_cv_ts, simcoe_k45_2006_w_cv_ts)
simcoe_spec_2006_w_sp_ts <- simcoe_spec_2006_w_sp_ts %>% select(Species, w_pop_cv_2006, Species_ID, Station_ID)
simcoe_spec_2006_w_sp1_ts <- simcoe_spec_2006_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_2006_sp1 = w_pop_cv_2006) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_2006_sp1, Species_ID, station_1)

simcoe_2006_cor_cv_ts <- merge(simcoe_2006_cor_ts, simcoe_spec_2006_w_sp1_ts, by = "species_1")

simcoe_spec_2006_w_sp2_ts <- simcoe_spec_2006_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_2006_sp2 = w_pop_cv_2006) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_2006_sp2, Species_ID, station_2)

simcoe_2006_cor_cv_ts <- merge(simcoe_2006_cor_cv_ts, simcoe_spec_2006_w_sp2_ts, by = "species_2")

simcoe_2006_cor_cv_ts_omit <- na.omit(simcoe_2006_cor_cv_ts)

simcoe_2006_ind_W_ts <- simcoe_2006_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_2006 = w_pop_cv_2006_sp1*w_pop_cv_2006_sp2,
         ind_W_2006 = (1 - corr)*(w_pop_cv_2006_sp1*w_pop_cv_2006_sp2),
         Year = 2006)

simcoe_2006_ind_W_ts$number_species_pairs <- nrow(simcoe_2006_ind_W_ts)

simcoe_mean_W_sp_pairs_2006 <- simcoe_2006_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_2006, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_2006_W_ts <- simcoe_2006_ind_W_ts %>% group_by(Year) %>%
  summarize(W_2006 = sum(ind_W_2006, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_2006_ind_W_ts_c9 <- filter(simcoe_2006_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_2006_ind_W_ts_k42 <- filter(simcoe_2006_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_2006_ind_W_ts_k45 <- filter(simcoe_2006_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_2006_c9 <- na.omit(simcoe_2006_ind_W_ts_c9)

simcoe_species_pairs_count_2006_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_2006_c9)

simcoe_mean_W_sp_pairs_2006_c9 <- simcoe_species_pairs_count_2006_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_2006, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2006_W_ts_c9 <- simcoe_2006_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_2006 = sum(ind_W_2006, na.rm = TRUE))

#K42
simcoe_species_pairs_count_2006_k42 <- na.omit(simcoe_2006_ind_W_ts_k42)

simcoe_species_pairs_count_2006_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_2006_k42)

simcoe_mean_W_sp_pairs_2006_k42 <- simcoe_species_pairs_count_2006_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_2006, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2006_W_ts_k42 <- simcoe_2006_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_2006 = sum(ind_W_2006, na.rm = TRUE))

#K45
simcoe_species_pairs_count_2006_k45 <- na.omit(simcoe_2006_ind_W_ts_k45)

simcoe_species_pairs_count_2006_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_2006_k45)

simcoe_mean_W_sp_pairs_2006_k45 <- simcoe_species_pairs_count_2006_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_2006, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2006_W_ts_k45 <- simcoe_2006_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_2006 = sum(ind_W_2006, na.rm = TRUE))

#Step 4
simcoe_2006_ind_W_as1 <- rbind(simcoe_2006_W_ts_c9, simcoe_2006_W_ts_k42, simcoe_2006_W_ts_k45)

simcoe_2006_W_as1 <- simcoe_2006_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_2006 = sum(W_2006, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_2006 <- rbind(simcoe_mean_W_sp_pairs_2006_c9, simcoe_mean_W_sp_pairs_2006_k42, simcoe_mean_W_sp_pairs_2006_k45)

simcoe_lc_sp_pairs_2006 <- simcoe_mean_lc_sp_pairs_2006 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_2006_ind_W_as4 <- filter(simcoe_2006_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_2006_mp  <- na.omit(simcoe_2006_ind_W_as4)

simcoe_species_pairs_count_2006_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_2006_mp)

simcoe_mean_mp_sp_pairs_2006  <- simcoe_species_pairs_count_2006_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_2006, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2006_sum_W_as4 <- simcoe_2006_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_2006 = sum(ind_W_2006, na.rm = TRUE))


#Step 3
simcoe_2006_sum_W_as4 <- simcoe_2006_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_2006 = sum(as4_2006, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_2006_ind_W_as5 <- filter(simcoe_2006_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_2006_cc <- na.omit(simcoe_2006_ind_W_as5)

simcoe_species_pairs_count_2006_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_2006_cc)

simcoe_mean_cc_sp_pairs_2006  <- simcoe_species_pairs_count_2006_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_2006, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2006_sum_W_as5 <- simcoe_2006_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_2006 = sum(ind_W_2006, na.rm = TRUE))


#Step 3
simcoe_2006_sum_W_as5 <- simcoe_2006_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_2006 = sum(as5_2006, na.rm = TRUE))

simcoe_2006_W <- merge(simcoe_2006_sum_W_as4, simcoe_2006_W_ts, by = c('Year'))
simcoe_2006_W <- merge(simcoe_2006_W, simcoe_2006_W_as1, by = c('Year'))
simcoe_2006_W <- merge(simcoe_2006_W, simcoe_2006_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_2006 <- merge(simcoe_mean_W_sp_pairs_2006, simcoe_lc_sp_pairs_2006, by=c("Year"))
simcoe_sp_pair_ave_stab_2006 <- merge(simcoe_sp_pair_ave_stab_2006, simcoe_mean_mp_sp_pairs_2006, by=c("Year"))
simcoe_sp_pair_ave_stab_2006 <- merge(simcoe_sp_pair_ave_stab_2006, simcoe_mean_cc_sp_pairs_2006, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_2006_final <- merge(simcoe_2006_W, simcoe_2006_wa_pop_var_2006, by = c('Year'))
simcoe_2006_final <- simcoe_2006_final %>% group_by(Year) %>%
  mutate(gcv_2006 = lcv - W_2006) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_2006, 
            lcv = lcv,
            W = W_2006,
            lc_stab = as2_2006,
            mp_stab = as4_2006,
            cc_stab = as5_2006,
            mc_cv = sqrt(gcv_2006),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_2006/lcv + as4_2006/lcv + as5_2006/lcv,
            lc_asynchrony = as2_2006/lcv,
            mp_asynchrony = as4_2006/lcv,
            cc_asynchrony= as5_2006/lcv)

simcoe_2006_cv <- simcoe_2006_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_2006_final$mc_sum_mean_density <- simcoe_2006_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 2006 #####
simcoe_mean_corr_2006_c9 <- simcoe_species_pairs_count_2006_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2006, ind_W_2006)
simcoe_mean_corr_2006_k42 <- simcoe_species_pairs_count_2006_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2006, ind_W_2006)
simcoe_mean_corr_2006_k45 <- simcoe_species_pairs_count_2006_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2006, ind_W_2006)
simcoe_mean_corr_2006_lc <- rbind(simcoe_mean_corr_2006_c9, simcoe_mean_corr_2006_k42, simcoe_mean_corr_2006_k45)
simcoe_mean_corr_2006_mp <- simcoe_species_pairs_count_2006_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_2006, ind_W_2006)
simcoe_mean_corr_2006_cc <- simcoe_species_pairs_count_2006_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_2006, ind_W_2006)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_2006_lc[sample(1:nrow(simcoe_mean_corr_2006_lc),  nrow(simcoe_mean_corr_2006_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_2006),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_2006),
              mean_lc_stab = mean(ind_W_2006),
              sd_lc_stab = sd(ind_W_2006))
  
  lc_corr_2006 <- as.data.frame(mean_lc)
}
lc_corr_2006

mp_corr_2006 <- simcoe_mean_corr_2006_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_2006),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_2006),
                                                       mean_mp_stab = mean(ind_W_2006),
                                                       sd_mp_stab = sd(ind_W_2006))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_2006_cc[sample(1:nrow(simcoe_mean_corr_2006_cc), nrow(simcoe_mean_corr_2006_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_2006),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_2006),
              mean_cc_stab = mean(ind_W_2006),
              sd_cc_stab = sd(ind_W_2006))
  
  cc_corr_2006 <- as.data.frame(mean_cc)
}
cc_corr_2006


corr_2006 <- cbind(lc_corr_2006, mp_corr_2006)
corr_2006 <- cbind(corr_2006, cc_corr_2006)
corr_2006 <- corr_2006 %>% mutate(Year = 2006)





##### 2007 #####
simcoe_2007 <- subset(simcoe, Year == '2007')


#Transpose dataframe to add 0s within years
simcoe_2007_t1 <- dcast(simcoe_2007, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2007_m <- melt(simcoe_2007_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_2007 <- subset(simcoe_2007_m, Station_ID == 'C9')
simcoe_c9_2007$Species <- paste("C9", simcoe_c9_2007$Species, sep=" ")
simcoe_k42_2007 <- subset(simcoe_2007_m, Station_ID == 'K42')
simcoe_k42_2007$Species <- paste("K42", simcoe_k42_2007$Species, sep=" ")
simcoe_k45_2007 <- subset(simcoe_2007_m, Station_ID == 'K45')
simcoe_k45_2007$Species <- paste("K45", simcoe_k45_2007$Species, sep=" ")

#recombine dataframes
simcoe_spec_2007 <- rbind(simcoe_c9_2007, simcoe_k42_2007, simcoe_k45_2007)
simcoe_spec_2007 <- simcoe_spec_2007 %>% select(- Station_ID)

#transpose dataframe
simcoe_2007_t2 <- dcast(simcoe_spec_2007, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2007_t2 <- log(simcoe_2007_t2 + 1)
simcoe_2007_t2 <- simcoe_2007_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_2007_cv <- simcoe_2007_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_2007 = mean(Density),
            pop_sd_2007 = sd(Density)) %>%
  mutate(uw_pop_cv_2007 = pop_sd_2007/pop_mean_2007)

simcoe_2007_cv <- simcoe_2007_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_2007),
         mean_pop_cv = mean(uw_pop_cv_2007, na.rm = T),
         mean_pop_density = sum(pop_mean_2007),
         mean_pop_variance = sum(pop_sd_2007))

simcoe_2007_cv <- simcoe_2007_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_2007/mc_sum_mean_density)%>%
  mutate(w_pop_cv_2007 = (pop_mean_2007/mc_sum_mean_density)*uw_pop_cv_2007) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_2007_wa_pop_var_2007 <- simcoe_2007_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_2007 = sum(w_pop_cv_2007, na.rm = T),
            lcv = sum(w_pop_cv_2007, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_2007_cor_list_ts <- as.dist(round(cor(simcoe_2007_t2[]),2))
simcoe_2007_cor_ts <- stack(simcoe_2007_cor_list_ts, dim.names = TRUE)
simcoe_2007_cor_ts <- simcoe_2007_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_2007_cor - simcoe_2007_cv
simcoe_c9_2007_w_cv_ts <- subset(simcoe_2007_cv, Station_ID == 'C9')
simcoe_c9_2007_w_cv_ts$Species <- paste("C9", simcoe_c9_2007_w_cv_ts$Species, sep=" ")
simcoe_c9_2007_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_2007_w_cv_ts)
simcoe_k42_2007_w_cv_ts <- subset(simcoe_2007_cv, Station_ID == 'K42')
simcoe_k42_2007_w_cv_ts$Species <- paste("K42", simcoe_k42_2007_w_cv_ts$Species, sep=" ")
simcoe_k42_2007_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_2007_w_cv_ts)
simcoe_k45_2007_w_cv_ts <- subset(simcoe_2007_cv, Station_ID == 'K45')
simcoe_k45_2007_w_cv_ts$Species <- paste("K45", simcoe_k45_2007_w_cv_ts$Species, sep=" ")
simcoe_k45_2007_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_2007_w_cv_ts)

simcoe_spec_2007_w_sp_ts <- rbind(simcoe_c9_2007_w_cv_ts, simcoe_k42_2007_w_cv_ts, simcoe_k45_2007_w_cv_ts)
simcoe_spec_2007_w_sp_ts <- simcoe_spec_2007_w_sp_ts %>% select(Species, w_pop_cv_2007, Species_ID, Station_ID)
simcoe_spec_2007_w_sp1_ts <- simcoe_spec_2007_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_2007_sp1 = w_pop_cv_2007) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_2007_sp1, Species_ID, station_1)

simcoe_2007_cor_cv_ts <- merge(simcoe_2007_cor_ts, simcoe_spec_2007_w_sp1_ts, by = "species_1")

simcoe_spec_2007_w_sp2_ts <- simcoe_spec_2007_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_2007_sp2 = w_pop_cv_2007) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_2007_sp2, Species_ID, station_2)

simcoe_2007_cor_cv_ts <- merge(simcoe_2007_cor_cv_ts, simcoe_spec_2007_w_sp2_ts, by = "species_2")

simcoe_2007_cor_cv_ts_omit <- na.omit(simcoe_2007_cor_cv_ts)

simcoe_2007_ind_W_ts <- simcoe_2007_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_2007 = w_pop_cv_2007_sp1*w_pop_cv_2007_sp2,
         ind_W_2007 = (1 - corr)*(w_pop_cv_2007_sp1*w_pop_cv_2007_sp2),
         Year = 2007)

simcoe_2007_ind_W_ts$number_species_pairs <- nrow(simcoe_2007_ind_W_ts)

simcoe_mean_W_sp_pairs_2007 <- simcoe_2007_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_2007, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_2007_W_ts <- simcoe_2007_ind_W_ts %>% group_by(Year) %>%
  summarize(W_2007 = sum(ind_W_2007, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_2007_ind_W_ts_c9 <- filter(simcoe_2007_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_2007_ind_W_ts_k42 <- filter(simcoe_2007_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_2007_ind_W_ts_k45 <- filter(simcoe_2007_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_2007_c9 <- na.omit(simcoe_2007_ind_W_ts_c9)

simcoe_species_pairs_count_2007_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_2007_c9)

simcoe_mean_W_sp_pairs_2007_c9 <- simcoe_species_pairs_count_2007_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_2007, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2007_W_ts_c9 <- simcoe_2007_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_2007 = sum(ind_W_2007, na.rm = TRUE))

#K42
simcoe_species_pairs_count_2007_k42 <- na.omit(simcoe_2007_ind_W_ts_k42)

simcoe_species_pairs_count_2007_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_2007_k42)

simcoe_mean_W_sp_pairs_2007_k42 <- simcoe_species_pairs_count_2007_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_2007, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2007_W_ts_k42 <- simcoe_2007_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_2007 = sum(ind_W_2007, na.rm = TRUE))

#K45
simcoe_species_pairs_count_2007_k45 <- na.omit(simcoe_2007_ind_W_ts_k45)

simcoe_species_pairs_count_2007_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_2007_k45)

simcoe_mean_W_sp_pairs_2007_k45 <- simcoe_species_pairs_count_2007_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_2007, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2007_W_ts_k45 <- simcoe_2007_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_2007 = sum(ind_W_2007, na.rm = TRUE))

#Step 4
simcoe_2007_ind_W_as1 <- rbind(simcoe_2007_W_ts_c9, simcoe_2007_W_ts_k42, simcoe_2007_W_ts_k45)

simcoe_2007_W_as1 <- simcoe_2007_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_2007 = sum(W_2007, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_2007 <- rbind(simcoe_mean_W_sp_pairs_2007_c9, simcoe_mean_W_sp_pairs_2007_k42, simcoe_mean_W_sp_pairs_2007_k45)

simcoe_lc_sp_pairs_2007 <- simcoe_mean_lc_sp_pairs_2007 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_2007_ind_W_as4 <- filter(simcoe_2007_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_2007_mp  <- na.omit(simcoe_2007_ind_W_as4)

simcoe_species_pairs_count_2007_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_2007_mp)

simcoe_mean_mp_sp_pairs_2007  <- simcoe_species_pairs_count_2007_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_2007, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2007_sum_W_as4 <- simcoe_2007_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_2007 = sum(ind_W_2007, na.rm = TRUE))


#Step 3
simcoe_2007_sum_W_as4 <- simcoe_2007_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_2007 = sum(as4_2007, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_2007_ind_W_as5 <- filter(simcoe_2007_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_2007_cc <- na.omit(simcoe_2007_ind_W_as5)

simcoe_species_pairs_count_2007_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_2007_cc)

simcoe_mean_cc_sp_pairs_2007  <- simcoe_species_pairs_count_2007_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_2007, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2007_sum_W_as5 <- simcoe_2007_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_2007 = sum(ind_W_2007, na.rm = TRUE))


#Step 3
simcoe_2007_sum_W_as5 <- simcoe_2007_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_2007 = sum(as5_2007, na.rm = TRUE))

simcoe_2007_W <- merge(simcoe_2007_sum_W_as4, simcoe_2007_W_ts, by = c('Year'))
simcoe_2007_W <- merge(simcoe_2007_W, simcoe_2007_W_as1, by = c('Year'))
simcoe_2007_W <- merge(simcoe_2007_W, simcoe_2007_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_2007 <- merge(simcoe_mean_W_sp_pairs_2007, simcoe_lc_sp_pairs_2007, by=c("Year"))
simcoe_sp_pair_ave_stab_2007 <- merge(simcoe_sp_pair_ave_stab_2007, simcoe_mean_mp_sp_pairs_2007, by=c("Year"))
simcoe_sp_pair_ave_stab_2007 <- merge(simcoe_sp_pair_ave_stab_2007, simcoe_mean_cc_sp_pairs_2007, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_2007_final <- merge(simcoe_2007_W, simcoe_2007_wa_pop_var_2007, by = c('Year'))
simcoe_2007_final <- simcoe_2007_final %>% group_by(Year) %>%
  mutate(gcv_2007 = lcv - W_2007) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_2007, 
            lcv = lcv,
            W = W_2007,
            lc_stab = as2_2007,
            mp_stab = as4_2007,
            cc_stab = as5_2007,
            mc_cv = sqrt(gcv_2007),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_2007/lcv + as4_2007/lcv + as5_2007/lcv,
            lc_asynchrony = as2_2007/lcv,
            mp_asynchrony = as4_2007/lcv,
            cc_asynchrony= as5_2007/lcv)

simcoe_2007_cv <- simcoe_2007_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_2007_final$mc_sum_mean_density <- simcoe_2007_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 2007 #####
simcoe_mean_corr_2007_c9 <- simcoe_species_pairs_count_2007_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2007, ind_W_2007)
simcoe_mean_corr_2007_k42 <- simcoe_species_pairs_count_2007_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2007, ind_W_2007)
simcoe_mean_corr_2007_k45 <- simcoe_species_pairs_count_2007_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2007, ind_W_2007)
simcoe_mean_corr_2007_lc <- rbind(simcoe_mean_corr_2007_c9, simcoe_mean_corr_2007_k42, simcoe_mean_corr_2007_k45)
simcoe_mean_corr_2007_mp <- simcoe_species_pairs_count_2007_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_2007, ind_W_2007)
simcoe_mean_corr_2007_cc <- simcoe_species_pairs_count_2007_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_2007, ind_W_2007)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_2007_lc[sample(1:nrow(simcoe_mean_corr_2007_lc),  nrow(simcoe_mean_corr_2007_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_2007),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_2007),
              mean_lc_stab = mean(ind_W_2007),
              sd_lc_stab = sd(ind_W_2007))
  
  lc_corr_2007 <- as.data.frame(mean_lc)
}
lc_corr_2007

mp_corr_2007 <- simcoe_mean_corr_2007_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_2007),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_2007),
                                                       mean_mp_stab = mean(ind_W_2007),
                                                       sd_mp_stab = sd(ind_W_2007))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_2007_cc[sample(1:nrow(simcoe_mean_corr_2007_cc), nrow(simcoe_mean_corr_2007_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_2007),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_2007),
              mean_cc_stab = mean(ind_W_2007),
              sd_cc_stab = sd(ind_W_2007))
  
  cc_corr_2007 <- as.data.frame(mean_cc)
}
cc_corr_2007


corr_2007 <- cbind(lc_corr_2007, mp_corr_2007)
corr_2007 <- cbind(corr_2007, cc_corr_2007)
corr_2007 <- corr_2007 %>% mutate(Year = 2007)






##### 2008 #####
simcoe_2008 <- subset(simcoe, Year == '2008')


#Transpose dataframe to add 0s within years
simcoe_2008_t1 <- dcast(simcoe_2008, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2008_m <- melt(simcoe_2008_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_2008 <- subset(simcoe_2008_m, Station_ID == 'C9')
simcoe_c9_2008$Species <- paste("C9", simcoe_c9_2008$Species, sep=" ")
simcoe_k42_2008 <- subset(simcoe_2008_m, Station_ID == 'K42')
simcoe_k42_2008$Species <- paste("K42", simcoe_k42_2008$Species, sep=" ")
simcoe_k45_2008 <- subset(simcoe_2008_m, Station_ID == 'K45')
simcoe_k45_2008$Species <- paste("K45", simcoe_k45_2008$Species, sep=" ")

#recombine dataframes
simcoe_spec_2008 <- rbind(simcoe_c9_2008, simcoe_k42_2008, simcoe_k45_2008)
simcoe_spec_2008 <- simcoe_spec_2008 %>% select(- Station_ID)

#transpose dataframe
simcoe_2008_t2 <- dcast(simcoe_spec_2008, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2008_t2 <- log(simcoe_2008_t2 + 1)
simcoe_2008_t2 <- simcoe_2008_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_2008_cv <- simcoe_2008_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_2008 = mean(Density),
            pop_sd_2008 = sd(Density)) %>%
  mutate(uw_pop_cv_2008 = pop_sd_2008/pop_mean_2008)

simcoe_2008_cv <- simcoe_2008_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_2008),
         mean_pop_cv = mean(uw_pop_cv_2008, na.rm = T),
         mean_pop_density = sum(pop_mean_2008),
         mean_pop_variance = sum(pop_sd_2008))

simcoe_2008_cv <- simcoe_2008_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_2008/mc_sum_mean_density)%>%
  mutate(w_pop_cv_2008 = (pop_mean_2008/mc_sum_mean_density)*uw_pop_cv_2008) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_2008_wa_pop_var_2008 <- simcoe_2008_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_2008 = sum(w_pop_cv_2008, na.rm = T),
            lcv = sum(w_pop_cv_2008, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_2008_cor_list_ts <- as.dist(round(cor(simcoe_2008_t2[]),2))
simcoe_2008_cor_ts <- stack(simcoe_2008_cor_list_ts, dim.names = TRUE)
simcoe_2008_cor_ts <- simcoe_2008_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_2008_cor - simcoe_2008_cv
simcoe_c9_2008_w_cv_ts <- subset(simcoe_2008_cv, Station_ID == 'C9')
simcoe_c9_2008_w_cv_ts$Species <- paste("C9", simcoe_c9_2008_w_cv_ts$Species, sep=" ")
simcoe_c9_2008_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_2008_w_cv_ts)
simcoe_k42_2008_w_cv_ts <- subset(simcoe_2008_cv, Station_ID == 'K42')
simcoe_k42_2008_w_cv_ts$Species <- paste("K42", simcoe_k42_2008_w_cv_ts$Species, sep=" ")
simcoe_k42_2008_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_2008_w_cv_ts)
simcoe_k45_2008_w_cv_ts <- subset(simcoe_2008_cv, Station_ID == 'K45')
simcoe_k45_2008_w_cv_ts$Species <- paste("K45", simcoe_k45_2008_w_cv_ts$Species, sep=" ")
simcoe_k45_2008_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_2008_w_cv_ts)

simcoe_spec_2008_w_sp_ts <- rbind(simcoe_c9_2008_w_cv_ts, simcoe_k42_2008_w_cv_ts, simcoe_k45_2008_w_cv_ts)
simcoe_spec_2008_w_sp_ts <- simcoe_spec_2008_w_sp_ts %>% select(Species, w_pop_cv_2008, Species_ID, Station_ID)
simcoe_spec_2008_w_sp1_ts <- simcoe_spec_2008_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_2008_sp1 = w_pop_cv_2008) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_2008_sp1, Species_ID, station_1)

simcoe_2008_cor_cv_ts <- merge(simcoe_2008_cor_ts, simcoe_spec_2008_w_sp1_ts, by = "species_1")

simcoe_spec_2008_w_sp2_ts <- simcoe_spec_2008_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_2008_sp2 = w_pop_cv_2008) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_2008_sp2, Species_ID, station_2)

simcoe_2008_cor_cv_ts <- merge(simcoe_2008_cor_cv_ts, simcoe_spec_2008_w_sp2_ts, by = "species_2")

simcoe_2008_cor_cv_ts_omit <- na.omit(simcoe_2008_cor_cv_ts)

simcoe_2008_ind_W_ts <- simcoe_2008_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_2008 = w_pop_cv_2008_sp1*w_pop_cv_2008_sp2,
         ind_W_2008 = (1 - corr)*(w_pop_cv_2008_sp1*w_pop_cv_2008_sp2),
         Year = 2008)

simcoe_2008_ind_W_ts$number_species_pairs <- nrow(simcoe_2008_ind_W_ts)

simcoe_mean_W_sp_pairs_2008 <- simcoe_2008_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_2008, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_2008_W_ts <- simcoe_2008_ind_W_ts %>% group_by(Year) %>%
  summarize(W_2008 = sum(ind_W_2008, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_2008_ind_W_ts_c9 <- filter(simcoe_2008_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_2008_ind_W_ts_k42 <- filter(simcoe_2008_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_2008_ind_W_ts_k45 <- filter(simcoe_2008_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_2008_c9 <- na.omit(simcoe_2008_ind_W_ts_c9)

simcoe_species_pairs_count_2008_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_2008_c9)

simcoe_mean_W_sp_pairs_2008_c9 <- simcoe_species_pairs_count_2008_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_2008, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2008_W_ts_c9 <- simcoe_2008_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_2008 = sum(ind_W_2008, na.rm = TRUE))

#K42
simcoe_species_pairs_count_2008_k42 <- na.omit(simcoe_2008_ind_W_ts_k42)

simcoe_species_pairs_count_2008_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_2008_k42)

simcoe_mean_W_sp_pairs_2008_k42 <- simcoe_species_pairs_count_2008_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_2008, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2008_W_ts_k42 <- simcoe_2008_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_2008 = sum(ind_W_2008, na.rm = TRUE))

#K45
simcoe_species_pairs_count_2008_k45 <- na.omit(simcoe_2008_ind_W_ts_k45)

simcoe_species_pairs_count_2008_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_2008_k45)

simcoe_mean_W_sp_pairs_2008_k45 <- simcoe_species_pairs_count_2008_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_2008, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2008_W_ts_k45 <- simcoe_2008_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_2008 = sum(ind_W_2008, na.rm = TRUE))

#Step 4
simcoe_2008_ind_W_as1 <- rbind(simcoe_2008_W_ts_c9, simcoe_2008_W_ts_k42, simcoe_2008_W_ts_k45)

simcoe_2008_W_as1 <- simcoe_2008_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_2008 = sum(W_2008, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_2008 <- rbind(simcoe_mean_W_sp_pairs_2008_c9, simcoe_mean_W_sp_pairs_2008_k42, simcoe_mean_W_sp_pairs_2008_k45)

simcoe_lc_sp_pairs_2008 <- simcoe_mean_lc_sp_pairs_2008 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_2008_ind_W_as4 <- filter(simcoe_2008_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_2008_mp  <- na.omit(simcoe_2008_ind_W_as4)

simcoe_species_pairs_count_2008_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_2008_mp)

simcoe_mean_mp_sp_pairs_2008  <- simcoe_species_pairs_count_2008_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_2008, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2008_sum_W_as4 <- simcoe_2008_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_2008 = sum(ind_W_2008, na.rm = TRUE))


#Step 3
simcoe_2008_sum_W_as4 <- simcoe_2008_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_2008 = sum(as4_2008, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_2008_ind_W_as5 <- filter(simcoe_2008_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_2008_cc <- na.omit(simcoe_2008_ind_W_as5)

simcoe_species_pairs_count_2008_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_2008_cc)

simcoe_mean_cc_sp_pairs_2008  <- simcoe_species_pairs_count_2008_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_2008, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2008_sum_W_as5 <- simcoe_2008_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_2008 = sum(ind_W_2008, na.rm = TRUE))


#Step 3
simcoe_2008_sum_W_as5 <- simcoe_2008_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_2008 = sum(as5_2008, na.rm = TRUE))

simcoe_2008_W <- merge(simcoe_2008_sum_W_as4, simcoe_2008_W_ts, by = c('Year'))
simcoe_2008_W <- merge(simcoe_2008_W, simcoe_2008_W_as1, by = c('Year'))
simcoe_2008_W <- merge(simcoe_2008_W, simcoe_2008_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_2008 <- merge(simcoe_mean_W_sp_pairs_2008, simcoe_lc_sp_pairs_2008, by=c("Year"))
simcoe_sp_pair_ave_stab_2008 <- merge(simcoe_sp_pair_ave_stab_2008, simcoe_mean_mp_sp_pairs_2008, by=c("Year"))
simcoe_sp_pair_ave_stab_2008 <- merge(simcoe_sp_pair_ave_stab_2008, simcoe_mean_cc_sp_pairs_2008, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_2008_final <- merge(simcoe_2008_W, simcoe_2008_wa_pop_var_2008, by = c('Year'))
simcoe_2008_final <- simcoe_2008_final %>% group_by(Year) %>%
  mutate(gcv_2008 = lcv - W_2008) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_2008, 
            lcv = lcv,
            W = W_2008,
            lc_stab = as2_2008,
            mp_stab = as4_2008,
            cc_stab = as5_2008,
            mc_cv = sqrt(gcv_2008),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_2008/lcv + as4_2008/lcv + as5_2008/lcv,
            lc_asynchrony = as2_2008/lcv,
            mp_asynchrony = as4_2008/lcv,
            cc_asynchrony= as5_2008/lcv)

simcoe_2008_cv <- simcoe_2008_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_2008_final$mc_sum_mean_density <- simcoe_2008_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 2008 #####
simcoe_mean_corr_2008_c9 <- simcoe_species_pairs_count_2008_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2008, ind_W_2008)
simcoe_mean_corr_2008_k42 <- simcoe_species_pairs_count_2008_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2008, ind_W_2008)
simcoe_mean_corr_2008_k45 <- simcoe_species_pairs_count_2008_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2008, ind_W_2008)
simcoe_mean_corr_2008_lc <- rbind(simcoe_mean_corr_2008_c9, simcoe_mean_corr_2008_k42, simcoe_mean_corr_2008_k45)
simcoe_mean_corr_2008_mp <- simcoe_species_pairs_count_2008_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_2008, ind_W_2008)
simcoe_mean_corr_2008_cc <- simcoe_species_pairs_count_2008_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_2008, ind_W_2008)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_2008_lc[sample(1:nrow(simcoe_mean_corr_2008_lc),  nrow(simcoe_mean_corr_2008_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_2008),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_2008),
              mean_lc_stab = mean(ind_W_2008),
              sd_lc_stab = sd(ind_W_2008))
  
  lc_corr_2008 <- as.data.frame(mean_lc)
}
lc_corr_2008

mp_corr_2008 <- simcoe_mean_corr_2008_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_2008),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_2008),
                                                       mean_mp_stab = mean(ind_W_2008),
                                                       sd_mp_stab = sd(ind_W_2008))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_2008_cc[sample(1:nrow(simcoe_mean_corr_2008_cc), nrow(simcoe_mean_corr_2008_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_2008),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_2008),
              mean_cc_stab = mean(ind_W_2008),
              sd_cc_stab = sd(ind_W_2008))
  
  cc_corr_2008 <- as.data.frame(mean_cc)
}
cc_corr_2008


corr_2008 <- cbind(lc_corr_2008, mp_corr_2008)
corr_2008 <- cbind(corr_2008, cc_corr_2008)
corr_2008 <- corr_2008 %>% mutate(Year = 2008)




##### 2009 #####
simcoe_2009 <- subset(simcoe, Year == '2009')


#Transpose dataframe to add 0s within years
simcoe_2009_t1 <- dcast(simcoe_2009, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2009_m <- melt(simcoe_2009_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_2009 <- subset(simcoe_2009_m, Station_ID == 'C9')
simcoe_c9_2009$Species <- paste("C9", simcoe_c9_2009$Species, sep=" ")
simcoe_k42_2009 <- subset(simcoe_2009_m, Station_ID == 'K42')
simcoe_k42_2009$Species <- paste("K42", simcoe_k42_2009$Species, sep=" ")
simcoe_k45_2009 <- subset(simcoe_2009_m, Station_ID == 'K45')
simcoe_k45_2009$Species <- paste("K45", simcoe_k45_2009$Species, sep=" ")

#recombine dataframes
simcoe_spec_2009 <- rbind(simcoe_c9_2009, simcoe_k42_2009, simcoe_k45_2009)
simcoe_spec_2009 <- simcoe_spec_2009 %>% select(- Station_ID)

#transpose dataframe
simcoe_2009_t2 <- dcast(simcoe_spec_2009, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2009_t2 <- log(simcoe_2009_t2 + 1)
simcoe_2009_t2 <- simcoe_2009_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_2009_cv <- simcoe_2009_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_2009 = mean(Density),
            pop_sd_2009 = sd(Density)) %>%
  mutate(uw_pop_cv_2009 = pop_sd_2009/pop_mean_2009)

simcoe_2009_cv <- simcoe_2009_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_2009),
         mean_pop_cv = mean(uw_pop_cv_2009, na.rm = T),
         mean_pop_density = sum(pop_mean_2009),
         mean_pop_variance = sum(pop_sd_2009))

simcoe_2009_cv <- simcoe_2009_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_2009/mc_sum_mean_density)%>%
  mutate(w_pop_cv_2009 = (pop_mean_2009/mc_sum_mean_density)*uw_pop_cv_2009) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_2009_wa_pop_var_2009 <- simcoe_2009_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_2009 = sum(w_pop_cv_2009, na.rm = T),
            lcv = sum(w_pop_cv_2009, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_2009_cor_list_ts <- as.dist(round(cor(simcoe_2009_t2[]),2))
simcoe_2009_cor_ts <- stack(simcoe_2009_cor_list_ts, dim.names = TRUE)
simcoe_2009_cor_ts <- simcoe_2009_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_2009_cor - simcoe_2009_cv
simcoe_c9_2009_w_cv_ts <- subset(simcoe_2009_cv, Station_ID == 'C9')
simcoe_c9_2009_w_cv_ts$Species <- paste("C9", simcoe_c9_2009_w_cv_ts$Species, sep=" ")
simcoe_c9_2009_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_2009_w_cv_ts)
simcoe_k42_2009_w_cv_ts <- subset(simcoe_2009_cv, Station_ID == 'K42')
simcoe_k42_2009_w_cv_ts$Species <- paste("K42", simcoe_k42_2009_w_cv_ts$Species, sep=" ")
simcoe_k42_2009_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_2009_w_cv_ts)
simcoe_k45_2009_w_cv_ts <- subset(simcoe_2009_cv, Station_ID == 'K45')
simcoe_k45_2009_w_cv_ts$Species <- paste("K45", simcoe_k45_2009_w_cv_ts$Species, sep=" ")
simcoe_k45_2009_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_2009_w_cv_ts)

simcoe_spec_2009_w_sp_ts <- rbind(simcoe_c9_2009_w_cv_ts, simcoe_k42_2009_w_cv_ts, simcoe_k45_2009_w_cv_ts)
simcoe_spec_2009_w_sp_ts <- simcoe_spec_2009_w_sp_ts %>% select(Species, w_pop_cv_2009, Species_ID, Station_ID)
simcoe_spec_2009_w_sp1_ts <- simcoe_spec_2009_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_2009_sp1 = w_pop_cv_2009) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_2009_sp1, Species_ID, station_1)

simcoe_2009_cor_cv_ts <- merge(simcoe_2009_cor_ts, simcoe_spec_2009_w_sp1_ts, by = "species_1")

simcoe_spec_2009_w_sp2_ts <- simcoe_spec_2009_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_2009_sp2 = w_pop_cv_2009) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_2009_sp2, Species_ID, station_2)

simcoe_2009_cor_cv_ts <- merge(simcoe_2009_cor_cv_ts, simcoe_spec_2009_w_sp2_ts, by = "species_2")

simcoe_2009_cor_cv_ts_omit <- na.omit(simcoe_2009_cor_cv_ts)

simcoe_2009_ind_W_ts <- simcoe_2009_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_2009 = w_pop_cv_2009_sp1*w_pop_cv_2009_sp2,
         ind_W_2009 = (1 - corr)*(w_pop_cv_2009_sp1*w_pop_cv_2009_sp2),
         Year = 2009)

simcoe_2009_ind_W_ts$number_species_pairs <- nrow(simcoe_2009_ind_W_ts)

simcoe_mean_W_sp_pairs_2009 <- simcoe_2009_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_2009, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_2009_W_ts <- simcoe_2009_ind_W_ts %>% group_by(Year) %>%
  summarize(W_2009 = sum(ind_W_2009, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_2009_ind_W_ts_c9 <- filter(simcoe_2009_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_2009_ind_W_ts_k42 <- filter(simcoe_2009_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_2009_ind_W_ts_k45 <- filter(simcoe_2009_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_2009_c9 <- na.omit(simcoe_2009_ind_W_ts_c9)

simcoe_species_pairs_count_2009_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_2009_c9)

simcoe_mean_W_sp_pairs_2009_c9 <- simcoe_species_pairs_count_2009_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_2009, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2009_W_ts_c9 <- simcoe_2009_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_2009 = sum(ind_W_2009, na.rm = TRUE))

#K42
simcoe_species_pairs_count_2009_k42 <- na.omit(simcoe_2009_ind_W_ts_k42)

simcoe_species_pairs_count_2009_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_2009_k42)

simcoe_mean_W_sp_pairs_2009_k42 <- simcoe_species_pairs_count_2009_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_2009, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2009_W_ts_k42 <- simcoe_2009_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_2009 = sum(ind_W_2009, na.rm = TRUE))

#K45
simcoe_species_pairs_count_2009_k45 <- na.omit(simcoe_2009_ind_W_ts_k45)

simcoe_species_pairs_count_2009_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_2009_k45)

simcoe_mean_W_sp_pairs_2009_k45 <- simcoe_species_pairs_count_2009_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_2009, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2009_W_ts_k45 <- simcoe_2009_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_2009 = sum(ind_W_2009, na.rm = TRUE))

#Step 4
simcoe_2009_ind_W_as1 <- rbind(simcoe_2009_W_ts_c9, simcoe_2009_W_ts_k42, simcoe_2009_W_ts_k45)

simcoe_2009_W_as1 <- simcoe_2009_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_2009 = sum(W_2009, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_2009 <- rbind(simcoe_mean_W_sp_pairs_2009_c9, simcoe_mean_W_sp_pairs_2009_k42, simcoe_mean_W_sp_pairs_2009_k45)

simcoe_lc_sp_pairs_2009 <- simcoe_mean_lc_sp_pairs_2009 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_2009_ind_W_as4 <- filter(simcoe_2009_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_2009_mp  <- na.omit(simcoe_2009_ind_W_as4)

simcoe_species_pairs_count_2009_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_2009_mp)

simcoe_mean_mp_sp_pairs_2009  <- simcoe_species_pairs_count_2009_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_2009, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2009_sum_W_as4 <- simcoe_2009_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_2009 = sum(ind_W_2009, na.rm = TRUE))


#Step 3
simcoe_2009_sum_W_as4 <- simcoe_2009_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_2009 = sum(as4_2009, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_2009_ind_W_as5 <- filter(simcoe_2009_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_2009_cc <- na.omit(simcoe_2009_ind_W_as5)

simcoe_species_pairs_count_2009_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_2009_cc)

simcoe_mean_cc_sp_pairs_2009  <- simcoe_species_pairs_count_2009_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_2009, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2009_sum_W_as5 <- simcoe_2009_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_2009 = sum(ind_W_2009, na.rm = TRUE))


#Step 3
simcoe_2009_sum_W_as5 <- simcoe_2009_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_2009 = sum(as5_2009, na.rm = TRUE))

simcoe_2009_W <- merge(simcoe_2009_sum_W_as4, simcoe_2009_W_ts, by = c('Year'))
simcoe_2009_W <- merge(simcoe_2009_W, simcoe_2009_W_as1, by = c('Year'))
simcoe_2009_W <- merge(simcoe_2009_W, simcoe_2009_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_2009 <- merge(simcoe_mean_W_sp_pairs_2009, simcoe_lc_sp_pairs_2009, by=c("Year"))
simcoe_sp_pair_ave_stab_2009 <- merge(simcoe_sp_pair_ave_stab_2009, simcoe_mean_mp_sp_pairs_2009, by=c("Year"))
simcoe_sp_pair_ave_stab_2009 <- merge(simcoe_sp_pair_ave_stab_2009, simcoe_mean_cc_sp_pairs_2009, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_2009_final <- merge(simcoe_2009_W, simcoe_2009_wa_pop_var_2009, by = c('Year'))
simcoe_2009_final <- simcoe_2009_final %>% group_by(Year) %>%
  mutate(gcv_2009 = lcv - W_2009) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_2009, 
            lcv = lcv,
            W = W_2009,
            lc_stab = as2_2009,
            mp_stab = as4_2009,
            cc_stab = as5_2009,
            mc_cv = sqrt(gcv_2009),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_2009/lcv + as4_2009/lcv + as5_2009/lcv,
            lc_asynchrony = as2_2009/lcv,
            mp_asynchrony = as4_2009/lcv,
            cc_asynchrony= as5_2009/lcv)

simcoe_2009_cv <- simcoe_2009_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_2009_final$mc_sum_mean_density <- simcoe_2009_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 2009 #####
simcoe_mean_corr_2009_c9 <- simcoe_species_pairs_count_2009_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2009, ind_W_2009)
simcoe_mean_corr_2009_k42 <- simcoe_species_pairs_count_2009_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2009, ind_W_2009)
simcoe_mean_corr_2009_k45 <- simcoe_species_pairs_count_2009_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2009, ind_W_2009)
simcoe_mean_corr_2009_lc <- rbind(simcoe_mean_corr_2009_c9, simcoe_mean_corr_2009_k42, simcoe_mean_corr_2009_k45)
simcoe_mean_corr_2009_mp <- simcoe_species_pairs_count_2009_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_2009, ind_W_2009)
simcoe_mean_corr_2009_cc <- simcoe_species_pairs_count_2009_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_2009, ind_W_2009)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_2009_lc[sample(1:nrow(simcoe_mean_corr_2009_lc),  nrow(simcoe_mean_corr_2009_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_2009),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_2009),
              mean_lc_stab = mean(ind_W_2009),
              sd_lc_stab = sd(ind_W_2009))
  
  lc_corr_2009 <- as.data.frame(mean_lc)
}
lc_corr_2009

mp_corr_2009 <- simcoe_mean_corr_2009_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_2009),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_2009),
                                                       mean_mp_stab = mean(ind_W_2009),
                                                       sd_mp_stab = sd(ind_W_2009))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_2009_cc[sample(1:nrow(simcoe_mean_corr_2009_cc), nrow(simcoe_mean_corr_2009_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_2009),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_2009),
              mean_cc_stab = mean(ind_W_2009),
              sd_cc_stab = sd(ind_W_2009))
  
  cc_corr_2009 <- as.data.frame(mean_cc)
}
cc_corr_2009


corr_2009 <- cbind(lc_corr_2009, mp_corr_2009)
corr_2009 <- cbind(corr_2009, cc_corr_2009)
corr_2009 <- corr_2009 %>% mutate(Year = 2009)





##### 2010 #####
simcoe_2010 <- subset(simcoe, Year == '2010')


#Transpose dataframe to add 0s within years
simcoe_2010_t1 <- dcast(simcoe_2010, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2010_m <- melt(simcoe_2010_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_2010 <- subset(simcoe_2010_m, Station_ID == 'C9')
simcoe_c9_2010$Species <- paste("C9", simcoe_c9_2010$Species, sep=" ")
simcoe_k42_2010 <- subset(simcoe_2010_m, Station_ID == 'K42')
simcoe_k42_2010$Species <- paste("K42", simcoe_k42_2010$Species, sep=" ")
simcoe_k45_2010 <- subset(simcoe_2010_m, Station_ID == 'K45')
simcoe_k45_2010$Species <- paste("K45", simcoe_k45_2010$Species, sep=" ")

#recombine dataframes
simcoe_spec_2010 <- rbind(simcoe_c9_2010, simcoe_k42_2010, simcoe_k45_2010)
simcoe_spec_2010 <- simcoe_spec_2010 %>% select(- Station_ID)

#transpose dataframe
simcoe_2010_t2 <- dcast(simcoe_spec_2010, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2010_t2 <- log(simcoe_2010_t2 + 1)
simcoe_2010_t2 <- simcoe_2010_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_2010_cv <- simcoe_2010_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_2010 = mean(Density),
            pop_sd_2010 = sd(Density)) %>%
  mutate(uw_pop_cv_2010 = pop_sd_2010/pop_mean_2010)

simcoe_2010_cv <- simcoe_2010_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_2010),
         mean_pop_cv = mean(uw_pop_cv_2010, na.rm = T),
         mean_pop_density = sum(pop_mean_2010),
         mean_pop_variance = sum(pop_sd_2010))

simcoe_2010_cv <- simcoe_2010_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_2010/mc_sum_mean_density)%>%
  mutate(w_pop_cv_2010 = (pop_mean_2010/mc_sum_mean_density)*uw_pop_cv_2010) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_2010_wa_pop_var_2010 <- simcoe_2010_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_2010 = sum(w_pop_cv_2010, na.rm = T),
            lcv = sum(w_pop_cv_2010, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_2010_cor_list_ts <- as.dist(round(cor(simcoe_2010_t2[]),2))
simcoe_2010_cor_ts <- stack(simcoe_2010_cor_list_ts, dim.names = TRUE)
simcoe_2010_cor_ts <- simcoe_2010_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_2010_cor - simcoe_2010_cv
simcoe_c9_2010_w_cv_ts <- subset(simcoe_2010_cv, Station_ID == 'C9')
simcoe_c9_2010_w_cv_ts$Species <- paste("C9", simcoe_c9_2010_w_cv_ts$Species, sep=" ")
simcoe_c9_2010_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_2010_w_cv_ts)
simcoe_k42_2010_w_cv_ts <- subset(simcoe_2010_cv, Station_ID == 'K42')
simcoe_k42_2010_w_cv_ts$Species <- paste("K42", simcoe_k42_2010_w_cv_ts$Species, sep=" ")
simcoe_k42_2010_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_2010_w_cv_ts)
simcoe_k45_2010_w_cv_ts <- subset(simcoe_2010_cv, Station_ID == 'K45')
simcoe_k45_2010_w_cv_ts$Species <- paste("K45", simcoe_k45_2010_w_cv_ts$Species, sep=" ")
simcoe_k45_2010_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_2010_w_cv_ts)

simcoe_spec_2010_w_sp_ts <- rbind(simcoe_c9_2010_w_cv_ts, simcoe_k42_2010_w_cv_ts, simcoe_k45_2010_w_cv_ts)
simcoe_spec_2010_w_sp_ts <- simcoe_spec_2010_w_sp_ts %>% select(Species, w_pop_cv_2010, Species_ID, Station_ID)
simcoe_spec_2010_w_sp1_ts <- simcoe_spec_2010_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_2010_sp1 = w_pop_cv_2010) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_2010_sp1, Species_ID, station_1)

simcoe_2010_cor_cv_ts <- merge(simcoe_2010_cor_ts, simcoe_spec_2010_w_sp1_ts, by = "species_1")

simcoe_spec_2010_w_sp2_ts <- simcoe_spec_2010_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_2010_sp2 = w_pop_cv_2010) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_2010_sp2, Species_ID, station_2)

simcoe_2010_cor_cv_ts <- merge(simcoe_2010_cor_cv_ts, simcoe_spec_2010_w_sp2_ts, by = "species_2")

simcoe_2010_cor_cv_ts_omit <- na.omit(simcoe_2010_cor_cv_ts)

simcoe_2010_ind_W_ts <- simcoe_2010_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_2010 = w_pop_cv_2010_sp1*w_pop_cv_2010_sp2,
         ind_W_2010 = (1 - corr)*(w_pop_cv_2010_sp1*w_pop_cv_2010_sp2),
         Year = 2010)

simcoe_2010_ind_W_ts$number_species_pairs <- nrow(simcoe_2010_ind_W_ts)

simcoe_mean_W_sp_pairs_2010 <- simcoe_2010_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_2010, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_2010_W_ts <- simcoe_2010_ind_W_ts %>% group_by(Year) %>%
  summarize(W_2010 = sum(ind_W_2010, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_2010_ind_W_ts_c9 <- filter(simcoe_2010_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_2010_ind_W_ts_k42 <- filter(simcoe_2010_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_2010_ind_W_ts_k45 <- filter(simcoe_2010_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_2010_c9 <- na.omit(simcoe_2010_ind_W_ts_c9)

simcoe_species_pairs_count_2010_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_2010_c9)

simcoe_mean_W_sp_pairs_2010_c9 <- simcoe_species_pairs_count_2010_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_2010, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2010_W_ts_c9 <- simcoe_2010_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_2010 = sum(ind_W_2010, na.rm = TRUE))

#K42
simcoe_species_pairs_count_2010_k42 <- na.omit(simcoe_2010_ind_W_ts_k42)

simcoe_species_pairs_count_2010_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_2010_k42)

simcoe_mean_W_sp_pairs_2010_k42 <- simcoe_species_pairs_count_2010_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_2010, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2010_W_ts_k42 <- simcoe_2010_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_2010 = sum(ind_W_2010, na.rm = TRUE))

#K45
simcoe_species_pairs_count_2010_k45 <- na.omit(simcoe_2010_ind_W_ts_k45)

simcoe_species_pairs_count_2010_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_2010_k45)

simcoe_mean_W_sp_pairs_2010_k45 <- simcoe_species_pairs_count_2010_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_2010, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2010_W_ts_k45 <- simcoe_2010_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_2010 = sum(ind_W_2010, na.rm = TRUE))

#Step 4
simcoe_2010_ind_W_as1 <- rbind(simcoe_2010_W_ts_c9, simcoe_2010_W_ts_k42, simcoe_2010_W_ts_k45)

simcoe_2010_W_as1 <- simcoe_2010_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_2010 = sum(W_2010, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_2010 <- rbind(simcoe_mean_W_sp_pairs_2010_c9, simcoe_mean_W_sp_pairs_2010_k42, simcoe_mean_W_sp_pairs_2010_k45)

simcoe_lc_sp_pairs_2010 <- simcoe_mean_lc_sp_pairs_2010 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_2010_ind_W_as4 <- filter(simcoe_2010_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_2010_mp  <- na.omit(simcoe_2010_ind_W_as4)

simcoe_species_pairs_count_2010_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_2010_mp)

simcoe_mean_mp_sp_pairs_2010  <- simcoe_species_pairs_count_2010_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_2010, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2010_sum_W_as4 <- simcoe_2010_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_2010 = sum(ind_W_2010, na.rm = TRUE))


#Step 3
simcoe_2010_sum_W_as4 <- simcoe_2010_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_2010 = sum(as4_2010, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_2010_ind_W_as5 <- filter(simcoe_2010_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_2010_cc <- na.omit(simcoe_2010_ind_W_as5)

simcoe_species_pairs_count_2010_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_2010_cc)

simcoe_mean_cc_sp_pairs_2010  <- simcoe_species_pairs_count_2010_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_2010, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2010_sum_W_as5 <- simcoe_2010_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_2010 = sum(ind_W_2010, na.rm = TRUE))


#Step 3
simcoe_2010_sum_W_as5 <- simcoe_2010_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_2010 = sum(as5_2010, na.rm = TRUE))

simcoe_2010_W <- merge(simcoe_2010_sum_W_as4, simcoe_2010_W_ts, by = c('Year'))
simcoe_2010_W <- merge(simcoe_2010_W, simcoe_2010_W_as1, by = c('Year'))
simcoe_2010_W <- merge(simcoe_2010_W, simcoe_2010_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_2010 <- merge(simcoe_mean_W_sp_pairs_2010, simcoe_lc_sp_pairs_2010, by=c("Year"))
simcoe_sp_pair_ave_stab_2010 <- merge(simcoe_sp_pair_ave_stab_2010, simcoe_mean_mp_sp_pairs_2010, by=c("Year"))
simcoe_sp_pair_ave_stab_2010 <- merge(simcoe_sp_pair_ave_stab_2010, simcoe_mean_cc_sp_pairs_2010, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_2010_final <- merge(simcoe_2010_W, simcoe_2010_wa_pop_var_2010, by = c('Year'))
simcoe_2010_final <- simcoe_2010_final %>% group_by(Year) %>%
  mutate(gcv_2010 = lcv - W_2010) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_2010, 
            lcv = lcv,
            W = W_2010,
            lc_stab = as2_2010,
            mp_stab = as4_2010,
            cc_stab = as5_2010,
            mc_cv = sqrt(gcv_2010),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_2010/lcv + as4_2010/lcv + as5_2010/lcv,
            lc_asynchrony = as2_2010/lcv,
            mp_asynchrony = as4_2010/lcv,
            cc_asynchrony= as5_2010/lcv)

simcoe_2010_cv <- simcoe_2010_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_2010_final$mc_sum_mean_density <- simcoe_2010_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 2010 #####
simcoe_mean_corr_2010_c9 <- simcoe_species_pairs_count_2010_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2010, ind_W_2010)
simcoe_mean_corr_2010_k42 <- simcoe_species_pairs_count_2010_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2010, ind_W_2010)
simcoe_mean_corr_2010_k45 <- simcoe_species_pairs_count_2010_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2010, ind_W_2010)
simcoe_mean_corr_2010_lc <- rbind(simcoe_mean_corr_2010_c9, simcoe_mean_corr_2010_k42, simcoe_mean_corr_2010_k45)
simcoe_mean_corr_2010_mp <- simcoe_species_pairs_count_2010_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_2010, ind_W_2010)
simcoe_mean_corr_2010_cc <- simcoe_species_pairs_count_2010_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_2010, ind_W_2010)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_2010_lc[sample(1:nrow(simcoe_mean_corr_2010_lc),  nrow(simcoe_mean_corr_2010_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_2010),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_2010),
              mean_lc_stab = mean(ind_W_2010),
              sd_lc_stab = sd(ind_W_2010))
  
  lc_corr_2010 <- as.data.frame(mean_lc)
}
lc_corr_2010

mp_corr_2010 <- simcoe_mean_corr_2010_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_2010),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_2010),
                                                       mean_mp_stab = mean(ind_W_2010),
                                                       sd_mp_stab = sd(ind_W_2010))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_2010_cc[sample(1:nrow(simcoe_mean_corr_2010_cc), nrow(simcoe_mean_corr_2010_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_2010),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_2010),
              mean_cc_stab = mean(ind_W_2010),
              sd_cc_stab = sd(ind_W_2010))
  
  cc_corr_2010 <- as.data.frame(mean_cc)
}
cc_corr_2010


corr_2010 <- cbind(lc_corr_2010, mp_corr_2010)
corr_2010 <- cbind(corr_2010, cc_corr_2010)
corr_2010 <- corr_2010 %>% mutate(Year = 2010)





##### 2011 #####
simcoe_2011 <- subset(simcoe, Year == '2011')


#Transpose dataframe to add 0s within years
simcoe_2011_t1 <- dcast(simcoe_2011, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2011_m <- melt(simcoe_2011_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_2011 <- subset(simcoe_2011_m, Station_ID == 'C9')
simcoe_c9_2011$Species <- paste("C9", simcoe_c9_2011$Species, sep=" ")
simcoe_k42_2011 <- subset(simcoe_2011_m, Station_ID == 'K42')
simcoe_k42_2011$Species <- paste("K42", simcoe_k42_2011$Species, sep=" ")
simcoe_k45_2011 <- subset(simcoe_2011_m, Station_ID == 'K45')
simcoe_k45_2011$Species <- paste("K45", simcoe_k45_2011$Species, sep=" ")

#recombine dataframes
simcoe_spec_2011 <- rbind(simcoe_c9_2011, simcoe_k42_2011, simcoe_k45_2011)
simcoe_spec_2011 <- simcoe_spec_2011 %>% select(- Station_ID)

#transpose dataframe
simcoe_2011_t2 <- dcast(simcoe_spec_2011, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2011_t2 <- log(simcoe_2011_t2 + 1)
simcoe_2011_t2 <- simcoe_2011_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_2011_cv <- simcoe_2011_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_2011 = mean(Density),
            pop_sd_2011 = sd(Density)) %>%
  mutate(uw_pop_cv_2011 = pop_sd_2011/pop_mean_2011)

simcoe_2011_cv <- simcoe_2011_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_2011),
         mean_pop_cv = mean(uw_pop_cv_2011, na.rm = T),
         mean_pop_density = sum(pop_mean_2011),
         mean_pop_variance = sum(pop_sd_2011))

simcoe_2011_cv <- simcoe_2011_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_2011/mc_sum_mean_density)%>%
  mutate(w_pop_cv_2011 = (pop_mean_2011/mc_sum_mean_density)*uw_pop_cv_2011) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_2011_wa_pop_var_2011 <- simcoe_2011_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_2011 = sum(w_pop_cv_2011, na.rm = T),
            lcv = sum(w_pop_cv_2011, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_2011_cor_list_ts <- as.dist(round(cor(simcoe_2011_t2[]),2))
simcoe_2011_cor_ts <- stack(simcoe_2011_cor_list_ts, dim.names = TRUE)
simcoe_2011_cor_ts <- simcoe_2011_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_2011_cor - simcoe_2011_cv
simcoe_c9_2011_w_cv_ts <- subset(simcoe_2011_cv, Station_ID == 'C9')
simcoe_c9_2011_w_cv_ts$Species <- paste("C9", simcoe_c9_2011_w_cv_ts$Species, sep=" ")
simcoe_c9_2011_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_2011_w_cv_ts)
simcoe_k42_2011_w_cv_ts <- subset(simcoe_2011_cv, Station_ID == 'K42')
simcoe_k42_2011_w_cv_ts$Species <- paste("K42", simcoe_k42_2011_w_cv_ts$Species, sep=" ")
simcoe_k42_2011_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_2011_w_cv_ts)
simcoe_k45_2011_w_cv_ts <- subset(simcoe_2011_cv, Station_ID == 'K45')
simcoe_k45_2011_w_cv_ts$Species <- paste("K45", simcoe_k45_2011_w_cv_ts$Species, sep=" ")
simcoe_k45_2011_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_2011_w_cv_ts)

simcoe_spec_2011_w_sp_ts <- rbind(simcoe_c9_2011_w_cv_ts, simcoe_k42_2011_w_cv_ts, simcoe_k45_2011_w_cv_ts)
simcoe_spec_2011_w_sp_ts <- simcoe_spec_2011_w_sp_ts %>% select(Species, w_pop_cv_2011, Species_ID, Station_ID)
simcoe_spec_2011_w_sp1_ts <- simcoe_spec_2011_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_2011_sp1 = w_pop_cv_2011) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_2011_sp1, Species_ID, station_1)

simcoe_2011_cor_cv_ts <- merge(simcoe_2011_cor_ts, simcoe_spec_2011_w_sp1_ts, by = "species_1")

simcoe_spec_2011_w_sp2_ts <- simcoe_spec_2011_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_2011_sp2 = w_pop_cv_2011) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_2011_sp2, Species_ID, station_2)

simcoe_2011_cor_cv_ts <- merge(simcoe_2011_cor_cv_ts, simcoe_spec_2011_w_sp2_ts, by = "species_2")

simcoe_2011_cor_cv_ts_omit <- na.omit(simcoe_2011_cor_cv_ts)

simcoe_2011_ind_W_ts <- simcoe_2011_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_2011 = w_pop_cv_2011_sp1*w_pop_cv_2011_sp2,
         ind_W_2011 = (1 - corr)*(w_pop_cv_2011_sp1*w_pop_cv_2011_sp2),
         Year = 2011)

simcoe_2011_ind_W_ts$number_species_pairs <- nrow(simcoe_2011_ind_W_ts)

simcoe_mean_W_sp_pairs_2011 <- simcoe_2011_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_2011, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_2011_W_ts <- simcoe_2011_ind_W_ts %>% group_by(Year) %>%
  summarize(W_2011 = sum(ind_W_2011, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_2011_ind_W_ts_c9 <- filter(simcoe_2011_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_2011_ind_W_ts_k42 <- filter(simcoe_2011_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_2011_ind_W_ts_k45 <- filter(simcoe_2011_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_2011_c9 <- na.omit(simcoe_2011_ind_W_ts_c9)

simcoe_species_pairs_count_2011_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_2011_c9)

simcoe_mean_W_sp_pairs_2011_c9 <- simcoe_species_pairs_count_2011_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_2011, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2011_W_ts_c9 <- simcoe_2011_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_2011 = sum(ind_W_2011, na.rm = TRUE))

#K42
simcoe_species_pairs_count_2011_k42 <- na.omit(simcoe_2011_ind_W_ts_k42)

simcoe_species_pairs_count_2011_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_2011_k42)

simcoe_mean_W_sp_pairs_2011_k42 <- simcoe_species_pairs_count_2011_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_2011, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2011_W_ts_k42 <- simcoe_2011_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_2011 = sum(ind_W_2011, na.rm = TRUE))

#K45
simcoe_species_pairs_count_2011_k45 <- na.omit(simcoe_2011_ind_W_ts_k45)

simcoe_species_pairs_count_2011_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_2011_k45)

simcoe_mean_W_sp_pairs_2011_k45 <- simcoe_species_pairs_count_2011_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_2011, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2011_W_ts_k45 <- simcoe_2011_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_2011 = sum(ind_W_2011, na.rm = TRUE))

#Step 4
simcoe_2011_ind_W_as1 <- rbind(simcoe_2011_W_ts_c9, simcoe_2011_W_ts_k42, simcoe_2011_W_ts_k45)

simcoe_2011_W_as1 <- simcoe_2011_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_2011 = sum(W_2011, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_2011 <- rbind(simcoe_mean_W_sp_pairs_2011_c9, simcoe_mean_W_sp_pairs_2011_k42, simcoe_mean_W_sp_pairs_2011_k45)

simcoe_lc_sp_pairs_2011 <- simcoe_mean_lc_sp_pairs_2011 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_2011_ind_W_as4 <- filter(simcoe_2011_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_2011_mp  <- na.omit(simcoe_2011_ind_W_as4)

simcoe_species_pairs_count_2011_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_2011_mp)

simcoe_mean_mp_sp_pairs_2011  <- simcoe_species_pairs_count_2011_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_2011, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2011_sum_W_as4 <- simcoe_2011_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_2011 = sum(ind_W_2011, na.rm = TRUE))


#Step 3
simcoe_2011_sum_W_as4 <- simcoe_2011_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_2011 = sum(as4_2011, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_2011_ind_W_as5 <- filter(simcoe_2011_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_2011_cc <- na.omit(simcoe_2011_ind_W_as5)

simcoe_species_pairs_count_2011_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_2011_cc)

simcoe_mean_cc_sp_pairs_2011  <- simcoe_species_pairs_count_2011_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_2011, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2011_sum_W_as5 <- simcoe_2011_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_2011 = sum(ind_W_2011, na.rm = TRUE))


#Step 3
simcoe_2011_sum_W_as5 <- simcoe_2011_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_2011 = sum(as5_2011, na.rm = TRUE))

simcoe_2011_W <- merge(simcoe_2011_sum_W_as4, simcoe_2011_W_ts, by = c('Year'))
simcoe_2011_W <- merge(simcoe_2011_W, simcoe_2011_W_as1, by = c('Year'))
simcoe_2011_W <- merge(simcoe_2011_W, simcoe_2011_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_2011 <- merge(simcoe_mean_W_sp_pairs_2011, simcoe_lc_sp_pairs_2011, by=c("Year"))
simcoe_sp_pair_ave_stab_2011 <- merge(simcoe_sp_pair_ave_stab_2011, simcoe_mean_mp_sp_pairs_2011, by=c("Year"))
simcoe_sp_pair_ave_stab_2011 <- merge(simcoe_sp_pair_ave_stab_2011, simcoe_mean_cc_sp_pairs_2011, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_2011_final <- merge(simcoe_2011_W, simcoe_2011_wa_pop_var_2011, by = c('Year'))
simcoe_2011_final <- simcoe_2011_final %>% group_by(Year) %>%
  mutate(gcv_2011 = lcv - W_2011) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_2011, 
            lcv = lcv,
            W = W_2011,
            lc_stab = as2_2011,
            mp_stab = as4_2011,
            cc_stab = as5_2011,
            mc_cv = sqrt(gcv_2011),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_2011/lcv + as4_2011/lcv + as5_2011/lcv,
            lc_asynchrony = as2_2011/lcv,
            mp_asynchrony = as4_2011/lcv,
            cc_asynchrony= as5_2011/lcv)

simcoe_2011_cv <- simcoe_2011_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_2011_final$mc_sum_mean_density <- simcoe_2011_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 2011 #####
simcoe_mean_corr_2011_c9 <- simcoe_species_pairs_count_2011_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2011, ind_W_2011)
simcoe_mean_corr_2011_k42 <- simcoe_species_pairs_count_2011_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2011, ind_W_2011)
simcoe_mean_corr_2011_k45 <- simcoe_species_pairs_count_2011_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2011, ind_W_2011)
simcoe_mean_corr_2011_lc <- rbind(simcoe_mean_corr_2011_c9, simcoe_mean_corr_2011_k42, simcoe_mean_corr_2011_k45)
simcoe_mean_corr_2011_mp <- simcoe_species_pairs_count_2011_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_2011, ind_W_2011)
simcoe_mean_corr_2011_cc <- simcoe_species_pairs_count_2011_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_2011, ind_W_2011)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_2011_lc[sample(1:nrow(simcoe_mean_corr_2011_lc),  nrow(simcoe_mean_corr_2011_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_2011),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_2011),
              mean_lc_stab = mean(ind_W_2011),
              sd_lc_stab = sd(ind_W_2011))
  
  lc_corr_2011 <- as.data.frame(mean_lc)
}
lc_corr_2011

mp_corr_2011 <- simcoe_mean_corr_2011_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_2011),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_2011),
                                                       mean_mp_stab = mean(ind_W_2011),
                                                       sd_mp_stab = sd(ind_W_2011))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_2011_cc[sample(1:nrow(simcoe_mean_corr_2011_cc), nrow(simcoe_mean_corr_2011_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_2011),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_2011),
              mean_cc_stab = mean(ind_W_2011),
              sd_cc_stab = sd(ind_W_2011))
  
  cc_corr_2011 <- as.data.frame(mean_cc)
}
cc_corr_2011


corr_2011 <- cbind(lc_corr_2011, mp_corr_2011)
corr_2011 <- cbind(corr_2011, cc_corr_2011)
corr_2011 <- corr_2011 %>% mutate(Year = 2011)





##### 2012 #####
simcoe_2012 <- subset(simcoe, Year == '2012')


#Transpose dataframe to add 0s within years
simcoe_2012_t1 <- dcast(simcoe_2012, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2012_m <- melt(simcoe_2012_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_2012 <- subset(simcoe_2012_m, Station_ID == 'C9')
simcoe_c9_2012$Species <- paste("C9", simcoe_c9_2012$Species, sep=" ")
simcoe_k42_2012 <- subset(simcoe_2012_m, Station_ID == 'K42')
simcoe_k42_2012$Species <- paste("K42", simcoe_k42_2012$Species, sep=" ")
simcoe_k45_2012 <- subset(simcoe_2012_m, Station_ID == 'K45')
simcoe_k45_2012$Species <- paste("K45", simcoe_k45_2012$Species, sep=" ")

#recombine dataframes
simcoe_spec_2012 <- rbind(simcoe_c9_2012, simcoe_k42_2012, simcoe_k45_2012)
simcoe_spec_2012 <- simcoe_spec_2012 %>% select(- Station_ID)

#transpose dataframe
simcoe_2012_t2 <- dcast(simcoe_spec_2012, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2012_t2 <- log(simcoe_2012_t2 + 1)
simcoe_2012_t2 <- simcoe_2012_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_2012_cv <- simcoe_2012_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_2012 = mean(Density),
            pop_sd_2012 = sd(Density)) %>%
  mutate(uw_pop_cv_2012 = pop_sd_2012/pop_mean_2012)

simcoe_2012_cv <- simcoe_2012_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_2012),
         mean_pop_cv = mean(uw_pop_cv_2012, na.rm = T),
         mean_pop_density = sum(pop_mean_2012),
         mean_pop_variance = sum(pop_sd_2012))

simcoe_2012_cv <- simcoe_2012_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_2012/mc_sum_mean_density)%>%
  mutate(w_pop_cv_2012 = (pop_mean_2012/mc_sum_mean_density)*uw_pop_cv_2012) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_2012_wa_pop_var_2012 <- simcoe_2012_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_2012 = sum(w_pop_cv_2012, na.rm = T),
            lcv = sum(w_pop_cv_2012, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_2012_cor_list_ts <- as.dist(round(cor(simcoe_2012_t2[]),2))
simcoe_2012_cor_ts <- stack(simcoe_2012_cor_list_ts, dim.names = TRUE)
simcoe_2012_cor_ts <- simcoe_2012_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_2012_cor - simcoe_2012_cv
simcoe_c9_2012_w_cv_ts <- subset(simcoe_2012_cv, Station_ID == 'C9')
simcoe_c9_2012_w_cv_ts$Species <- paste("C9", simcoe_c9_2012_w_cv_ts$Species, sep=" ")
simcoe_c9_2012_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_2012_w_cv_ts)
simcoe_k42_2012_w_cv_ts <- subset(simcoe_2012_cv, Station_ID == 'K42')
simcoe_k42_2012_w_cv_ts$Species <- paste("K42", simcoe_k42_2012_w_cv_ts$Species, sep=" ")
simcoe_k42_2012_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_2012_w_cv_ts)
simcoe_k45_2012_w_cv_ts <- subset(simcoe_2012_cv, Station_ID == 'K45')
simcoe_k45_2012_w_cv_ts$Species <- paste("K45", simcoe_k45_2012_w_cv_ts$Species, sep=" ")
simcoe_k45_2012_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_2012_w_cv_ts)

simcoe_spec_2012_w_sp_ts <- rbind(simcoe_c9_2012_w_cv_ts, simcoe_k42_2012_w_cv_ts, simcoe_k45_2012_w_cv_ts)
simcoe_spec_2012_w_sp_ts <- simcoe_spec_2012_w_sp_ts %>% select(Species, w_pop_cv_2012, Species_ID, Station_ID)
simcoe_spec_2012_w_sp1_ts <- simcoe_spec_2012_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_2012_sp1 = w_pop_cv_2012) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_2012_sp1, Species_ID, station_1)

simcoe_2012_cor_cv_ts <- merge(simcoe_2012_cor_ts, simcoe_spec_2012_w_sp1_ts, by = "species_1")

simcoe_spec_2012_w_sp2_ts <- simcoe_spec_2012_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_2012_sp2 = w_pop_cv_2012) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_2012_sp2, Species_ID, station_2)

simcoe_2012_cor_cv_ts <- merge(simcoe_2012_cor_cv_ts, simcoe_spec_2012_w_sp2_ts, by = "species_2")

simcoe_2012_cor_cv_ts_omit <- na.omit(simcoe_2012_cor_cv_ts)

simcoe_2012_ind_W_ts <- simcoe_2012_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_2012 = w_pop_cv_2012_sp1*w_pop_cv_2012_sp2,
         ind_W_2012 = (1 - corr)*(w_pop_cv_2012_sp1*w_pop_cv_2012_sp2),
         Year = 2012)

simcoe_2012_ind_W_ts$number_species_pairs <- nrow(simcoe_2012_ind_W_ts)

simcoe_mean_W_sp_pairs_2012 <- simcoe_2012_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_2012, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_2012_W_ts <- simcoe_2012_ind_W_ts %>% group_by(Year) %>%
  summarize(W_2012 = sum(ind_W_2012, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_2012_ind_W_ts_c9 <- filter(simcoe_2012_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_2012_ind_W_ts_k42 <- filter(simcoe_2012_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_2012_ind_W_ts_k45 <- filter(simcoe_2012_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_2012_c9 <- na.omit(simcoe_2012_ind_W_ts_c9)

simcoe_species_pairs_count_2012_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_2012_c9)

simcoe_mean_W_sp_pairs_2012_c9 <- simcoe_species_pairs_count_2012_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_2012, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2012_W_ts_c9 <- simcoe_2012_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_2012 = sum(ind_W_2012, na.rm = TRUE))

#K42
simcoe_species_pairs_count_2012_k42 <- na.omit(simcoe_2012_ind_W_ts_k42)

simcoe_species_pairs_count_2012_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_2012_k42)

simcoe_mean_W_sp_pairs_2012_k42 <- simcoe_species_pairs_count_2012_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_2012, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2012_W_ts_k42 <- simcoe_2012_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_2012 = sum(ind_W_2012, na.rm = TRUE))

#K45
simcoe_species_pairs_count_2012_k45 <- na.omit(simcoe_2012_ind_W_ts_k45)

simcoe_species_pairs_count_2012_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_2012_k45)

simcoe_mean_W_sp_pairs_2012_k45 <- simcoe_species_pairs_count_2012_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_2012, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2012_W_ts_k45 <- simcoe_2012_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_2012 = sum(ind_W_2012, na.rm = TRUE))

#Step 4
simcoe_2012_ind_W_as1 <- rbind(simcoe_2012_W_ts_c9, simcoe_2012_W_ts_k42, simcoe_2012_W_ts_k45)

simcoe_2012_W_as1 <- simcoe_2012_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_2012 = sum(W_2012, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_2012 <- rbind(simcoe_mean_W_sp_pairs_2012_c9, simcoe_mean_W_sp_pairs_2012_k42, simcoe_mean_W_sp_pairs_2012_k45)

simcoe_lc_sp_pairs_2012 <- simcoe_mean_lc_sp_pairs_2012 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_2012_ind_W_as4 <- filter(simcoe_2012_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_2012_mp  <- na.omit(simcoe_2012_ind_W_as4)

simcoe_species_pairs_count_2012_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_2012_mp)

simcoe_mean_mp_sp_pairs_2012  <- simcoe_species_pairs_count_2012_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_2012, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2012_sum_W_as4 <- simcoe_2012_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_2012 = sum(ind_W_2012, na.rm = TRUE))


#Step 3
simcoe_2012_sum_W_as4 <- simcoe_2012_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_2012 = sum(as4_2012, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_2012_ind_W_as5 <- filter(simcoe_2012_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_2012_cc <- na.omit(simcoe_2012_ind_W_as5)

simcoe_species_pairs_count_2012_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_2012_cc)

simcoe_mean_cc_sp_pairs_2012  <- simcoe_species_pairs_count_2012_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_2012, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2012_sum_W_as5 <- simcoe_2012_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_2012 = sum(ind_W_2012, na.rm = TRUE))


#Step 3
simcoe_2012_sum_W_as5 <- simcoe_2012_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_2012 = sum(as5_2012, na.rm = TRUE))

simcoe_2012_W <- merge(simcoe_2012_sum_W_as4, simcoe_2012_W_ts, by = c('Year'))
simcoe_2012_W <- merge(simcoe_2012_W, simcoe_2012_W_as1, by = c('Year'))
simcoe_2012_W <- merge(simcoe_2012_W, simcoe_2012_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_2012 <- merge(simcoe_mean_W_sp_pairs_2012, simcoe_lc_sp_pairs_2012, by=c("Year"))
simcoe_sp_pair_ave_stab_2012 <- merge(simcoe_sp_pair_ave_stab_2012, simcoe_mean_mp_sp_pairs_2012, by=c("Year"))
simcoe_sp_pair_ave_stab_2012 <- merge(simcoe_sp_pair_ave_stab_2012, simcoe_mean_cc_sp_pairs_2012, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_2012_final <- merge(simcoe_2012_W, simcoe_2012_wa_pop_var_2012, by = c('Year'))
simcoe_2012_final <- simcoe_2012_final %>% group_by(Year) %>%
  mutate(gcv_2012 = lcv - W_2012) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_2012, 
            lcv = lcv,
            W = W_2012,
            lc_stab = as2_2012,
            mp_stab = as4_2012,
            cc_stab = as5_2012,
            mc_cv = sqrt(gcv_2012),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_2012/lcv + as4_2012/lcv + as5_2012/lcv,
            lc_asynchrony = as2_2012/lcv,
            mp_asynchrony = as4_2012/lcv,
            cc_asynchrony= as5_2012/lcv)

simcoe_2012_cv <- simcoe_2012_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_2012_final$mc_sum_mean_density <- simcoe_2012_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 2012 #####
simcoe_mean_corr_2012_c9 <- simcoe_species_pairs_count_2012_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2012, ind_W_2012)
simcoe_mean_corr_2012_k42 <- simcoe_species_pairs_count_2012_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2012, ind_W_2012)
simcoe_mean_corr_2012_k45 <- simcoe_species_pairs_count_2012_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2012, ind_W_2012)
simcoe_mean_corr_2012_lc <- rbind(simcoe_mean_corr_2012_c9, simcoe_mean_corr_2012_k42, simcoe_mean_corr_2012_k45)
simcoe_mean_corr_2012_mp <- simcoe_species_pairs_count_2012_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_2012, ind_W_2012)
simcoe_mean_corr_2012_cc <- simcoe_species_pairs_count_2012_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_2012, ind_W_2012)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_2012_lc[sample(1:nrow(simcoe_mean_corr_2012_lc),  nrow(simcoe_mean_corr_2012_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_2012),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_2012),
              mean_lc_stab = mean(ind_W_2012),
              sd_lc_stab = sd(ind_W_2012))
  
  lc_corr_2012 <- as.data.frame(mean_lc)
}
lc_corr_2012

mp_corr_2012 <- simcoe_mean_corr_2012_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_2012),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_2012),
                                                       mean_mp_stab = mean(ind_W_2012),
                                                       sd_mp_stab = sd(ind_W_2012))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_2012_cc[sample(1:nrow(simcoe_mean_corr_2012_cc), nrow(simcoe_mean_corr_2012_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_2012),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_2012),
              mean_cc_stab = mean(ind_W_2012),
              sd_cc_stab = sd(ind_W_2012))
  
  cc_corr_2012 <- as.data.frame(mean_cc)
}
cc_corr_2012


corr_2012 <- cbind(lc_corr_2012, mp_corr_2012)
corr_2012 <- cbind(corr_2012, cc_corr_2012)
corr_2012 <- corr_2012 %>% mutate(Year = 2012)





##### 2013 #####
simcoe_2013 <- subset(simcoe, Year == '2013')


#Transpose dataframe to add 0s within years
simcoe_2013_t1 <- dcast(simcoe_2013, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2013_m <- melt(simcoe_2013_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_2013 <- subset(simcoe_2013_m, Station_ID == 'C9')
simcoe_c9_2013$Species <- paste("C9", simcoe_c9_2013$Species, sep=" ")
simcoe_k42_2013 <- subset(simcoe_2013_m, Station_ID == 'K42')
simcoe_k42_2013$Species <- paste("K42", simcoe_k42_2013$Species, sep=" ")
simcoe_k45_2013 <- subset(simcoe_2013_m, Station_ID == 'K45')
simcoe_k45_2013$Species <- paste("K45", simcoe_k45_2013$Species, sep=" ")

#recombine dataframes
simcoe_spec_2013 <- rbind(simcoe_c9_2013, simcoe_k42_2013, simcoe_k45_2013)
simcoe_spec_2013 <- simcoe_spec_2013 %>% select(- Station_ID)

#transpose dataframe
simcoe_2013_t2 <- dcast(simcoe_spec_2013, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2013_t2 <- log(simcoe_2013_t2 + 1)
simcoe_2013_t2 <- simcoe_2013_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_2013_cv <- simcoe_2013_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_2013 = mean(Density),
            pop_sd_2013 = sd(Density)) %>%
  mutate(uw_pop_cv_2013 = pop_sd_2013/pop_mean_2013)

simcoe_2013_cv <- simcoe_2013_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_2013),
         mean_pop_cv = mean(uw_pop_cv_2013, na.rm = T),
         mean_pop_density = sum(pop_mean_2013),
         mean_pop_variance = sum(pop_sd_2013))

simcoe_2013_cv <- simcoe_2013_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_2013/mc_sum_mean_density)%>%
  mutate(w_pop_cv_2013 = (pop_mean_2013/mc_sum_mean_density)*uw_pop_cv_2013) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_2013_wa_pop_var_2013 <- simcoe_2013_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_2013 = sum(w_pop_cv_2013, na.rm = T),
            lcv = sum(w_pop_cv_2013, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_2013_cor_list_ts <- as.dist(round(cor(simcoe_2013_t2[]),2))
simcoe_2013_cor_ts <- stack(simcoe_2013_cor_list_ts, dim.names = TRUE)
simcoe_2013_cor_ts <- simcoe_2013_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_2013_cor - simcoe_2013_cv
simcoe_c9_2013_w_cv_ts <- subset(simcoe_2013_cv, Station_ID == 'C9')
simcoe_c9_2013_w_cv_ts$Species <- paste("C9", simcoe_c9_2013_w_cv_ts$Species, sep=" ")
simcoe_c9_2013_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_2013_w_cv_ts)
simcoe_k42_2013_w_cv_ts <- subset(simcoe_2013_cv, Station_ID == 'K42')
simcoe_k42_2013_w_cv_ts$Species <- paste("K42", simcoe_k42_2013_w_cv_ts$Species, sep=" ")
simcoe_k42_2013_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_2013_w_cv_ts)
simcoe_k45_2013_w_cv_ts <- subset(simcoe_2013_cv, Station_ID == 'K45')
simcoe_k45_2013_w_cv_ts$Species <- paste("K45", simcoe_k45_2013_w_cv_ts$Species, sep=" ")
simcoe_k45_2013_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_2013_w_cv_ts)

simcoe_spec_2013_w_sp_ts <- rbind(simcoe_c9_2013_w_cv_ts, simcoe_k42_2013_w_cv_ts, simcoe_k45_2013_w_cv_ts)
simcoe_spec_2013_w_sp_ts <- simcoe_spec_2013_w_sp_ts %>% select(Species, w_pop_cv_2013, Species_ID, Station_ID)
simcoe_spec_2013_w_sp1_ts <- simcoe_spec_2013_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_2013_sp1 = w_pop_cv_2013) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_2013_sp1, Species_ID, station_1)

simcoe_2013_cor_cv_ts <- merge(simcoe_2013_cor_ts, simcoe_spec_2013_w_sp1_ts, by = "species_1")

simcoe_spec_2013_w_sp2_ts <- simcoe_spec_2013_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_2013_sp2 = w_pop_cv_2013) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_2013_sp2, Species_ID, station_2)

simcoe_2013_cor_cv_ts <- merge(simcoe_2013_cor_cv_ts, simcoe_spec_2013_w_sp2_ts, by = "species_2")

simcoe_2013_cor_cv_ts_omit <- na.omit(simcoe_2013_cor_cv_ts)

simcoe_2013_ind_W_ts <- simcoe_2013_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_2013 = w_pop_cv_2013_sp1*w_pop_cv_2013_sp2,
         ind_W_2013 = (1 - corr)*(w_pop_cv_2013_sp1*w_pop_cv_2013_sp2),
         Year = 2013)

simcoe_2013_ind_W_ts$number_species_pairs <- nrow(simcoe_2013_ind_W_ts)

simcoe_mean_W_sp_pairs_2013 <- simcoe_2013_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_2013, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_2013_W_ts <- simcoe_2013_ind_W_ts %>% group_by(Year) %>%
  summarize(W_2013 = sum(ind_W_2013, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_2013_ind_W_ts_c9 <- filter(simcoe_2013_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_2013_ind_W_ts_k42 <- filter(simcoe_2013_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_2013_ind_W_ts_k45 <- filter(simcoe_2013_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_2013_c9 <- na.omit(simcoe_2013_ind_W_ts_c9)

simcoe_species_pairs_count_2013_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_2013_c9)

simcoe_mean_W_sp_pairs_2013_c9 <- simcoe_species_pairs_count_2013_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_2013, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2013_W_ts_c9 <- simcoe_2013_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_2013 = sum(ind_W_2013, na.rm = TRUE))

#K42
simcoe_species_pairs_count_2013_k42 <- na.omit(simcoe_2013_ind_W_ts_k42)

simcoe_species_pairs_count_2013_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_2013_k42)

simcoe_mean_W_sp_pairs_2013_k42 <- simcoe_species_pairs_count_2013_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_2013, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2013_W_ts_k42 <- simcoe_2013_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_2013 = sum(ind_W_2013, na.rm = TRUE))

#K45
simcoe_species_pairs_count_2013_k45 <- na.omit(simcoe_2013_ind_W_ts_k45)

simcoe_species_pairs_count_2013_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_2013_k45)

simcoe_mean_W_sp_pairs_2013_k45 <- simcoe_species_pairs_count_2013_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_2013, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2013_W_ts_k45 <- simcoe_2013_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_2013 = sum(ind_W_2013, na.rm = TRUE))

#Step 4
simcoe_2013_ind_W_as1 <- rbind(simcoe_2013_W_ts_c9, simcoe_2013_W_ts_k42, simcoe_2013_W_ts_k45)

simcoe_2013_W_as1 <- simcoe_2013_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_2013 = sum(W_2013, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_2013 <- rbind(simcoe_mean_W_sp_pairs_2013_c9, simcoe_mean_W_sp_pairs_2013_k42, simcoe_mean_W_sp_pairs_2013_k45)

simcoe_lc_sp_pairs_2013 <- simcoe_mean_lc_sp_pairs_2013 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_2013_ind_W_as4 <- filter(simcoe_2013_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_2013_mp  <- na.omit(simcoe_2013_ind_W_as4)

simcoe_species_pairs_count_2013_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_2013_mp)

simcoe_mean_mp_sp_pairs_2013  <- simcoe_species_pairs_count_2013_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_2013, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2013_sum_W_as4 <- simcoe_2013_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_2013 = sum(ind_W_2013, na.rm = TRUE))


#Step 3
simcoe_2013_sum_W_as4 <- simcoe_2013_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_2013 = sum(as4_2013, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_2013_ind_W_as5 <- filter(simcoe_2013_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_2013_cc <- na.omit(simcoe_2013_ind_W_as5)

simcoe_species_pairs_count_2013_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_2013_cc)

simcoe_mean_cc_sp_pairs_2013  <- simcoe_species_pairs_count_2013_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_2013, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2013_sum_W_as5 <- simcoe_2013_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_2013 = sum(ind_W_2013, na.rm = TRUE))


#Step 3
simcoe_2013_sum_W_as5 <- simcoe_2013_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_2013 = sum(as5_2013, na.rm = TRUE))

simcoe_2013_W <- merge(simcoe_2013_sum_W_as4, simcoe_2013_W_ts, by = c('Year'))
simcoe_2013_W <- merge(simcoe_2013_W, simcoe_2013_W_as1, by = c('Year'))
simcoe_2013_W <- merge(simcoe_2013_W, simcoe_2013_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_2013 <- merge(simcoe_mean_W_sp_pairs_2013, simcoe_lc_sp_pairs_2013, by=c("Year"))
simcoe_sp_pair_ave_stab_2013 <- merge(simcoe_sp_pair_ave_stab_2013, simcoe_mean_mp_sp_pairs_2013, by=c("Year"))
simcoe_sp_pair_ave_stab_2013 <- merge(simcoe_sp_pair_ave_stab_2013, simcoe_mean_cc_sp_pairs_2013, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_2013_final <- merge(simcoe_2013_W, simcoe_2013_wa_pop_var_2013, by = c('Year'))
simcoe_2013_final <- simcoe_2013_final %>% group_by(Year) %>%
  mutate(gcv_2013 = lcv - W_2013) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_2013, 
            lcv = lcv,
            W = W_2013,
            lc_stab = as2_2013,
            mp_stab = as4_2013,
            cc_stab = as5_2013,
            mc_cv = sqrt(gcv_2013),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_2013/lcv + as4_2013/lcv + as5_2013/lcv,
            lc_asynchrony = as2_2013/lcv,
            mp_asynchrony = as4_2013/lcv,
            cc_asynchrony= as5_2013/lcv)

simcoe_2013_cv <- simcoe_2013_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_2013_final$mc_sum_mean_density <- simcoe_2013_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 2013 #####
simcoe_mean_corr_2013_c9 <- simcoe_species_pairs_count_2013_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2013, ind_W_2013)
simcoe_mean_corr_2013_k42 <- simcoe_species_pairs_count_2013_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2013, ind_W_2013)
simcoe_mean_corr_2013_k45 <- simcoe_species_pairs_count_2013_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2013, ind_W_2013)
simcoe_mean_corr_2013_lc <- rbind(simcoe_mean_corr_2013_c9, simcoe_mean_corr_2013_k42, simcoe_mean_corr_2013_k45)
simcoe_mean_corr_2013_mp <- simcoe_species_pairs_count_2013_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_2013, ind_W_2013)
simcoe_mean_corr_2013_cc <- simcoe_species_pairs_count_2013_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_2013, ind_W_2013)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_2013_lc[sample(1:nrow(simcoe_mean_corr_2013_lc),  nrow(simcoe_mean_corr_2013_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_2013),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_2013),
              mean_lc_stab = mean(ind_W_2013),
              sd_lc_stab = sd(ind_W_2013))
  
  lc_corr_2013 <- as.data.frame(mean_lc)
}
lc_corr_2013

mp_corr_2013 <- simcoe_mean_corr_2013_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_2013),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_2013),
                                                       mean_mp_stab = mean(ind_W_2013),
                                                       sd_mp_stab = sd(ind_W_2013))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_2013_cc[sample(1:nrow(simcoe_mean_corr_2013_cc), nrow(simcoe_mean_corr_2013_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_2013),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_2013),
              mean_cc_stab = mean(ind_W_2013),
              sd_cc_stab = sd(ind_W_2013))
  
  cc_corr_2013 <- as.data.frame(mean_cc)
}
cc_corr_2013


corr_2013 <- cbind(lc_corr_2013, mp_corr_2013)
corr_2013 <- cbind(corr_2013, cc_corr_2013)
corr_2013 <- corr_2013 %>% mutate(Year = 2013)





##### 2014 #####
simcoe_2014 <- subset(simcoe, Year == '2014')


#Transpose dataframe to add 0s within years
simcoe_2014_t1 <- dcast(simcoe_2014, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2014_m <- melt(simcoe_2014_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_2014 <- subset(simcoe_2014_m, Station_ID == 'C9')
simcoe_c9_2014$Species <- paste("C9", simcoe_c9_2014$Species, sep=" ")
simcoe_k42_2014 <- subset(simcoe_2014_m, Station_ID == 'K42')
simcoe_k42_2014$Species <- paste("K42", simcoe_k42_2014$Species, sep=" ")
simcoe_k45_2014 <- subset(simcoe_2014_m, Station_ID == 'K45')
simcoe_k45_2014$Species <- paste("K45", simcoe_k45_2014$Species, sep=" ")

#recombine dataframes
simcoe_spec_2014 <- rbind(simcoe_c9_2014, simcoe_k42_2014, simcoe_k45_2014)
simcoe_spec_2014 <- simcoe_spec_2014 %>% select(- Station_ID)

#transpose dataframe
simcoe_2014_t2 <- dcast(simcoe_spec_2014, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2014_t2 <- log(simcoe_2014_t2 + 1)
simcoe_2014_t2 <- simcoe_2014_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_2014_cv <- simcoe_2014_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_2014 = mean(Density),
            pop_sd_2014 = sd(Density)) %>%
  mutate(uw_pop_cv_2014 = pop_sd_2014/pop_mean_2014)

simcoe_2014_cv <- simcoe_2014_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_2014),
         mean_pop_cv = mean(uw_pop_cv_2014, na.rm = T),
         mean_pop_density = sum(pop_mean_2014),
         mean_pop_variance = sum(pop_sd_2014))

simcoe_2014_cv <- simcoe_2014_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_2014/mc_sum_mean_density)%>%
  mutate(w_pop_cv_2014 = (pop_mean_2014/mc_sum_mean_density)*uw_pop_cv_2014) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_2014_wa_pop_var_2014 <- simcoe_2014_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_2014 = sum(w_pop_cv_2014, na.rm = T),
            lcv = sum(w_pop_cv_2014, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_2014_cor_list_ts <- as.dist(round(cor(simcoe_2014_t2[]),2))
simcoe_2014_cor_ts <- stack(simcoe_2014_cor_list_ts, dim.names = TRUE)
simcoe_2014_cor_ts <- simcoe_2014_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_2014_cor - simcoe_2014_cv
simcoe_c9_2014_w_cv_ts <- subset(simcoe_2014_cv, Station_ID == 'C9')
simcoe_c9_2014_w_cv_ts$Species <- paste("C9", simcoe_c9_2014_w_cv_ts$Species, sep=" ")
simcoe_c9_2014_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_2014_w_cv_ts)
simcoe_k42_2014_w_cv_ts <- subset(simcoe_2014_cv, Station_ID == 'K42')
simcoe_k42_2014_w_cv_ts$Species <- paste("K42", simcoe_k42_2014_w_cv_ts$Species, sep=" ")
simcoe_k42_2014_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_2014_w_cv_ts)
simcoe_k45_2014_w_cv_ts <- subset(simcoe_2014_cv, Station_ID == 'K45')
simcoe_k45_2014_w_cv_ts$Species <- paste("K45", simcoe_k45_2014_w_cv_ts$Species, sep=" ")
simcoe_k45_2014_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_2014_w_cv_ts)

simcoe_spec_2014_w_sp_ts <- rbind(simcoe_c9_2014_w_cv_ts, simcoe_k42_2014_w_cv_ts, simcoe_k45_2014_w_cv_ts)
simcoe_spec_2014_w_sp_ts <- simcoe_spec_2014_w_sp_ts %>% select(Species, w_pop_cv_2014, Species_ID, Station_ID)
simcoe_spec_2014_w_sp1_ts <- simcoe_spec_2014_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_2014_sp1 = w_pop_cv_2014) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_2014_sp1, Species_ID, station_1)

simcoe_2014_cor_cv_ts <- merge(simcoe_2014_cor_ts, simcoe_spec_2014_w_sp1_ts, by = "species_1")

simcoe_spec_2014_w_sp2_ts <- simcoe_spec_2014_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_2014_sp2 = w_pop_cv_2014) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_2014_sp2, Species_ID, station_2)

simcoe_2014_cor_cv_ts <- merge(simcoe_2014_cor_cv_ts, simcoe_spec_2014_w_sp2_ts, by = "species_2")

simcoe_2014_cor_cv_ts_omit <- na.omit(simcoe_2014_cor_cv_ts)

simcoe_2014_ind_W_ts <- simcoe_2014_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_2014 = w_pop_cv_2014_sp1*w_pop_cv_2014_sp2,
         ind_W_2014 = (1 - corr)*(w_pop_cv_2014_sp1*w_pop_cv_2014_sp2),
         Year = 2014)

simcoe_2014_ind_W_ts$number_species_pairs <- nrow(simcoe_2014_ind_W_ts)

simcoe_mean_W_sp_pairs_2014 <- simcoe_2014_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_2014, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_2014_W_ts <- simcoe_2014_ind_W_ts %>% group_by(Year) %>%
  summarize(W_2014 = sum(ind_W_2014, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_2014_ind_W_ts_c9 <- filter(simcoe_2014_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_2014_ind_W_ts_k42 <- filter(simcoe_2014_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_2014_ind_W_ts_k45 <- filter(simcoe_2014_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_2014_c9 <- na.omit(simcoe_2014_ind_W_ts_c9)

simcoe_species_pairs_count_2014_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_2014_c9)

simcoe_mean_W_sp_pairs_2014_c9 <- simcoe_species_pairs_count_2014_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_2014, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2014_W_ts_c9 <- simcoe_2014_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_2014 = sum(ind_W_2014, na.rm = TRUE))

#K42
simcoe_species_pairs_count_2014_k42 <- na.omit(simcoe_2014_ind_W_ts_k42)

simcoe_species_pairs_count_2014_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_2014_k42)

simcoe_mean_W_sp_pairs_2014_k42 <- simcoe_species_pairs_count_2014_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_2014, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2014_W_ts_k42 <- simcoe_2014_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_2014 = sum(ind_W_2014, na.rm = TRUE))

#K45
simcoe_species_pairs_count_2014_k45 <- na.omit(simcoe_2014_ind_W_ts_k45)

simcoe_species_pairs_count_2014_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_2014_k45)

simcoe_mean_W_sp_pairs_2014_k45 <- simcoe_species_pairs_count_2014_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_2014, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2014_W_ts_k45 <- simcoe_2014_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_2014 = sum(ind_W_2014, na.rm = TRUE))

#Step 4
simcoe_2014_ind_W_as1 <- rbind(simcoe_2014_W_ts_c9, simcoe_2014_W_ts_k42, simcoe_2014_W_ts_k45)

simcoe_2014_W_as1 <- simcoe_2014_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_2014 = sum(W_2014, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_2014 <- rbind(simcoe_mean_W_sp_pairs_2014_c9, simcoe_mean_W_sp_pairs_2014_k42, simcoe_mean_W_sp_pairs_2014_k45)

simcoe_lc_sp_pairs_2014 <- simcoe_mean_lc_sp_pairs_2014 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_2014_ind_W_as4 <- filter(simcoe_2014_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_2014_mp  <- na.omit(simcoe_2014_ind_W_as4)

simcoe_species_pairs_count_2014_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_2014_mp)

simcoe_mean_mp_sp_pairs_2014  <- simcoe_species_pairs_count_2014_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_2014, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2014_sum_W_as4 <- simcoe_2014_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_2014 = sum(ind_W_2014, na.rm = TRUE))


#Step 3
simcoe_2014_sum_W_as4 <- simcoe_2014_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_2014 = sum(as4_2014, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_2014_ind_W_as5 <- filter(simcoe_2014_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_2014_cc <- na.omit(simcoe_2014_ind_W_as5)

simcoe_species_pairs_count_2014_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_2014_cc)

simcoe_mean_cc_sp_pairs_2014  <- simcoe_species_pairs_count_2014_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_2014, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2014_sum_W_as5 <- simcoe_2014_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_2014 = sum(ind_W_2014, na.rm = TRUE))


#Step 3
simcoe_2014_sum_W_as5 <- simcoe_2014_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_2014 = sum(as5_2014, na.rm = TRUE))

simcoe_2014_W <- merge(simcoe_2014_sum_W_as4, simcoe_2014_W_ts, by = c('Year'))
simcoe_2014_W <- merge(simcoe_2014_W, simcoe_2014_W_as1, by = c('Year'))
simcoe_2014_W <- merge(simcoe_2014_W, simcoe_2014_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_2014 <- merge(simcoe_mean_W_sp_pairs_2014, simcoe_lc_sp_pairs_2014, by=c("Year"))
simcoe_sp_pair_ave_stab_2014 <- merge(simcoe_sp_pair_ave_stab_2014, simcoe_mean_mp_sp_pairs_2014, by=c("Year"))
simcoe_sp_pair_ave_stab_2014 <- merge(simcoe_sp_pair_ave_stab_2014, simcoe_mean_cc_sp_pairs_2014, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_2014_final <- merge(simcoe_2014_W, simcoe_2014_wa_pop_var_2014, by = c('Year'))
simcoe_2014_final <- simcoe_2014_final %>% group_by(Year) %>%
  mutate(gcv_2014 = lcv - W_2014) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_2014, 
            lcv = lcv,
            W = W_2014,
            lc_stab = as2_2014,
            mp_stab = as4_2014,
            cc_stab = as5_2014,
            mc_cv = sqrt(gcv_2014),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_2014/lcv + as4_2014/lcv + as5_2014/lcv,
            lc_asynchrony = as2_2014/lcv,
            mp_asynchrony = as4_2014/lcv,
            cc_asynchrony= as5_2014/lcv)

simcoe_2014_cv <- simcoe_2014_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_2014_final$mc_sum_mean_density <- simcoe_2014_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 2014 #####
simcoe_mean_corr_2014_c9 <- simcoe_species_pairs_count_2014_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2014, ind_W_2014)
simcoe_mean_corr_2014_k42 <- simcoe_species_pairs_count_2014_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2014, ind_W_2014)
simcoe_mean_corr_2014_k45 <- simcoe_species_pairs_count_2014_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2014, ind_W_2014)
simcoe_mean_corr_2014_lc <- rbind(simcoe_mean_corr_2014_c9, simcoe_mean_corr_2014_k42, simcoe_mean_corr_2014_k45)
simcoe_mean_corr_2014_mp <- simcoe_species_pairs_count_2014_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_2014, ind_W_2014)
simcoe_mean_corr_2014_cc <- simcoe_species_pairs_count_2014_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_2014, ind_W_2014)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_2014_lc[sample(1:nrow(simcoe_mean_corr_2014_lc),  nrow(simcoe_mean_corr_2014_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_2014),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_2014),
              mean_lc_stab = mean(ind_W_2014),
              sd_lc_stab = sd(ind_W_2014))
  
  lc_corr_2014 <- as.data.frame(mean_lc)
}
lc_corr_2014

mp_corr_2014 <- simcoe_mean_corr_2014_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_2014),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_2014),
                                                       mean_mp_stab = mean(ind_W_2014),
                                                       sd_mp_stab = sd(ind_W_2014))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_2014_cc[sample(1:nrow(simcoe_mean_corr_2014_cc), nrow(simcoe_mean_corr_2014_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_2014),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_2014),
              mean_cc_stab = mean(ind_W_2014),
              sd_cc_stab = sd(ind_W_2014))
  
  cc_corr_2014 <- as.data.frame(mean_cc)
}
cc_corr_2014


corr_2014 <- cbind(lc_corr_2014, mp_corr_2014)
corr_2014 <- cbind(corr_2014, cc_corr_2014)
corr_2014 <- corr_2014 %>% mutate(Year = 2014)





##### 2015 #####
simcoe_2015 <- subset(simcoe, Year == '2015')


#Transpose dataframe to add 0s within years
simcoe_2015_t1 <- dcast(simcoe_2015, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2015_m <- melt(simcoe_2015_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_2015 <- subset(simcoe_2015_m, Station_ID == 'C9')
simcoe_c9_2015$Species <- paste("C9", simcoe_c9_2015$Species, sep=" ")
simcoe_k42_2015 <- subset(simcoe_2015_m, Station_ID == 'K42')
simcoe_k42_2015$Species <- paste("K42", simcoe_k42_2015$Species, sep=" ")
simcoe_k45_2015 <- subset(simcoe_2015_m, Station_ID == 'K45')
simcoe_k45_2015$Species <- paste("K45", simcoe_k45_2015$Species, sep=" ")

#recombine dataframes
simcoe_spec_2015 <- rbind(simcoe_c9_2015, simcoe_k42_2015, simcoe_k45_2015)
simcoe_spec_2015 <- simcoe_spec_2015 %>% select(- Station_ID)

#transpose dataframe
simcoe_2015_t2 <- dcast(simcoe_spec_2015, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2015_t2 <- log(simcoe_2015_t2 + 1)
simcoe_2015_t2 <- simcoe_2015_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_2015_cv <- simcoe_2015_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_2015 = mean(Density),
            pop_sd_2015 = sd(Density)) %>%
  mutate(uw_pop_cv_2015 = pop_sd_2015/pop_mean_2015)

simcoe_2015_cv <- simcoe_2015_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_2015),
         mean_pop_cv = mean(uw_pop_cv_2015, na.rm = T),
         mean_pop_density = sum(pop_mean_2015),
         mean_pop_variance = sum(pop_sd_2015))

simcoe_2015_cv <- simcoe_2015_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_2015/mc_sum_mean_density)%>%
  mutate(w_pop_cv_2015 = (pop_mean_2015/mc_sum_mean_density)*uw_pop_cv_2015) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_2015_wa_pop_var_2015 <- simcoe_2015_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_2015 = sum(w_pop_cv_2015, na.rm = T),
            lcv = sum(w_pop_cv_2015, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_2015_cor_list_ts <- as.dist(round(cor(simcoe_2015_t2[]),2))
simcoe_2015_cor_ts <- stack(simcoe_2015_cor_list_ts, dim.names = TRUE)
simcoe_2015_cor_ts <- simcoe_2015_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_2015_cor - simcoe_2015_cv
simcoe_c9_2015_w_cv_ts <- subset(simcoe_2015_cv, Station_ID == 'C9')
simcoe_c9_2015_w_cv_ts$Species <- paste("C9", simcoe_c9_2015_w_cv_ts$Species, sep=" ")
simcoe_c9_2015_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_2015_w_cv_ts)
simcoe_k42_2015_w_cv_ts <- subset(simcoe_2015_cv, Station_ID == 'K42')
simcoe_k42_2015_w_cv_ts$Species <- paste("K42", simcoe_k42_2015_w_cv_ts$Species, sep=" ")
simcoe_k42_2015_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_2015_w_cv_ts)
simcoe_k45_2015_w_cv_ts <- subset(simcoe_2015_cv, Station_ID == 'K45')
simcoe_k45_2015_w_cv_ts$Species <- paste("K45", simcoe_k45_2015_w_cv_ts$Species, sep=" ")
simcoe_k45_2015_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_2015_w_cv_ts)

simcoe_spec_2015_w_sp_ts <- rbind(simcoe_c9_2015_w_cv_ts, simcoe_k42_2015_w_cv_ts, simcoe_k45_2015_w_cv_ts)
simcoe_spec_2015_w_sp_ts <- simcoe_spec_2015_w_sp_ts %>% select(Species, w_pop_cv_2015, Species_ID, Station_ID)
simcoe_spec_2015_w_sp1_ts <- simcoe_spec_2015_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_2015_sp1 = w_pop_cv_2015) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_2015_sp1, Species_ID, station_1)

simcoe_2015_cor_cv_ts <- merge(simcoe_2015_cor_ts, simcoe_spec_2015_w_sp1_ts, by = "species_1")

simcoe_spec_2015_w_sp2_ts <- simcoe_spec_2015_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_2015_sp2 = w_pop_cv_2015) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_2015_sp2, Species_ID, station_2)

simcoe_2015_cor_cv_ts <- merge(simcoe_2015_cor_cv_ts, simcoe_spec_2015_w_sp2_ts, by = "species_2")

simcoe_2015_cor_cv_ts_omit <- na.omit(simcoe_2015_cor_cv_ts)

simcoe_2015_ind_W_ts <- simcoe_2015_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_2015 = w_pop_cv_2015_sp1*w_pop_cv_2015_sp2,
         ind_W_2015 = (1 - corr)*(w_pop_cv_2015_sp1*w_pop_cv_2015_sp2),
         Year = 2015)

simcoe_2015_ind_W_ts$number_species_pairs <- nrow(simcoe_2015_ind_W_ts)

simcoe_mean_W_sp_pairs_2015 <- simcoe_2015_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_2015, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_2015_W_ts <- simcoe_2015_ind_W_ts %>% group_by(Year) %>%
  summarize(W_2015 = sum(ind_W_2015, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_2015_ind_W_ts_c9 <- filter(simcoe_2015_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_2015_ind_W_ts_k42 <- filter(simcoe_2015_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_2015_ind_W_ts_k45 <- filter(simcoe_2015_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_2015_c9 <- na.omit(simcoe_2015_ind_W_ts_c9)

simcoe_species_pairs_count_2015_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_2015_c9)

simcoe_mean_W_sp_pairs_2015_c9 <- simcoe_species_pairs_count_2015_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_2015, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2015_W_ts_c9 <- simcoe_2015_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_2015 = sum(ind_W_2015, na.rm = TRUE))

#K42
simcoe_species_pairs_count_2015_k42 <- na.omit(simcoe_2015_ind_W_ts_k42)

simcoe_species_pairs_count_2015_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_2015_k42)

simcoe_mean_W_sp_pairs_2015_k42 <- simcoe_species_pairs_count_2015_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_2015, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2015_W_ts_k42 <- simcoe_2015_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_2015 = sum(ind_W_2015, na.rm = TRUE))

#K45
simcoe_species_pairs_count_2015_k45 <- na.omit(simcoe_2015_ind_W_ts_k45)

simcoe_species_pairs_count_2015_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_2015_k45)

simcoe_mean_W_sp_pairs_2015_k45 <- simcoe_species_pairs_count_2015_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_2015, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2015_W_ts_k45 <- simcoe_2015_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_2015 = sum(ind_W_2015, na.rm = TRUE))

#Step 4
simcoe_2015_ind_W_as1 <- rbind(simcoe_2015_W_ts_c9, simcoe_2015_W_ts_k42, simcoe_2015_W_ts_k45)

simcoe_2015_W_as1 <- simcoe_2015_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_2015 = sum(W_2015, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_2015 <- rbind(simcoe_mean_W_sp_pairs_2015_c9, simcoe_mean_W_sp_pairs_2015_k42, simcoe_mean_W_sp_pairs_2015_k45)

simcoe_lc_sp_pairs_2015 <- simcoe_mean_lc_sp_pairs_2015 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_2015_ind_W_as4 <- filter(simcoe_2015_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_2015_mp  <- na.omit(simcoe_2015_ind_W_as4)

simcoe_species_pairs_count_2015_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_2015_mp)

simcoe_mean_mp_sp_pairs_2015  <- simcoe_species_pairs_count_2015_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_2015, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2015_sum_W_as4 <- simcoe_2015_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_2015 = sum(ind_W_2015, na.rm = TRUE))


#Step 3
simcoe_2015_sum_W_as4 <- simcoe_2015_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_2015 = sum(as4_2015, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_2015_ind_W_as5 <- filter(simcoe_2015_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_2015_cc <- na.omit(simcoe_2015_ind_W_as5)

simcoe_species_pairs_count_2015_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_2015_cc)

simcoe_mean_cc_sp_pairs_2015  <- simcoe_species_pairs_count_2015_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_2015, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2015_sum_W_as5 <- simcoe_2015_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_2015 = sum(ind_W_2015, na.rm = TRUE))


#Step 3
simcoe_2015_sum_W_as5 <- simcoe_2015_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_2015 = sum(as5_2015, na.rm = TRUE))

simcoe_2015_W <- merge(simcoe_2015_sum_W_as4, simcoe_2015_W_ts, by = c('Year'))
simcoe_2015_W <- merge(simcoe_2015_W, simcoe_2015_W_as1, by = c('Year'))
simcoe_2015_W <- merge(simcoe_2015_W, simcoe_2015_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_2015 <- merge(simcoe_mean_W_sp_pairs_2015, simcoe_lc_sp_pairs_2015, by=c("Year"))
simcoe_sp_pair_ave_stab_2015 <- merge(simcoe_sp_pair_ave_stab_2015, simcoe_mean_mp_sp_pairs_2015, by=c("Year"))
simcoe_sp_pair_ave_stab_2015 <- merge(simcoe_sp_pair_ave_stab_2015, simcoe_mean_cc_sp_pairs_2015, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_2015_final <- merge(simcoe_2015_W, simcoe_2015_wa_pop_var_2015, by = c('Year'))
simcoe_2015_final <- simcoe_2015_final %>% group_by(Year) %>%
  mutate(gcv_2015 = lcv - W_2015) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_2015, 
            lcv = lcv,
            W = W_2015,
            lc_stab = as2_2015,
            mp_stab = as4_2015,
            cc_stab = as5_2015,
            mc_cv = sqrt(gcv_2015),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_2015/lcv + as4_2015/lcv + as5_2015/lcv,
            lc_asynchrony = as2_2015/lcv,
            mp_asynchrony = as4_2015/lcv,
            cc_asynchrony= as5_2015/lcv)

simcoe_2015_cv <- simcoe_2015_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_2015_final$mc_sum_mean_density <- simcoe_2015_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 2015 #####
simcoe_mean_corr_2015_c9 <- simcoe_species_pairs_count_2015_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2015, ind_W_2015)
simcoe_mean_corr_2015_k42 <- simcoe_species_pairs_count_2015_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2015, ind_W_2015)
simcoe_mean_corr_2015_k45 <- simcoe_species_pairs_count_2015_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2015, ind_W_2015)
simcoe_mean_corr_2015_lc <- rbind(simcoe_mean_corr_2015_c9, simcoe_mean_corr_2015_k42, simcoe_mean_corr_2015_k45)
simcoe_mean_corr_2015_mp <- simcoe_species_pairs_count_2015_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_2015, ind_W_2015)
simcoe_mean_corr_2015_cc <- simcoe_species_pairs_count_2015_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_2015, ind_W_2015)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_2015_lc[sample(1:nrow(simcoe_mean_corr_2015_lc),  nrow(simcoe_mean_corr_2015_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_2015),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_2015),
              mean_lc_stab = mean(ind_W_2015),
              sd_lc_stab = sd(ind_W_2015))
  
  lc_corr_2015 <- as.data.frame(mean_lc)
}
lc_corr_2015

mp_corr_2015 <- simcoe_mean_corr_2015_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_2015),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_2015),
                                                       mean_mp_stab = mean(ind_W_2015),
                                                       sd_mp_stab = sd(ind_W_2015))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_2015_cc[sample(1:nrow(simcoe_mean_corr_2015_cc), nrow(simcoe_mean_corr_2015_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_2015),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_2015),
              mean_cc_stab = mean(ind_W_2015),
              sd_cc_stab = sd(ind_W_2015))
  
  cc_corr_2015 <- as.data.frame(mean_cc)
}
cc_corr_2015


corr_2015 <- cbind(lc_corr_2015, mp_corr_2015)
corr_2015 <- cbind(corr_2015, cc_corr_2015)
corr_2015 <- corr_2015 %>% mutate(Year = 2015)




##### 2016 #####
simcoe_2016 <- subset(simcoe, Year == '2016')


#Transpose dataframe to add 0s within years
simcoe_2016_t1 <- dcast(simcoe_2016, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2016_m <- melt(simcoe_2016_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_2016 <- subset(simcoe_2016_m, Station_ID == 'C9')
simcoe_c9_2016$Species <- paste("C9", simcoe_c9_2016$Species, sep=" ")
simcoe_k42_2016 <- subset(simcoe_2016_m, Station_ID == 'K42')
simcoe_k42_2016$Species <- paste("K42", simcoe_k42_2016$Species, sep=" ")
simcoe_k45_2016 <- subset(simcoe_2016_m, Station_ID == 'K45')
simcoe_k45_2016$Species <- paste("K45", simcoe_k45_2016$Species, sep=" ")

#recombine dataframes
simcoe_spec_2016 <- rbind(simcoe_c9_2016, simcoe_k42_2016, simcoe_k45_2016)
simcoe_spec_2016 <- simcoe_spec_2016 %>% select(- Station_ID)

#transpose dataframe
simcoe_2016_t2 <- dcast(simcoe_spec_2016, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2016_t2 <- log(simcoe_2016_t2 + 1)
simcoe_2016_t2 <- simcoe_2016_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_2016_cv <- simcoe_2016_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_2016 = mean(Density),
            pop_sd_2016 = sd(Density)) %>%
  mutate(uw_pop_cv_2016 = pop_sd_2016/pop_mean_2016)

simcoe_2016_cv <- simcoe_2016_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_2016),
         mean_pop_cv = mean(uw_pop_cv_2016, na.rm = T),
         mean_pop_density = sum(pop_mean_2016),
         mean_pop_variance = sum(pop_sd_2016))

simcoe_2016_cv <- simcoe_2016_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_2016/mc_sum_mean_density)%>%
  mutate(w_pop_cv_2016 = (pop_mean_2016/mc_sum_mean_density)*uw_pop_cv_2016) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_2016_wa_pop_var_2016 <- simcoe_2016_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_2016 = sum(w_pop_cv_2016, na.rm = T),
            lcv = sum(w_pop_cv_2016, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_2016_cor_list_ts <- as.dist(round(cor(simcoe_2016_t2[]),2))
simcoe_2016_cor_ts <- stack(simcoe_2016_cor_list_ts, dim.names = TRUE)
simcoe_2016_cor_ts <- simcoe_2016_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_2016_cor - simcoe_2016_cv
simcoe_c9_2016_w_cv_ts <- subset(simcoe_2016_cv, Station_ID == 'C9')
simcoe_c9_2016_w_cv_ts$Species <- paste("C9", simcoe_c9_2016_w_cv_ts$Species, sep=" ")
simcoe_c9_2016_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_2016_w_cv_ts)
simcoe_k42_2016_w_cv_ts <- subset(simcoe_2016_cv, Station_ID == 'K42')
simcoe_k42_2016_w_cv_ts$Species <- paste("K42", simcoe_k42_2016_w_cv_ts$Species, sep=" ")
simcoe_k42_2016_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_2016_w_cv_ts)
simcoe_k45_2016_w_cv_ts <- subset(simcoe_2016_cv, Station_ID == 'K45')
simcoe_k45_2016_w_cv_ts$Species <- paste("K45", simcoe_k45_2016_w_cv_ts$Species, sep=" ")
simcoe_k45_2016_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_2016_w_cv_ts)

simcoe_spec_2016_w_sp_ts <- rbind(simcoe_c9_2016_w_cv_ts, simcoe_k42_2016_w_cv_ts, simcoe_k45_2016_w_cv_ts)
simcoe_spec_2016_w_sp_ts <- simcoe_spec_2016_w_sp_ts %>% select(Species, w_pop_cv_2016, Species_ID, Station_ID)
simcoe_spec_2016_w_sp1_ts <- simcoe_spec_2016_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_2016_sp1 = w_pop_cv_2016) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_2016_sp1, Species_ID, station_1)

simcoe_2016_cor_cv_ts <- merge(simcoe_2016_cor_ts, simcoe_spec_2016_w_sp1_ts, by = "species_1")

simcoe_spec_2016_w_sp2_ts <- simcoe_spec_2016_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_2016_sp2 = w_pop_cv_2016) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_2016_sp2, Species_ID, station_2)

simcoe_2016_cor_cv_ts <- merge(simcoe_2016_cor_cv_ts, simcoe_spec_2016_w_sp2_ts, by = "species_2")

simcoe_2016_cor_cv_ts_omit <- na.omit(simcoe_2016_cor_cv_ts)

simcoe_2016_ind_W_ts <- simcoe_2016_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_2016 = w_pop_cv_2016_sp1*w_pop_cv_2016_sp2,
         ind_W_2016 = (1 - corr)*(w_pop_cv_2016_sp1*w_pop_cv_2016_sp2),
         Year = 2016)

simcoe_2016_ind_W_ts$number_species_pairs <- nrow(simcoe_2016_ind_W_ts)

simcoe_mean_W_sp_pairs_2016 <- simcoe_2016_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_2016, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_2016_W_ts <- simcoe_2016_ind_W_ts %>% group_by(Year) %>%
  summarize(W_2016 = sum(ind_W_2016, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_2016_ind_W_ts_c9 <- filter(simcoe_2016_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_2016_ind_W_ts_k42 <- filter(simcoe_2016_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_2016_ind_W_ts_k45 <- filter(simcoe_2016_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_2016_c9 <- na.omit(simcoe_2016_ind_W_ts_c9)

simcoe_species_pairs_count_2016_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_2016_c9)

simcoe_mean_W_sp_pairs_2016_c9 <- simcoe_species_pairs_count_2016_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_2016, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2016_W_ts_c9 <- simcoe_2016_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_2016 = sum(ind_W_2016, na.rm = TRUE))

#K42
simcoe_species_pairs_count_2016_k42 <- na.omit(simcoe_2016_ind_W_ts_k42)

simcoe_species_pairs_count_2016_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_2016_k42)

simcoe_mean_W_sp_pairs_2016_k42 <- simcoe_species_pairs_count_2016_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_2016, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2016_W_ts_k42 <- simcoe_2016_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_2016 = sum(ind_W_2016, na.rm = TRUE))

#K45
simcoe_species_pairs_count_2016_k45 <- na.omit(simcoe_2016_ind_W_ts_k45)

simcoe_species_pairs_count_2016_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_2016_k45)

simcoe_mean_W_sp_pairs_2016_k45 <- simcoe_species_pairs_count_2016_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_2016, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2016_W_ts_k45 <- simcoe_2016_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_2016 = sum(ind_W_2016, na.rm = TRUE))

#Step 4
simcoe_2016_ind_W_as1 <- rbind(simcoe_2016_W_ts_c9, simcoe_2016_W_ts_k42, simcoe_2016_W_ts_k45)

simcoe_2016_W_as1 <- simcoe_2016_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_2016 = sum(W_2016, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_2016 <- rbind(simcoe_mean_W_sp_pairs_2016_c9, simcoe_mean_W_sp_pairs_2016_k42, simcoe_mean_W_sp_pairs_2016_k45)

simcoe_lc_sp_pairs_2016 <- simcoe_mean_lc_sp_pairs_2016 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_2016_ind_W_as4 <- filter(simcoe_2016_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_2016_mp  <- na.omit(simcoe_2016_ind_W_as4)

simcoe_species_pairs_count_2016_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_2016_mp)

simcoe_mean_mp_sp_pairs_2016  <- simcoe_species_pairs_count_2016_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_2016, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2016_sum_W_as4 <- simcoe_2016_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_2016 = sum(ind_W_2016, na.rm = TRUE))


#Step 3
simcoe_2016_sum_W_as4 <- simcoe_2016_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_2016 = sum(as4_2016, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_2016_ind_W_as5 <- filter(simcoe_2016_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_2016_cc <- na.omit(simcoe_2016_ind_W_as5)

simcoe_species_pairs_count_2016_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_2016_cc)

simcoe_mean_cc_sp_pairs_2016  <- simcoe_species_pairs_count_2016_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_2016, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2016_sum_W_as5 <- simcoe_2016_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_2016 = sum(ind_W_2016, na.rm = TRUE))


#Step 3
simcoe_2016_sum_W_as5 <- simcoe_2016_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_2016 = sum(as5_2016, na.rm = TRUE))

simcoe_2016_W <- merge(simcoe_2016_sum_W_as4, simcoe_2016_W_ts, by = c('Year'))
simcoe_2016_W <- merge(simcoe_2016_W, simcoe_2016_W_as1, by = c('Year'))
simcoe_2016_W <- merge(simcoe_2016_W, simcoe_2016_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_2016 <- merge(simcoe_mean_W_sp_pairs_2016, simcoe_lc_sp_pairs_2016, by=c("Year"))
simcoe_sp_pair_ave_stab_2016 <- merge(simcoe_sp_pair_ave_stab_2016, simcoe_mean_mp_sp_pairs_2016, by=c("Year"))
simcoe_sp_pair_ave_stab_2016 <- merge(simcoe_sp_pair_ave_stab_2016, simcoe_mean_cc_sp_pairs_2016, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_2016_final <- merge(simcoe_2016_W, simcoe_2016_wa_pop_var_2016, by = c('Year'))
simcoe_2016_final <- simcoe_2016_final %>% group_by(Year) %>%
  mutate(gcv_2016 = lcv - W_2016) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_2016, 
            lcv = lcv,
            W = W_2016,
            lc_stab = as2_2016,
            mp_stab = as4_2016,
            cc_stab = as5_2016,
            mc_cv = sqrt(gcv_2016),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_2016/lcv + as4_2016/lcv + as5_2016/lcv,
            lc_asynchrony = as2_2016/lcv,
            mp_asynchrony = as4_2016/lcv,
            cc_asynchrony= as5_2016/lcv)

simcoe_2016_cv <- simcoe_2016_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_2016_final$mc_sum_mean_density <- simcoe_2016_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 2016 #####
simcoe_mean_corr_2016_c9 <- simcoe_species_pairs_count_2016_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2016, ind_W_2016)
simcoe_mean_corr_2016_k42 <- simcoe_species_pairs_count_2016_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2016, ind_W_2016)
simcoe_mean_corr_2016_k45 <- simcoe_species_pairs_count_2016_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2016, ind_W_2016)
simcoe_mean_corr_2016_lc <- rbind(simcoe_mean_corr_2016_c9, simcoe_mean_corr_2016_k42, simcoe_mean_corr_2016_k45)
simcoe_mean_corr_2016_mp <- simcoe_species_pairs_count_2016_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_2016, ind_W_2016)
simcoe_mean_corr_2016_cc <- simcoe_species_pairs_count_2016_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_2016, ind_W_2016)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_2016_lc[sample(1:nrow(simcoe_mean_corr_2016_lc),  nrow(simcoe_mean_corr_2016_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_2016),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_2016),
              mean_lc_stab = mean(ind_W_2016),
              sd_lc_stab = sd(ind_W_2016))
  
  lc_corr_2016 <- as.data.frame(mean_lc)
}
lc_corr_2016

mp_corr_2016 <- simcoe_mean_corr_2016_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_2016),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_2016),
                                                       mean_mp_stab = mean(ind_W_2016),
                                                       sd_mp_stab = sd(ind_W_2016))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_2016_cc[sample(1:nrow(simcoe_mean_corr_2016_cc), nrow(simcoe_mean_corr_2016_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_2016),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_2016),
              mean_cc_stab = mean(ind_W_2016),
              sd_cc_stab = sd(ind_W_2016))
  
  cc_corr_2016 <- as.data.frame(mean_cc)
}
cc_corr_2016


corr_2016 <- cbind(lc_corr_2016, mp_corr_2016)
corr_2016 <- cbind(corr_2016, cc_corr_2016)
corr_2016 <- corr_2016 %>% mutate(Year = 2016)





##### 2017 #####
simcoe_2017 <- subset(simcoe, Year == '2017')


#Transpose dataframe to add 0s within years
simcoe_2017_t1 <- dcast(simcoe_2017, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2017_m <- melt(simcoe_2017_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_2017 <- subset(simcoe_2017_m, Station_ID == 'C9')
simcoe_c9_2017$Species <- paste("C9", simcoe_c9_2017$Species, sep=" ")
simcoe_k42_2017 <- subset(simcoe_2017_m, Station_ID == 'K42')
simcoe_k42_2017$Species <- paste("K42", simcoe_k42_2017$Species, sep=" ")
simcoe_k45_2017 <- subset(simcoe_2017_m, Station_ID == 'K45')
simcoe_k45_2017$Species <- paste("K45", simcoe_k45_2017$Species, sep=" ")

#recombine dataframes
simcoe_spec_2017 <- rbind(simcoe_c9_2017, simcoe_k42_2017, simcoe_k45_2017)
simcoe_spec_2017 <- simcoe_spec_2017 %>% select(- Station_ID)

#transpose dataframe
simcoe_2017_t2 <- dcast(simcoe_spec_2017, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2017_t2 <- log(simcoe_2017_t2 + 1)
simcoe_2017_t2 <- simcoe_2017_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_2017_cv <- simcoe_2017_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_2017 = mean(Density),
            pop_sd_2017 = sd(Density)) %>%
  mutate(uw_pop_cv_2017 = pop_sd_2017/pop_mean_2017)

simcoe_2017_cv <- simcoe_2017_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_2017),
         mean_pop_cv = mean(uw_pop_cv_2017, na.rm = T),
         mean_pop_density = sum(pop_mean_2017),
         mean_pop_variance = sum(pop_sd_2017))

simcoe_2017_cv <- simcoe_2017_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_2017/mc_sum_mean_density)%>%
  mutate(w_pop_cv_2017 = (pop_mean_2017/mc_sum_mean_density)*uw_pop_cv_2017) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_2017_wa_pop_var_2017 <- simcoe_2017_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_2017 = sum(w_pop_cv_2017, na.rm = T),
            lcv = sum(w_pop_cv_2017, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_2017_cor_list_ts <- as.dist(round(cor(simcoe_2017_t2[]),2))
simcoe_2017_cor_ts <- stack(simcoe_2017_cor_list_ts, dim.names = TRUE)
simcoe_2017_cor_ts <- simcoe_2017_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_2017_cor - simcoe_2017_cv
simcoe_c9_2017_w_cv_ts <- subset(simcoe_2017_cv, Station_ID == 'C9')
simcoe_c9_2017_w_cv_ts$Species <- paste("C9", simcoe_c9_2017_w_cv_ts$Species, sep=" ")
simcoe_c9_2017_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_2017_w_cv_ts)
simcoe_k42_2017_w_cv_ts <- subset(simcoe_2017_cv, Station_ID == 'K42')
simcoe_k42_2017_w_cv_ts$Species <- paste("K42", simcoe_k42_2017_w_cv_ts$Species, sep=" ")
simcoe_k42_2017_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_2017_w_cv_ts)
simcoe_k45_2017_w_cv_ts <- subset(simcoe_2017_cv, Station_ID == 'K45')
simcoe_k45_2017_w_cv_ts$Species <- paste("K45", simcoe_k45_2017_w_cv_ts$Species, sep=" ")
simcoe_k45_2017_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_2017_w_cv_ts)

simcoe_spec_2017_w_sp_ts <- rbind(simcoe_c9_2017_w_cv_ts, simcoe_k42_2017_w_cv_ts, simcoe_k45_2017_w_cv_ts)
simcoe_spec_2017_w_sp_ts <- simcoe_spec_2017_w_sp_ts %>% select(Species, w_pop_cv_2017, Species_ID, Station_ID)
simcoe_spec_2017_w_sp1_ts <- simcoe_spec_2017_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_2017_sp1 = w_pop_cv_2017) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_2017_sp1, Species_ID, station_1)

simcoe_2017_cor_cv_ts <- merge(simcoe_2017_cor_ts, simcoe_spec_2017_w_sp1_ts, by = "species_1")

simcoe_spec_2017_w_sp2_ts <- simcoe_spec_2017_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_2017_sp2 = w_pop_cv_2017) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_2017_sp2, Species_ID, station_2)

simcoe_2017_cor_cv_ts <- merge(simcoe_2017_cor_cv_ts, simcoe_spec_2017_w_sp2_ts, by = "species_2")

simcoe_2017_cor_cv_ts_omit <- na.omit(simcoe_2017_cor_cv_ts)

simcoe_2017_ind_W_ts <- simcoe_2017_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_2017 = w_pop_cv_2017_sp1*w_pop_cv_2017_sp2,
         ind_W_2017 = (1 - corr)*(w_pop_cv_2017_sp1*w_pop_cv_2017_sp2),
         Year = 2017)

simcoe_2017_ind_W_ts$number_species_pairs <- nrow(simcoe_2017_ind_W_ts)

simcoe_mean_W_sp_pairs_2017 <- simcoe_2017_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_2017, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_2017_W_ts <- simcoe_2017_ind_W_ts %>% group_by(Year) %>%
  summarize(W_2017 = sum(ind_W_2017, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_2017_ind_W_ts_c9 <- filter(simcoe_2017_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_2017_ind_W_ts_k42 <- filter(simcoe_2017_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_2017_ind_W_ts_k45 <- filter(simcoe_2017_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_2017_c9 <- na.omit(simcoe_2017_ind_W_ts_c9)

simcoe_species_pairs_count_2017_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_2017_c9)

simcoe_mean_W_sp_pairs_2017_c9 <- simcoe_species_pairs_count_2017_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_2017, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2017_W_ts_c9 <- simcoe_2017_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_2017 = sum(ind_W_2017, na.rm = TRUE))

#K42
simcoe_species_pairs_count_2017_k42 <- na.omit(simcoe_2017_ind_W_ts_k42)

simcoe_species_pairs_count_2017_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_2017_k42)

simcoe_mean_W_sp_pairs_2017_k42 <- simcoe_species_pairs_count_2017_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_2017, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2017_W_ts_k42 <- simcoe_2017_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_2017 = sum(ind_W_2017, na.rm = TRUE))

#K45
simcoe_species_pairs_count_2017_k45 <- na.omit(simcoe_2017_ind_W_ts_k45)

simcoe_species_pairs_count_2017_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_2017_k45)

simcoe_mean_W_sp_pairs_2017_k45 <- simcoe_species_pairs_count_2017_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_2017, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2017_W_ts_k45 <- simcoe_2017_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_2017 = sum(ind_W_2017, na.rm = TRUE))

#Step 4
simcoe_2017_ind_W_as1 <- rbind(simcoe_2017_W_ts_c9, simcoe_2017_W_ts_k42, simcoe_2017_W_ts_k45)

simcoe_2017_W_as1 <- simcoe_2017_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_2017 = sum(W_2017, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_2017 <- rbind(simcoe_mean_W_sp_pairs_2017_c9, simcoe_mean_W_sp_pairs_2017_k42, simcoe_mean_W_sp_pairs_2017_k45)

simcoe_lc_sp_pairs_2017 <- simcoe_mean_lc_sp_pairs_2017 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_2017_ind_W_as4 <- filter(simcoe_2017_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_2017_mp  <- na.omit(simcoe_2017_ind_W_as4)

simcoe_species_pairs_count_2017_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_2017_mp)

simcoe_mean_mp_sp_pairs_2017  <- simcoe_species_pairs_count_2017_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_2017, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2017_sum_W_as4 <- simcoe_2017_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_2017 = sum(ind_W_2017, na.rm = TRUE))


#Step 3
simcoe_2017_sum_W_as4 <- simcoe_2017_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_2017 = sum(as4_2017, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_2017_ind_W_as5 <- filter(simcoe_2017_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_2017_cc <- na.omit(simcoe_2017_ind_W_as5)

simcoe_species_pairs_count_2017_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_2017_cc)

simcoe_mean_cc_sp_pairs_2017  <- simcoe_species_pairs_count_2017_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_2017, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2017_sum_W_as5 <- simcoe_2017_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_2017 = sum(ind_W_2017, na.rm = TRUE))


#Step 3
simcoe_2017_sum_W_as5 <- simcoe_2017_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_2017 = sum(as5_2017, na.rm = TRUE))

simcoe_2017_W <- merge(simcoe_2017_sum_W_as4, simcoe_2017_W_ts, by = c('Year'))
simcoe_2017_W <- merge(simcoe_2017_W, simcoe_2017_W_as1, by = c('Year'))
simcoe_2017_W <- merge(simcoe_2017_W, simcoe_2017_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_2017 <- merge(simcoe_mean_W_sp_pairs_2017, simcoe_lc_sp_pairs_2017, by=c("Year"))
simcoe_sp_pair_ave_stab_2017 <- merge(simcoe_sp_pair_ave_stab_2017, simcoe_mean_mp_sp_pairs_2017, by=c("Year"))
simcoe_sp_pair_ave_stab_2017 <- merge(simcoe_sp_pair_ave_stab_2017, simcoe_mean_cc_sp_pairs_2017, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_2017_final <- merge(simcoe_2017_W, simcoe_2017_wa_pop_var_2017, by = c('Year'))
simcoe_2017_final <- simcoe_2017_final %>% group_by(Year) %>%
  mutate(gcv_2017 = lcv - W_2017) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_2017, 
            lcv = lcv,
            W = W_2017,
            lc_stab = as2_2017,
            mp_stab = as4_2017,
            cc_stab = as5_2017,
            mc_cv = sqrt(gcv_2017),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_2017/lcv + as4_2017/lcv + as5_2017/lcv,
            lc_asynchrony = as2_2017/lcv,
            mp_asynchrony = as4_2017/lcv,
            cc_asynchrony= as5_2017/lcv)

simcoe_2017_cv <- simcoe_2017_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_2017_final$mc_sum_mean_density <- simcoe_2017_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 2017 #####
simcoe_mean_corr_2017_c9 <- simcoe_species_pairs_count_2017_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2017, ind_W_2017)
simcoe_mean_corr_2017_k42 <- simcoe_species_pairs_count_2017_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2017, ind_W_2017)
simcoe_mean_corr_2017_k45 <- simcoe_species_pairs_count_2017_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2017, ind_W_2017)
simcoe_mean_corr_2017_lc <- rbind(simcoe_mean_corr_2017_c9, simcoe_mean_corr_2017_k42, simcoe_mean_corr_2017_k45)
simcoe_mean_corr_2017_mp <- simcoe_species_pairs_count_2017_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_2017, ind_W_2017)
simcoe_mean_corr_2017_cc <- simcoe_species_pairs_count_2017_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_2017, ind_W_2017)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_2017_lc[sample(1:nrow(simcoe_mean_corr_2017_lc),  nrow(simcoe_mean_corr_2017_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_2017),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_2017),
              mean_lc_stab = mean(ind_W_2017),
              sd_lc_stab = sd(ind_W_2017))
  
  lc_corr_2017 <- as.data.frame(mean_lc)
}
lc_corr_2017

mp_corr_2017 <- simcoe_mean_corr_2017_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_2017),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_2017),
                                                       mean_mp_stab = mean(ind_W_2017),
                                                       sd_mp_stab = sd(ind_W_2017))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_2017_cc[sample(1:nrow(simcoe_mean_corr_2017_cc), nrow(simcoe_mean_corr_2017_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_2017),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_2017),
              mean_cc_stab = mean(ind_W_2017),
              sd_cc_stab = sd(ind_W_2017))
  
  cc_corr_2017 <- as.data.frame(mean_cc)
}
cc_corr_2017


corr_2017 <- cbind(lc_corr_2017, mp_corr_2017)
corr_2017 <- cbind(corr_2017, cc_corr_2017)
corr_2017 <- corr_2017 %>% mutate(Year = 2017)






##### 2018 #####
simcoe_2018 <- subset(simcoe, Year == '2018')


#Transpose dataframe to add 0s within years
simcoe_2018_t1 <- dcast(simcoe_2018, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2018_m <- melt(simcoe_2018_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_2018 <- subset(simcoe_2018_m, Station_ID == 'C9')
simcoe_c9_2018$Species <- paste("C9", simcoe_c9_2018$Species, sep=" ")
simcoe_k42_2018 <- subset(simcoe_2018_m, Station_ID == 'K42')
simcoe_k42_2018$Species <- paste("K42", simcoe_k42_2018$Species, sep=" ")
simcoe_k45_2018 <- subset(simcoe_2018_m, Station_ID == 'K45')
simcoe_k45_2018$Species <- paste("K45", simcoe_k45_2018$Species, sep=" ")

#recombine dataframes
simcoe_spec_2018 <- rbind(simcoe_c9_2018, simcoe_k42_2018, simcoe_k45_2018)
simcoe_spec_2018 <- simcoe_spec_2018 %>% select(- Station_ID)

#transpose dataframe
simcoe_2018_t2 <- dcast(simcoe_spec_2018, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2018_t2 <- log(simcoe_2018_t2 + 1)
simcoe_2018_t2 <- simcoe_2018_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_2018_cv <- simcoe_2018_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_2018 = mean(Density),
            pop_sd_2018 = sd(Density)) %>%
  mutate(uw_pop_cv_2018 = pop_sd_2018/pop_mean_2018)

simcoe_2018_cv <- simcoe_2018_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_2018),
         mean_pop_cv = mean(uw_pop_cv_2018, na.rm = T),
         mean_pop_density = sum(pop_mean_2018),
         mean_pop_variance = sum(pop_sd_2018))

simcoe_2018_cv <- simcoe_2018_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_2018/mc_sum_mean_density)%>%
  mutate(w_pop_cv_2018 = (pop_mean_2018/mc_sum_mean_density)*uw_pop_cv_2018) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_2018_wa_pop_var_2018 <- simcoe_2018_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_2018 = sum(w_pop_cv_2018, na.rm = T),
            lcv = sum(w_pop_cv_2018, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_2018_cor_list_ts <- as.dist(round(cor(simcoe_2018_t2[]),2))
simcoe_2018_cor_ts <- stack(simcoe_2018_cor_list_ts, dim.names = TRUE)
simcoe_2018_cor_ts <- simcoe_2018_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_2018_cor - simcoe_2018_cv
simcoe_c9_2018_w_cv_ts <- subset(simcoe_2018_cv, Station_ID == 'C9')
simcoe_c9_2018_w_cv_ts$Species <- paste("C9", simcoe_c9_2018_w_cv_ts$Species, sep=" ")
simcoe_c9_2018_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_2018_w_cv_ts)
simcoe_k42_2018_w_cv_ts <- subset(simcoe_2018_cv, Station_ID == 'K42')
simcoe_k42_2018_w_cv_ts$Species <- paste("K42", simcoe_k42_2018_w_cv_ts$Species, sep=" ")
simcoe_k42_2018_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_2018_w_cv_ts)
simcoe_k45_2018_w_cv_ts <- subset(simcoe_2018_cv, Station_ID == 'K45')
simcoe_k45_2018_w_cv_ts$Species <- paste("K45", simcoe_k45_2018_w_cv_ts$Species, sep=" ")
simcoe_k45_2018_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_2018_w_cv_ts)

simcoe_spec_2018_w_sp_ts <- rbind(simcoe_c9_2018_w_cv_ts, simcoe_k42_2018_w_cv_ts, simcoe_k45_2018_w_cv_ts)
simcoe_spec_2018_w_sp_ts <- simcoe_spec_2018_w_sp_ts %>% select(Species, w_pop_cv_2018, Species_ID, Station_ID)
simcoe_spec_2018_w_sp1_ts <- simcoe_spec_2018_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_2018_sp1 = w_pop_cv_2018) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_2018_sp1, Species_ID, station_1)

simcoe_2018_cor_cv_ts <- merge(simcoe_2018_cor_ts, simcoe_spec_2018_w_sp1_ts, by = "species_1")

simcoe_spec_2018_w_sp2_ts <- simcoe_spec_2018_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_2018_sp2 = w_pop_cv_2018) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_2018_sp2, Species_ID, station_2)

simcoe_2018_cor_cv_ts <- merge(simcoe_2018_cor_cv_ts, simcoe_spec_2018_w_sp2_ts, by = "species_2")

simcoe_2018_cor_cv_ts_omit <- na.omit(simcoe_2018_cor_cv_ts)

simcoe_2018_ind_W_ts <- simcoe_2018_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_2018 = w_pop_cv_2018_sp1*w_pop_cv_2018_sp2,
         ind_W_2018 = (1 - corr)*(w_pop_cv_2018_sp1*w_pop_cv_2018_sp2),
         Year = 2018)

simcoe_2018_ind_W_ts$number_species_pairs <- nrow(simcoe_2018_ind_W_ts)

simcoe_mean_W_sp_pairs_2018 <- simcoe_2018_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_2018, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_2018_W_ts <- simcoe_2018_ind_W_ts %>% group_by(Year) %>%
  summarize(W_2018 = sum(ind_W_2018, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_2018_ind_W_ts_c9 <- filter(simcoe_2018_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_2018_ind_W_ts_k42 <- filter(simcoe_2018_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_2018_ind_W_ts_k45 <- filter(simcoe_2018_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_2018_c9 <- na.omit(simcoe_2018_ind_W_ts_c9)

simcoe_species_pairs_count_2018_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_2018_c9)

simcoe_mean_W_sp_pairs_2018_c9 <- simcoe_species_pairs_count_2018_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_2018, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2018_W_ts_c9 <- simcoe_2018_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_2018 = sum(ind_W_2018, na.rm = TRUE))

#K42
simcoe_species_pairs_count_2018_k42 <- na.omit(simcoe_2018_ind_W_ts_k42)

simcoe_species_pairs_count_2018_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_2018_k42)

simcoe_mean_W_sp_pairs_2018_k42 <- simcoe_species_pairs_count_2018_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_2018, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2018_W_ts_k42 <- simcoe_2018_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_2018 = sum(ind_W_2018, na.rm = TRUE))

#K45
simcoe_species_pairs_count_2018_k45 <- na.omit(simcoe_2018_ind_W_ts_k45)

simcoe_species_pairs_count_2018_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_2018_k45)

simcoe_mean_W_sp_pairs_2018_k45 <- simcoe_species_pairs_count_2018_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_2018, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2018_W_ts_k45 <- simcoe_2018_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_2018 = sum(ind_W_2018, na.rm = TRUE))

#Step 4
simcoe_2018_ind_W_as1 <- rbind(simcoe_2018_W_ts_c9, simcoe_2018_W_ts_k42, simcoe_2018_W_ts_k45)

simcoe_2018_W_as1 <- simcoe_2018_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_2018 = sum(W_2018, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_2018 <- rbind(simcoe_mean_W_sp_pairs_2018_c9, simcoe_mean_W_sp_pairs_2018_k42, simcoe_mean_W_sp_pairs_2018_k45)

simcoe_lc_sp_pairs_2018 <- simcoe_mean_lc_sp_pairs_2018 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_2018_ind_W_as4 <- filter(simcoe_2018_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_2018_mp  <- na.omit(simcoe_2018_ind_W_as4)

simcoe_species_pairs_count_2018_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_2018_mp)

simcoe_mean_mp_sp_pairs_2018  <- simcoe_species_pairs_count_2018_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_2018, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2018_sum_W_as4 <- simcoe_2018_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_2018 = sum(ind_W_2018, na.rm = TRUE))


#Step 3
simcoe_2018_sum_W_as4 <- simcoe_2018_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_2018 = sum(as4_2018, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_2018_ind_W_as5 <- filter(simcoe_2018_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_2018_cc <- na.omit(simcoe_2018_ind_W_as5)

simcoe_species_pairs_count_2018_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_2018_cc)

simcoe_mean_cc_sp_pairs_2018  <- simcoe_species_pairs_count_2018_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_2018, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2018_sum_W_as5 <- simcoe_2018_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_2018 = sum(ind_W_2018, na.rm = TRUE))


#Step 3
simcoe_2018_sum_W_as5 <- simcoe_2018_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_2018 = sum(as5_2018, na.rm = TRUE))

simcoe_2018_W <- merge(simcoe_2018_sum_W_as4, simcoe_2018_W_ts, by = c('Year'))
simcoe_2018_W <- merge(simcoe_2018_W, simcoe_2018_W_as1, by = c('Year'))
simcoe_2018_W <- merge(simcoe_2018_W, simcoe_2018_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_2018 <- merge(simcoe_mean_W_sp_pairs_2018, simcoe_lc_sp_pairs_2018, by=c("Year"))
simcoe_sp_pair_ave_stab_2018 <- merge(simcoe_sp_pair_ave_stab_2018, simcoe_mean_mp_sp_pairs_2018, by=c("Year"))
simcoe_sp_pair_ave_stab_2018 <- merge(simcoe_sp_pair_ave_stab_2018, simcoe_mean_cc_sp_pairs_2018, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_2018_final <- merge(simcoe_2018_W, simcoe_2018_wa_pop_var_2018, by = c('Year'))
simcoe_2018_final <- simcoe_2018_final %>% group_by(Year) %>%
  mutate(gcv_2018 = lcv - W_2018) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_2018, 
            lcv = lcv,
            W = W_2018,
            lc_stab = as2_2018,
            mp_stab = as4_2018,
            cc_stab = as5_2018,
            mc_cv = sqrt(gcv_2018),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_2018/lcv + as4_2018/lcv + as5_2018/lcv,
            lc_asynchrony = as2_2018/lcv,
            mp_asynchrony = as4_2018/lcv,
            cc_asynchrony= as5_2018/lcv)

simcoe_2018_cv <- simcoe_2018_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_2018_final$mc_sum_mean_density <- simcoe_2018_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 2018 #####
simcoe_mean_corr_2018_c9 <- simcoe_species_pairs_count_2018_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2018, ind_W_2018)
simcoe_mean_corr_2018_k42 <- simcoe_species_pairs_count_2018_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2018, ind_W_2018)
simcoe_mean_corr_2018_k45 <- simcoe_species_pairs_count_2018_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2018, ind_W_2018)
simcoe_mean_corr_2018_lc <- rbind(simcoe_mean_corr_2018_c9, simcoe_mean_corr_2018_k42, simcoe_mean_corr_2018_k45)
simcoe_mean_corr_2018_mp <- simcoe_species_pairs_count_2018_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_2018, ind_W_2018)
simcoe_mean_corr_2018_cc <- simcoe_species_pairs_count_2018_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_2018, ind_W_2018)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_2018_lc[sample(1:nrow(simcoe_mean_corr_2018_lc),  nrow(simcoe_mean_corr_2018_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_2018),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_2018),
              mean_lc_stab = mean(ind_W_2018),
              sd_lc_stab = sd(ind_W_2018))
  
  lc_corr_2018 <- as.data.frame(mean_lc)
}
lc_corr_2018

mp_corr_2018 <- simcoe_mean_corr_2018_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_2018),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_2018),
                                                       mean_mp_stab = mean(ind_W_2018),
                                                       sd_mp_stab = sd(ind_W_2018))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_2018_cc[sample(1:nrow(simcoe_mean_corr_2018_cc), nrow(simcoe_mean_corr_2018_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_2018),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_2018),
              mean_cc_stab = mean(ind_W_2018),
              sd_cc_stab = sd(ind_W_2018))
  
  cc_corr_2018 <- as.data.frame(mean_cc)
}
cc_corr_2018


corr_2018 <- cbind(lc_corr_2018, mp_corr_2018)
corr_2018 <- cbind(corr_2018, cc_corr_2018)
corr_2018 <- corr_2018 %>% mutate(Year = 2018)




##### 2019 #####
simcoe_2019 <- subset(simcoe, Year == '2019')


#Transpose dataframe to add 0s within years
simcoe_2019_t1 <- dcast(simcoe_2019, Year + Week + Station_ID ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2019_m <- melt(simcoe_2019_t1, id.vars = c("Year", "Week", "Station_ID"),
                      variable.name = "Species",
                      value.name = "Density")

#subset dataframes by station to add a station tag to species id
simcoe_c9_2019 <- subset(simcoe_2019_m, Station_ID == 'C9')
simcoe_c9_2019$Species <- paste("C9", simcoe_c9_2019$Species, sep=" ")
simcoe_k42_2019 <- subset(simcoe_2019_m, Station_ID == 'K42')
simcoe_k42_2019$Species <- paste("K42", simcoe_k42_2019$Species, sep=" ")
simcoe_k45_2019 <- subset(simcoe_2019_m, Station_ID == 'K45')
simcoe_k45_2019$Species <- paste("K45", simcoe_k45_2019$Species, sep=" ")

#recombine dataframes
simcoe_spec_2019 <- rbind(simcoe_c9_2019, simcoe_k42_2019, simcoe_k45_2019)
simcoe_spec_2019 <- simcoe_spec_2019 %>% select(- Station_ID)

#transpose dataframe
simcoe_2019_t2 <- dcast(simcoe_spec_2019, Week ~ Species, fun.aggregate = sum, value.var = "Density")
simcoe_2019_t2 <- log(simcoe_2019_t2 + 1)
simcoe_2019_t2 <- simcoe_2019_t2[-1]

##### W - Total Stabilization - asynchrony among all local populations - Hammond et al 2020 #####
# Equation - W = Sum[(1 - between population Pearson correlation coefficient) * weighted CV population i @ site k * weighted CV population j @ site l]
#Steps
# 1 - Calculate weighted population CV for each population
# 1 A - Calculate weighted-average population variability - lcv = (sum(population mean/metacommunity mean)*population CV)^2
# 2 - Calculate population pair Pearson correlation coefficients
# 3 - Plug into equation for each population pair to calculate stabilization from each population pair
# 4 - Sum all to calculate Total Stabilization (W)


#Step 1
simcoe_2019_cv <- simcoe_2019_m %>% group_by(Year, Species, Station_ID)%>%
  summarize(pop_mean_2019 = mean(Density),
            pop_sd_2019 = sd(Density)) %>%
  mutate(uw_pop_cv_2019 = pop_sd_2019/pop_mean_2019)

simcoe_2019_cv <- simcoe_2019_cv %>% group_by(Year) %>%
  mutate(mc_sum_mean_density = sum(pop_mean_2019),
         mean_pop_cv = mean(uw_pop_cv_2019, na.rm = T),
         mean_pop_density = sum(pop_mean_2019),
         mean_pop_variance = sum(pop_sd_2019))

simcoe_2019_cv <- simcoe_2019_cv %>% group_by(Species, Station_ID)%>%
  mutate(weight = pop_mean_2019/mc_sum_mean_density)%>%
  mutate(w_pop_cv_2019 = (pop_mean_2019/mc_sum_mean_density)*uw_pop_cv_2019) %>%
  group_by(Year) %>%
  mutate(sum_weight = sum(weight))

### Step 1 A - lcv - Weighted average population variability ###
simcoe_2019_wa_pop_var_2019 <- simcoe_2019_cv %>% group_by(Year) %>%
  summarize(wa_pop_mean_cv = mean(mean_pop_cv, na.rm = T),
            wa_pop_mean_density = mean(mean_pop_density),
            wa_pop_mean_variance = mean(mean_pop_variance),
            sum_w_pop_cv_2019 = sum(w_pop_cv_2019, na.rm = T),
            lcv = sum(w_pop_cv_2019, na.rm = TRUE)^2)

#Step 2
library(mefa)
simcoe_2019_cor_list_ts <- as.dist(round(cor(simcoe_2019_t2[]),2))
simcoe_2019_cor_ts <- stack(simcoe_2019_cor_list_ts, dim.names = TRUE)
simcoe_2019_cor_ts <- simcoe_2019_cor_ts %>%
  summarize(species_1 = row,
            species_2 = col,
            corr = dist)
detach("package:mefa", unload=TRUE)

#Step 3 - simcoe_2019_cor - simcoe_2019_cv
simcoe_c9_2019_w_cv_ts <- subset(simcoe_2019_cv, Station_ID == 'C9')
simcoe_c9_2019_w_cv_ts$Species <- paste("C9", simcoe_c9_2019_w_cv_ts$Species, sep=" ")
simcoe_c9_2019_w_cv_ts$Species_ID <- 1:nrow(simcoe_c9_2019_w_cv_ts)
simcoe_k42_2019_w_cv_ts <- subset(simcoe_2019_cv, Station_ID == 'K42')
simcoe_k42_2019_w_cv_ts$Species <- paste("K42", simcoe_k42_2019_w_cv_ts$Species, sep=" ")
simcoe_k42_2019_w_cv_ts$Species_ID <- 1:nrow(simcoe_k42_2019_w_cv_ts)
simcoe_k45_2019_w_cv_ts <- subset(simcoe_2019_cv, Station_ID == 'K45')
simcoe_k45_2019_w_cv_ts$Species <- paste("K45", simcoe_k45_2019_w_cv_ts$Species, sep=" ")
simcoe_k45_2019_w_cv_ts$Species_ID <- 1:nrow(simcoe_k45_2019_w_cv_ts)

simcoe_spec_2019_w_sp_ts <- rbind(simcoe_c9_2019_w_cv_ts, simcoe_k42_2019_w_cv_ts, simcoe_k45_2019_w_cv_ts)
simcoe_spec_2019_w_sp_ts <- simcoe_spec_2019_w_sp_ts %>% select(Species, w_pop_cv_2019, Species_ID, Station_ID)
simcoe_spec_2019_w_sp1_ts <- simcoe_spec_2019_w_sp_ts %>%
  mutate(species_1 = Species,
         station_1 = Station_ID,
         w_pop_cv_2019_sp1 = w_pop_cv_2019) %>%
  ungroup()%>%
  select(species_1, w_pop_cv_2019_sp1, Species_ID, station_1)

simcoe_2019_cor_cv_ts <- merge(simcoe_2019_cor_ts, simcoe_spec_2019_w_sp1_ts, by = "species_1")

simcoe_spec_2019_w_sp2_ts <- simcoe_spec_2019_w_sp_ts %>%
  mutate(species_2 = Species,
         station_2 = Station_ID,
         w_pop_cv_2019_sp2 = w_pop_cv_2019) %>%
  ungroup()%>%
  select(species_2, w_pop_cv_2019_sp2, Species_ID, station_2)

simcoe_2019_cor_cv_ts <- merge(simcoe_2019_cor_cv_ts, simcoe_spec_2019_w_sp2_ts, by = "species_2")

simcoe_2019_cor_cv_ts_omit <- na.omit(simcoe_2019_cor_cv_ts)

simcoe_2019_ind_W_ts <- simcoe_2019_cor_cv_ts_omit %>% group_by(species_1, species_2) %>%
  mutate(ind_pop_pair_cv_2019 = w_pop_cv_2019_sp1*w_pop_cv_2019_sp2,
         ind_W_2019 = (1 - corr)*(w_pop_cv_2019_sp1*w_pop_cv_2019_sp2),
         Year = 2019)

simcoe_2019_ind_W_ts$number_species_pairs <- nrow(simcoe_2019_ind_W_ts)

simcoe_mean_W_sp_pairs_2019 <- simcoe_2019_ind_W_ts %>% group_by(Year) %>%
  summarize(mean_W = mean(ind_W_2019, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

#Step 4
simcoe_2019_W_ts <- simcoe_2019_ind_W_ts %>% group_by(Year) %>%
  summarize(W_2019 = sum(ind_W_2019, na.rm = TRUE))



##### Type 2 Asynchrony - Local Stabilization - asynchrony among local species - Hammond et al 2020 #####
# Equation - delta = Sum[(1 - between population Pearson correlation coefficient @ site k) * weighted CV population i @ site k * weighted CV population j @ site k]
#Steps
# 1 - Create dataframes - populations separated by station
# 1 - Calculate population pair Pearson correlation coefficients for species by site
# 2 - Plug into equation for each population pair to calculate stabilization from each population pair & Calculate CV, SD, Mean for Local Communities
# 3 - Sum all to calculate Type 1 Asynchrony - Local Stabilization

#Step 1 
simcoe_2019_ind_W_ts_c9 <- filter(simcoe_2019_ind_W_ts, station_1 == "C9" & station_2 == "C9")
simcoe_2019_ind_W_ts_k42 <- filter(simcoe_2019_ind_W_ts, station_1 == "K42" & station_2 == "K42")
simcoe_2019_ind_W_ts_k45 <- filter(simcoe_2019_ind_W_ts, station_1 == "K45" & station_2 == "K45")


#Step 2 - Calculate Stabilization for each Station
#C9
simcoe_species_pairs_count_2019_c9 <- na.omit(simcoe_2019_ind_W_ts_c9)

simcoe_species_pairs_count_2019_c9$number_species_pairs <- nrow(simcoe_species_pairs_count_2019_c9)

simcoe_mean_W_sp_pairs_2019_c9 <- simcoe_species_pairs_count_2019_c9 %>% group_by(Year) %>%
  summarize(Station_ID = "C9",
            mean_lc_stab = mean(ind_W_2019, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2019_W_ts_c9 <- simcoe_2019_ind_W_ts_c9 %>% group_by(Year) %>%
  summarize(W_2019 = sum(ind_W_2019, na.rm = TRUE))

#K42
simcoe_species_pairs_count_2019_k42 <- na.omit(simcoe_2019_ind_W_ts_k42)

simcoe_species_pairs_count_2019_k42$number_species_pairs <- nrow(simcoe_species_pairs_count_2019_k42)

simcoe_mean_W_sp_pairs_2019_k42 <- simcoe_species_pairs_count_2019_k42 %>% group_by(Year) %>%
  summarize(Station_ID = "K42",
            mean_lc_stab = mean(ind_W_2019, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2019_W_ts_k42 <- simcoe_2019_ind_W_ts_k42 %>% group_by(Year) %>%
  summarize(W_2019 = sum(ind_W_2019, na.rm = TRUE))

#K45
simcoe_species_pairs_count_2019_k45 <- na.omit(simcoe_2019_ind_W_ts_k45)

simcoe_species_pairs_count_2019_k45$number_species_pairs <- nrow(simcoe_species_pairs_count_2019_k45)

simcoe_mean_W_sp_pairs_2019_k45 <- simcoe_species_pairs_count_2019_k45 %>% group_by(Year) %>%
  summarize(Station_ID = "K45",
            mean_lc_stab = mean(ind_W_2019, na.rm = T),
            total_species_pairs = mean(number_species_pairs))

simcoe_2019_W_ts_k45 <- simcoe_2019_ind_W_ts_k45 %>% group_by(Year) %>%
  summarize(W_2019 = sum(ind_W_2019, na.rm = TRUE))

#Step 4
simcoe_2019_ind_W_as1 <- rbind(simcoe_2019_W_ts_c9, simcoe_2019_W_ts_k42, simcoe_2019_W_ts_k45)

simcoe_2019_W_as1 <- simcoe_2019_ind_W_as1 %>% group_by(Year) %>%
  summarize(as2_2019 = sum(W_2019, na.rm = TRUE))

simcoe_mean_lc_sp_pairs_2019 <- rbind(simcoe_mean_W_sp_pairs_2019_c9, simcoe_mean_W_sp_pairs_2019_k42, simcoe_mean_W_sp_pairs_2019_k45)

simcoe_lc_sp_pairs_2019 <- simcoe_mean_lc_sp_pairs_2019 %>% group_by(Year) %>%
  summarize(mean_lc_stab = mean(mean_lc_stab),
            lc_sp_pairs = sum(total_species_pairs))

##### Type 4 Asynchrony - Metapopulation Stabilization - asynchrony of populations within metapopulations - Hammond et al 2020 #####
# Equation - Beta-mp = Sum[(1 - between population Pearson correlation coefficient of species i across sites) * weighted CV population i @ site k * weighted CV population i @ site l]
#Steps
# 1 - Create dataframes - populations separated by station
# 2 - Sum metapopulation asynchrony by species
# 3 - Sum all to calculate Type 4 Asynchrony - Metapopulation Stabilization

#Step 1
simcoe_2019_ind_W_as4 <- filter(simcoe_2019_ind_W_ts, Species_ID.x == Species_ID.y)

simcoe_species_pairs_count_2019_mp  <- na.omit(simcoe_2019_ind_W_as4)

simcoe_species_pairs_count_2019_mp$number_species_pairs <- nrow(simcoe_species_pairs_count_2019_mp)

simcoe_mean_mp_sp_pairs_2019  <- simcoe_species_pairs_count_2019_mp  %>% group_by(Year) %>%
  summarize(mean_mp_stab = mean(ind_W_2019, na.rm = T),
            mp_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2019_sum_W_as4 <- simcoe_2019_ind_W_as4 %>% group_by(Year, Species_ID.x) %>%
  summarize(as4_2019 = sum(ind_W_2019, na.rm = TRUE))


#Step 3
simcoe_2019_sum_W_as4 <- simcoe_2019_sum_W_as4 %>% group_by(Year) %>%
  summarize(as4_2019 = sum(as4_2019, na.rm = TRUE))



##### Type 5 Asynchrony - Cross-Community Stabilization - asynchrony of different species in different patches - Hammond et al 2020 #####
#Step 1
simcoe_2019_ind_W_as5 <- filter(simcoe_2019_ind_W_ts, station_1 != station_2 & Species_ID.x != Species_ID.y)


simcoe_species_pairs_count_2019_cc <- na.omit(simcoe_2019_ind_W_as5)

simcoe_species_pairs_count_2019_cc$number_species_pairs <- nrow(simcoe_species_pairs_count_2019_cc)

simcoe_mean_cc_sp_pairs_2019  <- simcoe_species_pairs_count_2019_cc  %>% group_by(Year) %>%
  summarize(mean_cc_stab = mean(ind_W_2019, na.rm = T),
            cc_species_pairs = mean(number_species_pairs))

#Step 2
simcoe_2019_sum_W_as5 <- simcoe_2019_ind_W_as5 %>% group_by(Year, Species_ID.x) %>%
  summarize(as5_2019 = sum(ind_W_2019, na.rm = TRUE))


#Step 3
simcoe_2019_sum_W_as5 <- simcoe_2019_sum_W_as5 %>% group_by(Year) %>%
  summarize(as5_2019 = sum(as5_2019, na.rm = TRUE))

simcoe_2019_W <- merge(simcoe_2019_sum_W_as4, simcoe_2019_W_ts, by = c('Year'))
simcoe_2019_W <- merge(simcoe_2019_W, simcoe_2019_W_as1, by = c('Year'))
simcoe_2019_W <- merge(simcoe_2019_W, simcoe_2019_sum_W_as5, by = c('Year'))


##### Combining Species Pairs Correlation & Stabilization Averages DFs #####
simcoe_sp_pair_ave_stab_2019 <- merge(simcoe_mean_W_sp_pairs_2019, simcoe_lc_sp_pairs_2019, by=c("Year"))
simcoe_sp_pair_ave_stab_2019 <- merge(simcoe_sp_pair_ave_stab_2019, simcoe_mean_mp_sp_pairs_2019, by=c("Year"))
simcoe_sp_pair_ave_stab_2019 <- merge(simcoe_sp_pair_ave_stab_2019, simcoe_mean_cc_sp_pairs_2019, by=c("Year"))



##### Create Final Dataframe with all Metrics #####
simcoe_2019_final <- merge(simcoe_2019_W, simcoe_2019_wa_pop_var_2019, by = c('Year'))
simcoe_2019_final <- simcoe_2019_final %>% group_by(Year) %>%
  mutate(gcv_2019 = lcv - W_2019) %>%
  group_by(Year) %>%
  summarize(gcv = gcv_2019, 
            lcv = lcv,
            W = W_2019,
            lc_stab = as2_2019,
            mp_stab = as4_2019,
            cc_stab = as5_2019,
            mc_cv = sqrt(gcv_2019),
            wa_pop_cv = sqrt(lcv),
            wa_pop_mean_cv = wa_pop_mean_cv,
            wa_pop_mean_density = wa_pop_mean_density,
            wa_pop_mean_variance = wa_pop_mean_variance,
            total_asynchrony = as2_2019/lcv + as4_2019/lcv + as5_2019/lcv,
            lc_asynchrony = as2_2019/lcv,
            mp_asynchrony = as4_2019/lcv,
            cc_asynchrony= as5_2019/lcv)

simcoe_2019_cv <- simcoe_2019_cv %>% group_by(Year) %>%
  summarize(mc_sum_mean_density = mean(mc_sum_mean_density))

simcoe_2019_final$mc_sum_mean_density <- simcoe_2019_cv$mc_sum_mean_density

##### Calculating Un-biased Mean Correlations by Scale for 2019 #####
simcoe_mean_corr_2019_c9 <- simcoe_species_pairs_count_2019_c9 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2019, ind_W_2019)
simcoe_mean_corr_2019_k42 <- simcoe_species_pairs_count_2019_k42 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2019, ind_W_2019)
simcoe_mean_corr_2019_k45 <- simcoe_species_pairs_count_2019_k45 %>% ungroup() %>% select(corr, ind_pop_pair_cv_2019, ind_W_2019)
simcoe_mean_corr_2019_lc <- rbind(simcoe_mean_corr_2019_c9, simcoe_mean_corr_2019_k42, simcoe_mean_corr_2019_k45)
simcoe_mean_corr_2019_mp <- simcoe_species_pairs_count_2019_mp %>% ungroup() %>% select(corr, ind_pop_pair_cv_2019, ind_W_2019)
simcoe_mean_corr_2019_cc <- simcoe_species_pairs_count_2019_cc %>% ungroup() %>% select(corr, ind_pop_pair_cv_2019, ind_W_2019)


rep <- 10000
for(irep in 1:rep) {
  
  mean_lc <- simcoe_mean_corr_2019_lc[sample(1:nrow(simcoe_mean_corr_2019_lc),  nrow(simcoe_mean_corr_2019_mp)), ]
  
  mean_lc <- mean_lc %>% 
    summarize(mean_lc_corr = mean(corr),
              sd_lc_corr = sd(corr),
              mean_lc_pop_pair_cv = mean(ind_pop_pair_cv_2019),
              sd_lc_pop_pair_cv = sd(ind_pop_pair_cv_2019),
              mean_lc_stab = mean(ind_W_2019),
              sd_lc_stab = sd(ind_W_2019))
  
  lc_corr_2019 <- as.data.frame(mean_lc)
}
lc_corr_2019

mp_corr_2019 <- simcoe_mean_corr_2019_mp %>% summarize(mean_mp_corr = mean(corr),
                                                       sd_mp_corr = sd(corr),
                                                       mean_mp_pop_pair_cv = mean(ind_pop_pair_cv_2019),
                                                       sd_mp_pop_pair_cv = sd(ind_pop_pair_cv_2019),
                                                       mean_mp_stab = mean(ind_W_2019),
                                                       sd_mp_stab = sd(ind_W_2019))

rep <- 10000
for(irep in 1:rep) {
  
  mean_cc <- simcoe_mean_corr_2019_cc[sample(1:nrow(simcoe_mean_corr_2019_cc), nrow(simcoe_mean_corr_2019_mp)), ]
  
  mean_cc <- mean_cc %>% 
    summarize(mean_cc_corr = mean(corr),
              sd_cc_corr = sd(corr),
              mean_cc_pop_pair_cv = mean(ind_pop_pair_cv_2019),
              sd_cc_pop_pair_cv = sd(ind_pop_pair_cv_2019),
              mean_cc_stab = mean(ind_W_2019),
              sd_cc_stab = sd(ind_W_2019))
  
  cc_corr_2019 <- as.data.frame(mean_cc)
}
cc_corr_2019


corr_2019 <- cbind(lc_corr_2019, mp_corr_2019)
corr_2019 <- cbind(corr_2019, cc_corr_2019)
corr_2019 <- corr_2019 %>% mutate(Year = 2019)




##### Compile all final dataframes #####
simcoe_hammond <- rbind(simcoe_1986_final, simcoe_1987_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_1988_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_1989_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_1990_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_1991_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_1992_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_1993_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_1994_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_1995_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_1996_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_1997_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_1998_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_1999_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_2000_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_2001_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_2002_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_2003_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_2004_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_2005_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_2006_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_2007_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_2008_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_2009_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_2010_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_2011_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_2012_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_2013_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_2014_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_2015_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_2016_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_2017_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_2018_final)
simcoe_hammond <- rbind(simcoe_hammond, simcoe_2019_final)



simcoe_hammond <- simcoe_hammond %>% group_by(Year) %>%
  mutate(mc_variance = mc_cv*mc_sum_mean_density,
         wa_pop_variance = wa_pop_cv*mc_sum_mean_density)

simcoe_final <- merge(simcoe_hammond, simcoe_diversity, by = c('Year'))


write.csv(simcoe_final, "simcoe_final_weeks_3.csv")




simcoe_corr <- rbind(corr_1986, corr_1987)
simcoe_corr <- rbind(simcoe_corr, corr_1988)
simcoe_corr <- rbind(simcoe_corr, corr_1989)
simcoe_corr <- rbind(simcoe_corr, corr_1990)
simcoe_corr <- rbind(simcoe_corr, corr_1991)
simcoe_corr <- rbind(simcoe_corr, corr_1992)
simcoe_corr <- rbind(simcoe_corr, corr_1993)
simcoe_corr <- rbind(simcoe_corr, corr_1994)
simcoe_corr <- rbind(simcoe_corr, corr_1995)
simcoe_corr <- rbind(simcoe_corr, corr_1996)
simcoe_corr <- rbind(simcoe_corr, corr_1997)
simcoe_corr <- rbind(simcoe_corr, corr_1998)
simcoe_corr <- rbind(simcoe_corr, corr_1999)
simcoe_corr <- rbind(simcoe_corr, corr_2000)
simcoe_corr <- rbind(simcoe_corr, corr_2001)
simcoe_corr <- rbind(simcoe_corr, corr_2002)
simcoe_corr <- rbind(simcoe_corr, corr_2003)
simcoe_corr <- rbind(simcoe_corr, corr_2004)
simcoe_corr <- rbind(simcoe_corr, corr_2005)
simcoe_corr <- rbind(simcoe_corr, corr_2006)
simcoe_corr <- rbind(simcoe_corr, corr_2007)
simcoe_corr <- rbind(simcoe_corr, corr_2008)
simcoe_corr <- rbind(simcoe_corr, corr_2009)
simcoe_corr <- rbind(simcoe_corr, corr_2010)
simcoe_corr <- rbind(simcoe_corr, corr_2011)
simcoe_corr <- rbind(simcoe_corr, corr_2012)
simcoe_corr <- rbind(simcoe_corr, corr_2013)
simcoe_corr <- rbind(simcoe_corr, corr_2014)
simcoe_corr <- rbind(simcoe_corr, corr_2015)
simcoe_corr <- rbind(simcoe_corr, corr_2016)
simcoe_corr <- rbind(simcoe_corr, corr_2017)
simcoe_corr <- rbind(simcoe_corr, corr_2018)
simcoe_corr <- rbind(simcoe_corr, corr_2019)



write.csv(simcoe_corr, "simcoe_final_corr_weeks_3.csv")



simcoe_pair_stab <- rbind(simcoe_sp_pair_ave_stab_1986, simcoe_sp_pair_ave_stab_1987)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_1988)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_1989)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_1990)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_1991)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_1992)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_1993)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_1994)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_1995)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_1996)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_1997)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_1998)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_1999)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_2000)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_2001)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_2002)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_2003)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_2004)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_2005)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_2006)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_2007)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_2008)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_2009)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_2010)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_2011)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_2012)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_2013)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_2014)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_2015)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_2016)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_2017)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_2018)
simcoe_pair_stab <- rbind(simcoe_pair_stab, simcoe_sp_pair_ave_stab_2019)


write.csv(simcoe_pair_stab, "simcoe_final_pair_stab_weeks_3.csv")












