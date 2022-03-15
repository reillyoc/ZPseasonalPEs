simcoe <- read.csv('simcoe_final_weeks_5.csv', header = T)
simcoe_div <- read.csv('simcoe_final_5_weeks_diversity', header = T)
simcoe_corr <- read.csv('simcoe_final_corr_weeks_5.csv', header = T)
simcoe_pairst <- read.csv('simcoe_final_pair_stab_weeks_5.csv', header = T)


library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

Year <- c(1986:1999)
Global_Change <- c('Pre')
simcoe_gc_pre <- data.frame(Year, Global_Change)

Year <- c(2008:2019)
Global_Change <- c('Post')
simcoe_gc_post <- data.frame(Year, Global_Change)

simcoe_gc <- rbind(simcoe_gc_pre, simcoe_gc_post)

simcoe <- merge(simcoe_gc, simcoe, by = c('Year'))
simcoe_corr <- merge(simcoe_gc, simcoe_corr, by = c('Year'))
simcoe_pairst <- merge(simcoe_gc, simcoe_pairst, by = c('Year'))

###### Asynchrony - Stabilization Figures #####
#Asynchrony vs Stabilization


lm_a1_stab <- lm(data = simcoe, W~total_asynchrony)
summary(lm_a1_stab)

lm_a2_stab <- lm(data = simcoe, lc_stab~lc_asynchrony)
summary(lm_a2_stab)

lm_a4_stab <- lm(data = simcoe, mp_stab~mp_asynchrony)
summary(lm_a4_stab)

lm_a5_stab <- lm(data = simcoe, cc_stab~cc_asynchrony)
summary(lm_a5_stab)


gg_simcoe_a1_stab <- simcoe %>% select(total_asynchrony, lcv, W, Shannon_Div_Local, Evenness_Local) %>%
  summarize(Asynchrony = total_asynchrony,
            lcv = lcv,
            Stabilization = W,
            Shannon_Div_Local = Shannon_Div_Local, 
            Evenness_Local = Evenness_Local,
            Color = 'Total')

gg_simcoe_a2_stab <- simcoe %>% select(lc_asynchrony, lcv, lc_stab, Shannon_Div_Local, Evenness_Local) %>%
  summarize(Asynchrony = lc_asynchrony,
            lcv = lcv,
            Stabilization = lc_stab,
            Shannon_Div_Local = Shannon_Div_Local, 
            Evenness_Local = Evenness_Local,
            Color = 'Local Community')

gg_simcoe_a4_stab <- simcoe %>% select(mp_asynchrony, lcv, mp_stab, Shannon_Div_Local, Evenness_Local) %>%
  summarize(Asynchrony = mp_asynchrony,
            lcv = lcv,
            Stabilization = mp_stab,
            Shannon_Div_Local = Shannon_Div_Local, 
            Evenness_Local = Evenness_Local,
            Color = 'Metapopulation')

gg_simcoe_a5_stab <- simcoe %>% select(cc_asynchrony, lcv, cc_stab, Shannon_Div_Local, Evenness_Local) %>%
  summarize(Asynchrony = cc_asynchrony,
            lcv = lcv,
            Stabilization = cc_stab,
            Shannon_Div_Local = Shannon_Div_Local, 
            Evenness_Local = Evenness_Local,
            Color = 'Cross-Community')

gg_simcoe_asynch_stab <- rbind(gg_simcoe_a1_stab , gg_simcoe_a2_stab)
gg_simcoe_asynch_stab <- rbind(gg_simcoe_asynch_stab , gg_simcoe_a4_stab)
gg_simcoe_asynch_stab <- rbind(gg_simcoe_asynch_stab , gg_simcoe_a5_stab)

gg_simcoe_asynch_stab$Color <- factor(gg_simcoe_asynch_stab$Color, 
                                      levels = c("Local Community",
                                                 "Metapopulation",
                                                 "Cross-Community",
                                                 "Total"))

#Asynchrony vs Stabilization
gg_simcoe_asynch_stab_plot <- ggplot(data = gg_simcoe_asynch_stab, 
                                     aes(x = Asynchrony, y = Stabilization, group = Color, color = Color)) +
  geom_point() +
  geom_smooth(se = F, method = "lm") +
  xlab('Asynchrony') +
  ylab('Stabilization') +
  theme_classic() +
  scale_color_brewer(palette = 'Dark2') +
  theme(legend.position = 'none')
gg_simcoe_asynch_stab_plot

#Weighted average local population variability vs Stabilization
lm_lcv1_stab <- lm(data = simcoe, W~lcv)
summary(lm_lcv1_stab)

lm_lcv2_stab <- lm(data = simcoe, lc_stab~lcv)
summary(lm_lcv2_stab)

lm_lcv4_stab <- lm(data = simcoe, mp_stab~lcv)
summary(lm_lcv4_stab)

lm_lcv5_stab <- lm(data = simcoe, cc_stab~lcv)
summary(lm_lcv5_stab)

#Weighted average local population variability  vs Stabilization
gg_simcoe_lcv_stab_plot <- ggplot(data = gg_simcoe_asynch_stab, 
                                  aes(x = lcv, y = Stabilization, group = Color, color = Color)) +
  geom_point() +
  geom_smooth(se = F, method = "lm") +
  xlab('Local Population Variability') +
  ylab('Stabilization') +
  theme_classic() +
  scale_color_brewer(palette = 'Dark2') +
  theme(legend.position = 'none')
gg_simcoe_lcv_stab_plot


#Weighted average local population variability  vs Asynchrony
lm_lcv1_stab <- lm(data = simcoe, lcv~total_asynchrony)
summary(lm_lcv1_stab)

lm_lcv2_stab <- lm(data = simcoe, lcv~lc_asynchrony)
summary(lm_lcv2_stab)

lm_lcv4_stab <- lm(data = simcoe, lcv~mp_asynchrony)
summary(lm_lcv4_stab)

lm_lcv5_stab <- lm(data = simcoe, lcv~cc_asynchrony)
summary(lm_lcv5_stab)

#Weighted average local population variability  vs Asynchrony
gg_simcoe_lcv_stab_plot <- ggplot(data = gg_simcoe_asynch_stab, 
                                  aes(x = lcv, y = Asynchrony, group = Color, color = Color)) +
  geom_point() +
  geom_smooth(se = F, method = "lm") +
  xlab('Local Population Variability') +
  ylab('Asynchrony') +
  theme_classic() +
  scale_color_brewer(palette = 'Dark2') +
  theme(legend.position = 'none')
gg_simcoe_lcv_stab_plot





simcoe_ave_asynch_app <- simcoe %>% select(Year, lc_asynchrony, mp_asynchrony, cc_asynchrony)

simcoe_ave_asynch_app$Species_Richness <- simcoe_div$Species_Richness

simcoe_ave_asynch_app <- simcoe_ave_asynch_app %>% group_by(Year) %>%
  summarize(lc_asynch_ppair = lc_asynchrony/5,
            mp_asynch_ppair = mp_asynchrony/Species_Richness,
            cc_asynch_ppair = cc_asynchrony/10)

###### Create new data frame subset with categories for Average Asynchrony per Population Pair for anova ##### 
simcoe_lc_asynch_app <- simcoe_ave_asynch_app %>% select(lc_asynch_ppair) %>%
  summarize(Stabilization = lc_asynch_ppair,
            Scale = 'Local Community')
simcoe_mp_asynch_app <- simcoe_ave_asynch_app %>% select(mp_asynch_ppair) %>%
  summarize(Stabilization = mp_asynch_ppair,
            Scale = 'Metapopulation')
simcoe_cc_asynch_app <- simcoe_ave_asynch_app %>% select(cc_asynch_ppair) %>%
  summarize(Stabilization = cc_asynch_ppair,
            Scale = 'Cross-Community')

simcoe_asynch_app <- rbind(simcoe_lc_asynch_app, simcoe_mp_asynch_app)
simcoe_asynch_app <- rbind(simcoe_asynch_app, simcoe_cc_asynch_app)

asynch_app_aov <- aov(Stabilization~Scale, data = simcoe_asynch_app)
summary(asynch_app_aov)
TukeyHSD(asynch_app_aov)



##### Create new data frame subset with categories for Average Asynchrony per Population Pair for gg barplot #####

gg_simcoe_lc_asynch_app <- simcoe_ave_asynch_app %>% select(lc_asynch_ppair) %>%
  summarize(CV = mean(lc_asynch_ppair),
            sd_stab = sd(lc_asynch_ppair),
            count = n(),
            se_stab = (sd_stab/(sqrt(count))),
            Scale2 = 'Stabilization',
            Scale = 'Local Community',
            Color = 'blue')
gg_simcoe_mp_asynch_app <- simcoe_ave_asynch_app %>% select(mp_asynch_ppair) %>%
  summarize(CV = mean(mp_asynch_ppair),
            sd_stab = sd(mp_asynch_ppair),
            count = n(),
            se_stab = (sd_stab/(sqrt(count))),
            Scale2 = 'Stabilization',
            Scale = 'Metapopulation',
            Color = 'green')
gg_simcoe_cc_asynch_app <- simcoe_ave_asynch_app %>% select(cc_asynch_ppair) %>%
  summarize(CV = mean(cc_asynch_ppair),
            sd_stab = sd(cc_asynch_ppair),
            count = n(),
            se_stab = (sd_stab/(sqrt(count))),
            Scale2 = 'Stabilization',
            Scale = 'Cross-Community',
            Color = 'red')

gg_simcoe_asynch_app <- rbind(gg_simcoe_mp_asynch_app, gg_simcoe_lc_asynch_app)
gg_simcoe_asynch_app <- rbind(gg_simcoe_asynch_app, gg_simcoe_cc_asynch_app)

gg_simcoe_asynch_app$Scale <- factor(gg_simcoe_asynch_app$Scale, levels = c("Local Community", 
                                                                            "Metapopulation", 
                                                                            "Cross-Community"))


gg_box_sim_asynch_app <- ggplot(data = gg_simcoe_asynch_app, aes(x = Scale2, y = CV, fill = Color)) +
  geom_bar(stat = 'identity',
           color = 'black',
           width = 0.5,
           position = position_dodge2()) +
  geom_errorbar(aes(ymin = CV - se_stab, ymax = CV + se_stab), 
                width = 0.5,
                position = position_dodge2()) +
  ylab("Average Annual Asynchrony per Unit") +
  xlab("Scale") +
  theme_classic() +
  theme(legend.position = 'none') +
  scale_fill_brewer(palette = 'Dark2')
gg_box_sim_asynch_app



###### Create new data frame subset with categories for Asynchrony anova #####
simcoe_ts <- simcoe %>% select(total_asynchrony) %>%
  summarize(Stabilization = total_asynchrony,
            Scale = 'Total')
simcoe_lc_asynchrony <- simcoe %>% select(lc_asynchrony) %>%
  summarize(Stabilization = lc_asynchrony,
            Scale = 'Local Community')
simcoe_mp_asynchrony <- simcoe %>% select(mp_asynchrony) %>%
  summarize(Stabilization = mp_asynchrony,
            Scale = 'Metapopulation')
simcoe_cc_asynchrony <- simcoe %>% select(cc_asynchrony) %>%
  summarize(Stabilization = cc_asynchrony,
            Scale = 'Cross-Community')

simcoe_asynchrony <- rbind(simcoe_lc_asynchrony, simcoe_mp_asynchrony, simcoe_cc_asynchrony)

stab_aov <- aov(Stabilization~Scale, data = simcoe_asynchrony)
summary(stab_aov)
TukeyHSD(stab_aov)


##### Create new data frame subset with categories for Asynchrony for gg barplot #####

gg_simcoe_mc_asynchrony <- simcoe %>% select(total_asynchrony) %>%
  summarize(CV = mean(total_asynchrony),
            sd_asynchrony = sd(total_asynchrony),
            count = n(),
            se_asynchrony = (sd_asynchrony/(sqrt(count))),
            Scale = 'Total',
            Color = 'blue')
gg_simcoe_lc_asynchrony <- simcoe %>% select(lc_asynchrony) %>%
  summarize(CV = mean(lc_asynchrony),
            sd_asynchrony = sd(lc_asynchrony),
            count = n(),
            se_asynchrony = (sd_asynchrony/(sqrt(count))),
            Scale2 = 'Asynchrony',
            Scale = 'Local Community',
            Color = 'blue')
gg_simcoe_mp_asynchrony <- simcoe %>% select(mp_asynchrony) %>%
  summarize(CV = mean(mp_asynchrony),
            sd_asynchrony = sd(mp_asynchrony),
            count = n(),
            se_asynchrony = (sd_asynchrony/(sqrt(count))),
            Scale2 = 'Asynchrony',
            Scale = 'Metapopulation',
            Color = 'green')
gg_simcoe_cc_asynchrony <- simcoe %>% select(cc_asynchrony) %>%
  summarize(CV = mean(cc_asynchrony),
            sd_asynchrony = sd(cc_asynchrony),
            count = n(),
            se_asynchrony = (sd_asynchrony/(sqrt(count))),
            Scale2 = 'Asynchrony',
            Scale = 'Cross-Community',
            Color = 'red')

gg_simcoe_asynchrony <- rbind(gg_simcoe_mp_asynchrony, gg_simcoe_lc_asynchrony)
gg_simcoe_asynchrony <- rbind(gg_simcoe_asynchrony, gg_simcoe_cc_asynchrony)

gg_simcoe_asynchrony$Scale <- factor(gg_simcoe_asynchrony$Scale, levels = c("Total", 
                                                                            "Local Community", 
                                                                            "Metapopulation", 
                                                                            "Cross-Community"))


gg_box_sim_asynchrony <- ggplot(data = gg_simcoe_asynchrony, aes(x = Scale2, y = CV, fill = Color)) +
  geom_bar(stat = 'identity',
           width = 0.3,
           color = 'black') +
  ylab("Average Annual Asynchrony") +
  xlab("Scale") +
  theme_classic() +
  theme(legend.position = 'none') +
  ylim(0, 1) +
  scale_fill_brewer(palette = 'Dark2')
gg_box_sim_asynchrony



#Asynchrony vs Evenness
lm_even_a1 <- lm(data = simcoe, total_asynchrony~Evenness_Local)
summary(lm_even_a1)

lm_even_a2 <- lm(data = simcoe, lc_asynchrony~Evenness_Local)
summary(lm_even_a2)

lm_even_a4 <- lm(data = simcoe, mp_asynchrony~Evenness_Local)
summary(lm_even_a4)

lm_even_a5 <- lm(data = simcoe, cc_asynchrony~Evenness_Local)
summary(lm_even_a5)


gg_simcoe_shann_even_a1 <- simcoe %>% select(total_asynchrony, Evenness_Local) %>%
  summarize(Asynchrony = total_asynchrony,
            Evenness_Local = Evenness_Local,
            Color = 'Total')

gg_simcoe_shann_even_a2 <- simcoe %>% select(lc_asynchrony, Evenness_Local) %>%
  summarize(Asynchrony = lc_asynchrony,
            Evenness_Local = Evenness_Local,
            Color = 'Local Community')

gg_simcoe_shann_even_a4 <- simcoe %>% select(mp_asynchrony, Evenness_Local) %>%
  summarize(Asynchrony = mp_asynchrony,
            Evenness_Local = Evenness_Local,
            Color = 'Metapopulation')

gg_simcoe_shann_even_a5 <- simcoe %>% select(cc_asynchrony, Evenness_Local) %>%
  summarize(Asynchrony = cc_asynchrony,
            Evenness_Local = Evenness_Local,
            Color = 'Cross-Community')

gg_simcoe_shan_even_asynch <- rbind(gg_simcoe_shann_even_a1 , gg_simcoe_shann_even_a2)
gg_simcoe_shan_even_asynch <- rbind(gg_simcoe_shan_even_asynch , gg_simcoe_shann_even_a4)
gg_simcoe_shan_even_asynch <- rbind(gg_simcoe_shan_even_asynch , gg_simcoe_shann_even_a5)

gg_simcoe_shan_even_asynch$Color <- factor(gg_simcoe_shan_even_asynch$Color, 
                                           levels = c("Local Community",
                                                      "Metapopulation",
                                                      "Cross-Community",
                                                      "Total"))

#Evenness vs Asynchrony
gg_sim_ham_26_even_asynch <- ggplot() +
  geom_point(data = simcoe, aes(x = Evenness_Local, y = total_asynchrony)) +
  geom_smooth(data = simcoe, aes(x = Evenness_Local, y = total_asynchrony, color = "purple"), se = F, method = "lm") +
  geom_point(data = simcoe, aes(x = Evenness_Local, y = mp_asynchrony)) +
  geom_smooth(data = simcoe, aes(x = Evenness_Local, y = mp_asynchrony, color = "orange"), se = F, method = "lm", linetype = 'dashed') +
  geom_point(data = simcoe, aes(x = Evenness_Local, y = lc_asynchrony)) +
  geom_smooth(data = simcoe, aes(x = Evenness_Local, y = lc_asynchrony, color = "green"), se = F, method = "lm") +
  geom_point(data = simcoe, aes(x = Evenness_Local, y = cc_asynchrony)) +
  geom_smooth(data = simcoe, aes(x = Evenness_Local, y = cc_asynchrony, color = "pink"), se = F, method = "lm") +
  xlab('Evenness') +
  ylab('Asynchrony') +
  theme_classic() +
  scale_color_brewer(palette = 'Dark2') +
  theme(legend.position = 'none')
gg_sim_ham_26_even_asynch


#Shannon Diversity vs Asynchrony
gg_sim_ham_26_shan_asynch <- ggplot() +
  geom_point(data = simcoe, aes(x = Shannon_Div_Local, y = total_asynchrony)) +
  geom_smooth(data = simcoe, aes(x = Shannon_Div_Local, y = total_asynchrony, color = "purple"), se = F, method = "lm") +
  geom_point(data = simcoe, aes(x = Shannon_Div_Local, y = mp_asynchrony)) +
  geom_smooth(data = simcoe, aes(x = Shannon_Div_Local, y = mp_asynchrony, color = "orange"), se = F, method = "lm", linetype = 'dashed') +
  geom_point(data = simcoe, aes(x = Shannon_Div_Local, y = lc_asynchrony)) +
  geom_smooth(data = simcoe, aes(x = Shannon_Div_Local, y = lc_asynchrony, color = "green"), se = F, method = "lm") +
  geom_point(data = simcoe, aes(x = Shannon_Div_Local, y = cc_asynchrony)) +
  geom_smooth(data = simcoe, aes(x = Shannon_Div_Local, y = cc_asynchrony, color = "pink"), se = F, method = "lm") +
  xlab('Shannon Diversity') +
  ylab('Asynchrony') +
  theme_classic() +
  scale_color_brewer(palette = 'Dark2') +
  theme(legend.position = 'none')
gg_sim_ham_26_shan_asynch



##### Create new data frame subset with categories for Asynchrony for gg barplot #####

gg_simcoe_lc_asynch <- simcoe %>% select(lc_asynchrony, Global_Change) %>%
  group_by(Global_Change) %>%
  summarize(CV = mean(lc_asynchrony),
            sd_asynchrony = sd(lc_asynchrony),
            count = n(),
            se_asynchrony = (sd_asynchrony/(sqrt(count))),
            Scale2 = 'Asynchrony',
            Scale = 'Local Community',
            Color = 'blue')
gg_simcoe_mp_asynch  <- simcoe %>% select(mp_asynchrony, Global_Change) %>%
  group_by(Global_Change) %>%
  summarize(CV = mean(mp_asynchrony),
            sd_asynchrony = sd(mp_asynchrony),
            count = n(),
            se_asynchrony = (sd_asynchrony/(sqrt(count))),
            Scale2 = 'Asynchrony',
            Scale = 'Metapopulation',
            Color = 'green')
gg_simcoe_cc_asynch <- simcoe %>% select(cc_asynchrony, Global_Change) %>%
  group_by(Global_Change) %>%
  summarize(CV = mean(cc_asynchrony),
            sd_asynchrony = sd(cc_asynchrony),
            count = n(),
            se_asynchrony = (sd_asynchrony/(sqrt(count))),
            Scale2 = 'Asynchrony',
            Scale = 'Cross-Community',
            Color = 'red')

gg_simcoe_asynch  <- rbind(gg_simcoe_mp_asynch , gg_simcoe_lc_asynch )
gg_simcoe_asynch  <- rbind(gg_simcoe_asynch, gg_simcoe_cc_asynch )

gg_simcoe_asynch$Scale <- factor(gg_simcoe_asynch$Scale, levels = c( 
  "Local Community", 
  "Metapopulation", 
  "Cross-Community"))

gg_simcoe_asynch$Global_Change <- factor(gg_simcoe_asynch$Global_Change, 
                                         levels = c("Pre", 
                                                    "Post"))

gg_box_sim_stab <- ggplot(data = gg_simcoe_asynch, aes(x = Global_Change, y = CV, fill = Color)) +
  geom_bar(stat = 'identity',
           color = 'black',
           width = 0.3) +
  ylab("Average Annual Asynchrony") +
  xlab("Scale") +
  theme_classic() +
  ylim(0, 1.4) +
  theme(legend.position = 'none') +
  scale_fill_brewer(palette = 'Dark2')
gg_box_sim_stab



##### Global Change T-Tests #####
t.test(total_asynchrony~Global_Change, data = simcoe)
t.test(lc_asynchrony~Global_Change, data = simcoe)
t.test(mp_asynchrony~Global_Change, data = simcoe)
t.test(cc_asynchrony~Global_Change, data = simcoe)


###### Create new data frame subset with categories for Average Stabilization per Population Pair for anova ##### 
simcoe_lc_stab_app <- simcoe_corr %>% select(mean_lc_stab) %>%
  summarize(Stabilization = mean_lc_stab,
            Scale = 'Local Community')
simcoe_mp_stab_app <- simcoe_corr %>% select(mean_mp_stab) %>%
  summarize(Stabilization = mean_mp_stab,
            Scale = 'Metapopulation')
simcoe_cc_stab_app <- simcoe_corr %>% select(mean_cc_stab) %>%
  summarize(Stabilization = mean_cc_stab,
            Scale = 'Cross-Community')

simcoe_stab_app <- rbind(simcoe_lc_stab_app, simcoe_mp_stab_app)
simcoe_stab_app <- rbind(simcoe_stab_app, simcoe_cc_stab_app)

stab_app_aov <- aov(Stabilization~Scale, data = simcoe_stab_app)
summary(stab_app_aov)
TukeyHSD(stab_app_aov)


##### Create new data frame subset with categories for Average Stabilization per Population Pair for gg barplot #####

gg_simcoe_lc_stab_app <- simcoe_corr %>% select(mean_lc_stab) %>%
  summarize(CV = mean(mean_lc_stab),
            sd_stab = sd(mean_lc_stab),
            count = n(),
            se_stab = (sd_stab/(sqrt(count))),
            Scale2 = 'Stabilization',
            Scale = 'Local Community',
            Color = 'blue')
gg_simcoe_mp_stab_app <- simcoe_corr %>% select(mean_mp_stab) %>%
  summarize(CV = mean(mean_mp_stab),
            sd_stab = sd(mean_mp_stab),
            count = n(),
            se_stab = (sd_stab/(sqrt(count))),
            Scale2 = 'Stabilization',
            Scale = 'Metapopulation',
            Color = 'green')
gg_simcoe_cc_stab_app <- simcoe_corr %>% select(mean_cc_stab) %>%
  summarize(CV = mean(mean_cc_stab),
            sd_stab = sd(mean_cc_stab),
            count = n(),
            se_stab = (sd_stab/(sqrt(count))),
            Scale2 = 'Stabilization',
            Scale = 'Cross-Community',
            Color = 'red')

gg_simcoe_stab_app <- rbind(gg_simcoe_mp_stab_app, gg_simcoe_lc_stab_app)
gg_simcoe_stab_app <- rbind(gg_simcoe_stab_app, gg_simcoe_cc_stab_app)

gg_simcoe_stab_app$Scale <- factor(gg_simcoe_stab_app$Scale, levels = c("Local Community", 
                                                                        "Metapopulation", 
                                                                        "Cross-Community"))


gg_box_sim_stab_app <- ggplot(data = gg_simcoe_stab_app, aes(x = Scale2, y = CV, fill = Color)) +
  geom_bar(stat = 'identity',
           color = 'black',
           width = 0.5,
           position = position_dodge2()) +
  geom_errorbar(aes(ymin = CV - se_stab, ymax = CV + se_stab), 
                width = 0.5,
                position = position_dodge2()) +
  ylab("Average Annual Stabilization per Population Pair") +
  xlab("Scale") +
  theme_classic() +
  theme(legend.position = 'none') +
  scale_fill_brewer(palette = 'Dark2')
gg_box_sim_stab_app


##### Species Pairs Barplots - simcoe_pairs #####
#Subset new dataframes
simcoe_ave_stab <- simcoe_pairst %>% select(Year, mean_W, mean_lc_stab, mean_mp_stab, mean_cc_stab)
simcoe_num_pairs <- simcoe_pairst %>% select(Year, lc_sp_pairs, mp_species_pairs, cc_species_pairs)


simcoe_num_pairs <- simcoe_num_pairs %>% group_by(Year) %>%
  mutate(log_lc_sp_pairs = log(lc_sp_pairs),
         log_mp_sp_pairs = log(mp_species_pairs),
         log_cc_sp_pairs = log(cc_species_pairs))

simcoe_ave_stab <- simcoe_ave_stab %>% group_by(Year) %>%
  mutate(log_total_sp_pairs = log(mean_W),
         log_lc_stab = log(mean_lc_stab),
         log_mp_stab = log(mean_mp_stab),
         log_cc_stab = log(mean_cc_stab))

###### Create new data frame subset with categories for log population pairs anova #####
simcoe_lc_pairs <- simcoe_num_pairs %>% select(log_lc_sp_pairs) %>%
  summarize(Species_Pairs = log_lc_sp_pairs,
            Scale = 'Local Community')
simcoe_mp_pairs <- simcoe_num_pairs %>% select(log_mp_sp_pairs) %>%
  summarize(Species_Pairs = log_mp_sp_pairs,
            Scale = 'Metapopulation')
simcoe_cc_pairs <- simcoe_num_pairs %>% select(log_cc_sp_pairs) %>%
  summarize(Species_Pairs = log_cc_sp_pairs,
            Scale = 'Cross-Community')


gg_simcoe_log_pairs <- rbind(simcoe_lc_pairs, simcoe_mp_pairs)
gg_simcoe_log_pairs <- rbind(gg_simcoe_log_pairs, simcoe_cc_pairs)

pairs_aov <- aov(Species_Pairs~Scale, data = gg_simcoe_log_pairs)
summary(pairs_aov)
TukeyHSD(pairs_aov)




##### Create new data frame subset with categories for Asynchrony for gg barplot #####

gg_simcoe_lc_pairs <- simcoe_num_pairs %>% ungroup() %>% 
  select(log_lc_sp_pairs) %>%
  summarize(CV = mean(log_lc_sp_pairs),
            sd_asynchrony = sd(log_lc_sp_pairs),
            count = n(),
            se_asynchrony = (sd_asynchrony/(sqrt(count))),
            Scale2 = 'Asynchrony',
            Scale = 'Local Community',
            Color = 'blue')
gg_simcoe_mp_pairs <- simcoe_num_pairs %>%  ungroup() %>% 
  select(log_mp_sp_pairs) %>%
  summarize(CV = mean(log_mp_sp_pairs),
            sd_asynchrony = sd(log_mp_sp_pairs),
            count = n(),
            se_asynchrony = (sd_asynchrony/(sqrt(count))),
            Scale2 = 'Asynchrony',
            Scale = 'Metapopulation',
            Color = 'green')
gg_simcoe_cc_pairs <- simcoe_num_pairs %>%  ungroup() %>% 
  select(log_cc_sp_pairs) %>%
  summarize(CV = mean(log_cc_sp_pairs),
            sd_asynchrony = sd(log_cc_sp_pairs),
            count = n(),
            se_asynchrony = (sd_asynchrony/(sqrt(count))),
            Scale2 = 'Asynchrony',
            Scale = 'Cross-Community',
            Color = 'red')

gg_simcoe_pairs <- rbind(gg_simcoe_lc_pairs, gg_simcoe_mp_pairs)
gg_simcoe_pairs <- rbind(gg_simcoe_pairs, gg_simcoe_cc_pairs)

gg_simcoe_pairs$Scale <- factor(gg_simcoe_pairs$Scale, levels = c("Local Community", 
                                                                  "Metapopulation", 
                                                                  "Cross-Community"))


gg_box_sim_pairs <- ggplot(data = gg_simcoe_pairs, aes(x = Scale2, y = CV, fill = Color)) +
  geom_bar(stat = 'identity',
           color = 'black',
           width = 0.5,
           position = position_dodge2()) +
  geom_errorbar(aes(ymin = CV - se_asynchrony, ymax = CV + se_asynchrony), 
                width = 0.5,
                position = position_dodge2()) +
  ylab("Average Number of Population Pairs (log)") +
  xlab("Scale") +
  theme_classic() +
  theme(legend.position = 'none') +
  coord_cartesian(ylim = c(4, 8.5)) +
  scale_fill_brewer(palette = 'Dark2')
gg_box_sim_pairs




###### Create new data frame subse t with categories for Mean Correlation for anova ##### 
simcoe_lc_corr <- simcoe_corr %>% select(mean_lc_corr) %>%
  summarize(Mean_Correlation = mean_lc_corr,
            Scale = 'Local Community')
simcoe_mp_corr <- simcoe_corr %>% select(mean_mp_corr) %>%
  summarize(Mean_Correlation = mean_mp_corr,
            Scale = 'Metapopulation')
simcoe_cc_corr <- simcoe_corr %>% select(mean_cc_corr) %>%
  summarize(Mean_Correlation = mean_cc_corr,
            Scale = 'Cross-Community')

simcoe_mean_corr <- rbind(simcoe_lc_corr, simcoe_mp_corr)
simcoe_mean_corr <- rbind(simcoe_mean_corr, simcoe_cc_corr)

corr_aov <- aov(Mean_Correlation~Scale, data = simcoe_mean_corr)
summary(corr_aov)
TukeyHSD(corr_aov)


##### Create new data frame subset with categories for Mean Correlation for gg barplot #####
gg_simcoe_lc_corr <- simcoe_corr %>% select(mean_lc_corr) %>%
  summarize(CV = mean(mean_lc_corr),
            sd_stab = sd(mean_lc_corr),
            count = n(),
            se_stab = (sd_stab/(sqrt(count))),
            Scale2 = 'Mean Correlation',
            Scale = 'Local Community',
            Color = 'blue')
gg_simcoe_mp_corr <- simcoe_corr %>% select(mean_mp_corr) %>%
  summarize(CV = 1 - mean(mean_mp_corr),
            sd_stab = sd(mean_mp_corr),
            count = n(),
            se_stab = (sd_stab/(sqrt(count))),
            Scale2 = 'Mean Correlation',
            Scale = 'Metapopulation',
            Color = 'green')
gg_simcoe_cc_corr <- simcoe_corr %>% select(mean_cc_corr) %>%
  summarize(CV = mean(mean_cc_corr),
            sd_stab = sd(mean_cc_corr),
            count = n(),
            se_stab = (sd_stab/(sqrt(count))),
            Scale2 = 'Mean Correlation',
            Scale = 'Cross-Community',
            Color = 'red')

gg_simcoe_corr <- rbind(gg_simcoe_lc_corr, gg_simcoe_mp_corr)
gg_simcoe_corr <- rbind(gg_simcoe_corr, gg_simcoe_cc_corr)

gg_simcoe_corr$Scale <- factor(gg_simcoe_corr$Scale, levels = c("Local Community", 
                                                                "Metapopulation", 
                                                                "Cross-Community"))


gg_box_sim_corr <- ggplot(data = gg_simcoe_corr, aes(x = Scale2, y = CV, fill = Color)) +
  geom_bar(stat = 'identity',
           color = 'black',
           width = 0.5,
           position = position_dodge2()) +
  geom_errorbar(aes(ymin = CV - se_stab, ymax = CV + se_stab), 
                width = 0.5,
                position = position_dodge2()) +
  ylab("Average Annual Mean Correlation P(ik,jl)") +
  xlab("Scale") +
  ylim(0,1) +
  theme_classic() +
  theme(legend.position = 'none') +
  scale_fill_brewer(palette = 'Dark2')
gg_box_sim_corr


###### Create new data frame subset with categories for Mean Pop Pair CV for anova ##### 
simcoe_lc_popp <- simcoe_corr %>% select(mean_lc_pop_pair_cv) %>%
  summarize(Mean_Correlation = mean_lc_pop_pair_cv,
            Scale = 'Local Community')
simcoe_mp_popp <- simcoe_corr %>% select(mean_mp_pop_pair_cv) %>%
  summarize(Mean_Correlation = mean_mp_pop_pair_cv,
            Scale = 'Metapopulation')
simcoe_cc_popp <- simcoe_corr %>% select(mean_cc_pop_pair_cv) %>%
  summarize(Mean_Correlation = mean_cc_pop_pair_cv,
            Scale = 'Cross-Community')

simcoe_mean_popp <- rbind(simcoe_lc_popp, simcoe_mp_popp)
simcoe_mean_popp <- rbind(simcoe_mean_popp, simcoe_cc_popp)

popp_aov <- aov(Mean_Correlation~Scale, data = simcoe_mean_popp)
summary(popp_aov)
TukeyHSD(popp_aov)


##### Create new data frame subset with categories for Mean Correlation for gg barplot #####
gg_simcoe_lc_popp <- simcoe_corr %>% select(mean_lc_pop_pair_cv) %>%
  summarize(CV = mean(mean_lc_pop_pair_cv),
            sd_stab = sd(mean_lc_pop_pair_cv),
            count = n(),
            se_stab = (sd_stab/(sqrt(count))),
            Scale2 = 'Mean Correlation',
            Scale = 'Local Community',
            Color = 'blue')
gg_simcoe_mp_popp <- simcoe_corr %>% select(mean_mp_pop_pair_cv) %>%
  summarize(CV = mean(mean_mp_pop_pair_cv),
            sd_stab = sd(mean_mp_pop_pair_cv),
            count = n(),
            se_stab = (sd_stab/(sqrt(count))),
            Scale2 = 'Mean Correlation',
            Scale = 'Metapopulation',
            Color = 'green')
gg_simcoe_cc_popp <- simcoe_corr %>% select(mean_cc_pop_pair_cv) %>%
  summarize(CV = mean(mean_cc_pop_pair_cv),
            sd_stab = sd(mean_cc_pop_pair_cv),
            count = n(),
            se_stab = (sd_stab/(sqrt(count))),
            Scale2 = 'Mean Correlation',
            Scale = 'Cross-Community',
            Color = 'red')

gg_simcoe_popp <- rbind(gg_simcoe_lc_popp, gg_simcoe_mp_popp)
gg_simcoe_popp <- rbind(gg_simcoe_popp, gg_simcoe_cc_popp)

gg_simcoe_popp$Scale <- factor(gg_simcoe_popp$Scale, levels = c("Local Community", 
                                                                "Metapopulation", 
                                                                "Cross-Community"))


gg_box_sim_popp <- ggplot(data = gg_simcoe_popp, aes(x = Scale2, y = CV, fill = Color)) +
  geom_bar(stat = 'identity',
           color = 'black',
           width = 0.5,
           position = position_dodge2()) +
  geom_errorbar(aes(ymin = CV - se_stab, ymax = CV + se_stab), 
                width = 0.5,
                position = position_dodge2()) +
  ylab("Average Annual Mean Population Pair Variability (CVik X CVjl))") +
  xlab("Scale") +
  theme_classic() +
  theme(legend.position = 'none') +
  scale_fill_brewer(palette = 'Dark2')
gg_box_sim_popp





##### Create Variability & Stabilization Panel ##### 
corr_stab_pop <- ggarrange(gg_box_sim_pairs, gg_box_sim_stab_app, gg_box_sim_corr, gg_box_sim_popp,
                           labels = c("A", "B","C","D"),
                           ncol = 2, nrow = 2)
corr_stab_pop




