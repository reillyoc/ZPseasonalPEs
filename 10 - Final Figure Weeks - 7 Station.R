simcoe <- read.csv('simcoe_final_weeks_7.csv', header = T)
simcoe_div <- read.csv('simcoe_final_7_weeks_diversity', header = T)
simcoe_corr <- read.csv('simcoe_final_corr_weeks_7.csv', header = T)
simcoe_pairst <- read.csv('simcoe_final_pair_stab_weeks_7.csv', header = T)


library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

Year <- c(1986:1999)
Global_Change <- c('Pre')
simcoe_gc_pre <- data.frame(Year, Global_Change)

Year <- c(2008:2013)
Global_Change <- c('Post')
simcoe_gc_post <- data.frame(Year, Global_Change)

simcoe_gc <- rbind(simcoe_gc_pre, simcoe_gc_post)

simcoe <- merge(simcoe_gc, simcoe, by = c('Year'))

###### Create new data frame subset with categories for Variability for t.test #####
simcoe_gcv <- simcoe %>% select(gcv) %>%
  summarize(Variability = gcv,
            Scale = 'Metacommunity')
simcoe_lcv <- simcoe %>% select(lcv) %>%
  summarize(Variability = lcv,
            Scale = 'Weighted Average Population')

simcoe_var <- rbind(simcoe_gcv, simcoe_lcv)

t.test(Variability~Scale, data = simcoe_var)


##### Create new data frame subset with categories for Variability for gg barplot #####

gg_simcoe_mc_variability <- simcoe %>% select(gcv) %>%
  summarize(CV = mean(gcv),
            sd_variability = sd(gcv),
            count = n(),
            se_variability = (sd_variability/(sqrt(count))),
            Scale = 'Metacommunity',
            Color = 'blue')
gg_simcoe_wap_variability <- simcoe %>% select(lcv) %>%
  summarize(CV = mean(lcv),
            sd_variability = sd(lcv),
            count = n(),
            se_variability = (sd_variability/(sqrt(count))),
            Scale = 'Local Population',
            Color = 'green')

gg_simcoe_variability <- rbind(gg_simcoe_mc_variability, gg_simcoe_wap_variability)
gg_simcoe_variability$Scale <- factor(gg_simcoe_variability$Scale, levels = c("Metacommunity", 
                                                                              "Local Population"))

gg_box_sim_variability <- ggplot(data = gg_simcoe_variability, aes(x = Scale, y = CV, fill = Color)) +
  geom_bar(stat = 'identity',
           color = 'black',
           width = 0.3,
           position = position_dodge()) +
  geom_errorbar(aes(ymin = CV - se_variability, ymax = CV + se_variability), width = 0.1) +
  ylab("Average Annual Variability (CV Squared)") +
  xlab("Scale") +
  theme_classic() +
  theme(legend.position = 'none') +
  ylim(0, 2.3) +
  scale_fill_brewer(palette = 'Set1')
gg_box_sim_variability


###### Create new data frame subset with categories for Stabilization for anova ##### 
simcoe_W <- simcoe %>% select(W) %>%
  summarize(Stabilization = W,
            Scale = 'Total')
simcoe_lc_stab <- simcoe %>% select(lc_stab) %>%
  summarize(Stabilization = lc_stab,
            Scale = 'Local Community')
simcoe_mp_stab <- simcoe %>% select(mp_stab) %>%
  summarize(Stabilization = mp_stab,
            Scale = 'Metapopulation')
simcoe_cc_stab <- simcoe %>% select(cc_stab) %>%
  summarize(Stabilization = cc_stab,
            Scale = 'Cross-Community')

simcoe_stab <- rbind(simcoe_lc_stab, simcoe_mp_stab, simcoe_cc_stab)

stab_aov <- aov(Stabilization~Scale, data = simcoe_stab)
summary(stab_aov)
TukeyHSD(stab_aov)


##### Create new data frame subset with categories for Stabilization for gg barplot #####

gg_simcoe_mc_stab <- simcoe %>% select(W) %>%
  summarize(CV = mean(W),
            sd_stab = sd(W),
            count = n(),
            se_stab = (sd_stab/(sqrt(count))),
            Scale = 'Total',
            Color = 'blue')
gg_simcoe_lc_stab <- simcoe %>% select(lc_stab) %>%
  summarize(CV = mean(lc_stab),
            sd_stab = sd(lc_stab),
            count = n(),
            se_stab = (sd_stab/(sqrt(count))),
            Scale2 = 'Stabilization',
            Scale = 'Local Community',
            Color = 'blue')
gg_simcoe_mp_stab <- simcoe %>% select(mp_stab) %>%
  summarize(CV = mean(mp_stab),
            sd_stab = sd(mp_stab),
            count = n(),
            se_stab = (sd_stab/(sqrt(count))),
            Scale2 = 'Stabilization',
            Scale = 'Metapopulation',
            Color = 'green')
gg_simcoe_cc_stab <- simcoe %>% select(cc_stab) %>%
  summarize(CV = mean(cc_stab),
            sd_stab = sd(cc_stab),
            count = n(),
            se_stab = (sd_stab/(sqrt(count))),
            Scale2 = 'Stabilization',
            Scale = 'Cross-Community',
            Color = 'red')

gg_simcoe_stab <- rbind(gg_simcoe_mp_stab, gg_simcoe_lc_stab)
gg_simcoe_stab <- rbind(gg_simcoe_stab, gg_simcoe_cc_stab)

gg_simcoe_stab$Scale <- factor(gg_simcoe_stab$Scale, levels = c("Total", 
                                                                "Local Community", 
                                                                "Metapopulation", 
                                                                "Cross-Community"))


gg_box_sim_stab <- ggplot(data = gg_simcoe_stab, aes(x = Scale2, y = CV, fill = Color)) +
  geom_bar(stat = 'identity',
           color = 'black',
           width = 0.3) +
  ylab("Average Annual Stabilization") +
  xlab("Scale") +
  theme_classic() +
  ylim(0, 2.3) +
  theme(legend.position = 'none') +
  scale_fill_brewer(palette = 'Dark2')
gg_box_sim_stab




##### Create Variability & Stabilization Panel ##### 
stab_asynch <- ggarrange(gg_box_sim_variability, gg_box_sim_stab,
                         labels = c("A", "B"),
                         ncol = 2, nrow = 1)
stab_asynch



###### Shannon Diversity vs Asynchrony Plots #####
#Asynchrony vs Shannon diversity
lm_shan_a1 <- lm(data = simcoe, total_asynchrony~Shannon_Div_Local)
summary(lm_shan_a1)

lm_shan_a2 <- lm(data = simcoe, lc_asynchrony~Shannon_Div_Local)
summary(lm_shan_a2)

lm_shan_a4 <- lm(data = simcoe, mp_asynchrony~Shannon_Div_Local)
summary(lm_shan_a4)

lm_shan_a5 <- lm(data = simcoe, cc_asynchrony~Shannon_Div_Local)
summary(lm_shan_a5)


#Shannon Diversity vs Asynchrony
gg_sim_ham_26_shan_asynch <- ggplot() +
  geom_point(data = simcoe, aes(x = Shannon_Div_Local, y = total_asynchrony, color = "purple")) +
  geom_smooth(data = simcoe, aes(x = Shannon_Div_Local, y = total_asynchrony, color = "purple"), se = F, method = "lm", linetype = 'dashed') +
  geom_point(data = simcoe, aes(x = Shannon_Div_Local, y = mp_asynchrony, color = "orange")) +
  geom_smooth(data = simcoe, aes(x = Shannon_Div_Local, y = mp_asynchrony, color = "orange"), se = F, method = "lm", linetype = 'dashed') +
  geom_point(data = simcoe, aes(x = Shannon_Div_Local, y = lc_asynchrony, color = "green")) +
  geom_smooth(data = simcoe, aes(x = Shannon_Div_Local, y = lc_asynchrony, color = "green"), se = F, method = "lm", linetype = 'dashed') +
  geom_point(data = simcoe, aes(x = Shannon_Div_Local, y = cc_asynchrony, color = "pink")) +
  geom_smooth(data = simcoe, aes(x = Shannon_Div_Local, y = cc_asynchrony, color = "pink"), se = F, method = "lm", linetype = 'dashed') +
  xlab('Shannon Diversity') +
  ylab('Asynchrony') +
  theme_classic() +
  scale_color_brewer(palette = 'Dark2') +
  theme(legend.position = 'none')
gg_sim_ham_26_shan_asynch


##### Global Change Panel #####


##### Create new data frame subset with categories for Variability for gg barplot #####

gg_simcoe_mc_variability_invasive <- simcoe %>% select(gcv, Global_Change) %>%
  group_by(Global_Change) %>%
  summarize(CV = mean(gcv),
            sd_variability = sd(gcv),
            count = n(),
            se_variability = (sd_variability/(sqrt(count))),
            Color = 'Metacommunity')

gg_simcoe_wap_variability_invasive <- simcoe %>% select(lcv, Global_Change) %>%
  group_by(Global_Change) %>%
  summarize(CV = mean(lcv),
            sd_variability = sd(lcv),
            count = n(),
            se_variability = (sd_variability/(sqrt(count))),
            Color = 'Local_Population')


gg_simcoe_var_stab_invasive <- rbind(gg_simcoe_mc_variability_invasive , gg_simcoe_wap_variability_invasive)


gg_simcoe_var_stab_invasive$Global_Change <- factor(gg_simcoe_var_stab_invasive$Global_Change, 
                                                    levels = c("Pre", 
                                                               "Post"))

gg_simcoe_var_stab_invasive$Color <- factor(gg_simcoe_var_stab_invasive$Color, 
                                            levels = c("Metacommunity",
                                                       "Local_Population"))

gg_box_gc_var <- ggplot(data = gg_simcoe_var_stab_invasive, aes(x = Global_Change, y = CV, fill = Color)) +
  geom_bar(stat = 'identity',
           color = 'black',
           width = 0.5,
           position = position_dodge2()) +
  geom_errorbar(aes(ymin = CV - se_variability, ymax = CV + se_variability), 
                width = 0.5,
                position = position_dodge2()) +
  ylab("Average Annual Variability (CV Squared)") +
  xlab("Global Change") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = 'Set1')
gg_box_gc_var


##### Create new data frame subset with categories for Stabilization for gg barplot #####
gg_simcoe_total_stab <- simcoe %>% select(W, Global_Change) %>%
  group_by(Global_Change) %>%
  summarize(CV = mean(W),
            sd_stab = sd(W),
            count = n(),
            se_stab = (sd_stab/(sqrt(count))),
            Scale2 = 'Stabilization',
            Scale = 'Total',
            Color = 'blue')
gg_simcoe_lc_stab <- simcoe %>% select(lc_stab, Global_Change) %>%
  group_by(Global_Change) %>%
  summarize(CV = mean(lc_stab),
            sd_stab = sd(lc_stab),
            count = n(),
            se_stab = (sd_stab/(sqrt(count))),
            Scale2 = 'Stabilization',
            Scale = 'Local Community',
            Color = 'blue')
gg_simcoe_mp_stab <- simcoe %>% select(mp_stab, Global_Change) %>%
  group_by(Global_Change) %>%
  summarize(CV = mean(mp_stab),
            sd_stab = sd(mp_stab),
            count = n(),
            se_stab = (sd_stab/(sqrt(count))),
            Scale2 = 'Stabilization',
            Scale = 'Metapopulation',
            Color = 'green')
gg_simcoe_cc_stab <- simcoe %>% select(cc_stab, Global_Change) %>%
  group_by(Global_Change) %>%
  summarize(CV = mean(cc_stab),
            sd_stab = sd(cc_stab),
            count = n(),
            se_stab = (sd_stab/(sqrt(count))),
            Scale2 = 'Stabilization',
            Scale = 'Cross-Community',
            Color = 'red')

gg_simcoe_stab <- rbind(gg_simcoe_mp_stab, gg_simcoe_lc_stab)
gg_simcoe_stab <- rbind(gg_simcoe_stab, gg_simcoe_cc_stab)

gg_simcoe_stab$Scale <- factor(gg_simcoe_stab$Scale, levels = c( 
  "Local Community", 
  "Metapopulation", 
  "Cross-Community"))

gg_simcoe_stab$Global_Change <- factor(gg_simcoe_stab$Global_Change, levels = c(
  "Pre", 
  "Post"))


gg_box_sim_stab <- ggplot(data = gg_simcoe_stab, aes(x = Global_Change, y = CV, fill = Color)) +
  geom_bar(stat = 'identity',
           color = 'black',
           width = 0.3) +
  ylab("Average Annual Stabilization") +
  xlab("Scale") +
  theme_classic() +
  ylim(0, 1.4) +
  theme(legend.position = 'none') +
  scale_fill_brewer(palette = 'Dark2')
gg_box_sim_stab


##### Create new data frame subset with categories for Shannon Diversity for gg barplot #####

gg_simcoe_shan <- simcoe %>% select(Shannon_Div_Local, Global_Change) %>%
  group_by(Global_Change) %>%
  summarize(CV = mean(Shannon_Div_Local),
            sd_stab = sd(Shannon_Div_Local),
            count = n(),
            se_shan = (sd_stab/(sqrt(count))),
            Color = 'blue')

gg_simcoe_shan$Global_Change <- factor(gg_simcoe_shan$Global_Change, levels = c(
  "Pre", 
  "Post"))

gg_box_gc_shan_asynch <- ggplot(data = gg_simcoe_shan, aes(x = Global_Change, y = CV, fill = Color)) +
  geom_bar(stat = 'identity',
           color = 'black',
           width = 0.2,
           position = position_dodge2()) +
  geom_errorbar(aes(ymin = CV - se_shan, ymax = CV + se_shan), 
                width = 0.2,
                position = position_dodge2()) +
  ylab("Average Annual Shannon Diversity") +
  xlab("Global Change") +
  theme_classic() +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(1.5, 2.3)) +
  scale_fill_brewer(palette = 'Paired')
gg_box_gc_shan_asynch

##### Create new data frame subset with categories for Metacommunity Mean Density for gg barplot #####

gg_simcoe_mc_mean_dens <- simcoe %>% select(mc_sum_mean_density, Global_Change) %>%
  group_by(Global_Change) %>%
  summarize(CV = mean(mc_sum_mean_density),
            sd_variability = sd(mc_sum_mean_density),
            count = n(),
            se_variability = (sd_variability/(sqrt(count))),
            Color = 'Metacommunity')

gg_simcoe_mc_mean_dens$Global_Change <- factor(gg_simcoe_mc_mean_dens$Global_Change, 
                                               levels = c("Pre", 
                                                          "Post"))


gg_box_mc_mean <- ggplot(data = gg_simcoe_mc_mean_dens, aes(x = Global_Change, y = CV, fill = Color)) +
  geom_bar(stat = 'identity',
           color = 'black',
           width = 0.2,
           position = position_dodge2()) +
  geom_errorbar(aes(ymin = CV - se_variability, ymax = CV + se_variability), 
                width = 0.2,
                position = position_dodge2()) +
  ylab("Average Annual Summed Mean Density") +
  xlab("Global Change") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = 'Paired')
gg_box_mc_mean



##### Create 4-panel Global Change ##### 
gc_var_asynch <- ggarrange(gg_box_gc_var, gg_box_sim_stab,
                           gg_box_mc_mean, gg_box_gc_shan_asynch,
                           labels = c("A", "B",
                                      "C", "D"),
                           ncol = 2, nrow = 2)
gc_var_asynch


##### Global Change T-Tests ##### 
#Variability vs Global Change
t.test(gcv~Global_Change, data = simcoe)
t.test(lcv~Global_Change, data = simcoe)

#Stabilization vs Global Change
t.test(W~Global_Change, data = simcoe)
t.test(lc_stab~Global_Change, data = simcoe)
t.test(mp_stab~Global_Change, data = simcoe)
t.test(cc_stab~Global_Change, data = simcoe)

#Shannon Diversity vs Global Change
t.test(Shannon_Div_Local~Global_Change, data = simcoe)

#Evenness vs Global Change
t.test(mc_sum_mean_density~Global_Change, data = simcoe)



