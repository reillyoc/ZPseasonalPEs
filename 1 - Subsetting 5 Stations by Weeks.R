sim_start <- read.csv("Simcoe_Zooplankton.csv", header = T)

library(tidyverse)
library(reshape2)


##### Preparing the Dataframe #####

#rename column headers
sim_start <- sim_start %>%
  summarize(Date = ZDATE,
            Station_ID = STN,
            Species = SPECIES,
            Density = DENSITY....m3.)

#Pull out year and month columns from date
sim_start$Year <- format(as.Date(sim_start$Date,format="%Y-%m-%d"),"%Y")
sim_start$Month <- format(as.Date(sim_start$Date,format="%Y-%m-%d"),"%m")


#Subset .csv to required columns
sim_start <- sim_start %>% select(Year, Month, Density, Date, Species, Station_ID)

#Species Re-Defined According to Joelle Species List
daph_recomb <- subset(sim_start, Species =='DAPHNIA (DAPHNIA) PULICARIA'|
                        Species == 'DAPHNIA (DAPHNIA) CATAWBA')

daph_recomb <- daph_recomb %>% group_by(Year, Month, Date, Station_ID) %>%
  summarize(Density = sum(Density))

sim_start <- sim_start %>%
  filter(! Species == 'DAPHNIA (DAPHNIA) PULICARIA'|
           Species == 'DAPHNIA (DAPHNIA) CATAWBA')

daph_recomb <- daph_recomb %>% mutate(Species = 'DAPHNIA (DAPHNIA) PULICARIA')

sim_start <- rbind(sim_start, daph_recomb)

sim_start <- subset(sim_start,
                    Station_ID == 'C6'|
                      Station_ID == 'C9'|
                      Station_ID == 'K42'|
                      Station_ID == 'K45'|
                      Station_ID == 'S15')

#Subset years - remove 2000 to 2007
sim_start <- filter(sim_start, Year < 2000 | Year > 2007)

#Create a column of weeks
library(aweek)
set_week_start("Sunday")
sim_start$Week <- date2week(sim_start$Date, numeric = T)

detach("package:aweek", unload=TRUE)

#Reformat weeks to 1:n per year as time steps and remove weeks that have not all stations sampled
sim_start <- sim_start %>% select(Year, Week, Density, Station_ID, Species)

sim_start <- sim_start %>% group_by(Year, Week) %>%
  mutate(station_week = length(unique(Station_ID)))

sim_start <- sim_start %>%
  filter(! Species == 'CALANOID COPEPODID')
sim_start <- sim_start %>%
  filter(! Species == 'CALANOID NAUPLIUS')
sim_start <- sim_start %>%
  filter(! Species == 'CYCLOPOID COPEPODID')
sim_start <- sim_start %>%
  filter(! Species == 'CYCLOPOID NAUPLIUS')
sim_start <- sim_start %>%
  filter(! Species == 'DREISSENA POLYMORPHA')
sim_start <- sim_start %>%
  filter(! Species == 'BYTHOTREPHES CEDERSTROEMI')
sim_start <- sim_start %>%
  filter(! Species == '#N/A')


write.csv(sim_start, "sim_start_weeks_manual_subset_5.csv")
