# Santa Rita Experimental Range
# Longterm transect and point data clean and join
# W. Austin Rutherford
# email:arutherford@email.arizona.edu
# 2020-12-22

# load packages
library(tidyverse)
library(ggpubr)
library(Metrics)

# read in data (cover and transects)
cover <- read_csv('data/ltcover_2018July10.csv')

trans_cols <- c("pasttran", "pasture", "transect", "utm_x", "utm_y")

trans <- read_csv('data/lttranutm.csv', skip = 1, col_names = trans_cols)

glimpse(trans)
glimpse(cover)


# select cover data for 2015
cover_2015 <- cover %>% 
  mutate(species = as.factor(SPECIES),
         ecosite = as.factor(ECOSITE),
         pasture = as.factor(PASTURE),
         transect = as.factor(TRANSECT)) %>% 
  # only pull out 2015 transects to match model
  select(pasture, transect, ecosite, species, YR2006, YR2009, YR2012, YR2015) %>%
  # filter to only woody species used in SRER long term data
  # https://cals.arizona.edu/SRER/species/species_code_translation.txt for species codes
  # Woody species match woody species classifications in "Figures for cover values on 100 transects through 2018" spreadsheet
  # https://cals.arizona.edu/SRER/data.html
  filter(species %in% c('ACAN',
                     'ACCO',
                     'ACGR',
                     'AGAV',
                     'ALWR',
                     'ANTH',
                     'APTE',
                     'ATCA',
                     'ATCO',
                     'BABR',
                     'BACC',
                     'BARA',
                     'BASA',
                     'CAAR',
                     'CAER',
                     'CEFL',
                     'CEMI',
                     'CEPA',
                     'ENFR',
                     'EPTR',
                     'ERIO',
                     'ERWR',
                     'FOSP',
                     'GUSA',
                     'HALA',
                     'HASP',
                     'KRGR',
                     'KRPA',
                     'LATR',
                     'LYAN',
                     'LYCI',
                     'LYPA',
                     'LYTO',
                     'MIBI',
                     'MIDY',
                     'OTHSHRUB',
                     'PRVE',
                     'YUEL',
                     'ZIOB',
                     'ZIPU')) %>%
  # change -9999 values to NAs
  mutate(YR2015 = na_if(YR2015, -9999)) %>% 
  mutate(YR2006 = na_if(YR2006, -9999)) %>%
  mutate(YR2009 = na_if(YR2009, -9999)) %>%
  mutate(YR2012 = na_if(YR2012, -9999)) %>%
  drop_na(YR2006, YR2009, YR2012, YR2015)

# create a transect ID to match with transect UTMs
cover_sum <- cover_2015 %>% 
  group_by(pasture, transect) %>% 
  summarise(cover15 = sum(YR2015),
            cover12 = sum(YR2012),
            cover09 = sum(YR2009),
            cover06 = sum(YR2006)) %>% 
  mutate(transect = str_remove(transect, "^0+")) %>% 
  unite("tranID",pasture:transect, sep = "_", remove = FALSE)

# create a transect ID to match with cover data
trans_clean <- trans %>% 
  select(-pasttran) %>% 
  unite("tranID",pasture:transect, sep = "_", remove = FALSE)

# join together transect and cover data to map in ArcPro
transect_cover <- left_join(cover_sum, trans_clean, by = "tranID") %>% 
  select(tranID, pasture.x, transect.x, cover15, cover12, cover09, cover06, utm_x, utm_y) %>% 
  mutate(Pasture = as.factor(pasture.x),
         Transect = as.factor(transect.x),
         tranID = as.factor(tranID)) %>% 
  select(tranID, Pasture, Transect, cover15, cover12, cover09, cover06, utm_x, utm_y) %>% 
  # cover calculations from .1 ft to meters: https://cals.arizona.edu/SRER/longterm/longtermstudy.txt
  mutate(cover15 = cover15/10) %>% # cover in 0.1ft to ft, along 100ft transect = % cover
  mutate(cover15 = cover15*0.3048) %>% # convert to meters
  mutate(cover15 = cover15/30.48)  %>%  # covert to model cover (not in percent) 
  mutate(cover12 = cover12/10) %>% # cover in 0.1ft to ft, along 100ft transect = % cover
  mutate(cover12 = cover12*0.3048) %>% # convert to meters
  mutate(cover12 = cover12/30.48) %>% # covert to model cover (not in percent)
  mutate(cover09 = cover09/10) %>% # cover in 0.1ft to ft, along 100ft transect = % cover
  mutate(cover09 = cover09*0.3048) %>% # convert to meters
  mutate(cover09 = cover09/30.48) %>% # covert to model cover (not in percent)
  mutate(cover06 = cover06/10) %>% # cover in 0.1ft to ft, along 100ft transect = % cover
  mutate(cover06 = cover06*0.3048) %>% # convert to meters
  mutate(cover06 = cover06/30.48) # covert to model cover (not in percent)


transect_cover_long <- transect_cover %>% pivot_longer(cover15:cover06, names_to = "year", values_to = "cover") %>% 
  mutate(year = as.factor(year)) %>% 
  dplyr::select(-utm_x, -utm_y)

# write final table to csv
#write_csv(transect_cover_long, file = "data/srer_transect_cover_long-term.csv")

# read in Random Forest model samples
rf_names <- c("point", "x", "y", "rf_cover")

rf_samps <- read_csv('data/2015_wc_rf_samples.dbf.csv', skip = 1, col_names = rf_names)

# make a row number ID
rf_row <- rf_samps %>% 
  rowid_to_column("rowID")

# make a row number ID to join predicted values to transect values
trans_row <- transect_cover %>% 
  rowid_to_column("rowID")

# join tables by row number
cover_comb <- left_join(trans_row, rf_row, by = "rowID") %>% 
  rowwise()

hist(cover_comb$cover15)
hist(cover_comb$rf_cover)

# graph of ground vs model wc values
cover_plot <- cover_comb %>%
  ggplot(aes(y = 100*rf_cover, x = 100*cover15)) +
  geom_point(size = 2, alpha = 0.5) +
  labs(x = 'Observed Woody Cover (%)', y = 'Predicted Woody Cover (%)') +
  ylim(c(0, 60))+
  xlim(c(0, 60))+
  scale_color_manual(values = c("gray80", "darkred"))+
  #geom_abline(intercept = 0, slope = 1, color = "blue", size = 1)+
  stat_smooth(method = "lm", formula = y ~ x, color = "red")+
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label =  paste(stat(rr.label), stat(p.value.label), sep = "*\", \"*")),
                        parse = TRUE)+
  theme_pubr()+
  labs_pubr(base_size = 18)

cover_plot

ggsave(filename = "2015_cover_compare.tiff",
       plot = cover_plot,
       dpi = 800,
       width = 11,
       height = 8,
       units = "in",
       compression = "lzw")

summary(lm(rf_cover~cover15, data = cover_comb))

rmse(cover_comb$cover15, cover_comb$rf_cover)

cor.test(cover_comb$cover15, cover_comb$rf_cover, method = c("pearson"))

mae(cover_comb$cover15, cover_comb$rf_cover)

# look at mean SRER transect cover (0.22)
summary(cover_comb$srer_mean)

# look at mean Random Forest cover (0.289)
summary(cover_comb$rf_cover)

# SRER transect cover is lower as expected comparing a 100ft x 1ft belt transect
# vs a 30m2 cell in the random forest