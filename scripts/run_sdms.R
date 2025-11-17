library(dplyr)
library(dismo)
library(mgcv)
library(glmnet)
library(caret)


atex_dat <- read.csv("data/clean/Ah_2000-2024OCT_annot.csv") %>%
  mutate(presence = as.integer(Ah.Status == "Present")) %>%
  mutate(Longitude = ifelse(Longitude > 0, (-1)*Longitude, Longitude))

library(ggplot2)
usa <- rnaturalearth::ne_states(country = "United States of America")
ggplot(atex_dat, aes(x = Longitude, y = Latitude, color = Ah.Status)) +
  geom_sf(data = usa, inherit.aes = F) +
  geom_point() +
  coord_sf(xlim = range(atex_dat$Longitude), ylim = range(atex_dat$Latitude))

atex_dat <- filter_all(atex_dat, ~!is.na(.x))

atex_x <- dplyr::select(atex_dat, evergreen_deciduous_needleleaf_trees:cation_exchange_capacity) %>%
  as.matrix()
atex_y <- atex_dat$presence

lasso <- glmnet::glmnet(x = atex_x, y = atex_y, family = "binomial")
maxent <- dismo::maxent(x = as.data.frame(atex_x), p = atex_y, silent = F)
