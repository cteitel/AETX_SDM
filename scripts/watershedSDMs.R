# Pull in data and run models?

######## load packages ############--------------------------------------------------
library(tidyverse)
library(dismo)
library(terra)
library(sf)
library(here)

######## load data ############--------------------------------------------------
pres_abs <- read_csv("data/clean/extracted_data.csv")
correlations <- read_csv("data/clean/envdata_correlations.csv")
bad_pairs <- filter(correlations, abs(value) >= 0.7)

#convert presence/absence data to spatial format
pres_abs_vect <- vect(pres_abs, geom = c("Longitude", "Latitude"), crs = "epsg:4326")

#read in Watershed Boundary Dataset
# https://www.usgs.gov/national-hydrography/access-national-hydrography-products
huc_dir <- "data/spatial/WBD"
huc10_files <- here(huc_dir, dir(huc_dir), "Shape/WBDHU10.shp")
huc10 <- map(huc10_files, vect)

huc10 <- vect(huc10)

######## convert to huc10 data ############--------------------------------------------------

#get HUC10s for presence/absence points
huc10s <- extract(huc10, pres_abs_vect)
pres_abs$huc10 <- huc10s$huc10
pres_abs$huc10_name <- huc10s$name

#filter HUC10 dataset to sampled watersheds only
huc10_samp <- huc10[huc10$huc10 %in% pres_abs$huc10,]
huc10_samp_df <- as.data.frame(huc10_samp)

# Summarize and plot some data
n_status <- pres_abs %>% group_by(huc10) %>%
  summarize(n_tot = n(),
            n_pres = sum(Ah.Status == "Present"), n_abs = sum(Ah.Status == "Absent")) %>%
  mutate(status = ifelse(n_pres > 0 & n_abs == 0, "Present",
                         ifelse(n_pres == 0 & n_abs > 0, "Absent", "Mixed")))
filter(n_status, n_pres > 0 & n_abs > 0)
hist(n_status$n_tot)

huc10_samp %>%
  st_as_sf() %>%
  left_join(n_status) %>%
  ggplot() +
  geom_sf(aes(fill = factor(n_tot)), color = NA) +
  # geom_point(data = pres_abs, aes(x = Longitude, y = Latitude)) +
  scale_fill_viridis_d("Total samples (HUC10)") +
  theme_bw() +
  NULL

huc10_samp %>%
  st_as_sf() %>%
  left_join(n_status) %>%
  ggplot() +
  geom_sf(aes(fill = status), color = NA) +
  scale_fill_manual("AETX Status (HUC10)", breaks = c("Present","Absent","Mixed"),
                       values = c("#D41159","#1A85FF","gray")) +
  theme_bw() +
  NULL

######## extract covariates at huc10s ############--------------------------------------------------

lc_avg <- rast("data/spatial/landcover_average.nc")
clim_avg <- rast("data/spatial/hydroclim_average+sum.nc")
soil_avg <- rast("data/spatial/soil_average.nc")
elev <- rast("data/spatial/elevation.nc")
slope <- rast("data/spatial/slope.nc")
var_names <- read.csv("data/spatial/earthenv_data_variables.csv")

names(lc_avg) <- var_names$variable_code[var_names$file_name == "landcover_average.nc"]
names(clim_avg) <- var_names$variable_code[var_names$file_name == "hydroclim_average+sum.nc"]
names(soil_avg) <- var_names$variable_code[var_names$file_name == "soil_average.nc"]
names(elev) <- var_names$variable_code[var_names$file_name == "elevation.nc"]
names(slope) <- var_names$variable_code[var_names$file_name == "slope.nc"]

#crop and reproject landcover raster
lc_avg <- crop(lc_avg, project(huc10_samp, lc_avg))
lc_avg <- project(lc_avg, "epsg:4326")
#crop all other rasters
ext <- ext(huc10_samp)
clim_avg <- crop(clim_avg, ext)
soil_avg <- crop(soil_avg, ext)
elev <- crop(elev, ext)
slope <- crop(slope, ext)

rast_all <- rast(list(clim_avg, soil_avg, elev, slope))
values <- extract(lc_avg, huc10_samp)
values2 <- extract(rast_all, huc10_samp)

#eliminate NAs (terrestrial cells)
#freq(rast_all, value = NA) #check that all layers have the same frequency of NAs
#freq(lc_avg, value = NA)
values <- filter(values, !is.na(lc_avg_01))
values2 <- filter(values2, !is.na(hydro_avg_01))

huc10_samp_df <- huc10_samp_df %>% left_join(n_status, by = "huc10")
values$huc10 <- huc10_samp_df$huc10[values$ID]
values2$huc10 <- huc10_samp_df$huc10[values2$ID]

# histograms of landcover vars
values %>%
  left_join(huc10_samp_df %>% select(huc10, status)) %>%
  pivot_longer(cols = names(lc_avg), names_to = "variable_code", values_to = "value") %>%
  left_join(var_names) %>%
  mutate(variable_explanation = str_wrap(variable_explanation, 15)) %>%
  ggplot(aes(x = value, fill = status)) +
  geom_histogram() +
  facet_wrap(~variable_explanation, strip.position = "bottom") +
  scale_fill_manual("AETX Status (HUC10)", breaks = c("Present","Absent","Mixed"),
                    values = c("#D41159","#1A85FF","gray")) +
  theme_bw() +
  theme(axis.title.x = element_blank(), strip.background = element_blank(),
        strip.placement = "outside", strip.text = element_text(hjust = 0))

# histograms of bilclim vars
values2 %>%
  left_join(huc10_samp_df %>% select(huc10, status)) %>%
  pivot_longer(cols = names(clim_avg), names_to = "variable_code", values_to = "value") %>%
  left_join(var_names) %>%
  mutate(variable_explanation = str_wrap(variable_explanation, 15)) %>%
  ggplot(aes(x = value, fill = status)) +
  geom_histogram() +
  facet_wrap(~variable_explanation, strip.position = "bottom",
             scales = "free") +
  scale_fill_manual("AETX Status (HUC10)", breaks = c("Present","Absent","Mixed"),
                    values = c("#D41159","#1A85FF","gray")) +
  theme_bw() +
  theme(axis.title.x = element_blank(), strip.background = element_blank(),
        strip.placement = "outside", strip.text = element_text(hjust = 0))

# histograms of soil vars
values2 %>%
  left_join(huc10_samp_df %>% select(huc10, status)) %>%
  pivot_longer(cols = names(soil_avg), names_to = "variable_code", values_to = "value") %>%
  left_join(var_names) %>%
  mutate(variable_explanation = str_wrap(variable_explanation, 15)) %>%
  ggplot(aes(x = value, fill = status)) +
  geom_histogram() +
  facet_wrap(~variable_explanation, strip.position = "bottom",
             scales = "free_x") +
  scale_fill_manual("AETX Status (HUC10)", breaks = c("Present","Absent","Mixed"),
                    values = c("#D41159","#1A85FF","gray")) +
  theme_bw() +
  theme(axis.title.x = element_blank(), strip.background = element_blank(),
        strip.placement = "outside", strip.text = element_text(hjust = 0))

# histograms of soil vars
values2 %>%
  left_join(huc10_samp_df %>% select(huc10, status)) %>%
  pivot_longer(cols = names(soil_avg), names_to = "variable_code", values_to = "value") %>%
  left_join(var_names) %>%
  mutate(variable_explanation = str_wrap(variable_explanation, 15)) %>%
  ggplot(aes(x = value, fill = status)) +
  geom_histogram() +
  facet_wrap(~variable_explanation, strip.position = "bottom",
             scales = "free_x") +
  scale_fill_manual("AETX Status (HUC10)", breaks = c("Present","Absent","Mixed"),
                    values = c("#D41159","#1A85FF","gray")) +
  theme_bw() +
  theme(axis.title.x = element_blank(), strip.background = element_blank(),
        strip.placement = "outside", strip.text = element_text(hjust = 0))

# histograms of elevation vars
values2 %>%
  left_join(huc10_samp_df %>% select(huc10, status)) %>%
  pivot_longer(cols = names(elev), names_to = "variable_code", values_to = "value") %>%
  left_join(var_names) %>%
  ggplot(aes(x = value, fill = status)) +
  geom_histogram() +
  facet_wrap(~variable_explanation, strip.position = "bottom",
             scales = "free_x") +
  scale_fill_manual("AETX Status (HUC10)", breaks = c("Present","Absent","Mixed"),
                    values = c("#D41159","#1A85FF","gray")) +
  theme_bw() +
  theme(axis.title.x = element_blank(), strip.background = element_blank(),
        strip.placement = "outside")

#histograms of slope vars
values2 %>%
  left_join(huc10_samp_df %>% select(huc10, status)) %>%
  pivot_longer(cols = names(slope), names_to = "variable_code", values_to = "value") %>%
  left_join(var_names) %>%
  ggplot(aes(x = value, fill = status)) +
  geom_histogram() +
  facet_wrap(~variable_explanation, strip.position = "bottom",
             scales = "free_x") +
  scale_fill_manual("AETX Status (HUC10)", breaks = c("Present","Absent","Mixed"),
                    values = c("#D41159","#1A85FF","gray")) +
  theme_bw() +
  theme(axis.title.x = element_blank(), strip.background = element_blank(),
        strip.placement = "outside")

######## calculate watershed-level means ############--------------------------------------------------

means <- values %>%
  left_join(huc10_samp_df %>% select(huc10, status))  %>%
  group_by(huc10, status) %>%
  summarize_all(mean)

means2 <- values2 %>%
  left_join(huc10_samp_df %>% select(huc10, status))  %>%
  group_by(huc10, status) %>%
  summarize_all(mean)


means <- left_join(means, means2)

ct <- crds(centroids(huc10_samp)) %>%
  as.data.frame() %>% select(long_cent = x, lat_cent = y) %>%
  mutate(huc10 = huc10_samp$huc10)

means <- left_join(means, ct)

var_names <- mutate(var_names, variable_exp = variable_explanation %>%
                      str_to_lower() %>%
                      str_replace("maximum", "max") %>%
                      str_replace("minimum", "min") %>%
                      str_replace("average", "avg") %>%
                      str_replace("probability", "prob") %>%
                      str_replace("content", "cont") %>%
                      str_replace("vegetation", "veg") %>%
                      str_sub(1, 15) %>% trimws() %>% str_replace_all(" |/", "_"))
# n_distinct(var_names$variable_exp[var_names$variable_code %in% names(means)])
# n_distinct(var_names$variable_explanation[var_names$variable_code %in% names(means)])
tt <- deframe(select(var_names, variable_exp, variable_code))
tt <- tt[tt %in% names(means)]
means <- means %>%
  rename(all_of(tt))

means <- mutate(means, present = ifelse(status %in% c("Present", "Mixed"), 1, 0))
glm(present ~ evergreen_decid + mixed_other_tre + cultivated_and + deciduous_broad +
      bioclim_1 + bioclim_13 +
      cation_exchange +
      slope_range + max_elevation,
    data = means %>% filter(lat_cent < 38), na.action = "na.fail") %>%
  summary()


