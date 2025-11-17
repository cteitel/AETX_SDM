library(terra); library(sf)
library(dplyr)
library(stringr)
library(tidyr); library(purrr)

#get files from readme
ee_files <- read.table("data/spatial/EarthEnv_ReadMe.txt", sep = "\n")
start <- which(str_detect(ee_files[,1], "netCDF file"))
end <- nrow(ee_files)-1
ee_files <- ee_files[start:end,]
ee_files <- ee_files[-2]
ee_files <- map(ee_files, ~str_split_fixed(.x, "[|]", 10) %>%
                     str_remove_all("[|]|\\[|\\]") %>%
                  trimws() %>% t() %>% as.data.frame())
ee_cols <- ee_files[[1]]
ee_files <- ee_files[-1]
ee_files <- bind_rows(ee_files)
names(ee_files) <- ee_cols
ee_files <- ee_files[,-1]


#layers from Gerrin report
gerrin_layers <- ee_files %>%
  filter(`Variable explanation` %in% c("Cation exchange capacity",
                                       "Bioclim 16",
                                       "Bioclim 3",
                                       "Bioclim 5",
                                       "Bioclim 2",
                                       "Soil organic carbon")) %>%
  filter(!str_detect(`netCDF file / Variable name`, "weighted|range|minimum|maximum"))

#download files
file_names <- unique(gerrin_layers$`netCDF file / Variable name`)
urls <- paste0("https://data.earthenv.org/streams/", file_names)
file_dests <- paste0("data/spatial/raw/",file_names) %>% str_replace("[+]","_")
dll_info_hydro <- bind_cols(url = urls, file_dest = file_dests)
dll_info_hydro <- filter(dll_info_hydro)
getOption('timeout')
options(timeout = 200)
for(i in 1:nrow(dll_info_hydro)){
  if(!file.exists(dll_info_hydro$file_dest[i])){
    tryCatch(download.file(dll_info_hydro$url[i], dll_info_hydro$file_dest[i], mode = "wb"),
             error = function(e) NULL)
  }
}

#land cover
lc_classes <- c(
  "1"="Evergreen/Deciduous Needleleaf Trees",
  "2"="Evergreen Broadleaf Trees",
  "3"="Deciduous Broadleaf Trees",
  "4"="Mixed/Other Trees",
  "5"="Shrubs",
  "6"="Herbaceous Vegetation",
  "7"="Cultivated and Managed Vegetation",
  "8"="Regularly Flooded Vegetation",
  "9"="Urban/Built-up",
  "10"="Snow/Ice",
  "11"="Barren",
  "12"="Open Water")
gerrin_classes <- c("Evergreen/Deciduous Needleleaf Trees","Evergreen Broadleaf Trees",
                    "Urban/Built-up","Deciduous Broadleaf Trees","Herbaceous Vegetation",
                    "Cultivated and Managed Vegetation","Mixed/Other Trees")
classes_used <- lc_classes[lc_classes %in% gerrin_classes]
urls_lc <- paste0("https://data.earthenv.org/consensus_landcover/with_DISCover/consensus_full_class_", names(classes_used),".tif")
file_dests_lc <- paste0("data/spatial/raw/landcover_",names(classes_used),".tif")
dll_info_lc <- bind_cols(url = urls_lc, file_dest = file_dests_lc)
for(i in 1:nrow(dll_info_lc)){
  if(!file.exists(dll_info_lc$file_dest[i])){
    tryCatch(download.file(dll_info_lc$url[i], dll_info_lc$file_dest[i], mode = "wb"),
             error = function(e) NULL)
  }
}

#crop to study area
atex_dat <- read.csv("data/raw/Ah_2000-2024OCT.csv")
atex_sf <- st_as_sf(atex_dat, coords = c("Longitude","Latitude"), crs = "epsg:4326")
atex_vect <- vect(atex_sf)
atex_ext <- c(range(atex_dat$Longitude), range(atex_dat$Latitude))
names(atex_ext) <- c("xmin","xmax","ymin","ymax")

rst_lc <- rast(lapply(dll_info_lc$file_dest, rast))
rst_hydro <- rast(lapply(dll_info_hydro$file_dest, rast))
gerrin_layers <- gerrin_layers %>% mutate(variable_code_lyr = str_replace(`Variable code`, "_0", "_variable=_"))
rst_hydro <- rst_hydro[[gerrin_layers$variable_code_lyr[gerrin_layers$variable_code_lyr %in% names(rst_hydro)]]]

rst_lc <- crop(rst_lc, atex_vect)
atex_vect_proj <- project(atex_vect, rst_hydro)
rst_hydro <- crop(rst_hydro, atex_vect_proj)

atex_vals_lc <- terra::extract(rst_lc, atex_vect)

#extract values from the closest non-NA pixel
atex_vect_proj_buff <- buffer(atex_vect_proj, 10000)
atex_vals_hydro <- terra::extract(rst_hydro, atex_vect_proj_buff, weights = T)
atex_vals_hydro <- atex_vals_hydro %>% filter_all(~!is.na(.x)) %>%
  group_by(ID) %>%
  summarize_at(vars(-weight), ~.x[which.max(weight)], .groups = "drop")

#another way to do it, very slow
# rst_hydro_cells <- as.polygons(rst_hydro[[1]], values = F, na.rm = T)
# dists <- distance(rst_hydro_cells, atex_vect_proj)
# # hist(apply(dists,2,min))
# ids <- apply(dists, 2, which.min)
# extract_cells <- rst_hydro_cells[ids,]
# extract_cells$ID <- 1:nrow(extract_cells)
# extract_cells <- centroids(extract_cells)
# atex_vals_hydro <- terra::extract(rst_hydro, extract_cells)


atex_vals <- left_join(atex_vals_lc, atex_vals_hydro)
new_names <- tibble::enframe(classes_used, name = "name", value = "new_name") %>%
  mutate(name = paste0("landcover_", name)) %>%
  bind_rows(select(gerrin_layers, name = variable_code_lyr, new_name = `Variable explanation`)) %>%
  mutate(new_name = str_to_lower(new_name) %>% str_replace_all("/|_| ", "_")) %>%
  bind_rows(bind_cols(name = "ID", new_name = "ID")) %>% tibble::deframe()
atex_vals <- setNames(atex_vals, unname(new_names[names(atex_vals)]))
atex_dat_annot <- bind_cols(atex_dat, atex_vals)

write.csv(atex_dat_annot, "data/clean/Ah_2000-2024OCT_annot.csv", row.names=F)




