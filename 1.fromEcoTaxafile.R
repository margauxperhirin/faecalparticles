#' Script to associate images to their respective plankton ID, get the images,
#' trim them and save the new ones in a new file, then save the
#' morphological features and transform them

remove(list = ls())
set.seed(121999)

# Libraries
library("tidyverse") # Data manipulation
library("morphr") # For morphological descriptors and images manipulation
library('plyr')

# Data frames
allplkt <- read.table("data/ecotaxa_export.tsv", header = TRUE, sep = "\t")
taxon_ref <- read_csv("data/taxon_ecotaxa.csv", col_names = T)

# 1. Check the dataset to get the good profiles ----
## Compute percent validated per sample
valid <- allplkt %>%
  group_by(sample_id) %>%
  dplyr::summarise(percent = sum(object_annotation_status == "validated")/n())
table(valid$percent) # -> 166 fully validated samples

## Keep only fully validated samples
fully_valid <- filter(valid, percent == 1)
plk <- allplkt %>%
  filter(sample_id %in% fully_valid$sample_id) %>%
  select(-object_annotation_status)
nrow(plk) # -> 1,233,938 validated objects
remove(valid)

## Check users that validated images
plk %>%
  group_by(object_annotation_person_name) %>%
  dplyr::summarise(validation = n()/nrow(plk)) %>%
  ungroup()

## Eliminate points to the south
plk <- plk %>%
  filter(object_lat > 60)

# 2. Define stations ----
## Check if all depth_min = depth_max, and create a new column
plk <- plk %>%
  mutate(depth = object_depth_max) %>% # - object_depth_min) %>%
  # select(-object_depth_max, -object_depth_min) # -> Yes, so no need to redo that, but keep only 1 column
  select(-object_depth_max, -object_depth_min)

## Reformat station name
plk <- plk %>%
  mutate(station = str_replace(sample_id, "ge_2016_", "")) %>%
  # and sort data to be clean
  arrange(station, depth)

## Compute stations
sta <- plk %>%
  distinct(sample_id, station, object_lat, object_lon, object_date, object_time) %>%
  mutate(date = parse_datetime(as.character(object_date), format = "%Y%m%d")) %>%
  select(-object_date)

## Check
summary(sta)
sta %>%
  ggplot() +
  geom_point(aes(object_lon, object_lat)) +
  coord_quickmap()
# -> check date (one in 1987...) : done

remove(sta)

# 3. Store plankton objects ----
## Rename and/or eliminate groups
plk_det_simplified <- plk %>%
  left_join(taxon_ref, by = c("object_annotation_category" = "taxon")) %>%
  filter(new_taxon != "NA")

plk_simplified <- plk_det_simplified %>%
  filter(new_taxon != "darksphere") %>%
  filter(new_taxon != "fiber") %>%
  filter(new_taxon != "detritus")
# write_csv(plk_simplified, file = "data/plankton.csv.gz") # Save it for plankton analyses

det_simplified <- plk_det_simplified %>%
  filter(new_taxon == "darksphere" | new_taxon == "fiber" | new_taxon == "detritus")
write_csv(det_simplified, file = "data/detritus.csv.gz") # Save it for detritus analyses
remove(det_simplified, plk_det_simplified)

