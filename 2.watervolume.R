#' This script computes the watervolume sampled per 10-meters depth bin
#' at each station and save the dataframe

remove(list = ls())
set.seed(121999)

# Libraries
library("tidyverse") # Data manipulation
library("morphr")
library('plyr')

# Data frames
Laure <- read_csv("data/zoo_env_greenedge_final.csv")

# Keep only the column we need + mutate in the same depth bins (10 m)
vLaure <- Laure %>%
  select(station, depth_bin, watervolume) %>%
  distinct() %>%
  mutate(depth_bin10 = round_any(depth_bin, 10, floor)) %>%
  select(-depth_bin)

# Create the dataset per station and depth bin
vol_d <- vLaure %>%
  group_by(station, depth_bin10) %>%
  dplyr::summarise(totvol_10m_m3 = sum(watervolume, na.rm = T)*0.001)

# Save the dataframe
write_csv(vol_d, file = "data/volper10mbin.csv.gz")

# Look at the distribution
ggplot(vol_d) +
  geom_boxplot(aes(x = 1, y = totvol_10m_m3, group = 1)) +
  scale_y_log10() +
  theme_classic(base_size = 25) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

mean(vol_d$totvol_10m_m3)*1000 / 1.15 / 10 # => In average, 15 pictures are taken per meter
