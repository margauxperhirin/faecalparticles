#' This script compares faecal pellets' size from the 3 clusters, the 18 UVP5
#' images, sediment traps and the copepods incubations.
#' It produces Figure 1F.

remove(list = ls())
set.seed(121999)

# Libraries
library("tidyverse") # Data manipulation
library("ggrepel") # For plot
library("RColorBrewer") # Colours in plot
library('cowplot') # For multiple plotting
library('plyr') # For round_any function
library('morphr') # To move images files

# Data
env <- read_csv("data/env_stations.csv.gz") # To get mean_sic, sic, owd
clusters <- read_csv("kmeans2/detritus_37-clusters-kmeans.csv.gz") %>%
  filter(cluster %in% c(13,32,35)) %>%
  mutate(Type = ifelse(cluster == 13, 'Cluster 13',
                       ifelse(cluster == 32, 'Cluster 32', 'Cluster 35')),
         Length_mm = major_mm) %>%
  dplyr::select(Type, Length_mm) # Faecal pellets-like particles from the 3 clusters

st <- read_delim("data/Arctic_Driftingtraps25m.csv",
                 delim = ";", escape_double = FALSE, trim_ws = TRUE) # Sediment traps metadata_mp
st_size <- read_csv("data/faecalpellets_inddata.csv") # Sediment traps
inc_size <- read_csv("data/faecalpellets_incubations_mm.csv") # Incubation

poop <- read_csv("data/poopingcopepods_mm.csv") # 18 UVP5 images

# 1. Prepare the data for comparison ----
## 1.a) In sediment traps data ----
st$Trap <- seq(1,12,1)
st_size <- left_join(st_size, dplyr::select(st, -mission, -date_time_close, -PI))
rm(st)
st_size <- st_size[,-12]

st_size <- st_size %>%
  mutate(Type = 'Sediment traps')
summary(dplyr::select(st_size, Length_mm, Type))

## 1.b) From incubations ----
inc_size <- inc_size %>%
  mutate(Type = 'Incubations') # To compensate for the cm to mm error
summary(dplyr::select(inc_size, Length_mm, Type))

## 1.c) From 18 images ----
poop <- poop %>%
  mutate(Length_mm = major_mm_ecotaxa,
         Type = "18 UVP5 images")

## 1.d) Merge them ----
allFP <- rbind(dplyr::select(clusters, Length_mm, Type),
               dplyr::select(st_size, Length_mm, Type),
               dplyr::select(inc_size, Length_mm, Type),
               dplyr::select(poop, Length_mm, Type))

# 2. Statistical test ----
# Are variances homogeneous ?
car::leveneTest(Length_mm ~ Type, data = allFP) # => Nope, so Kruskal-Wallis

# Kruskal-Wallis test
kruskal.test(Length_mm ~ Type, data = allFP) # Significantly different

# Which side ? (Dunn)
rstatix::dunn_test(data = allFP,
                   formula = Length_mm ~ Type,
                   p.adjust.method = 'bonferroni') # All of them are significantly different from another

# 3. Plot ----
allFP$Type <- factor(allFP$Type,
                     levels = c('Cluster 13','Cluster 32','Cluster 35',
                                '18 UVP5 images','Incubations','Sediment traps')) # Order the variable

majorp <- ggplot(allFP) +
  geom_boxplot(aes(x = Type, group = Type, y = Length_mm, colour = Type),
               outlier.size = 0.7) +
  scale_y_log10('Length (mm)') +
  geom_text(aes(x = 1, y = 20, label = 'a'),
            size = 6, colour = 'black') +
  geom_text(aes(x = 2, y = 20, label = 'b'),
            size = 6, colour = 'black') +
  geom_text(aes(x = 3, y = 20, label = 'c'),
            size = 6, colour = 'black') +
  geom_text(aes(x = 4, y = 20, label = 'cd'),
            size = 6, colour = 'black') +
  geom_text(aes(x = 5, y = 20, label = 'd'),
            size = 6, colour = 'black') +
  geom_text(aes(x = 6, y = 20, label = 'd'),
            size = 6, colour = 'black') + # Here letters derived from statistical tests above
  annotation_logticks(sides = 'l') +
  scale_colour_manual(values = c('#480355','#9448BC','pink2','#38761d','black','black')) +
  theme_light(base_size = 17) +
  theme(axis.title.x = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(colour = 'black')) # Figure 1F

plot_grid(majorp, nrow = 1, labels = c('F'), label_size = 25) # Add the label
ggsave("figures/Figure3F-bottom.jpeg",
       width = 10, height = 5, dpi = 600, limitsize = F) # Save it

