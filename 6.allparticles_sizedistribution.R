#' 7. Distribution of size from faecal pellets' size from the 3 clusters, the 18 UVP5
#' images, sediment traps and the copepods incubations.
#' It produces Figure 4.

remove(list = ls())
set.seed(121999)

# Libraries
library("tidyverse")
library("cowplot")
library("patchwork")

# Datasets
part <- readr::read_tsv("data/export_detailed_20240404_15_40_PAR_Aggregated.tsv",
                        locale = locale(encoding = "latin1")) # All particles counted (EcoPart)
det <- read_csv("data/detritus_good-descriptors.csv.gz") # All detritus imaged (EcoTaxa)
clusters <- read_csv("kmeans2/detritus_37-clusters-kmeans.csv.gz") %>%
  filter(cluster %in% c(13,32,35)) %>%
  mutate(Type = ifelse(cluster == 13, 'Cluster 13',
                       ifelse(cluster == 32, 'Cluster 32', 'Cluster 35'))) %>%
  mutate(Type = factor(Type, levels = c('Cluster 13', 'Cluster 32', 'Cluster 35'))) # Faecal pellets-like particles from the 3 clusters

st_size <- read_csv("data/faecalpellets_inddata.csv") # Sediment traps
inc_size <- read_csv("data/faecalpellets_incubations_mm.csv") # Incubation
poop <- read_csv("data/poopingcopepods_mm.csv") # From the 18 images
calanus <- read_csv("data/faecalpellets_calanus-sizerewiev.csv",
                    col_types = cols(Value = col_number())) # Literature data for Arctic copepods
groszoo <- read_csv("data/faecalpellets_macrozoo-sizerewiev.csv") # Literature data for zoo

# 1. Prepare the dataset ----
## 1.a) Particles counted from EcoPart ----
colnames(part)[7:51] <- c(c(1, 1.26, 1.59, 2, 2.52, 3.17, 4, 5.04, 6.35, 8, 10.1,
                          12.7, 16, 20.2, 25.4, 32, 40.3, 50.8, 64, 80.6, 102,
                          128, 161, 203, 256, 323, 406, 512, 645, 813)*0.001,
                          c(1.02, 1.29, 1.63, 2.05, 2.58, 3.25, 4.1, 5.16,
                            6.5, 8.19, 10.3, 13, 16.4, 20.6, 26))

part <- part %>%
  dplyr::select(`0.001`:`26`) %>%
  pivot_longer(cols = `0.001`:`26`, names_to = 'major_mm', values_to = 'n') %>%
  group_by(major_mm) %>%
  dplyr::summarise(sum = sum(n)) %>%
  mutate(major_mm = as.numeric(major_mm))

major_mm <- NA
for (i in 1:length(part$major_mm)){
  t <- rep(part$major_mm[i], times = part$sum[i])
  major_mm <- c(major_mm, t)
}
major_mm <- major_mm[-1]

particle <- as.data.frame(major_mm)
particle$type <- 'UVP5 - Counted particles'

## 1.b) Zoo literature ----
groszoo_mm <- groszoo %>%
  filter(Variable == 'Length') %>%
  mutate(major_mm = ifelse(Unit == 'µm', Value * 0.001, Value),
         Unit = 'mm',
         type = 'Euphausids, mysids and amphipods')

## 1.c) Calanus literature ----
calanus_mm <- calanus %>%
  filter(Variable == 'Length') %>%
  mutate(major_mm = ifelse(Unit == 'µm', Value * 0.001, Value),
         Unit = 'mm',
         type = 'Calanus spp.')

# 2. Merge the datasets ----
det <- det %>%
  dplyr::select(major_mm) %>%
  mutate(type = "UVP5 - Imaged particles")

clusters <- clusters %>%
  dplyr::select(major_mm, Type) %>%
  mutate(type = paste("UVP5 - ", Type, sep = '')) %>%
  dplyr::select(-Type)

st_size <- st_size %>%
  dplyr::select(Length_mm) %>%
  mutate(major_mm = Length_mm,
         type = "Sediment traps") %>%
  dplyr::select(-Length_mm)

inc_size <- inc_size %>%
  dplyr::select(Length_mm) %>%
  mutate(major_mm = Length_mm,
         type = "Incubations") %>%
  dplyr::select(-Length_mm)

all <- rbind(particle, det, clusters, st_size, inc_size)
all$type <- factor(all$type, levels = c('UVP5 - Counted particles',
                                        'UVP5 - Imaged particles',
                                        'UVP5 - Cluster 13',
                                        'UVP5 - Cluster 32',
                                        'UVP5 - Cluster 35',
                                        'Incubations','Sediment traps')) # Order the variable

# 3. Plot ----
p1 <- ggplot() +
  geom_boxplot(aes(x = major_mm_ecotaxa, y = 2.5),
             colour = '#38761d', fill = '#38761d', alpha = 0.1, linewidth = 1,
             data = poop) +

  geom_density(aes(x = major_mm, colour = type, fill = type, lty = type),
               linewidth = 1, bw = 0.03, alpha = 0.05,
               data = all %>% filter(! type %in% c('Incubations','Sediment traps'))) +
  geom_density(aes(x = major_mm, colour = type, fill = type, lty = type),
               linewidth = 1, bw = 0.03, alpha = 0.05,
               data = all %>% filter(type %in% c('Incubations','Sediment traps'))) +

  geom_text(aes(x = mean(poop$major_mm_ecotaxa), y = 3.25,
                label = 'UVP5 - 18 images'),
            size = 5, colour = '#38761d') +

  scale_colour_manual(values = c('#8c510a','#dfc27d','#480355','#9448BC','pink2','black','black')) +
  scale_fill_manual(values = c('#8c510a','#dfc27d','#480355','#9448BC','pink2','black','black')) +
  scale_linetype_manual(values = c('solid','solid','solid','solid','solid','solid','dashed')) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,4.1), position = 'bottom',
                     breaks = seq(0,5,0.25)) +
  xlab('Particles length (mm)') +
  scale_y_continuous(expand = c(0, 0), limits = c(0,10.5)) +
  ylab('Density') +
  theme_light(base_size = 17) +
  theme(legend.title = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0.85,0.87)) +
  guides(fill = guide_legend(label.position = 'left'),
         colour = guide_legend(label.position = 'left'))


p2 <- ggplot() +
  annotate("segment", x = 0.85, xend = 4.1, y = 5, yend = 5,
               colour = 'black', linewidth = 1, lty = 'dotted') +
  annotate("segment", x = 0.85, xend = 2, y = 5, yend = 5,
           colour = 'black', linewidth = 1.5, lty = 'solid') +
  geom_point(aes(x = major_mm, y = 1),
               colour = 'grey40', alpha = 0.25, size = 3,
               data = calanus_mm) +
  geom_point(aes(x = major_mm, y = 3),
               colour = 'grey40', alpha = 0.25, size = 3,
               data = groszoo_mm) +
  geom_text(aes(x = 1, y = 5.75,
                label = 'Size range visible through UVP5/6 '),
            size = 5, colour = 'black', hjust = 0) +
  geom_text(aes(x = quantile(calanus_mm$major_mm,0.25), y = 1.75,
                label = 'Calanus spp.'),
            size = 5, colour = 'grey40', hjust = 0, fontface = "italic") +
  geom_text(aes(x = quantile(calanus_mm$major_mm,0.25) + 0.55, y = 1.75,
                label = 'faecal pellets'),
            size = 5, colour = 'grey40', hjust = 0) +
  geom_text(aes(x = quantile(groszoo_mm$major_mm,0.25), y = 3.75,
                label = 'Euphausids, mysids and amphipods faecal pellets'),
            size = 5, colour = 'grey40', hjust = 0) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,4.1)) +
  theme_classic() +
  theme(axis.text = element_text(color = 'white'),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())

p1/p2 + plot_layout(ncol = 1, nrow = 2, heights = c(4, 1)) # Plot both parts

ggsave("figures/Figure4-kmeans37.jpg",
       dpi = 1200, height = 10, width = 10) # Save it

# 4. Try to compute intersection areas ----
part1à4 <- as.data.frame(density(particle$major_mm, from = 0, to = 4)[c("x", "y")], adjust = 0.000000001) %>%
  mutate(Length = x,
         Density = y,
         Type = 'UVP5 - Counted particles') %>%
  dplyr::select(Length, Density, Type)

det1à4 <- as.data.frame(density(det$major_mm, from = 0, to = 4)[c("x", "y")]) %>%
  mutate(Length = x,
         Density = y,
         Type = 'UVP5 - Imaged particles') %>%
  dplyr::select(Length, Density, Type)

cluster131à4 <- as.data.frame(density(dplyr::filter(clusters, cluster == 13)$major_mm, from = 0, to = 4)[c("x", "y")]) %>%
  mutate(Length = x,
         Density = y,
         Type = 'UVP5 - Cluster 13') %>%
  dplyr::select(Length, Density, Type)

cluster321à4 <- as.data.frame(density(dplyr::filter(clusters, cluster == 32)$major_mm, from = 0, to = 4)[c("x", "y")]) %>%
  mutate(Length = x,
         Density = y,
         Type = 'UVP5 - Cluster 32') %>%
  dplyr::select(Length, Density, Type)

cluster351à4 <- as.data.frame(density(dplyr::filter(clusters, cluster == 35)$major_mm, from = 0, to = 4)[c("x", "y")]) %>%
  mutate(Length = x,
         Density = y,
         Type = 'UVP5 - Cluster 35') %>%
  dplyr::select(Length, Density, Type)

st1à4 <- as.data.frame(density(st_size$major_mm, from = 0, to = 4)[c("x", "y")]) %>%
  mutate(Length = x,
         Density = y,
         Type = 'Sediment traps') %>%
  dplyr::select(Length, Density, Type)

inc1à4 <- as.data.frame(density(inc_size$major_mm, from = 0, to = 4)[c("x", "y")]) %>%
  mutate(Length = x,
         Density = y,
         Type = 'Incubations') %>%
  dplyr::select(Length, Density, Type)

all <- rbind(part1à4, det1à4,
             cluster131à4, cluster321à4, cluster351à4,
             st1à4, inc1à4)

all$Type <- factor(all$Type, levels = c('UVP5 - Counted particles',
                                        'UVP5 - Imaged particles',
                                        'UVP5 - Cluster 13',
                                        'UVP5 - Cluster 32',
                                        'UVP5 - Cluster 35',
                                        'Incubations','Sediment traps'))


IntDensity <- function(dx, dy){

  Xdensity <- approxfun(density(dx$major_mm, from = 0, to = 4)[c("x", "y")])
  Ydensity <- approxfun(density(dy$major_mm, from = 0, to = 4)[c("x", "y")])

  # XminusY <- function(x) (Xdensity(x) - Ydensity(x))
  #
  # ex <- uniroot(XminusY, c(0,4))$root

  dX <- as.data.frame(density(dx$major_mm, from = 0, to = 4)[c("x", "y")])
  dY <- as.data.frame(density(dy$major_mm, from = 0, to = 4)[c("x", "y")])

  dX <- dX %>%
    mutate(delta = dX$y - dY$y)
  ex <- NULL  # Store the intersection points
  ey <- NULL
  k = 0

  for (i in 2:length(dX$x)) {
    # Look for a change in sign of the difference in densities
    if (sign(dX$delta[i-1]) != sign(dX$delta[i])) {
      k = k + 1
      # Linearly interpolate
      ex[k] <- dX$x[i-1] +
        (dX$x[i]-dX$x[i-1]) * (0-dX$delta[i-1]) / (dX$delta[i]-dX$delta[i-1])
    }
  }

  area <- integrate(Xdensity, 0, ex)$value + integrate(Ydensity, ex, 4)$value

  return(ex)
}


dx = det
dy = dplyr::filter(clusters, cluster == 32)
IntDensity(st_size,dplyr::filter(clusters, cluster == 32))


ggplot() +
  geom_boxplot(aes(x = major_mm_ecotaxa, y = 2.5),
               colour = '#38761d', fill = '#38761d', alpha = 0.1, linewidth = 1,
               data = poop) +

  geom_path(aes(x = Length, y = Density, colour = Type, lty = Type),
               linewidth = 1, alpha = 1,
               data = all) +

  geom_text(aes(x = mean(poop$major_mm_ecotaxa), y = 3.25,
                label = 'UVP5 - 18 images'),
            size = 5, colour = '#38761d') +

  scale_colour_manual(values = c('#8c510a','#dfc27d','#480355','#9448BC','pink2','black','black')) +
  scale_fill_manual(values = c('#8c510a','#dfc27d','#480355','#9448BC','pink2','black','black')) +
  scale_linetype_manual(values = c('solid','solid','solid','solid','solid','solid','dashed')) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,4.1), position = 'bottom',
                     breaks = seq(0,5,0.25)) +
  xlab('Particles length (mm)') +
  scale_y_continuous(expand = c(0, 0), limits = c(0,10.5)) +
  ylab('Density') +
  theme_light(base_size = 17) +
  theme(legend.title = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0.85,0.87)) +
  guides(fill = guide_legend(label.position = 'left'),
         colour = guide_legend(label.position = 'left'))

