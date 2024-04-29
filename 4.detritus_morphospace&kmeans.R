#' This script computes detritus morphospace from 18 morphological descriptors,
#' hence identifies the main traits of variations. It also tries all the clustering
#' possibilities between 1 and 50.
#' The morphospace plot was inspired by JOI code on https://github.com/jiho/ptb_morphodiv

remove(list = ls())
set.seed(121999)

# Libraries
library("tidyverse") # Data manipulation
library("morphr")
library('plyr')
library('car')
library('corrplot')
library('ggrepel')
library("RColorBrewer") # Colours in plot
library('cowplot') # For multiple plotting
library('vegan')
library('FactoMineR')
library('factoextra')
library('ggimage') # For multiple images plot

# Data frames
det <- read_csv("data/detritus_good-descriptors.csv.gz") # Detritus with cleaned morphological descriptors
st <- read_delim("data/Arctic_Driftingtraps25m.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE) # Sediment traps metadata_mp
st_size <- read_csv("data/faecalpellets_inddata.csv") # Sediment traps
inc_size <- read_csv("data/faecalpellets_incubations_mm.csv") # Incubation
poop <- read_csv("data/poopingcopepods_mm.csv") # From the 18 images

# 1. Prepare data ----
det <- det %>%
  dplyr::select(object_id : depth, MajorAxis = major, MinorAxis = minor, Perimeter = perim,
                FeretDiameter = feret, Area = area, Elongation = elongation,
                Symmetry = symv, "Symmetry75%" = symvc, Circularity = circ,
                ThicknessRatio = thickr, MeanGrey = mean, MedianGrey = median,
                "Grey75%" = hist75, StdDevGrey = stddev, SkewnessGrey = skew, Fractal = fractal,
                "Perim/Feret" = perimferet, "Perim/Major" = perimmajor, sample_id : conc, station)

# 1. Morphospace through PCA ----
## 1.a) Perform the PCA ----
pca_var <- dplyr::select(det, MajorAxis:'Perim/Major') # Select relevant morphological features

pca <- FactoMineR::PCA(pca_var,
                       row.w = c(det$conc/sum(det$conc, na.rm = T)), graph = FALSE) # Each detritus is weighted by its concentration

## 1.b) Analyses the axes ----
# Eigenvalues
eig <- as.data.frame(pca$eig)
colnames(eig) <- c('eigenvalues', '%variance', 'cumulative%')
eig <- rownames_to_column(eig, 'Components')
eig %>% filter(eigenvalues > mean(eigenvalues))
remove(eig)

# -> 4 axis significant, 90.54%
# PC1 = 46.13%
# PC2 = 20.98%
# PC3 = 17.83%
# PC4 = 5.60%

## Check the contribution of each variables
corrplot(pca$var$cor[,1:4], tl.col = "black")

## 1.c) Save the coordinates ----
det_c <- bind_cols(det, as_tibble(pca$ind$coord[,1:4]) %>%
                     set_names(str_c("PC", 1:4))) %>%
  distinct()

det_c %>%
  write_csv("data/detritus_coordinates.csv.gz")

# 2. Prepare datasets for clustering ----
d_coord <- det_c %>%
  dplyr::select(PC1:PC2) # Non-living partivles coordinates

st$Trap <- seq(1,12,1)
st_size <- left_join(st_size, dplyr::select(st, -mission, -date_time_close, -PI))
rm(st)
st_size <- st_size[,-12]
st_size$Type <- 'Sediment traps'
inc_size <- inc_size %>%
  mutate(Type = 'Incubations') # Incubations size

poop <- poop %>%
  mutate(Length_mm = major_mm_ecotaxa,
         Radius_mm = minor_mm_ecotaxa,
         Type = 'UVP5 - 18 images') # 18 UVP5 images size

img_trimmed <- list.files("images/det_cut/", pattern = '.jpg') # Get the names of all the images files

wss <- 0 # To keep within sum of squares at each iteration
within <- as.data.frame(matrix(ncol = 2, nrow = 50)) # To save all of them
colnames(within) <- c('k','wss')
within$k <- seq(1,50,1)


# 3. Try k-means for k = 1 to k = 50 ----
for (i in 1:50) {
  km <- kmeans(d_coord, centers = i, iter.max = 1000, algorithm = "MacQueen") # K-means

  det_c$cluster <- km$cluster # Add cluster as a variable in the dataset

  for (j in 1:i){
    d <- det_c %>%
      dplyr::select('PC1','PC2','cluster') %>%
      filter(cluster == j)

    y_mean <- mean(d$PC1)
    tss_1 <- sum((d$PC1-y_mean)^2)
    y_mean <- mean(d$PC2)
    tss_2 <- sum((d$PC2-y_mean)^2)
    wss <- tss_1 + tss_2 + wss
  } # Compute the within sum of squares

  within$wss[i] <- wss # Save it in the dataset
  wss <- 0 # Back to zero for next iteration

  print(paste('Kmeans - ',i,'clusters', sep = ' '))

  det_c %>%
    write_csv(paste("kmeans4/detritus_",i,"-clusters-kmeans.csv.gz", sep = '')) # Save the dataset for later

  centroids <- as.data.frame(km$centers)
  centroids$cluster <- seq(1,i,1) # Save the centroids positions

  p1 <- imagesa +
    geom_point(aes(x = PC1, y = PC2, colour = as.factor(cluster)),
               alpha = 0.01, size = 1,
               data = det_c) +
    geom_point(aes(x = PC1, y = PC2, colour = as.factor(cluster)),
               alpha = 1, size = 3,
               data = centroids) +
    geom_text_repel(aes(x = PC1, y = PC2, colour = as.factor(cluster), label = cluster),
                    size = 10,
                    data = centroids) +
    theme_classic(base_size = 25) +
    theme(legend.position = 'none',
          legend.direction = 'horizontal',
          legend.title = element_blank()) # Plot clusters position in PC1PC2

  ggsave(paste("kmeans4/detritus_",i,"-clusters-kmeans_12.png", sep = ''),
         p1, width = 16, height = 16, limitsize = FALSE) # Save it

  filterPC <- det_c %>%
    dplyr::select('PC1','PC2','cluster') %>%
    distinct() %>%
    group_by(cluster) %>%
    dplyr::summarise(PC1m = mean(PC1),
                     PC2m = mean(PC2)) %>%
    filter(PC1m < 0 & PC2m < 0) # Filter clusters in the bottom left panel (where we look for potential faecal pellets)

  kcofinterest <- det_c %>%
    filter(cluster %in% filterPC$cluster) %>%
    dplyr::select('cluster','major_mm','minor_mm','object_id') %>%
    mutate(Type = cluster,
           Length_mm = major_mm,
           Radius_mm = minor_mm) %>%
    arrange(Type) # Get the sizes of particles in such clusters

  if (dim(kcofinterest)[1] == 0)
  {
    print(paste('No cluster with PC1 < 0 and PC2 < 0 for kmeans with',i,"clusters.", sep = " "))
  }
  else
  {
    print(paste("Among", length(unique(det_c$cluster)), "clusters,", length(unique(kcofinterest$cluster)), "are in PC1 and PC2 < 0."))

    allFP <- rbind(dplyr::select(kcofinterest, Length_mm, Radius_mm, Type), # Clusters, note that Type here = cluster
                   dplyr::select(poop, Length_mm, Radius_mm, Type), # 18 images
                   dplyr::select(st_size, Length_mm, Radius_mm, Type), # Sediment traps
                   dplyr::select(inc_size, Length_mm, Radius_mm, Type)) # Incubations

    allFP$Type <- factor(allFP$Type,
                         levels = c(unique(kcofinterest$Type),
                                    "UVP5 - 18 images","Incubations","Sediment traps")) # Order the "Type" variable

    boxplot <- ggplot(allFP) +
      geom_boxplot(aes(x = Type, y = Length_mm, group = Type, colour = Type),
                   outlier.size = 0.7) +
      scale_y_log10('Length (mm)') +
      annotation_logticks(sides = 'l') +
      ggtitle(paste(i, 'clusters, only PC1 < 0 & PC2 < 0', sep = " ")) +
      theme_light(base_size = 17) +
      theme(axis.title.x = element_blank(),
            legend.position = 'none',
            axis.text.x = element_text(colour = 'black')) # Boxplots of the lengths of particles from different origins

    ggsave(paste("kmeans4/detritus_",i,"-clusters-kmeans_comparison.png", sep = ''),
           boxplot, width = 12, height = 6, limitsize = FALSE) # Save it

    kcofinterest$img_path <- paste("images/det_cut/", kcofinterest$object_id, ".jpg", sep = '')
    kcofinterest$img_path <- ifelse(paste(kcofinterest$object_id, ".jpg", sep = '') %in% img_trimmed,
                                    kcofinterest$img_path, "images/det_cut/blank.png") # Get the file path for images in clusters

    for (j in 1:length(filterPC$cluster)) { # 1 cluster = 1 iteration
      truc <- kcofinterest %>%
        filter(Type == filterPC$cluster[j]) %>%
        sample_n(100, replace = TRUE) # Filter 100 random images in each cluster

      plot <- ggimg_grid(truc$img_path, scale = 0.003) +
        geom_rect(aes(xmin = -0.05, xmax = 1.05, ymin = -0.05, ymax = 1.05),
                  colour = 'black', fill = "white", alpha = 0, linewidth = 3) +
        ggtitle(paste('Cluster', filterPC$cluster[j], sep = ' ')) +
        theme_classic(base_size = 17) +
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              plot.title = element_text(hjust = 0.1, vjust = -6, colour = 'black', size = 30)) # Plot 100 random images from each cluster

      ggsave(paste("kmeans4/detritus_",i,"-clusters-kmeans_images-cluster",filterPC$cluster[j],".png", sep = ''),
             plot, width = 8, height = 8, limitsize = FALSE) # Save it
    }
  }

}

ggplot(within) +
  geom_point(aes(x = k, y = wss),
             shape = 16, size = 3) +
  geom_path(aes(x = k, y = wss),
            linewidth = 1) +
  scale_x_continuous(breaks = seq(0,51,1)) +
  theme_light(base_size = 17) # Plot the within sum of squares for each k

ggsave("kmeans3/wss.jpeg",
       width = 12, height = 6) # Save it

#' After looking at all the images from each cluster in the bottom left panel
#' and their size distribution, the best choice was for k = 37

# 4. Plot of the final clusters ----
## 4.a) 100 random images ----
clusters <- read_csv("kmeans2/detritus_37-clusters-kmeans.csv.gz") %>%
  filter(cluster %in% c(13,32,35)) %>%
  mutate(Type = ifelse(cluster == 13, 'Cluster 13',
                       ifelse(cluster == 32, 'Cluster 32', 'Cluster 35'))) # Final choice

p13 <- ggimg_grid(sample_n(filter(clusters, Type == 'Cluster 13'), 100, replace = T)$img_path2,
                  scale = 0.002) +
  geom_rect(aes(xmin = -0.05, xmax = 1.05, ymin = -0.05, ymax = 1.05),
            colour = 'pink2', fill = "white", alpha = 0, linewidth = 3) +
  theme_classic(base_size = 17) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())

p32 <- ggimg_grid(sample_n(filter(clusters, Type == 'Cluster 32'), 100, replace = T)$img_path2,
                  scale = 0.002) +
  geom_rect(aes(xmin = -0.05, xmax = 1.05, ymin = -0.05, ymax = 1.05),
            colour = '#9448BC', fill = "white", alpha = 0, linewidth = 3) +
  theme_classic(base_size = 17) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())

p35 <- ggimg_grid(sample_n(filter(clusters, Type == 'Cluster 35'), 100, replace = T)$img_path2,
                  scale = 0.002) +
  geom_rect(aes(xmin = -0.05, xmax = 1.05, ymin = -0.05, ymax = 1.05),
            colour = '#480355', fill = "white", alpha = 0, linewidth = 3) +
  theme_classic(base_size = 17) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())

plot_grid(p13, p32, p35, nrow = 1, align = 'hv',
          labels = c('A','B','C'), label_size = 25) # Get 100 images from each cluster side by side
ggsave("figures/Figure3-test11-top.jpeg",
       width = 10, height = 3, dpi = 600, limitsize = F) # Save it

## 4.b) Save 1000 images per cluster in a folder ----
nbcluster <- c(13,32,35) # Clusters number
for (i in nbcluster){
  subsample <- clusters %>%
    filter(Type == paste('Cluster ', i, sep = '')) %>%
    sample_n(1000, replace = T)

  FP <- paste("images/all_images/", subsample$sample_id, "/", subsample$object_id, ".jpg", sep = '')
  nFP <- paste("images/faecalparticles_",i,"/", sep = '')

  ok <- sapply((1:dim(subsample)[1]), function(j) {
    # make a small progress bar
    if (j %% 10 == 0) {cat(".")}
    picture <- str_split_fixed(FP[j], pattern = '/', n = 4)[4]
    img <- img_read(FP[j])
    path <- paste(nFP, picture, sep = '')
    img_write(img, path)
  })


}

## 4.c) PCA plot ----
det_c$img_path <- paste("images/det_cut/", det_c$object_id, ".jpg", sep = '')
img_trimmed <- list.files("images/det_cut/", pattern = '.jpg')
det_c$img_path2 <- ifelse(paste(det_c$object_id, ".jpg", sep = '') %in% img_trimmed,
                          det_c$img_path, "images/det_cut/blank.png")
which(det_c$img_path2 == "images/det_cut/blank.png") # Get images path, if images absent then white image

imagesa <- ggmorph_tile(pca, det_c$img_path2, dimensions = c(1,2),
                        steps = 40, n_imgs = 3, scale = 0.01, adjust_grey = F) # For PC1PC2
imagesb <- ggmorph_tile(pca, det_c$img_path2, dimensions = c(3,4),
                        steps = 40, n_imgs = 3, scale = 0.01, adjust_grey = F) # For PC3PC4

# Additional needs
circleFun <- function(center = c(0,0), diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0, 2*pi, length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
} # To draw a circle
circ <- circleFun(c(0,0), 2, npoints = 500)

# Change scaling of variables/columns from scaling 1 to 2
eig <- pca$eig[,1]
varscores <- as.data.frame(t(t(pca$var$coord) / sqrt(eig))) #de-scale
varscores2 <- as.data.frame(t(t(varscores) * sqrt(nrow(varscores) * eig))) #re-scale
pca_coord <- varscores2
pca_coord$Descriptors <- rownames(pca_coord)
pca_coord$type <- factor(c('Size','Size','Size','Size','Size',
                           'Shape','Shape','Shape','Shape','Shape',
                           'Transparency','Transparency','Transparency','Transparency','Transparency',
                           'Complexity','Complexity','Complexity'),
                         levels = c('Size','Transparency','Shape','Complexity')) # Type of morphological variables
pca_coord$x <- c(4,-0.2,4.8,4.8,2.45, # Size
                 1.75,4,1.95,-4,4.2, # Shape
                 3,2.5,4.1,-4,-4.6, # Transparency
                 4.3,1.7,3.1) # Complexity
pca_coord$y <- c(-1.75,3.35,0.1,-1.15,1.5,
                 -3.85,-2.55,-2.65,1.2,1.25,
                 0.85,0.2,0.65,-0.75,-0.45,
                 0.4,3.75,3.35)

## 4.d) Plots ----
pc12 <- imagesa +
  geom_path(aes(x = x*4, y = y*4),
            lty = 2, color = "grey50", alpha = 0.7,
            data = circ) + # Adapt scaling of circle to fit the arrows (here: 4 for PC12; 3 for PC34)
  geom_hline(yintercept = 0,
             lty = 2, color = "grey50", alpha = 0.9) +
  geom_vline(xintercept = 0,
             lty = 2, color = "grey50", alpha = 0.9) +
  geom_point(aes(x = PC1, y = PC2, colour = Type),
             data = clusters,
             shape = 16, alpha = 0.05, size = 1) + # Faecal particles points
  geom_segment(aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2, lty = type),
               colour = "black", linewidth = 1, alpha = 1,
               arrow = arrow(length = unit(0.65,"cm")),
               data = pca_coord) + # Descriptors arrows
  geom_text(aes(x = x, y = y, label = Descriptors),
            colour = "black", size = 6,
            fontface = 1, alpha = 1, show.legend = FALSE,
            data = pca_coord) + # Descriptors names
  scale_linetype_manual(values = c('dashed','solid','dotted','dotdash')) +
  scale_colour_manual(values = c('#480355','#9448BC','pink2')) +
  theme_classic(base_size = 17) +
  theme(legend.position = c(0.08,0.18),
        legend.spacing = unit(0, 'cm'),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.key.height = unit(1, "cm"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  xlab(paste("PC1 (", round(pca$eig[1,2], 1),"%)", sep = "")) +
  ylab(paste("PC2 (", round(pca$eig[2,2], 1),"%)", sep = "")) +
  guides(colour = guide_legend(override.aes = list(size = 6, alpha = 1))) # Figure 2A

ggsave("figures/Figure2-clusters-b.jpg",
       pc12, width = 12, height = 10, dpi = 600) # Save it

pc34 <- imagesb +
  geom_path(aes(x = x*3, y = y*3),
            lty = 2, color = "grey50", alpha = 0.7,
            data = circ) + # Adapt scaling of circle to fit the arrows (here: 4 for PC12; 3 for PC34)
  geom_hline(yintercept = 0,
             lty = 2, color = "grey50", alpha = 0.9) +
  geom_vline(xintercept = 0,
             lty = 2, color = "grey50", alpha = 0.9) +
  geom_point(aes(x = PC3, y = PC4, colour = Type),
             data = clusters,
             shape = 16, alpha = 0.05, size = 1) + # Faecal particles points
  geom_segment(aes(x = 0, y = 0, xend = Dim.3, yend = Dim.4, lty = type),
               colour = "black", linewidth = 1, alpha = 1,
               arrow = arrow(length = unit(0.65,"cm")),
               data = pca_coord) + # Descriptors arrows
  geom_text_repel(aes(x = Dim.3, y = Dim.4, label = Descriptors),
                  colour = "black", size = 6,
                  fontface = 1, alpha = 1, show.legend = FALSE,
                  data = pca_coord) + # Descriptors names
  scale_linetype_manual(values = c('dashed','solid','dotted','dotdash')) +
  scale_colour_manual(values = c('#480355','#9448BC','pink2')) +
  theme_classic(base_size = 17) +
  theme(legend.position = c(0.08,0.18),
        legend.spacing = unit(0, 'cm'),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.key.height = unit(1, "cm"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  xlab(paste("PC3 (", round(pca$eig[3,2], 1),"%)", sep ="")) +
  ylab(paste("PC4 (", round(pca$eig[4,2], 1),"%)", sep ="")) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) # Figure S2

ggsave("figures/FigureS2-clusters.jpg",
       pc34, width = 12, height = 10, dpi = 600) # Save it

remove(pc12,pc34)
