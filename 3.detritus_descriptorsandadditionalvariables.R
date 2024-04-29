#' Script to normalise and transform morphological descriptors using BoxCox
#' and compute additional variables such as the concentration for detritus,
#' and save it

remove(list = ls())
set.seed(121999)

# Libraries
library("tidyverse") # Data manipulation
library("morphr") # For morphological descriptors and images manipulation
library('plyr')
library('car')

# Data frames
detritus <- read_csv("data/detritus.csv.gz")
vol <- read_csv("data/volper10mbin.csv.gz", col_types = cols())

# 1. Re-format the file ----
## Add human readable names for features and remove unuseful variables
keep <- c('object_id', 'object_date', 'object_time', 'station', 'new_taxon', 'depth', 'object_area',
          'object_mean', 'object_stddev', 'object_perim.', 'object_major', 'object_minor',
          'object_circ.', 'object_feret', 'object_median', 'object_skew', 'object_fractal',
          'object_histcum1', 'object_symetriev', 'object_symetrievc', 'object_thickr',
          'object_elongation', 'object_perimferet', 'object_perimmajor', 'sample_id', # Yes until here
          'sample_bottomdepth')

pixel_mm <- unique(detritus$acq_pixel) # pixel size in mm 0.088

det_small <- detritus[keep]
remove(keep)

## Rename the descriptors into simpler names
names(det_small)[1:26] <- c("object_id", "date","time",'station',"category","depth",
                            "area","mean","stddev","perim","major","minor","circ",
                            "feret","median","skew","fractal","hist75","symv",
                            "symvc","thickr","elongation","perimferet","perimmajor",
                            'sample_id','sample_bottomdepth')

# 2. Transform and save images ----
# ## Get information relevant to images
# d_imgs <- list.files("images/images_total/", pattern = '.jpg', recursive = T) # Which images we have in the folders
# d_imgs_path <- paste("images/images_total/", d_imgs, sep = '')
# det <- read_csv("data/detritus_bad-descriptors.csv.gz")
# det_img <- det %>%
#   mutate(path = paste("images/images_total/", det$sample_id, "/", det$object_id, ".jpg", sep = '')) # Which images we want to trim
# d_to_trim <- intersect(d_imgs_path, det_img$path)
# # det_img$TF <- det_img$path %in% d_to_trim
#
# ## Trim images in a new file
# trim_imgs <- "images/detritus_trimmed2/" # New path for trimmed images
# i <- 0
#
# ok <- sapply(1:length(d_to_trim), function(i) {
#   i <- i + 1
#   # make a small progress bar
#   if (i %% 10000 == 0) {cat("=")}
#   picture <- strsplit(d_to_trim[[i]], split = '/', fixed = T)[[1]][4]
#   img <- img_read(d_to_trim[i]) %>%
#     img_extract_largest(quiet = TRUE)
#   path <- paste(trim_imgs, picture, sep = '')
#   img_write(img, path)
# })

# -> Done once, so good with that now
# Might need to do it several times (a lot of images)

# 3. Compute additional variables ----
### Major, minor in mm and then biovolume (in mm3)
det_small$major_mm <- det_small$major*pixel_mm
det_small$minor_mm <- det_small$minor*pixel_mm
det_small$biovolume <- (4/3)*pi*((det_small$minor_mm/2)^2)*(det_small$major_mm/2) # Ellipsoidal shape

### Depth bins and individual concentrations
det_small <- det_small %>%
  mutate(depth_bin10 = round_any(depth, 10, floor)) %>% # create bins
  left_join(distinct(vol), by = c("station","depth_bin10")) %>% # left join to have watervolume of the bin for each objects
  mutate(conc = 1/(1.15*0.000001)) # calculate individual concentration (1/watervolume)

# 4. Transform distribution of images features (from 1.prepare_zoo.R, https://github.com/jiho/ptb_morphodiv) ----
## Plot all histograms => Extreme values in hist75, median especially
gather(sample_frac(det_small, 0.1), key = "var", val = "val", area:perimmajor) %>%
  ggplot() +
  geom_histogram(aes(x = val), bins = 50) +
  facet_wrap(~var, scales = "free", ncol = 4)

# Trim extreme values from the descriptors
#' Trim extreme values
#'
#' Replace extreme values by NA
#'
#' @param x a vector
#' @param p proprotion of values to remove
#' @param side on which extreme to remove values (left=low, right=high, both=...both). Can be abbreviated
#'

trim <- function(x, p = 0.001, side = "right") {
  # check argument (and allow to abbreviate it)
  side <- match.arg(side, choices = c("both", "left", "right"))
  # compute quantiles
  q <- quantile(x, probs = c(p, 1-p), na.rm = TRUE)
  # mask extreme values
  if (side == "left" | side == "both") {
    x[x < q[1]] <- NA
  }
  if (side == "right" | side == "both") {
    x[x > q[2]] <- NA
  }
  return(x)
}

# Apply the function trim
det_t <- det_small %>% mutate(
  area = trim(area),
  mean = trim(mean, side = "left"),
  stddev = trim(stddev),
  perim = trim(perim),
  major = trim(major),
  minor = trim(minor),
  circ = trim(circ),
  feret = trim(feret),
  median = trim(median, side = "left"),
  skew = trim(skew, side = "both"),
  fractal = trim(fractal, side = "both"),
  hist75 = trim(hist75, side = "left"),
  symv = trim(symv),
  symvc = trim(symvc),
  thickr = trim(thickr, side = "both"),
  elongation = trim(elongation),
  perimferet = trim(perimferet, 0.005),
  perimmajor = trim(perimmajor)
)

# Transform with boxcox
# Switch the values below 1 to above 1
min.cols <- sapply(det_t[,7:24], function(x) min(x, na.rm = T))
idx <- which(min.cols < 1)

for (ii in idx) {
  det_t[,7:24][[ii]] <- det_t[,7:24][[ii]] - min.cols[ii] + 1
}
summary(det_t)
#-> no more values under 1 (problematic for powerTransform function of car package)

# Find the coefficients for the boxcox transformation
lambda <- c() # Calcultate lambda
for (ii in 7:24){
  x <- car::powerTransform(det_t[[ii]])
  lambda[ii] <- x$lambda
  cat("Columns",ii,"done, lambda is:",x$lambda,". ")
}
lambda <- lambda[7:24]
names(lambda) <- colnames(det_t)[7:24]

# Transform with boxcox
for (ii in colnames(det_t)[7:24]){
  det_t[[ii]] <- bcPower(det_t[[ii]], lambda[[ii]])
}

# Plot to check if data seem normalised
gather(sample_frac(det_t, 0.1), key = "var", val = "val", area:perimmajor) %>%
  ggplot() +
  geom_histogram(aes(x = val), bins = 50) +
  facet_wrap(~var, scales = "free", ncol = 4)

# Bind to have objects information
## Eliminate extreme inviduals = more than 5 features are NA
n_na <- select(det_t, area:perimmajor) %>%
  apply(1, function(x) {sum(is.na(x))})
det_na <- det_t[n_na <= 5,]
summary(det_na)
#-> 992,407 objects now


## Replace NAs by the mean of the column
for (col in names(select(det_p, area:perimmajor))) {
  det_p[[col]][is.na(det_p[[col]])] <- mean(det_p[[col]], na.rm = TRUE)
}

summary(det_p)
# -> no more NA

## Save dataframe with cleaned fetures
write_csv(det_p, "data/detritus_good-descriptors.csv.gz")

