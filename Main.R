# Script for processing AVIRIS/NEON AOP flightline trait data from 
# Phil Townsend Spectral Lab at UWisc
# Trait data based on Wang et al. (2020, New Phytologist)
# Code by Natasha Stavros, CU Boulder

library(rgdal)
library(raster)
library(gdalUtils)
library(sf)
library(ggplot2)
library(caTools)
library(mgcv)
#library(rdist)
library(RiemBase)
library(factoextra)
library(NbClust)
library(stringr)
source("/Users/natasha/Coding/RFunctions/Airborne_FunctionalTraits/Functions/Preprocess.R")
source("/Users/natasha/Coding/RFunctions/Airborne_FunctionalTraits/Initial_Parameters.R")
source("/Users/natasha/Coding/RFunctions/Airborne_FunctionalTraits/Functions/Diversity.R")


####### Prep Trait data into mosaics #######
AVIRISmosaics <- 
  list.files(AVIRISpath,pattern= paste0(outputfile_descriptor,"_mosaic_"))

if (length(AVIRISmosaics) == 0){
  # Mosaic Flightlines assuming all flightlines are saved to one path
  mosaic_flightlines(traits = traits, path = AVIRISpath,
                     bounding_vector = boundingKML, missingData = missingData,
                     datatype=datatype,
                     dataThresholds_upper = dataThresholds_upper,
                     dataThresholds_lower = dataThresholds_lower,
                     save_intermediate_output = FALSE, 
                     outputfile_descriptor = outputfile_descriptor) 
}

####### Do an Analysis #######
# Read in one trait mosaic to determine size
temp <- raster(paste0(AVIRISpath,AVIRISmosaics[1]),band=1)
mosaic_size <- dim(temp)
mosaic_size_1D <- mosaic_size[1]*mosaic_size[2]
mosaic_extent <- temp@extent
mosaic_crs <- temp@crs
rm(temp)

## Run K-means
# Prep your matrix
trait_datastack <- trait_matrix_prep(mosaic_size_1D, traits, traits_4_diversity, 
                                     path, AVIRISmosaics)
# Select number of clusters 
# Good Tutorial: https://statsandr.com/blog/clustering-analysis-k-means-and-hierarchical-clustering-by-hand-and-in-r/#kmeans-with-2-groups
# randomly select subset to test # of clusters
random_sample_indices <- sample(1:dim(na.omit(trait_datastack))[1],1000)

# Elbow method
fviz_nbclust(na.omit(trait_datastack)[random_sample_indices,], kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2) + # add line for better visualisation
  labs(subtitle = "Elbow method") # add subtitle

## Silhouette method
fviz_nbclust(na.omit(trait_datastack)[random_sample_indices,], kmeans, method = "silhouette") +
  labs(subtitle = "Silhouette method")

cluster_raster <- kmeans_trait_data(trait_datastack,k,AVIRISpath,outputfile_descriptor,
                              mosaic_extent,mosaic_crs)
# # Before removing NAs to run kmeans, store their indices
# not_nan_indices <- match(na.omit(trait_datastack)[,1],trait_datastack[,1])
# 
# # run your cluster
# functional_clusters <- kmeans(na.omit(trait_datastack),centers = 4)
# 
# # how good are our clusters?
# BSS <- functional_clusters$betweenss
# TSS <- functional_clusters$totss
# (percent_explained <- BSS / TSS * 100)

# # map clusters back to 2d
# cluster_matrix <- matrix(data=NA,nrow = dim(trait_datastack)[1], ncol = 1)
# cluster_matrix[not_nan_indices,1] <- functional_clusters$cluster
# cluster_2d <- matrix(cluster_matrix,nrow=mosaic_size[1], byrow = TRUE)
# 
# # export to geotiff
# cluster_raster <- raster(cluster_2d)
# extent(cluster_raster) <- mosaic_extent
# crs(cluster_raster) <- mosaic_crs
# writeRaster(cluster_raster,paste0(AVIRISpath,outputfile_descriptor,"_kmeans.tif"),overwrite=TRUE)








####### Calculate Functional Diversity #######
# Find new grid spacings at coarser res for diversity calculation
# Normalize trait rasters
# Loop through all coarse grid cells to calculate trait diversity matrix and then diversity metrics by grid cell
# put grid cells back together into raster


Trait_Diversity_Matrix <- trait_matrix_prep(mosaic_size_1D, traits, traits_4_diversity, 
                                                        AVIRISpath, AVIRISmosaics)

# Create Diversity Matrices
physio_divergence <- functional_divergence(Trait_Diversity_Matrix)

# Write Raster to Geotiff