library(rgdal)
library(raster)
library(gdalUtils)
library(sf)
library(ggplot2)
library(caTools)

raster2matrix <- function(RasterLayer){
  # Precondition: RasterLayer is a raster object
  
  # Postcondition: a matrix object of the RasterLayer data values
  
  if(dim(RasterLayer)[3] <= 1){
    A <- raster::as.matrix(RasterLayer)
    m <- t(A[nrow(A):1,])
  }else{
    m <- aperm(raster::as.array(flip(RasterLayer,direction=2)),c(2,1,3))
  }
  return(m)
}

normalize_traits <- function(RasterLayer, max_value = -9999){
  # Precondition: R is a raster object of trait data to be normalized; max_value
  #           is the upper bound of trait value where -9999 sets the max value to the 
  #           value within a single raster scene, otherwise it sets the max to whatever
  #           value you select across sites (e.g., max value across sites).
  
  # Postconditions: normalized trait raster to provide values between 0 and 1. 
  #           See Fabian D Schneider et al. (2017) Nature Communications
  #           DOI: 10.1038/s41467-017-01530-3.
  
  raster_matrix <- raster2matrix(RasterLayer)
  if (max_value == -9999){raster_max <- max(raster_matrix[!is.na(raster_matrix)])}
  else {raster_max <- max_value}
  raster_matrix_normalized <- 
    apply(raster_matrix, MARGIN = c(1,2), FUN = function(x) x/(raster_max))
  
  return (raster_matrix_normalized)
}

trait_matrix_prep <- function(mosaic_size_1D, traits, traits_4_diversity, 
                              path, AVIRISmosaics){
  # Preconditions: mosaic_size_1D is the total number of pixels in the mosaiced
  #               trait scene. traits is a list of strings state the traits available.
  #               traits_4_diversity are a list of strings with only the traits
  #               from the list available to consider for calculating diversity.
  #               Recommended 3 traits for diversity metric calculation. 
  #               path is the path to the geotiff raster mosaics. AVIRIS mosaics
  #               is a list of strings containing the names of each raster mosaic
  
  # Postcondition: a matrix set up with all pixels as rows and each trait needed
  #               to calculate the diversity metrics as a column
  
  trait_indices <- match(traits_4_diversity,traits)
  physiological_traits <- traits[trait_indices[!is.na(trait_indices)]]
  
  # Read in Files for trait mosaics for physiological diversity metrics and prepare
  # a trait matrix as input to diversity functions
  print("Initializing empty Trait Matrix...")
  Trait_Matrix <- matrix(,nrow= mosaic_size_1D, ncol=length(trait_indices))
  for (i in 1:length(trait_indices)){
    # Get Mosaiced Trait file names
    file_list <- Filter(function(x) grepl(physiological_traits[i],x), AVIRISmosaics)
    temp_file <- paste0(AVIRISpath,file_list[1]) # get ".tif" not .tif.aux.xml
    
    # Read in mosaic trait from Geotiff as Raster
    print(paste0("Reading in raster from Geotiff of: ",physiological_traits[i],"..."))
    trait_mosaic <- raster(temp_file,band=1)
    
    # Flatten trait matrix to 1D vector
    print("Flattening raster to 1-D vector...")
    trait_flat <- as.vector(trait_mosaic)
    
    # Prepare input matrix for Diversity metrics
    print("Adding to Trait Matrix...")
    Trait_Matrix[,i] <- trait_flat
    
    # clean up variables
    rm(file_list, temp_file,trait_mosaic,trait_flat)
  }
  return(Trait_Matrix)
}
