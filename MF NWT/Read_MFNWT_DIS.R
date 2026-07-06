read_mfnwt_dis <- function(data_path,
                           header = FALSE,
                           delimeter = ' ',
                           newline = '\t',
                           layer_text = c('top','botm'))
{
  #-------------------------------------------------------------------------------
  # read in the raw text values
  raw_txt <- read.delim(data_path,
                        header = header,
                        sep = newline)
  raw_txt <- as.vector(unlist(raw_txt))
  layer_indices <- c(grep(layer_text[1],raw_txt),
                     grep(layer_text[2],raw_txt))
  dims <- strsplit(raw_txt[2], delimeter)[[1]]
  dims <- dims[dims != '']
  nlay <- as.numeric(dims[1])
  nrow <- as.numeric(dims[2])
  ncol <- as.numeric(dims[3])
  end_dims <- grep('TR',raw_txt)[1]
  ncell <- nrow*ncol
  #-------------------------------------------------------------------------------

  
  #-------------------------------------------------------------------------------
  layer_defs <- list()
  for(i in 1:(length(layer_indices))){
    #-------------------------------------------------------------------------------
    if(i != length(layer_indices)){
      #-------------------------------------------------------------------------------
      elevations <- raw_txt[(layer_indices[i]+1):(layer_indices[i+1]-1)]
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # if its auto converted to a vector by R
      if(length(elevations) > nrow){
        elevations <- elevations[elevations != '']
        m <- matrix(data = as.numeric(elevations),
                    nrow = nrow,
                    ncol = ncol,
                    byrow = T)
      } else {
        elevations <- split(elevations, 1:length(elevations))
        elevations <- lapply(elevations, function(x){
          step1 <- strsplit(x, split = delimeter)[[1]]
          step2 <- step1[step1 != '']
          step2 <- as.numeric(step2)
        })
        names(elevations) <- NULL
        elevations <- do.call(rbind, elevations)
        m <- elevations
      }
      #-------------------------------------------------------------------------------
    }
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    if(i == length(layer_indices)){
      #-------------------------------------------------------------------------------
      elevations <- raw_txt[(layer_indices[i]+1):(end_dims-1)]
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # if its auto converted to a vector by R
      if(length(elevations) > nrow){
        elevations <- elevations[elevations != '']
        m <- matrix(data = as.numeric(elevations),
                    nrow = nrow,
                    ncol = ncol,
                    byrow = T)
      } else {
        elevations <- split(elevations, 1:length(elevations))
        elevations <- lapply(elevations, function(x){
          step1 <- strsplit(x, split = delimeter)[[1]]
          step2 <- step1[step1 != '']
          step2 <- as.numeric(step2)
        })
        names(elevations) <- NULL
        elevations <- do.call(rbind, elevations)
        m <- elevations
      }
      #-------------------------------------------------------------------------------
    }
    #-------------------------------------------------------------------------------
    
    layer_defs[[i]] <- m
  }
  names(layer_defs) <- c('model_top',
                         paste0('botm_layer_',c(0:(nlay-1))))
  #-------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  # calculate layer thickness
  layer_thick <- list()
  for(i in 1:(length(layer_defs)-1)){
    layer_thick[[i]] <- layer_defs[[i]] - layer_defs[[i+1]]
  }
  names(layer_thick) <- paste0('layer_thick_',c(0:(nlay-1)))
  #-------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  return(list(layer_defs,
              layer_thick))
  #-------------------------------------------------------------------------------
}
#-------------------------------------------------------------------------------