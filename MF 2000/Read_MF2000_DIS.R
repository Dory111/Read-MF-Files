read_mf2000_dis <- function(data_path,
                            header       = FALSE,
                            delimeter    = ' ',
                            newline      = '\t',
                            start_index  = 6)
{
  #-------------------------------------------------------------------------------
  # logging
  cat(paste0('\nReading MF2000 discretization file at:\n',
             data_path))
  #-------------------------------------------------------------------------------
  
  
  #-------------------------------------------------------------------------------
  if(file.exists(data_path) == FALSE){
    stop(paste0('\nFile at:\n',
                data_path,'\n',
                'Does not exist...exiting program'))
  }
  #-------------------------------------------------------------------------------
  
  
  #-------------------------------------------------------------------------------
  # read in the raw text values
  raw_txt <- read.delim(file   = data_path,
                        header = header,
                        sep    = newline)
  raw_txt <- as.vector(unlist(raw_txt))
  #-------------------------------------------------------------------------------

  
  
  
  #-------------------------------------------------------------------------------
  # get the dimensions of the grid
  # MF 2000 dis file organized like
  # NLAY NROW NCOL NPER ITMUNI LENUNI
  itmuni_def  <- data.frame(value = c(0:5),
                            units = c('undef', 'seconds', 'minutes', 'hours', 'days', 'years'))
  lenuni_def  <- data.frame(value = c(0:3),
                            units = c('undef',    'feet',  'meters', 'centimeters'))
  
  dims <- strsplit(raw_txt[2], delimeter)[[1]]
  dims <- dims[dims != '']
  nlay <- as.numeric(dims[1])
  nrow <- as.numeric(dims[2])
  ncol <- as.numeric(dims[3])
  nper <- as.numeric(dims[4])
  end_dims <- grep('TR',raw_txt)[1]
  ncell <- nrow*ncol
  #-------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  log_dims <- c(nlay,nrow,ncol,nper)
  for(i in 1:length(log_dims)){
    if(nchar(log_dims[i]) < 4){
      log_dims[i] <- paste0(log_dims[i], paste0(rep(' ', 4 - nchar(log_dims[i])), collapse = ''))
    }
  }
  #-------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  # logging
  cat(paste0('\nModel described as:\n',
             'NLAY NROW NCOL NPER\n',
             log_dims[1],' ', log_dims[2], ' ', log_dims[3], ' ', log_dims[4]))
  #-------------------------------------------------------------------------------
  
  
  #-------------------------------------------------------------------------------
  # units
  itmuni <- itmuni_def$units[which(itmuni_def$value == as.numeric(dims[5]))]
  lenuni <- lenuni_def$units[which(lenuni_def$value == as.numeric(dims[6]))]
  #-------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  # logging
  confining_layers <- strsplit(raw_txt[3], delimeter)[[1]]
  confining_layers <- confining_layers[confining_layers != '']
  if(all(confining_layers == '0')){
    cat('\nWith no confining layers')
  } else {
    ind <- which(confining_layers != '0')
    cat('\nWith no confining beds below layer(s):', paste0(ind, collapse = ','))
  }
  #-------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  # logging
  cat(paste0('\nOf length units: ', lenuni))
  cat(paste0('\nAnd time units: ', itmuni))
  #-------------------------------------------------------------------------------
  
  
  #-------------------------------------------------------------------------------
  ncell           <- 0
  layer_processed <- 0
  counter         <- 0
  max_layers      <- nlay+1 # account that first array is just the ground surface
  failed_max      <- 3
  failed_attempt  <- 0
  
  elev_list <- list()
  holdover  <- c()
  done      <- FALSE
  index     <- start_index
  lay_sep   <- c()
  lay_mult  <- 1
  while(done == FALSE){
    #-------------------------------------------------------------------------------
    counter <- counter + 1
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    item    <- strsplit(raw_txt[index], delimeter)[[1]]
    item    <- item[item != '']
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    # detect beginning of stress period information by character presence
    if(length(item) == 4 & is.na(as.numeric(item[4])) == TRUE){
      cat(paste0('\nStress period beginning at: ', index))
      cat(paste0('\nIf you believe this to be an error please check the .dis file'))
      lay_sep <- append(lay_sep, counter)
      done    <- TRUE
      next
    }
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    # detect layer control information
    if(length(item) == 3 & is.na(as.numeric(item[2])) == TRUE){
     lay_mult <- as.numeric(strsplit(item[2], '(', fixed = TRUE)[[1]][1])
     index   <- index + 1
     next
    }
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    item    <- as.numeric(item) * lay_mult
    ncell   <- ncell + length(item)
    index   <- index + 1
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    # was there an overflow on last row read?
    if(length(holdover) > 0){
      item     <- append(holdover, item)
      print(holdover)
      holdover <- c()
    }
    elev_list[[counter]] <- item
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    # detect end of file by NA presence
    if(failed_attempt >= failed_max){
      lay_sep <- append(lay_sep, counter - failed_max)
      done    <- TRUE
    }
    if(all(is.na(item)) == TRUE){
      failed_attempt <- failed_attempt + 1
      next
    }
    #-------------------------------------------------------------------------------
    
    
    
    #-------------------------------------------------------------------------------
    if(ncell >= nrow*ncol){
      if(ncell > nrow*ncol){
        overflow <- ncell - nrow*ncol
        item     <- elev_list[[counter]]
        holdover <- append(holdover, item[c((length(item) - overflow + 1):length(item))])
        item     <- item[-c((length(item) - overflow + 1):length(item))]
        elev_list[[counter]] <- item
        
        lay_sep         <- append(lay_sep, counter)
        layer_processed <- layer_processed + 1
        ncell           <- 0
        
      } else {
        lay_sep         <- append(lay_sep, counter)
        layer_processed <- layer_processed + 1
        ncell           <- 0
      }
      
      if(layer_processed == max_layers){
        done <- TRUE
      }
    }
    #-------------------------------------------------------------------------------
  }
  #-------------------------------------------------------------------------------
  

  #-------------------------------------------------------------------------------
  lay_sep    <- append(0, lay_sep)
  layer_defs <- list()
  for(i in 1:(length(lay_sep)-1)){
    blank      <- matrix(data  = as.vector(unlist(elev_list[c((lay_sep[i]+1):lay_sep[i+1])])),
                         nrow  = nrow,
                         ncol  = ncol,
                         byrow = TRUE)
    layer_defs[[i]] <- blank
  }
  names(layer_defs) <- paste0('bot_lay_', c(0:nlay))
  layer_thick <- list()
  for(i in 1:(length(layer_defs)-1)){
    layer_thick[[i]] <- matrix(data  = as.numeric(layer_defs[[i]]) - as.numeric(layer_defs[[i+1]]),
                               nrow  = nrow,
                               ncol  = ncol,
                               byrow = FALSE)
  }
  names(layer_thick) <- paste0('thick_lay_', c(1:nlay))
  #-------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  return(list(layer_defs,
              layer_thick))
  #-------------------------------------------------------------------------------
}
#-------------------------------------------------------------------------------