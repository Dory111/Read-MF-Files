read_mfnwt_bas <- function(data_path,
                           nlay,
                           ncol,
                           nrow,
                           header       = FALSE,
                           delimeter    = ' ',
                           newline      = '\t',
                           start_index  = 5)
{
  #-------------------------------------------------------------------------------
  # logging
  cat(paste0('\nReading MFNWT basic file at:\n',
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
  hnoflo    <- c()
  var_names <- c(paste0('ibound_lay_', 1:nlay), paste0('strt_lay_', 1:nlay))
  #-------------------------------------------------------------------------------

  #-------------------------------------------------------------------------------
  ncell           <- 0
  layer_processed <- 0
  counter         <- 0
  max_layers      <- length(var_names) + 1
  failed_max      <- 3
  failed_attempt  <- 0
  
  var_list <- list()
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
    if(layer_processed == nlay){
      index           <- index + 1
      hnoflo          <- as.numeric(item)
      layer_processed <- layer_processed + 1
      next
    }
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    # detect layer control information
    if(length(item) == 3 & is.na(as.numeric(item[2])) == TRUE){
      if(layer_processed < nlay+1){
        lay_mult <- 1
        index    <- index + 1
      } else {
        lay_mult <- as.numeric(strsplit(item[2], '(', fixed = TRUE)[[1]][1])
        index    <- index + 1
      }
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
      holdover <- c()
    }
    var_list[[counter]] <- item
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
        item     <- var_list[[counter]]
        holdover <- append(holdover, item[c((length(item) - overflow + 1):length(item))])
        item     <- item[-c((length(item) - overflow + 1):length(item))]
        var_list[[counter]] <- item
        
        lay_sep         <- append(lay_sep, counter)
        layer_processed <- layer_processed + 1
        ncell           <- overflow
        
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
  lay_sep    <- append(0, lay_sep)
  #-------------------------------------------------------------------------------
  

  #-------------------------------------------------------------------------------
  var_list_final <- list()
  for(i in 1:(length(lay_sep)-1)){
    blank      <- matrix(data  = as.vector(unlist(var_list[c((lay_sep[i]+1):lay_sep[i+1])])),
                         nrow  = nrow,
                         ncol  = ncol,
                         byrow = FALSE)
    var_list_final[[i]] <- blank
  }
  names(var_list_final) <- var_names
  #-------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  return(list(var_list_final,
              hnoflo))
  #-------------------------------------------------------------------------------
}
#-------------------------------------------------------------------------------