read_mfnwt_upw <- function(data_path,
                           nlay,
                           ncol,
                           nrow,
                           quasi_conf,
                           header       = FALSE,
                           delimeter    = ' ',
                           newline      = '\t',
                           start_index  = 8)
{
  #-------------------------------------------------------------------------------
  # logging
  cat(paste0('\nReading MFNWT upstream weighting file at:\n',
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
  # line 1
  preamble <- strsplit(raw_txt[2], delimeter)[[1]]
  preamble <- preamble[preamble != '']
  hdry     <- preamble[2]
  #-------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  # laydefs
  laydefs  <- data.frame(value = c(0:2),
                         def   = c('conf', 'unconf', 'mod-unconf'))
  laytyps  <- strsplit(raw_txt[3], delimeter)[[1]]
  laytyps  <- laytyps[laytyps != '']
  laytyps  <- vapply(laytyps, function(x){laydefs$def[which(laydefs$value == x)]}, character(1))
  nlay     <- length(laytyps)
  
  log_out  <- c()
  for(i in 1:length(laytyps)){
    log_out <- append(log_out,
                      paste0(i, paste0(rep(' ', nchar(laytyps[i]) - 1), collapse = ''), '|'))
  }
  log_out  <- paste0(log_out, collapse = '')
  cat('\nOf', length(laytyps), 'layers')
  cat('\nOf layer types: \n')
  cat(log_out, '\n')
  cat(paste0(paste0(laytyps, collapse = '|'),'|'), '\n')
  #-------------------------------------------------------------------------------
  
  
  #-------------------------------------------------------------------------------
  # chani
  chani <- strsplit(raw_txt[5], delimeter)[[1]]
  chani <- round(as.numeric(chani[chani != '']), 2)
  cat('\nOf horizontal anisotropy x-y of\n')
  chani_out <- vapply(chani, function(x){
    if(nchar(x) < nchar('chani ')){
      paste0(c(x, rep(' ', nchar('chani '))), collapse = '')
    } else {
      x
    }}, character(1))
  cat(paste0('chani ', 1:nlay, '|'), '\n')
  cat(paste0(paste0(chani_out, collapse = '| '),'|'), '\n')
  
  
  any_chani <- which(chani != 1)
  #-------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  # find all the variables that will be in the .upw file
  var_names <- c(paste0('HK_lay_', 1:nlay))
  if(length(any_chani) > 0){var_names <- append(var_names, paste0('CHANI_lay_', which(chani != 1)))}
  var_names <- append(var_names, paste0('VK_lay_', 1:nlay))
  var_names <- append(var_names, paste0('SS_lay_', 1:nlay))
  if(length(which(names(laytyps) != '0')) > 0){var_names <- append(var_names, paste0('SY_lay_', which(names(laytyps) != '0')))}
  if(length(which(quasi_conf != '0')) > 0){var_names <- append(var_names, paste0('VKCB_lay_', which(names(quasi_conf) != '0')))}
  
  cat('\nPredicted configuration of the upstream weighting file: \n')
  max_nchar <- c(0)
  for(i in 1:length(var_names)){
    if(nchar(var_names[i]) > max_nchar){max_nchar <- nchar(var_names[i])}
  }
  for(i in 1:length(var_names)){
    cat(paste0(c(var_names[i], rep(' ', max_nchar - nchar(var_names[i]))), collapse = ''), ':', nrow*ncol, 'entries\n')
  }
  #-------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  ncell           <- 0
  layer_processed <- 0
  counter         <- 0
  max_layers      <- length(var_names)
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
    # detect layer control information
    if(length(item) == 3 & is.na(as.numeric(item[2])) == TRUE){
     lay_mult <- as.numeric(strsplit(item[2], '(', fixed = TRUE)[[1]][1])
     index    <- index + 1
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
  return(var_list_final)
  #-------------------------------------------------------------------------------
}
#-------------------------------------------------------------------------------