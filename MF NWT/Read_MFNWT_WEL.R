read_mfnwt_wel <- function(data_path,
                           header = TRUE,
                           timestep_separator = 'Stress Period',
                           delimiter = ' ',
                           newline = '\t')
{
  #-------------------------------------------------------------------------------
  # read in raw text
  raw_txt <- read.delim(data_path,
                        header = header,
                        sep = newline)
  raw_txt <- as.vector(unlist(raw_txt))
  #-------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  # find which indices correspond to the start of a timestep
  timestep_indices <- grep(timestep_separator, raw_txt[1:length(raw_txt)])
  n_timesteps <- raw_txt[tail(timestep_indices,1)]
  n_timesteps <- strsplit(n_timesteps,' ')[[1]]
  n_timesteps <- n_timesteps[which(n_timesteps != '')]
  n_timesteps <- as.numeric(tail(n_timesteps,1))
  #-------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  # if header was incorrectly captured, redefine the text to only after the
  # first timestep
  if(timestep_indices[1] != 1){
    raw_txt <- raw_txt[timestep_indices[1]:length(raw_txt)]
    timestep_indices <- grep(timestep_separator, raw_txt[1:length(raw_txt)])
  } else {}
  #-------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  # split into list for lapply
  raw_txt <- split(raw_txt, 1:length(raw_txt))
  #-------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  # remove stress period text, split data into layer, row, column, flux
  txt_processed <- base::lapply(1:length(raw_txt), function(i){
    if(!i %in% timestep_indices){
      step1 <- strsplit(raw_txt[[i]], delimiter)
      step2 <- which(is.na(as.numeric(step1[[1]])) == FALSE)
      keyp <- step1[[1]][step2]
      processed <- c(step1[[1]][step2],
                     paste(c(keyp[1],keyp[2],keyp[3]),collapse = ','))
    }
  })
  #-------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  # make more user friendly
  names(txt_processed) <- NULL
  txt_processed <- do.call(rbind,txt_processed)
  txt_processed <- as.data.frame(txt_processed)
  colnames(txt_processed) <- c('Layer','Row','Column','Flux','Key')
  #-------------------------------------------------------------------------------

  #-------------------------------------------------------------------------------
  # getting coordinates for each unique well
  unique_keys <- unique(txt_processed$Key)
  # adding underscore after the layer to account for each well with a repeated key
  # the names that are repeated is repeated_names
  # the number of times they are repeated is n_shared
  tab <- table(txt_processed$Key)
  repeated_names <- names(which(tab != n_timesteps))
  if(length(repeated_names) > 0){
    n_shared <- as.numeric(tab[which(tab != n_timesteps)])/n_timesteps
    #-------------------------------------------------------------------------------
    for(i in 1:length(repeated_names)){
      replacement_inds <- which(txt_processed$Key == repeated_names[i])
      keys <- txt_processed$Key[replacement_inds]
      orig_key <- strsplit(repeated_names[i],',')[[1]]
      for(j in 1:n_shared[i]){
        skipped_sequence <- seq(from = j, to = length(keys), by = n_shared[i])
        keys[skipped_sequence] <- paste0(c(paste0(c(orig_key[1],'_',j), collapse = ''),
                                         orig_key[2],orig_key[3]), collapse = ',')
      }
      txt_processed$Key[replacement_inds] <- keys
    }
    #-------------------------------------------------------------------------------
  }
  #-------------------------------------------------------------------------------
  

  #-------------------------------------------------------------------------------
  # getting keys and splitting data to capture each well
  unique_keys <- unique(txt_processed$Key)
  
  skeys <- split(unique_keys, 1:length(unique_keys))
  names(skeys) <- NULL
  well_coords <- lapply(skeys, function(x){
    step1 <- strsplit(x,',')[[1]]
    step2 <- strsplit(step1[1],'_')[[1]]
    if(length(step2) == 1){
      step2 <- c(step2[1],NA)
    }
    c(step2,step1[2:length(step1)])
  })
  names(well_coords) <- NULL
  well_coords <- do.call(rbind, well_coords)
  #-------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  # making lookup table out of keys
  m <- matrix(nrow = length(unique_keys),
              ncol = 5)
  m[,1] <- as.numeric(c(1:length(unique_keys)))
  m[,2] <- as.numeric(well_coords[,1])
  m[,3] <- as.numeric(well_coords[,2])
  m[,4] <- as.numeric(well_coords[,3])
  m[,5] <- as.numeric(well_coords[,4])
  m <- as.data.frame(m)
  colnames(m) <- c('W','Layer','N','Row','Column')
  well_lookup <- m
  #-------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  m <- matrix(nrow = length(unique_keys),
              ncol = as.numeric(table(txt_processed$Key)[1]))
  for(i in 1:length(unique_keys)){
    m[i,] <- as.numeric(txt_processed$Flux[txt_processed$Key == unique_keys[i]])
  }
  m <- as.data.frame(m)
  colnames(m) <- paste0('T',1:ncol(m))
  wide_pumping <- m
  #-------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  return(list(well_lookup,
              wide_pumping))
  #-------------------------------------------------------------------------------
}
#-------------------------------------------------------------------------------