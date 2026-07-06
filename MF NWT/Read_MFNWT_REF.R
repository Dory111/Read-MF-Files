read_mfnwt_ref <- function(data_path,
                           header = FALSE,
                           delimeter = ' ',
                           newline = '\t')
{
  #-------------------------------------------------------------------------------
  # load files and get their names
  files <- list.files(file.path(data_path),pattern = '*.ref$', full.names = TRUE)
  file_names <- list.files(file.path(data_path),pattern = '*.ref$', full.names = FALSE)
  file_names <- split(file_names, 1:length(file_names))
  file_names <- lapply(file_names, function(x){
    strsplit(x, '.ref')[[1]][1]
  })
  file_names <- as.vector(unlist(file_names))
  #-------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  data_list <- list()
  for(i in 1:length(files)){
    #-------------------------------------------------------------------------------
    # read in the raw text values
    raw_txt <- read.delim(files[i],
                          header = header,
                          sep = newline)
    raw_txt <- as.vector(unlist(raw_txt))
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    # get rid of the delimiting spaces
    txt_processed <- split(raw_txt,
                           1:length(raw_txt))
    txt_processed <- lapply(txt_processed, function(x){
      step1 <- strsplit(x, split = delimeter)[[1]]
      step2 <- step1[which(step1 != '')]
      step2 <- as.numeric(step2)
    })
    txt_processed <- do.call(rbind, txt_processed)
    data_list[[i]] <- txt_processed
    #-------------------------------------------------------------------------------
  }
  names(data_list) <- file_names
  #-------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  return(data_list)
  #-------------------------------------------------------------------------------
}
#-------------------------------------------------------------------------------