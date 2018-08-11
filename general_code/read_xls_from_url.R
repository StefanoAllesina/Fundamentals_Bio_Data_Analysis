library(readxl)
library(httr)

read_xlsx_from_url <- function(my_url){
  GET(my_url, write_disk(tf <- tempfile(fileext = ".xlsx"))) # download into temp file
  df <- read_excel(tf) # read xlsx
  unlink(tf) # delete temporary file
  return(df)
}

read_xls_from_url <- function(my_url){
  GET(my_url, write_disk(tf <- tempfile(fileext = ".xls"))) # download into temp file
  df <- read_excel(tf) # read xls
  unlink(tf) # delete temporary file
  return(df)
}
