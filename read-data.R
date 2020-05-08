
default_ukb_fields <- function(){
  
  SRcolumns<-c("20001","20002","20004","20003")
  SRdatecolumns <- c("20006","20008","20010")
  Othercolumns <- c("53","40000","40001","40002") 
  
  c(SRcolumns,SRdatecolumns,Othercolumns)
  
}

# Creates a variable name from the field description.
#
# @param data Field-to-description table from html file
#
description_to_name <-  function(Vct) {
  #https://github.com/kenhanscombe/ukbtools/blob/master/R/dataset.R
  name <- tolower(Vct) %>%
    gsub(" - ", "_", x = .) %>%
    gsub(" ", "_", x = .) %>%
    gsub("uses.*data.*coding.*simple_list.$", "", x = .) %>%
    gsub("uses.*data.*coding.*hierarchical_tree.", "", x = .) %>%
    gsub("uses.*data.*coding_[0-9]*", "", x = .) %>%
    gsub("[^[:alnum:][:space:]_]", "", x = .) %>%
    gsub("__*", "_", x = .)
  
  return(name)
}

ReplaceNAWithNearestNonNAOnTheLeft <- function(Vct) {
  N <- length(Vct)
  na.pos <- which(is.na(Vct))
  if (length(na.pos) %in% c(0, N)) {
    print("0 or N na's")
    return(Vct)
  }
  non.na.pos <- which(!is.na(Vct))
  intervals  <- findInterval(na.pos, non.na.pos,
                             all.inside = TRUE)
  left.pos   <- non.na.pos[pmax(1, intervals)]
  right.pos  <- non.na.pos[pmin(N, intervals+1)]
  left.dist  <- na.pos - left.pos
  right.dist <- right.pos - na.pos
  
  # Vct[na.pos] <- ifelse(left.dist <= right.dist,
  #                       Vct[left.pos], Vct[right.pos])
  Vct[na.pos] <- Vct[left.pos]
  return(Vct)
}

parse_html_tables <- function(x){
  df<-data.frame(x,stringsAsFactors = F)
  names(df) <- gsub(pattern = "NULL.","",names(df))
  return(df)
}

read_ukb_metadata <- function(fhtml="/Volumes/data/ukb/ukb38326.html") {
  
  # fileset="ukb38326"
  # Column types as described by UKB
  # http://biobank.ctsu.ox.ac.uk/crystal/help.cgi?cd=value_type
  col_type <- c(
    "Sequence" = "integer",
    "Integer" = "integer",
    "Categorical (single)" = "character",
    "Categorical (multiple)" = "character",
    "Continuous" = "double",
    "Text" = "character",
    "Date" = "character",
    "Time" = "character",
    "Compound" = "character",
    "Binary object" = "character",
    "Records" = "character",
    "Curve" = "character"
  )
  print("reading .html")
  lst_html_tables <- XML::readHTMLTable(fhtml)
  #lst_html_tables <- lapply(lst_html_tables,function(x) parse_html_tables(x) )
  
  df_meta <- lst_html_tables[[2]]
  df_meta$Type <- ReplaceNAWithNearestNonNAOnTheLeft(df_meta$Type)
  df_meta$Description <- ReplaceNAWithNearestNonNAOnTheLeft(df_meta$Description)
  df_meta$Description <- description_to_name(df_meta$Description)
  
  lookup.reference <- tibble::tibble(
    field.number = df_meta$Column,
    field.count = as.numeric(df_meta$Count),
    field.showcase = gsub("-.*$", "", df_meta[, "UDI"]),
    field.html = df_meta[, "UDI"],
    field.tab = paste("f.", gsub("-", ".", df_meta$UDI), sep = ""),
    field.description = df_meta[, "Description"],
    col.type = df_meta[, "Type"],
    col.name = ifelse(
      field.showcase == "eid",
      "eid",
      stringr::str_c(
        df_meta$Description, "_f",
        stringr::str_replace_all(field.html, c("-" = "_", "\\." = "_"))
      )
    ),
    fread_column_type = col_type[as.character(df_meta$Type)]
    
  )
  
  return(lookup.reference)
}





read_ukb_data <- function(fukb, 
                          dfhtml,
                          fields_to_keep = default_ukb_fields()) {
  
  if (!exists("n_threads")){n_threads=1}
  
  if(!any(fields_to_keep %in% "eid" )){
    fields_to_keep <- c("eid", fields_to_keep)
  }
  
  fields_to_keep.tab <- dfhtml[dfhtml$field.showcase %in% fields_to_keep & dfhtml$field.count!=0,]$field.tab
  fields_to_keep.classes <- dfhtml[dfhtml$field.showcase %in% fields_to_keep  & dfhtml$field.count!=0,]$fread_column_type
  
  
  if(length(fields_to_keep.classes) != length(fields_to_keep.tab)){
    print("error, fields_to_keep.classes doesnt match fields_to_keep")
    break
  }
  
  dfhtml.sub <- dfhtml[dfhtml$field.tab %in% fields_to_keep.tab,]
  
  freadcolclasses <- rep("NULL",nrow(dfhtml))
  freadcolclasses[which(dfhtml$field.tab %in% fields_to_keep.tab)] <- dfhtml[which(dfhtml$field.tab %in% fields_to_keep.tab),]$fread_column_type
  
  tic("fread data")
  df <- fread( paste0(fukb), header=T,
               colClasses = freadcolclasses,
               sep = "\t",
               showProgress = TRUE)
               # nThread = if(n_threads == "max") {
               #   max(1,parallel::detectCores()-1)
               # } else  if (n_threads == "dt") {
               #   data.table::getDTthreads()
               # } else if (is.numeric(n_threads)) {
               #   min(n_threads, parallel::detectCores())
               # }
  #,select=cols.to_keep,colClasses = col.classes.to_keep)
  print(format(object.size(df), units = "Gb"))
  toc()
  
  # library(vroom)
  # 
  # which(freadcolclasses !="NULL")
  # 
  # tic("vroom")
  # df <- vroom(fukb,delim = "\t",col_types="", col_select =which(freadcolclasses !="NULL")  , progress=F) # which(freadcolclasses !="NULL")  #fields_to_keep.tab
  # toc()
  
  
  return(df)
}
#col.classes[col.classes %in% "integer64"] <- "character" #integer64 not supported, unsupported in disk.frame? but not in data.table



