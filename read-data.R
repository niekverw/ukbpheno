
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
  
  lst_html_tables <- XML::readHTMLTable(html_file)
  #lst_html_tables <- lapply(lst_html_tables,function(x) parse_html_tables(x) )
  
  df_meta <- lst_html_tables[[2]]
  df_meta$Type <- ReplaceNAWithNearestNonNAOnTheLeft(df_meta$Type)
  df_meta$Description <- ReplaceNAWithNearestNonNAOnTheLeft(df_meta$Description)
  df_meta$Description <- description_to_name(df_meta$Description)
  df_meta$field.tab <- paste("f.", gsub("-", ".", df_meta$UDI), sep = "")
  df_meta$field.showcase <- gsub("-.*$", "", df_meta[, "UDI"])
  
  lookup.reference <- tibble::tibble(
      field.showcase = gsub("-.*$", "", df_meta[, "UDI"]),
      field.html = df_meta[, "UDI"],
      field.tab = df_meta[,"field.tab"],
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
      
    )
    
    return(lookup.reference)
}





read_ukb_data <- function(f, 
                          fields_to_keep = c("20002","20008:numeric", #), #"20009", #non cancer codees, interpolated year, age
                                             "20001","20006:numeric", #"20007", # cancer codes , interpolated year, 
                                             "20004", "20010:numeric", #"20011"	
                                             "20003")
){
  
  ## ASSUMING INTEGER IF CLASS NOT PROVIDED
  fields_to_keep.classes =unlist( lapply( strsplit(fields_to_keep,split=":"),function(x) {if(length(x)==1){x=c(x,"integer")};return(x[[2]])} ))
  fields_to_keep = unlist(lapply(strsplit(fields_to_keep,split=":"),function(x) x[[1]]))
  
  if(!any(fields_to_keep %in% "53" )){
    fields_to_keep <- c("53", fields_to_keep)
    fields_to_keep.classes = c("character",fields_to_keep.classes)
  }
  
  if(length(fields_to_keep.classes) != length(fields_to_keep)){
    print("error, fields_to_keep.classes doesnt match fields_to_keep")
    break
  }
  
  
  df_header <- names(fread( paste0('head -1 ',f ))) # read header to find positions and match colclasses.
  df_header.colclasses <- rep("NULL",length(df_header))
  for (i in 1:length(fields_to_keep)){
    col <- df_header[grepl(paste(paste0("[^0-9]",fields_to_keep[i],"([^0-9])"),collapse="|"), df_header )]
    df_header.colclasses[df_header %in% col] <- rep(fields_to_keep.classes[i],length(col))
  }
  
  df_header.colclasses[df_header %in% "f.eid"] <- "character" # R
  df_header.colclasses[df_header %in% "n_eid"] <- "character" # Stata
  
  
  df <- fread( paste0(f ),header=T,colClasses = df_header.colclasses)#,select=cols.to_keep,colClasses = col.classes.to_keep)
  print(format(object.size(df), units = "Gb"))
  
  return(df)
}
#col.classes[col.classes %in% "integer64"] <- "character" #integer64 not supported, unsupported in disk.frame? but not in data.table



read_ukb_data_2 <- function(f,  fields_to_keep = c("20002","20008")){
  
  ## ASSUMING INTEGER IF CLASS NOT PROVIDED
  fields_to_keep.classes =unlist( lapply( strsplit(fields_to_keep,split=":"),function(x) {if(length(x)==1){x=c(x,"integer")};return(x[[2]])} ))
  fields_to_keep = unlist(lapply(strsplit(fields_to_keep,split=":"),function(x) x[[1]]))
  
  if(!any(fields_to_keep %in% "53" )){
    fields_to_keep <- c("53", fields_to_keep)
    fields_to_keep.classes = c("character",fields_to_keep.classes)
  }
  
  if(length(fields_to_keep.classes) != length(fields_to_keep)){
    print("error, fields_to_keep.classes doesnt match fields_to_keep")
    break
  }
  
  
  df_header <- names(fread( paste0('head -1 ',f ))) # read header to find positions and match colclasses.
  df_header.colclasses <- rep("NULL",length(df_header))
  for (i in 1:length(fields_to_keep)){
    col <- df_header[grepl(paste(paste0("[^0-9]",fields_to_keep[i],"([^0-9])"),collapse="|"), df_header )]
    df_header.colclasses[df_header %in% col] <- rep("character",length(col))
  }
  
  df_header.colclasses[df_header %in% "f.eid"] <- "character" # R
  df_header.colclasses[df_header %in% "n_eid"] <- "character" # Stata
  
  
  df <- fread( paste0(f ),header=T,colClasses = df_header.colclasses)#,select=cols.to_keep,colClasses = col.classes.to_keep)
  print(format(object.size(df), units = "Gb"))
  
  return(df)
}
#col.classes[col.classes %in% "integer64"] <- "character" #integer64 not supported, unsupported in disk.frame? but not in data.table

