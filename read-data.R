library(matrixStats)

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
  print(format(object.size(df), units = "Gb"))
  toc()
  
  # library(vroom)
  # which(freadcolclasses !="NULL")
  # tic("vroom")
  # df <- vroom(fukb,delim = "\t",col_types="", col_select =which(freadcolclasses !="NULL")  , progress=F) # which(freadcolclasses !="NULL")  #fields_to_keep.tab
  # toc()
  
  
  return(df)
}
#col.classes[col.classes %in% "integer64"] <- "character" #integer64 not supported, unsupported in disk.frame? but not in data.table





#### NOT WORKING ON MACBOOK... ON CLUSTER THIS WORKED; BUT THEN I CAN'T SEEM TO LOAD THE DATA.. 
convert_ukb_to_diskframe <- function(fukbtab,fhtml,outdir="diskframe/",rows_to_read=20000,ram_size=8) {
  
  library(disk.frame)
  dfhtml <- read_ukb_metadata(fhtml)
  rows_to_read=20000
  ram_size=8
  dfukb <- csv_to_disk.frame(fukbtab,
                             nchunks = recommend_nchunks(sum(file.size(fukb)),ram_size=ram_size),
                             in_chunk_size = rows_to_read,
                             outdir = outdir,
                             colClasses = as.vector(dfhtml$fread_column_type),
                             col.names = dfhtml$field.tab, sep="\t")
}
# 
# fukb="/data/pg-exp_cardio/UKBIO_database/12010_2019_10_31-38326/ukb38326.tab"
# fhtml="/data/pg-exp_cardio/UKBIO_database/12010_2019_10_31-38326/ukb38326.html"
# outdir='/data/pg-exp_cardio/UKBIO_database/12010_2019_10_31-38326/diskframe' #'/Volumes/data/ukb/diskframe'
# #dfukb <- csv_to_disk.frame(fukb, in_chunk_size = 100,colClasses = freadcolclasses)
# 
# convert_ukb_to_diskframe(fukb,fhtml,outdir)
# 
# dfhtml[dfhtml$field.tab %in% 'f.4258.0.1',]


read_hesin_data <- function(fhesin, fhesin_diag,fhesin_oper){
  
  
  # read hesin
  print("read hesin")
  dfhesin <- (fread(fhesin,header=T,sep="\t", stringsAsFactors=FALSE, na.strings=""))
  
  # read diag
  print("read diag")
  dfdiag <- (fread(fhesin_diag,header=T,sep="\t", stringsAsFactors=FALSE, na.strings=""))
  
  # read oper
  print("read oper")
  dfoper <- (fread(fhesin_oper,header=T,sep="\t", stringsAsFactors=FALSE, na.strings=""))
  dfoper$opdate <- as.integer(format(as.Date(as.character(dfoper$opdate),format="%d/%m/%Y"), "%Y%m%d")) # same format as dfhesin
  
  print("merging hesin + diagnosis")
  dfhesin_diag <- merge(dfhesin,dfdiag,by = c("eid","ins_index"),all=T)
  dfhesin_diag$eventdate <- dfhesin_diag$epistart
  dfhesin_diag[is.na(dfhesin_diag$eventdate),"eventdate"] <- dfhesin_diag[is.na(dfhesin_diag$eventdate),"admidate"]
  dfhesin_diag[is.na(dfhesin_diag$epiend),"epiend"] <- dfhesin_diag[is.na(dfhesin_diag$epiend),"disdate"]
  dfhesin_diag$epidur <- as.numeric(dfhesin_diag$eventdate - dfhesin_diag$epistart)
  dfhesin_diag$epidur[dfhesin_diag$epidur <0,] <- NA
  
  dfhesin_diag$event <- 1 
  dfhesin_diag[is.na(dfhesin_diag$eventdate)]$event <- 0
  dfhesin_diag <- dfhesin_diag[, event:=as.integer(event)]
  
  
  print("merging hesin + operation") # for duration, take episode duration. 
  dfhesin_oper <- merge(dfhesin,dfoper,by = c("eid","ins_index"),all=T)
  dfhesin_oper$eventdate <- dfhesin_oper$opdate
  dfhesin_diag[is.na(dfhesin_diag$epistart),"epistart"] <- dfhesin_diag[is.na(dfhesin_diag$epistart),"admidate"]
  dfhesin_oper[is.na(dfhesin_oper$eventdate),"eventdate"] <- dfhesin_oper[is.na(dfhesin_oper$eventdate),"epistart"]
  dfhesin_oper[is.na(dfhesin_oper$eventdate),"eventdate"] <- dfhesin_oper[is.na(dfhesin_oper$eventdate),"admidate"]
  dfhesin_oper[is.na(dfhesin_oper$epiend),"epiend"] <- dfhesin_oper[is.na(dfhesin_oper$epiend),"disdate"]
  dfhesin_oper <- dfhesin_oper[!is.na(dfhesin_oper$eventdate),]
  
  dfhesin_oper$eventdate <- as.Date(as.character(dfhesin_oper$eventdate),format="%Y%m%d")
  dfhesin_oper$epistart <- as.Date(as.character(dfhesin_oper$epistart),format="%Y%m%d")
  dfhesin_oper$epiend <- as.Date(as.character(dfhesin_oper$epiend),format="%Y%m%d")
  dfhesin_oper$epidur <- as.numeric(dfhesin_oper$epiend-dfhesin_oper$epistart)
  
  dfhesin_oper[dfhesin_oper$epidur <0,]$epidur <- NA
  
  dfhesin_oper$event <- 1 
  dfhesin_oper[is.na(dfhesin_oper$eventdate)]$event <- 0
  dfhesin_oper <- dfhesin_oper[, event:=as.integer(event)]
  
  
  tte.oper3.primary <- dfhesin_oper %>% filter(level==1 & !is.na(oper3))  %>% select(eid,epistart,epidur,oper3,event)  %>% rename(f.eid=eid,eventdate = epistart,epidur=epidur,code = oper3,event=event)  %>% as.data.table()
  tte.oper4.primary <- dfhesin_oper %>% filter(level==1 & !is.na(oper4))  %>% select(eid,epistart,epidur,oper4,event)  %>% rename(f.eid=eid,eventdate = epistart,epidur=epidur,code = oper4,event=event)  %>% as.data.table()
  tte.icd10.primary <- dfhesin_diag %>% filter( level==1 & !is.na(diag_icd10))  %>% select(eid,epistart,epidur,diag_icd10,event)  %>% rename(f.eid=eid,eventdate = epistart,epidur=epidur,code = diag_icd10,event=event)  %>% as.data.table()
  tte.icd9.primary <- dfhesin_diag %>% filter( level==1 & !is.na(diag_icd9))  %>% select(eid,epistart,epidur,diag_icd9,event)  %>% rename(f.eid=eid,eventdate = epistart,epidur=epidur,code = diag_icd9,event=event)  %>% as.data.table()
  
  tte.oper3.secondary <- dfhesin_oper %>% filter(level==2 & !is.na(oper3))  %>% select(eid,epistart,epidur,oper3,event)  %>% rename(f.eid=eid,eventdate = epistart,epidur=epidur,code = oper3,event=event)  %>% as.data.table()
  tte.oper4.secondary <- dfhesin_oper %>% filter(level==2 & !is.na(oper4))  %>% select(eid,epistart,epidur,oper4,event)  %>% rename(f.eid=eid,eventdate = epistart,epidur=epidur,code = oper4,event=event)  %>% as.data.table()
  tte.icd10.secondary <- dfhesin_diag %>% filter( level==2 & !is.na(diag_icd10))  %>% select(eid,epistart,epidur,diag_icd10,event)  %>% rename(f.eid=eid,eventdate = epistart,epidur=epidur,code = diag_icd10,event=event)  %>% as.data.table()
  tte.icd9.secondary <- dfhesin_diag %>% filter( level==2 & !is.na(diag_icd9))  %>% select(eid,epistart,epidur,diag_icd9,event)  %>% rename(f.eid=eid,eventdate = epistart,epidur=epidur,code = diag_icd9,event=event)  %>% as.data.table()
  
  
  
  
  setkey(tte.oper3.primary,f.eid)    
  setkey(tte.oper4.primary,f.eid)    
  setkey(tte.icd10.primary,f.eid)    
  setkey(tte.icd9.primary,f.eid)    
  
  setkey(tte.oper3.secondary,f.eid)    
  setkey(tte.oper4.secondary,f.eid)    
  setkey(tte.icd10.secondary,f.eid)    
  setkey(tte.icd9.secondary,f.eid)    
  
  lst <- list(tte.oper3.primary = tte.oper3.primary, 
              tte.oper4.primary = tte.oper4.primary,
              tte.icd10.primary = tte.icd10.primary, 
              tte.icd9.primary = tte.icd9.primary,
              
              tte.oper3.secondary = tte.oper3.secondary, 
              tte.oper4.secondary = tte.oper4.secondary,
              tte.icd10.secondary = tte.icd10.secondary, 
              tte.icd9.secondary = tte.icd9.secondary)
  return(lst)
  
}


