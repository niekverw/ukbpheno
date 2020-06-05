# function to load data. 

convert_year_to_date <- function(year){
  #https://stackoverflow.com/questions/29697436/how-to-convert-decimal-date-format-e-g-2011-580-to-normal-date-format
  
  inotna<-which(!is.na(year))
  out <- rep(NA,length(year))
  
  year <- year[!is.na(year)]
  start <- as.POSIXct(paste0(trunc(year),  "/01/01"), tz="UTC") 
  end   <- as.POSIXct(paste0(trunc(year)+1,"/01/01"), tz="UTC") 
  date <- start + (difftime(end, start, units="secs") * (year - trunc(year)))
  date <- format(date, format='%Y-%m-%d')
  out[inotna] <- date
  return(out)
}

convert_col_to_integer <- function(col){
  n <- length(which(is.na(col)))
  col_int <-suppressWarnings(as.integer(col)) ## warning 
  if (length(which(is.na(col_int))) > n){
    return(col)
  } else{
    return(col_int)
  }
}

convert_nurseinterview_to_episodedata <- function(df,field_sr_diagnosis = "20002",field_sr_date = "20008",field_sr_date_type="interpolated_year",qc_treshold_year=10){
  # TODO: IF NA, then list first visit answered yes. 
  # 
  # df = lst$df_sr # data.table
  # field_sr_diagnosis = "20002"
  # field_sr_date = NULL # "20008"# "20008" # interpolated year
  # field_sr_date_type="interpolated_year"
  # field_visit_date="53"
  # qc_treshold_year=10
  
  # 
  tic(paste("convert_nurseinterview_to_episodedata: ",field_sr_diagnosis))
  
  if(!is.null(field_sr_date)) { if(field_sr_date==""){field_sr_date=NULL}} 
  if(is.null(field_sr_date))  { print("field_sr_date == NULL; qc_treshold_year and field_sr_date_type will not be used.") }
  daysinyear=365.25
  field_visit_date="53"
  identifierfield = names(df)[grepl("eid", names(df))]
  visitdatefields = names(df)[grepl(paste0("[^0-9]",field_visit_date,"[^0-9]"), names(df))]
  srdiagnosisfields = names(df)[grepl( paste0("[^0-9]",field_sr_diagnosis,"[^0-9]"), names(df))]
  if(!is.null(field_sr_date)){srdiagnosisdatefields = names(df)[grepl(field_sr_date, names(df))]} else srdiagnosisdatefields=NULL
  
  visits = length(visitdatefields) #sum(grepl("53_", names(df)))
  
  # only need n_eid, visit dates and diag-codes + age-of-diag
  columns_to_keep = c(identifierfield,
                      visitdatefields,
                      srdiagnosisfields, 
                      srdiagnosisdatefields
  )
  df_ <- df[,columns_to_keep,with=FALSE]
  df_$dummy <- NA
  df_out <-  matrix(ncol=5, nrow=0) # initiate output 
  
  for (v in 0:(visits-1)){ # for each visit, 
    print(paste0("querying visit ",v))
    diagfields = names(df_)[grepl(paste0("[^0-9]",field_sr_diagnosis,"[^0-9]",v),names(df_))]
    if(length(diagfields)==0){print(paste0("no data on visit ",v));next}
    if(!is.null(field_sr_date)){diagdatefields = names(df_)[grepl(paste0("[^0-9]",field_sr_date,"[^0-9]",v),names(df_))]}
    visitdatefield = visitdatefields[v+1]
    for (i in 1:length(diagfields)) { # for each occurence of diagfield, find the corresponding age and convert it to  date - code and rbind() to df_out. 
      #agefield = paste0("age_",v)
      diagfield = diagfields[i]
      if(!is.null(field_sr_date)){ 
        diagdatefield = diagdatefields[i]
        if( length(diagdatefields) == 1) { diagdatefield = diagdatefields } # in case of death (40001/2)
      } else {
          diagdatefield = "dummy" # e.g. incase of medication
          }

      #print(paste0((diagfield), " - ", diagdatefield))
      if(all(is.na(df_[,diagfield,with=FALSE]))){next}
      df_sub <- df_[!is.na(get(diagfield) ),c(identifierfield,diagfield,diagdatefield,visitdatefield),with=FALSE]
      df_sub$visit <- v
      names(df_sub) <- c(identifierfield,"code","eventdate","visitdate","visit")
      df_out <- rbind(df_out,as.matrix(df_sub[,c(identifierfield,"code","eventdate","visit","visitdate"),with=FALSE]))
    }
    
  }
  
  print("convert to dataframe")
  df_out <- data.table(df_out,stringsAsFactors=F)
  df_out <- df_out[, visitdate:=as.Date(visitdate)]
  #df_out <- df_out[, code:=as.character(code)] #convert_col_to_integer(df_out$code) # df_out[, code:=as.integer(code)]
  df_out[, code:=lapply(.SD, trimws), .SDcols = "code"]
  df_out <- df_out[!code %in% "99999"]
  df_out$code <- convert_col_to_integer(df_out$code)
  
  
  if(!is.null(field_sr_date) & (field_sr_date_type=="interpolated_year"| field_sr_date_type=="interpolated_age" |  field_sr_date_type=="date" )) {

    if(field_sr_date_type=="interpolated_year") {
      df_out <- df_out[, eventdate:=as.numeric(eventdate)] ## as number. interpolated year 
      df_out[eventdate <0,'eventdate']<-NA
      df_out$eventdate <- as.Date(convert_year_to_date(df_out$eventdate))
    } else if (field_sr_date_type=="interpolated_age"){
      df_out <- df_out[, eventdate:=as.numeric(eventdate)] ## as number. interpolated age. 
      df_out[eventdate <0,'eventdate']<-NA
      df_out$eventdate = df_out[,"visitdate"] - (df_out[,"eventdate"]*daysinyear)
    } else if (field_sr_date_type=="date"){
      df_out = df_out[, eventdate:=as.Date(eventdate)]
    }
    # deduplicate, min/max/mean/sd <- not very efficient? 
    print("deduplicate")
    dfout_extrastats<- df_out %>% group_by(!!as.name(identifierfield),code) %>%
      mutate(mindt = min(eventdate, na.rm = TRUE),maxdt = max(eventdate, na.rm = TRUE),meandt = mean(eventdate, na.rm = TRUE))
    
    dfout_extrastats$diffdt <- (dfout_extrastats$maxdt - dfout_extrastats$mindt)/daysinyear
    dfout_extrastats[dfout_extrastats$diffdt>qc_treshold_year ,"meandt"] <- NA
    dfout_extrastats[dfout_extrastats$diffdt > 0,]
    df_out$eventdate <- dfout_extrastats$meandt
    df_out <- df_out[!duplicated(df_out[,c(identifierfield,"code","eventdate"),with=FALSE]),] #sorted on visit, so first occurence is always first visit. 
    
  }  else {
    df_out <- df_out[, eventdate:=as.Date(eventdate)]
  }
  
  # record which can be set as an event or not (when no event_date is reported, only visit)
  df_out$event <- 1
  df_out <- df_out[, event:=as.integer(event)]
  df_out[is.na(df_out$eventdate)]$event <- 0
  df_out[is.na(df_out$eventdate)]$eventdate <- df_out[is.na(df_out$eventdate)]$visitdate
  df_out <- df_out[,c(identifierfield,"code","eventdate","event"),with=FALSE]
  setkey(df_out,f.eid)    
  
  #head(df_out)
  gc()
  print(format(object.size(df_out), units = "Mb"))
  
  toc()
  return(df_out)
}

