
#' Convert year in decimal unit to date of format yyyy-mm-dd
#'
#' This function takes numerical year and converts into a date format.
#' @param year
#' @return  character vector of date in format yyyy-mm-dd
#' @keywords auxiliary
#' @export
#' @examples
#' convert_year_to_date(1999.25)
#' convert_year_to_date(c(1999.25,2020))
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



#' Convert numbers in character to integer
#'
#' This function takes numbers as string and converts into type integer, digits after decimal point is discarded
#' @param col
#' @return  character vector of number as integer
#' @keywords auxiliary
#' @export
#' @examples
#'convert_col_to_integer("2020")
#'convert_col_to_integer(c("2047","2020"))
convert_col_to_integer <- function(col){
  n <- length(which(is.na(col)))
  col_int <-suppressWarnings(as.integer(col)) ## warning
  if (length(which(is.na(col_int))) > n){
    return(col)
  } else{
    return(col_int)
  }
}


#' Convert self-reported illness codes from verbal interview with corresponding time of diagnosis to episodes format
#'
#' This function takes data fields containing illness codes/time of diagnosis distributed in multiple with each row representing one individual. It processed the data and return in episodes of event for all individuals
#' @param df dataframe containing the fields
#' @param field_sr_diagnosis data field number for illness code, default: 20002
#' @param field_sr_date data field number for corresponding time of diagnosis, default: 20008
#' @param field_sr_date_type type of time of diagnosis - interpolated_year/interpolated_year/interpolated_year, default:interpolated_year
#' @param qc_threshold_year  in case of multiple episodes for same illness code, if time difference is larger than qc threshold between the oldest and newest episodes, eventdate will be discarded i.e. set to NA, default:10
#' @param event_code code number to indicate if the episodes is a true event 1, self-reported/interpolated event 2 , for which a mean date is taken as event date in case of multiple occurence, or non-event 0 , for which date does not correspond to occurence of illness code. These codes are used in time-to-event analysis, default: 2
#' @return  a data.table object with all episodes
#' @keywords self-reported data
#' @export
#' @examples
#' convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20001",field_sr_date = "20006",qc_treshold_year = 10)
#' convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20002",field_sr_date = "20008",qc_treshold_year = 10)
convert_nurseinterview_to_episodedata <- function(df,field_sr_diagnosis = "20002",field_sr_date = "20008",field_sr_date_type="interpolated_year",qc_threshold_year=10,event_code=2,codetype="numeric"){
  # TODO: IF NA, then list first visit answered yes.
  #
  # df = lst$df_sr # data.table
  # field_sr_diagnosis = "20002"
  # field_sr_date = "20008" #NULL # "20008"# "20008" # interpolated year
  # field_sr_date_type="interpolated_year"
  # field_visit_date="53"
  # qc_threshold_year=10

  # df=dfukb
  tictoc::tic(paste("convert_nurseinterview_to_episodedata: ",field_sr_diagnosis))

  if(!is.null(field_sr_date)) { if(field_sr_date==""){field_sr_date=NULL}}
  if(is.null(field_sr_date))  { print("field_sr_date == NULL; qc_threshold_year and field_sr_date_type will not be used.") }
  daysinyear=365.25
  field_visit_date="53"
  # vector with name of the identifier col
  # identifierfield = identifier #names(df)[grepl("eid", names(df))]
  #  vector with names of all visit cols : "f.53.0.0" "f.53.1.0" "f.53.2.0" "f.53.3.0"
  visitdatefields = names(df)[grepl(paste0("[^0-9]",field_visit_date,"[^0-9]"), names(df))]

  # default:20002, Non-cancer illness code, self-reported
  srdiagnosisfields = names(df)[grepl( paste0("[^0-9]",field_sr_diagnosis,"[^0-9]"), names(df))]
  # default:20008 , 	Interpolated Year when non-cancer illness first diagnosed
  if(!is.null(field_sr_date)){srdiagnosisdatefields = names(df)[grepl(field_sr_date, names(df))]} else srdiagnosisdatefields=NULL
  visits = length(visitdatefields) #sum(grepl("53_", names(df)))

  # need for calculating diagdate from age of diagnosis
  field_birth_year ="34"
  field_birth_month="52"
  birthyearfield = names(df)[grepl(paste0("[^0-9]",field_birth_year,"[^0-9]"), names(df))]
  birthmonthfield = names(df)[grepl(paste0("[^0-9]",field_birth_month,"[^0-9]"), names(df))]
  # identifierfield = names(df)[grepl("eid", names(df))]
  identifierfield = names(df)[grepl("identifier", names(df))]



  # only need n_eid, visit dates and diag-codes + age-of-diag
  columns_to_keep = c(identifierfield,
                      visitdatefields,
                      srdiagnosisfields,
                      srdiagnosisdatefields,
                      birthyearfield,
                      birthmonthfield
  )
  # data.table - Setting with = FALSE disables the ability to refer to columns as if they are variables
  df_ <- df[,columns_to_keep,with=FALSE]
  df_$dummy <- NA
  df_$birthdt = as.Date(as.character(paste(df_[[birthyearfield]],df_[[birthmonthfield]], 1, sep = "-")),format="%Y-%m-%d")

  df_out <-  matrix(ncol=6, nrow=0) # initiate output

  for (v in 0:(visits-1)){ # for each visit,
    message(paste0("querying visit ",v))
    diagfields = names(df_)[grepl(paste0("[^0-9]",field_sr_diagnosis,"[^0-9]",v),names(df_))]
    if(length(diagfields)==0){print(paste0("no data on visit ",v));next}

    if(!is.null(field_sr_date)){diagdatefields = names(df_)[grepl(paste0("[^0-9]",field_sr_date,"[^0-9]",v),names(df_))]}
    # f.53.v.0
    visitdatefield = visitdatefields[v+1]
    for (i in 1:length(diagfields)) { # for each occurence of diagfield, find the corresponding age and convert it to  date - code and rbind() to df_out.
      #agefield = paste0("age_",v)
      diagfield = diagfields[i]

      if(!is.null(field_sr_date)){
        # should this line go under the if else block to avoid index out of bound error?
        diagdatefield = diagdatefields[i]
        if( length(diagdatefields) == 1) { diagdatefield = diagdatefields } # in case of death (40001/2)
      } else {
          diagdatefield = "dummy" # e.g. incase of medication
          }

      # print(paste0((diagfield), " - ", diagdatefield))
      #  empty diagnosis column
      if(all(is.na(df_[,diagfield,with=FALSE]))){next}
      # for rows with non-empty current diagfield, select identifier,diagfield,diagdatefield,visitdatefield
      df_sub <- df_[!is.na(get(diagfield) ),c("identifier",diagfield,diagdatefield,visitdatefield,"birthdt"),with=FALSE]
      df_sub$visit <- v
      names(df_sub) <- c("identifier","code","eventdate","visitdate","birthyearmonth","visit")
      df_out <- rbind(df_out,as.matrix(df_sub[,c("identifier","code","eventdate","visit","visitdate","birthyearmonth"),with=FALSE]))
    }

  }
  message("convert to dataframe")
  # df_out contains all visits , each row in df_out is a event
  df_out <- data.table::data.table(df_out,stringsAsFactors=F)
  df_out <- df_out[, visitdate:=as.Date(visitdate)]
  #df_out <- df_out[, code:=as.character(code)] #convert_col_to_integer(df_out$code) # df_out[, code:=as.integer(code)]
  # remove leading and trailing whitespace in column code
  df_out[, code:=lapply(.SD, trimws), .SDcols = "code"]
  # 99999 unclassifiable for ukb codings
  df_out <- df_out[!code %in% "99999"]

  # code type either numeric or character
  if (codetype == "numeric"){
    #  char -> integer
    df_out$code <- convert_col_to_integer(df_out$code)
  }



  # field_sr_date <-> field_sr_date_type, default :20008 <->	Interpolated Year
  if(!is.null(field_sr_date) & (field_sr_date_type=="interpolated_year"| field_sr_date_type=="interpolated_age" |  field_sr_date_type=="date" )) {

    if(field_sr_date_type=="interpolated_year") {
      df_out <- df_out[, eventdate:=as.numeric(eventdate)] ## as number. interpolated year
      # negative time = not meaningful. Refer coding 13
      df_out[eventdate <0,'eventdate']<-NA
      df_out$eventdate <- as.Date(convert_year_to_date(df_out$eventdate))

    } else if (field_sr_date_type=="interpolated_age"){
      df_out <- df_out[, eventdate:=as.numeric(eventdate)] ## as number. interpolated age.
      df_out[eventdate <0,'eventdate']<-NA
      # interpolate the event date as birth + age of diagnosis
      df_out[, birthyearmonth := as.Date(birthyearmonth)]
      df_out$eventdate = df_out[,"birthyearmonth"] + (df_out[,"eventdate"]*daysinyear)
    } else if (field_sr_date_type=="date"){
      df_out = df_out[, eventdate:=as.Date(eventdate)]
    }
    # remove rounding error from interpolation if self reported (can't self report after baseline.)
    if(event_code==2){
      df_out[df_out$eventdate > df_out$visitdate,'eventdate'] <- df_out[df_out$eventdate > df_out$visitdate,]$visitdate
    }
    # deduplicate, min/max/mean/sd <- not very efficient?!!
    message("deduplicate")
    # for each code in the same participant, compute min(oldest record)/max(newest record)/mean date
    #### slow: # dfout_extrastats <- df_out %>% group_by(identifier,code) %>% mutate(mindt = min(eventdate, na.rm = TRUE),maxdt = max(eventdate, na.rm = TRUE),meandt = mean(eventdate, na.rm = TRUE))
    data.table::setkey(df_out,identifier,code)
    dfout_extrastats <- suppressWarnings(df_out[, .(mindt= min(eventdate,na.rm = T),maxdt= max(eventdate,na.rm = T),meandt= mean(eventdate,na.rm=T) ), keyby=list(identifier,code)])
    dfout_extrastats <- merge(df_out[,c('identifier','code','eventdate')] ,dfout_extrastats,by=c('identifier','code'))

    # time between oldest and newest record in unit of year
    dfout_extrastats$diffdt <- (dfout_extrastats$maxdt - dfout_extrastats$mindt)/daysinyear
    # if this time difference is larger than qc threshold , mark NA in meandt
    dfout_extrastats[dfout_extrastats$diffdt>qc_threshold_year ,"meandt"] <- NA
    dfout_extrastats[dfout_extrastats$diffdt > 0,]
    #  take meandt as the event date , i.e. duplicate records with time difference > qc threshold will be changed to NA
    df_out$eventdate <- dfout_extrastats$meandt
    df_out <- df_out[order(df_out$visitdate),]
    df_out <- df_out[!duplicated(df_out[,c("identifier","code","eventdate"),with=FALSE]),] #sorted on visit, so first occurence is always first visit.

  }  else {
    df_out <- df_out[, eventdate:=as.Date(eventdate)]
  }

  # record which can be set as an event or not (when no event_date is reported, only visit)
  df_out$event <- event_code
  df_out <- df_out[, event:=as.integer(event)]
  # mark record without valid event date with 0
  df_out[is.na(df_out$eventdate)]$event <- 0
  # take visitdate as event date
  df_out[is.na(df_out$eventdate)]$eventdate <- df_out[is.na(df_out$eventdate)]$visitdate

  # add all visit dates as event=0 dates
  if(event_code==2){
    df_out_visit <- df_out
    df_out_visit$event<-0
    df_out_visit$eventdate <- df_out_visit$visitdate
    df_out<- unique(rbind(df_out,df_out_visit))
  }
  df_out <- df_out[,c("identifier","code","eventdate","event"),with=FALSE]

  data.table::setkey(df_out,code)
  df_out[, ('identifier') := lapply(.SD, as.character), .SDcols = 'identifier']

  gc()
  print(format(object.size(df_out), units = "Mb"))
  tictoc::toc()
  return(df_out)
}

