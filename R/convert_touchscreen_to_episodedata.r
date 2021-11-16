# library(tic) # library(tictoc)?

# # !!!!!!!backward imputation is problematic!!!!
# harmonize_agediag_bycols<-function (df,tsdiagnosisdatefields,qc_treshold_year=10){
#   # Age of diagnosis should not change across visits
#   # take mean and remove records with discrepancy exceeding threshold
#
#   df_extrastats<-df[,tsdiagnosisdatefields,with=FALSE]
#
#   # set negative values -1 (not know) /-3 (prefer not answer) to NA
#   for(j in tsdiagnosisdatefields){
#     set(df_extrastats, i= which(df_extrastats[[j]]<0), j= j, value=NA)
#   }
#   df_extrastats<-data.table(df_extrastats,stringsAsFactors=F)
#   # for each code in the same participant, compute min(oldest record)/max(newest record)/mean age
#   # https://stackoverflow.com/questions/31258547/data-table-row-wise-sum-mean-min-max-like-dplyr
#
#   df_extrastats[, `:=`(agemin = rowMins(as.matrix(.SD), na.rm=T),
#                        agemax = rowMaxs(as.matrix(.SD), na.rm=T),
#                        agemean = rowMeans(.SD, na.rm=T)), .SDcols=tsdiagnosisdatefields]
#   # time between oldest and newest record in unit of year
#   df_extrastats$agediff <- (df_extrastats$agemax - df_extrastats$agemin)
#   # if this time difference is larger than qc threshold , mark NA in meandt
#   df_extrastats[df_extrastats$agediff>qc_treshold_year ,"agemean"] <- NA
#   df_extrastats[is.nan(df_extrastats$agemean)]$agemean <- NA
#   unique(df_extrastats$agemean)
#   # show records with discrepancies
#   # df_extrastats[df_extrastats$agediff > 0,]
#   for (col in tsdiagnosisdatefields){
#     set(df, j = col, value = df_extrastats[["agemean"]])
#   }
#   return (df)
#
# }





#' Convert self-reported illness codes from touchscreen with optionally corresponding time of diagnosis to episodes format
#'
#' This function takes data fields containing illness codes/time of diagnosis distributed in multiple with each row representing one individual. It processed the data and return in episodes of event for all individuals
#' @param df dataframe containing the fields
#' @param ts_conditions touchscreen data fields for illness codes, each followed by the data field for corresponding time of diagnosis in square brackets
#' @param qc_threshold_year  in case of multiple episodes for same illness code, if time difference is larger than qc threshold between the oldest and newest episodes, eventdate will be discarded i.e. set to NA, default:10
#' @param event_code code number to indicate if the episodes is a true event 1, self-reported/interpolated event 2 , for which a mean date is taken as event date in case of multiple occurence, or non-event 0 , for which date does not correspond to occurence of illness code. These codes are used in time-to-event analysis, default: 2
#' @return  a data.table object with all episodes
#' @keywords self-reported data
#' @export
#' @examples
#' convert_touchscreen_to_episodedata(dfukb,ts_conditions = "6150[3894]")
#' convert_touchscreen_to_episodedata(dfukb,ts_conditions = dfDefinitions_processed$TS)
convert_touchscreen_to_episodedata<- function(df,ts_conditions=dfDefinitions_processed$TS,qc_treshold_year=10,event_code=2){
  tictoc::tic()
  ts_conditions <- unique(c(na.omit(unlist(strsplit(ts_conditions,",")))))
  print(paste("Input fields:",ts_conditions,sep=" ") )

  daysinyear=365.25
  # visit date code
  field_visit_date="53"
  # vector with name of the identifier col
  identifierfield = "identifier"
  #  vector with names of all visit cols : "f.53.0.0" "f.53.1.0" "f.53.2.0" "f.53.3.0"
  visitdatefields = names(df)[grepl(paste0("[^0-9]",field_visit_date,"[^0-9]"), names(df))]
  visits = length(visitdatefields) #sum(grepl("53_", names(df)))
  # need for calculating diagdate from age of diagnosis
  field_birth_year ="34"
  field_birth_month="52"
  birthyearfield = names(df)[grepl(paste0("[^0-9]",field_birth_year,"[^0-9]"), names(df))]
  birthmonthfield = names(df)[grepl(paste0("[^0-9]",field_birth_month,"[^0-9]"), names(df))]


  df_out <-  matrix(ncol=6, nrow=0) # initiate output

  # for each field listed in ts
  for (col in ts_conditions) {
    # col<-"3581≥0[3581]"
    message(glue::glue("process touchscreen data for {col}"))
    # parse the field and condition
    cdn<-stringr::str_extract(col,"[=|<|>|≥|≤|!][=]*\\d+")
    field_ts_diagnosis<-stringr::str_extract(col,"\\d+")
    tsdiagnosisfields = names(df)[grepl( paste0("[^0-9]",field_ts_diagnosis,"[^0-9]"), names(df))]


    # optional col for age of diagnosis specified in bracket []
    # https://stackoverflow.com/questions/52061753/r-capturing-string-inside-brackets
    # regmatches returns character(0) if there is no match
    age_diagnosis_col<-regmatches(col, regexpr("\\[\\K[^][]*", col,perl = TRUE))
    if(length(age_diagnosis_col)>0){tsdiagnosisdatefields = names(df)[grepl(age_diagnosis_col, names(df))]} else tsdiagnosisdatefields=NULL


  # only need n_eid, visit dates and diag-codes + age-of-diag
  columns_to_keep = c(identifierfield,
                      visitdatefields,
                      tsdiagnosisfields,
                      tsdiagnosisdatefields,
                      birthyearfield,
                      birthmonthfield
                     )
  # data.table - Setting with = FALSE disables the ability to refer to columns as if they are variables
  df_ <- df[,columns_to_keep,with=FALSE]
  df_$dummy <- NA
  # birthdate ,day set to first of the month
  df_$birthdt = as.Date(as.character(paste(df_[[birthyearfield]],df_[[birthmonthfield]], 1, sep = "-")),format="%Y-%m-%d")
  # COULD CREATE EVENTDATE > VISITDATE ! harmonize the age of diagnosis fields
  # if(!is.null(tsdiagnosisdatefields)){ df_<-harmonize_agediag_bycols (df_,tsdiagnosisdatefields)}


  for (v in 0:(visits-1)){ # for each visit,
    # v=0
    message(paste0("querying visit ",v))
    # f.xxxxx.v.0-9
    diagfields =unique(names(df_)[grepl(paste0("[^0-9]",field_ts_diagnosis,"[^0-9]",v),names(df_))])


    if(length(diagfields)==0){print(paste0("no data on visit ",v));next}
    if(!is.null(tsdiagnosisdatefields)){diagdatefields = unique(names(df_)[grepl(paste0("[^0-9]",age_diagnosis_col,"[^0-9]",v),names(df_))])}
    # f.53.v.0
    visitdatefield = visitdatefields[v+1]
    # for each occurence of diagfield, find the corresponding age and convert it to  date - code and rbind() to df_out.
    for (i in 1:length(diagfields)) {
      diagfield = diagfields[i]

      if(length(age_diagnosis_col)>0 ){
        diagdatefield <- diagdatefields
        } else {
        diagdatefield = "dummy"} # no known date of diag, fill with visit date below
      #  empty diagnosis column
      if(all(is.na(df_[,diagfield,with=FALSE]))){next}
      # diagfield example f.xxxxx.v.i
      # for rows with non-empty current diagfield, select identifier,diagfield,diagdatefield,visitdatefield
      df_sub <- df_[!is.na(get(diagfield) ),c(identifierfield,diagfield,diagdatefield,visitdatefield,"birthdt"),with=FALSE]
      # in case diagfield == diagdatefield
      colnames(df_sub)[3] <- paste(diagdatefield,"_",sep="")

      # find rows that fulfil the condition
      cdn_exp <-paste(diagfield,cdn,sep="") #"f.xxxxx.v.i ==1"

      ########################################################################
      #  if no expression is given , get all non-NA rows
      # this allows handling of field like 5441, but is it a good idea to allow this? e.g. minus values in certain fields => missing
      if (is.na(cdn)){

        df_sub<- df_sub %>% dplyr::filter(!is.na(diagfield))
      }else{
        df_sub<- df_sub %>% dplyr::filter(eval((parse(text=cdn_exp))))
      }
      ########################################################################
      # if no rows fulfil the condition
      if (nrow(df_sub)==0){next}
      # replace the diagfield content with the condition
      df_sub[[diagfield]]<-col #paste(field_ts_diagnosis,cdn,sep="")
      # add visit instance
      df_sub$visit <- v
      names(df_sub) <- c("identifier","code","eventdate","visitdate","birthyearmonth","visit")
      df_out <- rbind(df_out,as.matrix(df_sub[,c("identifier","code","eventdate","visit","visitdate","birthyearmonth"),with=FALSE]))
    }

  }
  }

  # after loop through all fields listed in ts
  message("convert to dataframe")
  # df_out contains all visits , each row in df_out is a event
  df_out <- data.table::data.table(df_out,stringsAsFactors=F)
  df_out <- df_out[, visitdate:=as.Date(visitdate)]
  df_out <- df_out[, birthyearmonth:=as.Date(birthyearmonth)]


  # compute the event date from visitdate and age of diagnosis
  df_out <- df_out[, eventdate:=as.numeric(eventdate)]
  # negative age not meaningful
  df_out[eventdate <0,'eventdate']<-NA
  # interpolate the event date as birth + age of diagnosis
  df_out$eventdate = df_out[,"birthyearmonth"] + (df_out[,"eventdate"]*daysinyear)


  # remove rounding error from interpolation
  df_out[df_out$eventdate > df_out$visitdate,'eventdate'] <- df_out[df_out$eventdate > df_out$visitdate,]$visitdate

  # deduplicate
  message("Deduplicate")
  data.table::setkey(df_out,identifier,code)
  dfout_extrastats <- suppressWarnings(df_out[, .(mindt= min(eventdate,na.rm = T),maxdt= max(eventdate,na.rm = T),meandt= mean(eventdate,na.rm=T) ), keyby=list(identifier,code)])
  dfout_extrastats <- merge(df_out[,c('identifier','code','eventdate')] ,dfout_extrastats,by=c('identifier','code'))

  # time between oldest and newest record in unit of year
  dfout_extrastats$diffdt <- (dfout_extrastats$maxdt - dfout_extrastats$mindt)/daysinyear
  # if this time difference is larger than qc threshold , mark NA in meandt
  dfout_extrastats[dfout_extrastats$diffdt>qc_treshold_year ,"meandt"] <- NA
  dfout_extrastats[dfout_extrastats$diffdt > 0,]
  #  take meandt as the event date , i.e. duplicate records with time difference > qc threshold will be changed to NA
  df_out$eventdate <- dfout_extrastats$meandt
  df_out <- df_out[order(df_out$visitdate),]
  df_out <- df_out[!duplicated(df_out[,c("identifier","code","eventdate"),with=FALSE]),] #sorted on visit, so first occurence is always first visit.


  # record which can be set as an event or not (when no event_date is reported, only visit)
  df_out$event <- event_code
  df_out <- df_out[, event:=as.integer(event)]
  # mark record without valid event date with 0
  df_out[is.na(df_out$eventdate)]$event <- 0
  # take visitdate as event date
  df_out[is.na(df_out$eventdate)]$eventdate <- df_out[is.na(df_out$eventdate)]$visitdate

  # add all visit dates as event=0 dates
  df_out_visit <- df_out
  df_out_visit$event<-0
  df_out_visit$eventdate <- df_out_visit$visitdate
  df_out<- unique(rbind(df_out,df_out_visit))

  df_out <- df_out[,c("identifier","code","eventdate","event"),with=FALSE]
  setkey(df_out,code)
  df_out[, ('identifier') := lapply(.SD, as.character), .SDcols = 'identifier']

  gc()
  print(format(object.size(df_out), units = "Mb"))

  tictoc::toc()
  return(df_out)


}




