#' Get data for phenotype incidence and prevalence
#'
#' Given a data.table with all events for a phenotype and reference dates per individual, compute the time to event data for each individual. If no reference date is given then the date of first available event will be taken as reference date for each individual. ...
#' @param all_event_dt data table containing all events
#' @param df.data.settings data frame containing data settings
#' @param reference_date reference dates for each individuals in the whole cohort as a named vector
#' @param include_secondary_recurrence
#' @param window_ref_days_include number of days around the reference date (visit) that should be used to indicates if individuals had the event on the reference date. Relevant if you want to know if participant took medication on the visit
#' @param window_fu_days_mask number of days that future events should not be counted
#' @return  a data table with individual with time to event information
#' @keywords time-to-event
#' @export
#' @examples
#' all_event_dt <- get_all_events(dfDefinitions_processed_expanded[1,],lst.data,df.data.settings)
#' get_incidence_prevalence(all_event_dt,df.data.settings, reference_date = setNames(as.Date(as.character(dfukb$f.53.0.0),format="%Y-%m-%d"),dfukb$identifier))
### get prevalence (move this.. )

get_incidence_prevalence <- function(all_event_dt,
                                     df.data.settings,
                                     df_reference_date=NULL,
                                     include_secondary_recurrence=FALSE,
                                     #return_dates=FALSE, # TODO: not working yet.
                                     window_ref_days_include=0,##  indicate number of days around the reference date (visit) that should be used to indicates if individuals had the event on the reference date. For example, relevant if you want to know if participant took medication on the visit
                                     window_fu_days_mask=0, ## indicates number of days that future events should not be counted; e.g. you could only count events after 10 days from the reference visit to avoid events related to the reference date. e.g. you could also use it to only count events after X years, in order to avoid assessment bias.
                                     verbose=TRUE
                                     ) {
  # event==0, event cannot be used forr age - of - diagnosiis or new events.
  # event==1, event can be used for any type of future events, as it is based on ICD10 type of data
  # event==2, event can be used for any type of future events and age of diagnoses but only if there is no evidence in history, as it is self reported.

  # define window arround reference_dates? -- e.g. if you're diagnosed with XX, count follow-up events only after 20 days of no-events.
  # define window to score the diagnosis on refereence date.
  # reference_dates # should be a vector of dates, with the identifier as name.
  #reference_date = setNames(as.Date(as.character(dfukb$f.53.0.0),format="%Y-%m-%d"),dfukb$identifier)

  ### Future: 1) count any new event==1 code (icd10 etc)
  ###         2) count any primary event==1 if Hx==1 (recurrence) if 'include_secondary_recurrence'==TRUE
  ###         3) countany new event==2 (self report) age of diagnosis in future

  # If a reference date is given, it only returns individuals that match non-missing reference dates.
  #print(glue::glue("window_ref_days_include: {window_ref_days_include}"))
  #print(glue::glue("window_fu_days_mask: {window_fu_days_mask}"))

  # reference_date <- reference_date[!is.na(reference_date)]
  # change class of the column to date
  # scenario 1: df_reference_date is not NULL
  if(!is.null(df_reference_date)){
    # 1.1 not NULL but empty -> same treatment as of NULL
    if(nrow(df_reference_date)==0){
      # message("empyty reference_date table given, taking the first available event as reference")
      df_reference_date <- suppressWarnings(all_event_dt[,.(reference_date= min(eventdate,na.rm = T)),by=identifier])
    } else{
    # 1.2 not NULL and not empty
    # referencedate <-data.table::data.table(reference_date)
    # referencedate$identifier <- names(reference_date)
    names(df_reference_date)<-c("identifier","reference_date")
    df_reference_date$reference_date<- as.Date(df_reference_date$reference)
    # comment this line out to allow everyone to have no reference date
    # df_reference_date<-df_reference_date[!is.na(df_reference_date$reference_date)]
    if(verbose){
      message(glue::glue("Number of individuals in df_reference_date: {nrow(df_reference_date)}"))
    }
    }
    #scenario 2: df_reference_date is NULL
  } else{
    # message("no reference_date given, taking the first available event as reference")
    df_reference_date <- suppressWarnings(all_event_dt[,.(reference_date= min(eventdate,na.rm = T)),by=identifier])
  }


  if(include_secondary_recurrence){
    sources_recurrence_events <- df.data.settings$datasource
  } else {
    sources_recurrence_events <- df.data.settings %>%  dplyr::filter(diagnosis==1) %>% dplyr::pull(datasource)
  }

  ###########################################################
  # in case of empty rows , the entire column will be cast to default of class NA i.e. logical
  # this could create error for downstream functions "Error in bmerge....Incompatible join types"
  #  move identifier column to type numeric
  col_num<-c("identifier","count","sum.epidur","median.epidur","max.epidur", "survival_days", "Death_primary","Death_any" , "Hx_days", "Fu_days" ,"Hx" , "Fu","Ref" , "first_diagnosis_days","Any")
  # col_chr<-"identifier"
  col_date<-"reference_date"
  # catch if there is no event, return an empty table which can be merged
  if (nrow(all_event_dt)==0){
    message("No event found.")
    results_cols<-c("identifier","count","sum.epidur","median.epidur","max.epidur", "survival_days", "Death_primary","Death_any" , "Hx_days", "Fu_days" ,"Hx" , "Fu","Ref" , "first_diagnosis_days", "reference_date","Any")
    all_event_dt.summary <- setNames(data.table(matrix(nrow = 0, ncol = 16)),results_cols )
    all_event_dt.summary[, (col_num) := lapply(.SD, as.numeric), .SDcols = col_num]
    # all_event_dt.summary[, (col_chr) := lapply(.SD, as.character), .SDcols = col_chr]
    all_event_dt.summary[, (col_date) := lapply(.SD, as.Date), .SDcols = col_date]
    return(all_event_dt.summary)
  }


  df <- merge(all_event_dt,df_reference_date,by = 'identifier') %>% dplyr::arrange(eventdate) %>% data.table::as.data.table()



  # df <- df %>% filter(!is.na(reference_date)) # comment out if missing identifiers is fixed.
  df$days <- df$eventdate - df$reference_date

  data.table::setkey(df,days) # i don't know why, but setkey was already on identifier and cannot refresh..
  data.table::setkey(df,identifier)


  #################################################################################
  # death
  dfDth<-df[df$.id %in% df.data.settings[df.data.settings$death,]$datasource,]

  if (nrow(dfDth)>0){
  # get death records
  # dfDth<-df[(df$.id %in% df.data.settings[df.data.settings$death,]$datasource),]
  # ### flag primary death records
  # dfDth$death.primary<-ifelse((df.data.settings[match(dfDth$.id ,df.data.settings$datasource),]$diagnosis==1),2,NA)
  # ### flag secondary death records
  # dfDth$death.secondary<-ifelse((df.data.settings[match(dfDth$.id ,df.data.settings$datasource),]$diagnosis==2),2,NA)
  ### this seems faster
  dfDth$Death_primary <- NA
  dfDth$Death_primary[df.data.settings[match(dfDth$.id ,df.data.settings$datasource),]$diagnosis==1] <- 2

  dfDth$death_secondary<- NA
  dfDth$death_secondary[df.data.settings[match(dfDth$.id ,df.data.settings$datasource),]$diagnosis==2] <- 2
  dfDth$death_secondary<-data.table::fcoalesce(dfDth$Death_primary,dfDth$death_secondary)
  # dfDth
  dfDth<-dfDth[,c("identifier","days","Death_primary","death_secondary")]
  colnames(dfDth)<-c("identifier","survival_days","Death_primary","Death_any")
  # in the case of duplicate death records,one per primary /secondary
  dfDth_extrastats <- suppressWarnings(dfDth[, .(mindy= min(survival_days,na.rm = T),maxdy= max(survival_days,na.rm = T),meandy= mean(survival_days,na.rm=T) ), keyby=list(identifier)])
  id.diff.deathdt<-unique(dfDth_extrastats[dfDth_extrastats$mindy!=dfDth_extrastats$meandy|dfDth_extrastats$maxdy!=dfDth_extrastats$meandy,]$identifier)

  if (length(id.diff.deathdt)>0){
  message(glue::glue("Inconsistent death records for these individuals: {glue::glue_collapse(id.diff.deathdt,sep = ',')}"))
  message("Mean survival day will be taken.")
  }
  # in case of multiple death records with different death dates , take mean
  dfDth<-stats::aggregate(x=dfDth[,!(names(dfDth) %in% c("identifier")),with=FALSE], by=list(identifier=dfDth$identifier), mean, na.rm = TRUE)
  # Nan to NA
  dfDth[is.na(dfDth)]<-NA

  }else{
    # no records, take identifier for the merge later
    dfDth<-unique(df[,"identifier"])
    dfDth[,c("survival_days","Death_primary","Death_any")]<-as.numeric(NA)
  }
  ###############################################################################################################################

  ### History
  dfHx <- df[days<=0]
  Hx_days <- suppressWarnings(unique(dfHx[event>0][, .(Hx_days=min(days,na.rm=T)), keyby=list(identifier)] ))
  dfHx[,Hx:=2]

  ### Future
  # window_fu_days_mask
  dfFu <- df[  ((days>(0+window_fu_days_mask) & event==1) & (!identifier %in% dfHx$identifier)) |
               (days>(0+window_fu_days_mask) & event==1 & (identifier %in% dfHx$identifier) & .id %in% sources_recurrence_events ) |
               (days>(0+window_fu_days_mask) & event==2  & (!identifier %in% dfHx$identifier)) ]
  dfFu[,Fu:=2] #unique(dfFu$identifier)
  Fu_days <- suppressWarnings( unique(dfFu[,.(Fu_days= min(days,na.rm=T)), keyby=list(identifier)] ))
  ### age of diagnosis
  #system.time({ df %>% filter(event>0) %>% group_by(identifier) %>% summarise(first_diagnosis_days=min(days)) }) # <- slow..
  #system.time({ df[df$event>0][,.(first_diagnosis_days=min(days,na.rm=T) ), by=identifier] }) # <- fast..
  #df[df$event>0][df[, .I[which.max(days)], by=identifier]$V1] # <- aanother way..
  #
  # first_diagnosis_days <- df[  ,.(first_diagnosis_days=min(days,na.rm=T)), by=identifier]
  #
  #
  # first_diagnosis_days <- df[  ,.SD[which.min(days)], by=identifier]

  ### Data if participant had event/med  reference date (visit);+/- x day
  #TODO I think episodes for visitdate (event=0 rows) is messing up with this flag! i.e. almost all of them have this flag on
  # MWY: NOW added df$event!=0 ,though there are still a lot Ref==2 using window of 0 day!
  dfRef <- df[df$eventdate>=(df$reference_date-window_ref_days_include) & df$eventdate<=(df$reference_date+window_ref_days_include) &df$event!=0,c("identifier")]

  dfRef[,Ref:=2]


  ## some other stats, and to include event==0 individuals:
  # allow abesence of HES data , sr/death does not have epidur
  if(!all(unique(all_event_dt$event==0)) & !all(is.na(all_event_dt$epidur))){
    stats <- suppressWarnings( all_event_dt[, .(count = .N,sum.epidur= sum(epidur,na.rm = T),median.epidur= median(epidur,na.rm = T),max.epidur= max(epidur,na.rm=T)), by = identifier] )
    stats[is.infinite(stats$max.epidur),]$max.epidur <-NA
  } else{
    stats <- suppressWarnings(all_event_dt[, .(count = .N,sum.epidur= NA,median.epidur= NA,max.epidur=NA), by = identifier])
  }



  # test <- Reduce(function(...) merge(..., all = TRUE,by='identifier'), list(df,
  #                                                                      Hx_days,
  #                                                                      Fu_days,
  #                                                                      unique(dfHx[,c("identifier","Hx")]),
  #                                                                      unique(dfFu[,c("identifier","Fu")]),
  #                                                                      unique(dfRef[,c("identifier","Ref")]),
  #                                                                      first_diagnosis_days,
  #                                                                      df.referencedate
  #                                                                      stats)
  # )
  # View(test)

  # stats
  # Hx_days
  # Fu_days
  # dfDth
  # dfHx[,c("identifier","Hx")]
  # dfHx
  # dfRef[,c("identifier","Ref")]

  all_event_dt.summary <- Reduce(function(...) merge(..., all = TRUE,by='identifier'), list(
      stats,
      unique(dfDth[ , c("identifier","survival_days","Death_primary","Death_any")]),
      Hx_days,
      Fu_days,
      unique(dfHx[,c("identifier","Hx")]),
      unique(dfFu[,c("identifier","Fu")]),
      unique(dfRef[,c("identifier","Ref")])
  ))


  all_event_dt.summary$first_diagnosis_days <- pmin(all_event_dt.summary$Hx_days,all_event_dt.summary$Fu_days,na.rm = T)
  all_event_dt.summary[ Hx==2 & is.na(Hx_days),'first_diagnosis_days'] <- NA

  all_event_dt.summary <- merge(all_event_dt.summary,df_reference_date,by="identifier")
  all_event_dt.summary[,Any:=2]
  all_event_dt.summary <- data.table::data.table(all_event_dt.summary)

  ####################################type check######################################################
  all_event_dt.summary[, (col_num) := lapply(.SD, as.numeric), .SDcols = col_num]
  # all_event_dt.summary[, (col_chr) := lapply(.SD, as.character), .SDcols = col_chr]
  all_event_dt.summary[, (col_date) := lapply(.SD, as.Date), .SDcols = col_date]
  #####################################################################################################
  return(all_event_dt.summary)


}


#' Get case for a phenotype
#'
#' Given a phenotype, a list of episode data and reference dates per individual, identify cases /exclude cases fulfilling exclusion criteria and compute the time to event data for these individuals.  If no reference date is given then the date of first available event will be taken as reference date for each individual.
#' @param definitions phenotype/trait specified in definition table (a row in the table)
#' @param lst.data list of data table with all episode data
#' @param df.data.settings data frame containing data settings
#' @param reference_date reference dates for each individuals in the whole cohort as a named vector
#' @param window_fu_days_mask number of days that future events should not be counted. e.g. you could only count events after 10 days from the reference visit to avoid events related to the reference date.

#' @return  a list of 2 data tables : all events for valid cases and an event summary containing time to event information for these individuals.
#' @keywords time-to-event
#' @export
#' @examples
#' get_cases(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="Nicm"), lst.data,df.data.settings, reference_date=setNames(as.Date(as.character(dfukb$f.53.0.0),format="%Y-%m-%d"),dfukb$identifier))
get_cases <- function(definitions,
                       lst.data,
                       df.data.settings,
                      df_reference_dt=NULL,
                      include_secondary_recurrence=FALSE,
                      verbose=TRUE,
                      window_fu_days_mask=0,
                       ...
                         ) {

  # define cases
  if(length(unique(definitions$TRAIT))>1){
    message("more than 1 TRAIT in definitions")
  }
  if(length(definitions$TRAIT)==0){
    message("No TRAIT in definitions.Stop.")
    return(0)
  }
  if(verbose){
  message("..Identify case status")
  }
  all_event_dt.Include_in_cases <- get_all_events(definitions %>% dplyr::filter(Definitions =="Include_in_cases"),lst.data,df.data.settings,verbose = verbose)
  # if(missing(window_fu_days_mask)) {
  #   all_event_dt.Include_in_cases.summary <- get_incidence_prevalence(all_event_dt = all_event_dt.Include_in_cases,df.data.settings,
  #                                                                     df_reference_date = df_reference_dt,verbose = verbose)
  #                                                                     #...)
  # } else {
    all_event_dt.Include_in_cases.summary <- get_incidence_prevalence(all_event_dt = all_event_dt.Include_in_cases,df.data.settings,
                                                                      df_reference_date = df_reference_dt,include_secondary_recurrence=include_secondary_recurrence,verbose = verbose,window_fu_days_mask=window_fu_days_mask)
    # }


  if(verbose){
  message(glue::glue("Including {nrow(all_event_dt.Include_in_cases.summary)} cases with reference dates..."))
  }
  all_event_dt.Exclude_from_cases <- get_all_events(definitions %>% dplyr::filter(Definitions =="Exclude_from_cases"),lst.data,df.data.settings)
  if(!is.null(all_event_dt.Exclude_from_cases)){
    exclude=all_event_dt.Include_in_cases.summary$identifier %in% unique(all_event_dt.Exclude_from_cases$identifier)
    if(verbose){
    message(glue::glue("Excluding {sum(exclude)} cases..."))
    }
    set_to_na <- names(all_event_dt.Include_in_cases.summary)[!names(all_event_dt.Include_in_cases.summary) %in% c("identifier","reference_date")]
    # the .summary needed to be flagged for the get_case_control() functoin
    all_event_dt.Include_in_cases.summary[exclude,(set_to_na):=-2]
    # MWY:added first_diagnosis_days as well because negative value is possible in this col
    all_event_dt.Include_in_cases.summary[exclude,("first_diagnosis_days"):=NA]
    # remove these in all_event_dt
    all_event_dt.Include_in_cases<-all_event_dt.Include_in_cases[! (all_event_dt.Include_in_cases$identifier %in% unique(all_event_dt.Exclude_from_cases$identifier)),]
    if(verbose){
    message(glue::glue("{nrow(all_event_dt.Include_in_cases)} events fulfiling criteria in include_in_cases"))
    }
  }
  return(list(all_event_dt.Include_in_cases=all_event_dt.Include_in_cases,
              all_event_dt.Include_in_cases.summary=all_event_dt.Include_in_cases.summary)
  )
}




#' Get case and controls for a phenotype
#'
#' Given a phenotype, a list of episode data and reference dates per individual, identify valid cases/controls and compute the time to event data for these individuals respectively.  If no reference date is given then the date of first available event will be taken as reference date for each individual.
#' @param definitions phenotype/trait specified in definition table (a row in the table)
#' @param lst.data list of data table with all episode data
#' @param df.data.settings data frame containing data settings
#' @param df_reference_date dataframe where first column is the identifier and second column the reference dates
#' @param vct.identifiers character vector listing the identifiers in the cohort. This is used to define controls if reference_date is not given
#' @param window_fu_days_mask number of days that future events should not be counted. e.g. you could only count events after 10 days from the reference visit to avoid events related to the reference date.
#' @return  a list of 2 data tables : all events for valid cases and an event summary containing time to event information for these individuals.
#' @keywords time-to-event
#' @export
#' @examples
#' get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="Nicm"), lst.data,df.data.settings,  reference_date=setNames(as.Date(as.character(dfukb$f.53.0.0),format="%Y-%m-%d"),dfukb$identifier))
get_cases_controls <-function(definitions,
                                lst.data,
                                df.data.settings,
                                df_reference_date=NULL,
                              include_secondary_recurrence=FALSE,
                                vct.identifiers=NULL,
                                verbose=TRUE,
                              window_fu_days_mask=0
                                ) {
  # vct.identifiers <-Used to define controls if reference_date is not given (NULL)
  if(nrow(definitions)==0){
    message("No definition is provided.Stop.")
    return(0)
  }
  if (length(definitions$TRAIT)==0){
    message("No TRAIT in definitions.Stop.")
    return(0)
  }

  #  scenario 1 : if df_reference_date is provided
  if (!is.null(df_reference_date)){
  # reference_date <- reference_date[!is.na(reference_date)]
    # reference_date <- reference_date[!is.na(names(reference_date))]
    names(df_reference_date)<-c("identifier","reference_date")
    df_reference_date <- df_reference_date[!is.na(df_reference_date$identifier),]
    # df_reference_date <- df_reference_date[!is.na(df_reference_date$reference_date),]
    if (nrow(df_reference_date)==0){
      message("Please verify there is individual in df_reference_date")
    }
    # if there are multiple reference dates for one individual, the calculation will mess up
    df_reference_date<-unique(df_reference_date,by="identifier")

  }else {
    # scenario 2: if df_reference_date is empty
    if(verbose){
    message("df_reference_date=NULL, read vct.identifiers")
    }
    if(length(vct.identifiers)<=1){
      # scenario 2.1 if no vct.identifiers
      if(verbose){
      message("...vct.identifiers is empty, create the vct.identifiers from lst.data")
      }
      # subset the identifier columns from list of data tables ->concatenate -> keep unique identifiers
      df_reference_date<-unique(data.table::rbindlist(lapply(lst.data,function(x) subset(x,select="identifier"))))
      # df_reference_date[, reference_date := as.Date(NA)]
    # scenario 2.2 if there is vct,identifiers
    }else{
      # reference_date = setNames(as.Date(rep(NA,length(vct.identifiers))),vct.identifiers)
      df_reference_date <- as.data.table(vct.identifiers)
      names(df_reference_date)<-c("identifier")
      # df_reference_date[, reference_date := as.Date(NA)]
      # names(df_reference_date)<-c("identifier","reference_date")
    }
    if(verbose){
    message(glue::glue("...No reference date information provided, take first event date as reference date"))
    }
    all_event_dt.first_occur <- get_all_events(definitions %>% dplyr::filter(Definitions =="Include_in_cases"),lst.data,df.data.settings,verbose = FALSE)
    # only keep date with real event
    all_event_dt.first_occur[all_event_dt.first_occur$event==0,]$eventdate <-NA
    # get everyone , take the date of relevant events respectively
    all_event_dt.first_occur.summary <- get_incidence_prevalence(all_event_dt = all_event_dt.first_occur,df.data.settings,df_reference_date = NULL,include_secondary_recurrence=include_secondary_recurrence,window_fu_days_mask = window_fu_days_mask,verbose = verbose)
    #  merge, NA for those without an event
    df_reference_date<-merge(df_reference_date,all_event_dt.first_occur.summary[, c("identifier","reference_date")],by="identifier",all.x=TRUE)
    #  should be already be in Date format but just to be sure
    df_reference_date$reference_date<-as.Date(as.character(df_reference_date$reference_date),format="%Y-%m-%d")
    # message(glue::glue("Population: {length(unique(df_reference_date$identifier))} individuals "))
   }
  # until here a df_reference_date with identifier and reference dates is prepared
  # now to the definitions


  # define population
  all_event_dt.population <- get_all_events(definitions %>% dplyr::filter(Definitions =="Study_population"),lst.data,df.data.settings,verbose = FALSE)
    # scenario 3: if there is study population stated in the definition, df_reference_date is replaced
    if(! is.null(all_event_dt.population)) {
    # Study population field of the definition is not empty , derive the df_reference_date from the study population trait
      if(verbose){
      message(glue::glue("...Study population column not empty, take first event date of this trait as reference date instead"))
      }
    # only consider event with real event date in study population, set those with visitdate to NA
    all_event_dt.population[all_event_dt.population$event==0,]$eventdate <-NA
    # get everyone , take the date of relevant events respectively
    all_event_dt.population.summary <- get_incidence_prevalence(all_event_dt = all_event_dt.population,df.data.settings,df_reference_date = NULL,include_secondary_recurrence=include_secondary_recurrence,window_fu_days_mask = window_fu_days_mask,verbose = verbose)
    # update the dates in the original df_reference_date
    # reference_date = setNames(as.Date(as.character(all_event_dt.population.summary$reference_date),format="%Y-%m-%d"),all_event_dt.population.summary$identifier)
    # df_reference_date <- all_event_dt.population.summary[, c("identifier","reference_date")]
    ###############################################
    ###############################################
    df_reference_date<-merge(df_reference_date[,c("identifier")],all_event_dt.population.summary[, c("identifier","reference_date")],by="identifier",all.y=TRUE)
    df_reference_date$reference_date<-as.Date(as.character(df_reference_date$reference_date),format="%Y-%m-%d")

    if(verbose){
      message(glue::glue("Population: {length(unique(df_reference_date$identifier))} individuals "))
      }

   } else {
    # Study population field of the definition is empty, keep reference date as it is
     if(verbose){
       message(glue::glue("Population: total non-missing reference date = {nrow(df_reference_date[!is.na(df_reference_date$reference_date),])}, total missing reference date= {nrow(df_reference_date[is.na(df_reference_date$reference_date),])} "))
     }
     }

  # define cases
  cases <- get_cases(definitions,lst.data,df.data.settings,df_reference_date,include_secondary_recurrence=include_secondary_recurrence,window_ref_days_include=0,window_fu_days_mask=window_fu_days_mask,verbose=verbose)

  all_event_dt.Include_in_cases.summary <- cases$all_event_dt.Include_in_cases.summary
  all_event_dt.Include_in_cases <- cases$all_event_dt.Include_in_cases
  # define exclude controls
  all_event_dt.Exclude_from_controls <- get_all_events(definitions %>% dplyr::filter(Definitions =="Exclude_from_controls"),lst.data,df.data.settings,verbose = FALSE)

  ### define case & control
  df.casecontrol <- df_reference_date


  # exclude id in case_include & case_exclude from summary => potential control
  df.casecontrol <- df.casecontrol[!df.casecontrol$identifier %in% all_event_dt.Include_in_cases.summary$identifier,]
  df.casecontrol$reference_date <- as.Date(as.character(df.casecontrol$reference_date),format="%Y-%m-%d")


  # case without eventdate
  # this line is not correct because it will also get ppl who do not pass the instance filter! i.e. they will have event but not appearing in case.summary!
  # potential_control<-df.casecontrol$identifier
  # non_proper_case<-dplyr::intersect(potential_control,all_event_dt.Include_in_cases$identifier)
  # hence replace with these 2 lines
  w_dt<-dplyr::union(all_event_dt.Include_in_cases[all_event_dt.Include_in_cases$event==1,]$identifier,all_event_dt.Include_in_cases[all_event_dt.Include_in_cases$event==2,]$identifier)
  non_proper_case<-dplyr::setdiff(all_event_dt.Include_in_cases[all_event_dt.Include_in_cases$event==0,]$identifier,w_dt)
  # restrict to only those included in df_reference_date,
  # because get_all_events() get all events in the the whole data
  non_proper_case<-dplyr::intersect(non_proper_case,df_reference_date$identifier)
  # merge it with the cases
  df.casecontrol <- merge(df.casecontrol,all_event_dt.Include_in_cases.summary,by=c("identifier","reference_date"),all=T)


  if(length(non_proper_case)>0) {
    if(verbose){
      non_proper_case_norefdt<- unique(df_reference_date[identifier %in% non_proper_case][is.na(reference_date)]$identifier)
    message(glue::glue("!!NOTE: {length(non_proper_case)} cases without valid event dates, of which {length(non_proper_case_norefdt)} without reference dates..."))
    }
    # Hx Fu not known only Any ==2
    df.casecontrol[df.casecontrol$identifier %in% non_proper_case,(c("Any")):=2]
  }


  if(!is.null(all_event_dt.Exclude_from_controls)) {
    # !is.na(df.casecontrol$Any) for cases i.e. from all_event_dt.Include_in_cases.summary
    # here to exclude = 1)not a case  2) in exclude_from_control
    exclude=(is.na(df.casecontrol$Any) & (df.casecontrol$identifier %in% all_event_dt.Exclude_from_controls$identifier))
    if(verbose){
    message(glue::glue("excluding {sum(exclude)} controls"))
    }
    # cols to be set to na
    set_to_na <- names(df.casecontrol)[!names(df.casecontrol) %in% c("identifier","reference_date")]
    df.casecontrol[exclude,(set_to_na):=-1] #mark the non-control as -1
    # MWY:added first_diagnosis_days as well because negative value is possible in this col
    df.casecontrol[exclude,("first_diagnosis_days"):=NA]
  }

  # 2 => case -2 => non-case  -1 => non-control , set control to 1
  df.casecontrol[is.na(df.casecontrol$Any),]$Any <- 1
  df.casecontrol[is.na(df.casecontrol$Death_primary),]$Death_primary <- 1
  df.casecontrol[is.na(df.casecontrol$Death_any),]$Death_any <- 1

  # as control has no event the event count is NA after merge, set it to 1
  df.casecontrol[is.na(df.casecontrol$count),]$count <- 1
  df.casecontrol[(is.na(df.casecontrol$Hx)&(df.casecontrol$Any == 1)),]$Hx <- 1
  df.casecontrol[(is.na(df.casecontrol$Fu)&(df.casecontrol$Any == 1)),]$Fu <- 1
  df.casecontrol[(is.na(df.casecontrol$Ref)&(df.casecontrol$Any == 1)),]$Ref <- 1
  # At this point NA is only possible in Hx/Fu/Ref for case_included i.e. Any==2
  # first consider Hx, cases either have event in past ==2 / future ==NA / unknown==NA [no ref date]
  # case that don't have event in the past but in future
  df.casecontrol[(is.na(df.casecontrol$Hx)&(df.casecontrol$Any == 2)&is.na(df.casecontrol$Hx_days)& (!is.na(df.casecontrol$Fu_days))),]$Hx <-  1
  # similar for Ref , cases either have event within window ==2 / no ==NA / unknown==NA [no ref date]
  # case that don't have event within window i.e. eventdate is known
  df.casecontrol[(is.na(df.casecontrol$Ref)&(df.casecontrol$Any == 2)&((!is.na(df.casecontrol$Hx_days)) | (!is.na(df.casecontrol$Fu_days)))),]$Ref <-  1
  # for Fu , cases either have event in past ==NA / future ==2/ unknown==NA [no ref date]
  # ******************but this is also affected by the sources_recurrence_events option in get_incidence_rate()?***************
  df.casecontrol[(is.na(df.casecontrol$Fu)&(df.casecontrol$Any == 2)& (!is.na(df.casecontrol$Hx_days))& (is.na(df.casecontrol$Fu_days))),]$Fu <-  1
  ############################################################################################################
  if(verbose){
      print(table(df.casecontrol$Any))
  }


  return(list(df.casecontrol=df.casecontrol,
              all_event_dt.Include_in_cases=all_event_dt.Include_in_cases[all_event_dt.Include_in_cases$identifier %in% df.casecontrol$identifier,],
              all_event_dt.Include_in_cases.summary=all_event_dt.Include_in_cases.summary)
  )
}



