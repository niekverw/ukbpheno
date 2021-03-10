
#' Get case-control status by data sources
#'
#' Given a phenotype, a list of episode data and reference dates per individual, identify cases and control status as of *reference date" by source.
#' @param source_lst list of data sources (as vector) each containing the repsective individual sources with names matched to datasource from lst.data.setting
#' @param definitions phenotype/trait specified in definition table (a row in the table)
#' @param lst.data list of data table with all episode data
#' @param lst.data.settings data frame containing data settings
#' @param reference_date reference dates for each individuals in the whole cohort as a named vector

#' @return  a data table : case control status according to different sources and one column based on any of the sources.
#' @keywords time-to-event
#' @export
#' @examples
#' cancer_source<-c("tte.cancer.icd10" ,"tte.cancer.icd9" )
#' get_case_status_by_source(cancer_source,definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="Nicm"), lst.data,lst.data.settings,lst.identifiers=dfukb$f.eid))
get_case_status_by_source <- function(source_lst,definitions,
                      lst.data,
                      lst.data.settings,
                      lst.identifiers,
                      reference_date
                      ) {
  lst.case_control <- get_cases_controls(definitions, lst.data,lst.data.settings,reference_date,lst.identifiers)
  ###########################################################################
  # NOTE: this function depends on the get_cases_controls! modify accordingly
  ##########################################################################
  # extract only relevant columns from df_casecontrol
  case_status<-lst.case_control$df.casecontrol[,c("f.eid","reference_date","Death_any","survival_days","Any","Hx")]

  # ##########################################################################
  # exclude future events
  # otherwise the comparison between sources is not fair
  # ##########################################################################
  event_dt.Include_in_cases_beforeOrOn_refdate<-merge(lst.case_control$all_event_dt.Include_in_cases,case_status[,c("f.eid","reference_date")],by="f.eid")
  # discard future events
  message("Remove future events at reference date")
  event_dt.Include_in_cases_beforeOrOn_refdate<-event_dt.Include_in_cases_beforeOrOn_refdate[event_dt.Include_in_cases_beforeOrOn_refdate$eventdate<=event_dt.Include_in_cases_beforeOrOn_refdate$reference_date]


  # create new columns for each source
  for (source_name in names(source_lst)) {
    varname<-paste('Hx',source_name,sep='.')
    # create new columns for each source
    case_status[, (varname)]<-as.numeric(NA)
    # discard non-event in which the event date is explicitly set to visit date? BUT some events simply have not event dates.
    # lookup source from df with all episodes i.e. all_event_dt.Include_in_cases
    # if rows in include_in_case which originated from the target source, look up the eid and set to 2
    # case_status[[varname]][case_status$f.eid %in% lst.case_control$all_event_dt.Include_in_cases[lst.case_control$all_event_dt.Include_in_cases$.id %in% source_lst[[source_name]],]$f.eid] <- 2
    case_status[[varname]][case_status$f.eid %in% event_dt.Include_in_cases_beforeOrOn_refdate[event_dt.Include_in_cases_beforeOrOn_refdate$.id %in% source_lst[[source_name]],]$f.eid] <-2
  }
  for (j in names(case_status)[grepl("Hx",names(case_status))]){
    # for all Hx.x columns
    # if Hx <0 , copy Hx (excluded)
    set(case_status,which(case_status$Hx<=0),j,case_status[case_status$Hx<=0]$Hx)
    #  if the cell value was NA -> set to 1 -> control
    set(case_status,which(is.na(case_status[[j]])),j,1)
    message(glue::glue("Control/Case by {j}"))

    print(dplyr::count(case_status,.data[[j]]))
  }
  return(case_status)
}
case_status<-get_case_status_by_source(source_lst=source_lst,definitions,lst.data,lst.data.settings,lst.identifiers=vct.identifiers,reference_date=setNames(as.Date(as.character(dfukb$f.53.0.0),format="%Y-%m-%d"),dfukb$f.eid))




# head(case_status)
#
# count(case_status,Any.sr)
# count(case_status,Any.hes)


#TODO
# a function to show (number of case) unique to each source

#add radarplots

