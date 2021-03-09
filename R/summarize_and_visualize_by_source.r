
#' Get case-control status by data sources
#'
#' Given a phenotype, a list of episode data and reference dates per individual, identify cases and control status by source. 
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
                      reference_date=NULL
                      ) {
  ###########################################################################
  # NOTE: this function depends on the get_cases_controls! modify accordingly
  ##########################################################################
  lst.case_control <- get_cases_controls(definitions, lst.data,lst.data.settings,reference_date,lst.identifiers)
  
  # extract only relevant column from df_casecontrol 
  case_status<-lst.case_control$df.casecontrol[,c("f.eid","reference_date","Any")]
  # create new columns for each source
  for (source_name in names(source_lst)) {
    varname<-paste('Any',source_name,sep='.')
    # create new columns for each source
    case_status[, (varname)]<-as.numeric(NA)
    print(head(case_status))
    
    # lookup source from df with all episodes i.e. all_event_dt.Include_in_cases 
    # if rows in include_in_case which originated from the target source, look up the eid and set to 2
    case_status[[varname]][case_status$f.eid %in% lst.case_control$all_event_dt.Include_in_cases[lst.case_control$all_event_dt.Include_in_cases$.id %in% source_lst[[source_name]],]$f.eid] <- 2
    print(head(case_status))
    
    # case_status$Any.sr[case_status$f.eid %in% lst.case_control$all_event_dt.Include_in_cases[lst.case_control$all_event_dt.Include_in_cases$.id %in% sr_source,]$f.eid]
  }
  
  # for all Any.x columns
  for (j in   names(case_status)[grepl("Any",names(case_status))]){
    # if Any <0 , copy Any (excluded)
    set(case_status,which(case_status$Any<=0),j,case_status[case_status$Any<=0]$Any)
    #  if the cell value was NA -> set to 1 -> control 
    set(case_status,which(is.na(case_status[[j]])),j,1)

  }
  return(case_status)
}

# case_status<-get_case_status_by_source(source_lst,definitions,lst.data,lst.data.settings,reference_date=setNames(as.Date(as.character(dfukb$f.53.0.0),format="%Y-%m-%d"),dfukb$f.eid))
# head(case_status)
# 
# count(case_status,Any.sr)
# count(case_status,Any.hes)


#TODO
# a function to show (number of case) unique to each source 

#add radarplots

