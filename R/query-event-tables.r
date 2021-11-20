#' @export
to_datatype <- function(vct=c(),datatype){
  if(datatype=='numeric'){return(as.numeric(vct)) }
  if(datatype=='character'){return(as.character(vct)) }
  return(vct)
}

#' Get all episodes for a phenotype
#'
#' Given a phenotype and a list of episode data , extract events for this phenotype from all data sources
#' @param definition phenotype/trait specified in definition table (a row in the table)
#' @param lst.data list of data table with all episode data collapsed to 1 datatable
#' @param df.data.settings data frame containing data settings
#' @return  a data table with all events
#' @keywords time-to-event
#' @export
#' @examples
#' get_all_events(dfDefinitions_processed_expanded[14,],lst.data,df.data.settings)
get_all_events <- function (definition,lst.data=lst.data,df.data.settings){

  # look up for all dataframes in list
  if(is.null(definition) | nrow(definition)==0 ){
    #message("no input given")
    return(NULL)
  }
  if(nrow(definition)>1){
    message("ERROR: Provide one definition")
    return(NULL)
  }

  message(paste("querying the following classifications: " ,paste(names(definition)[names(definition) %in% df.data.settings$classification],collapse=", ")))

  all_event_lst<-lapply(names(lst.data), function(x) {
    classification=df.data.settings %>% dplyr::filter(datasource == x) %>% dplyr::pull(classification)
    datatype=df.data.settings %>% dplyr::filter(datasource == x) %>% dplyr::pull(datatype)
    codes <- to_datatype(strsplit(definition[,classification],split = ",")[[1]],datatype)
    lst.data[[x]][.(codes),nomatch=NULL] # nomatch is important, otherwise it will return row with NA if it didnt find the code (but which on the other hand may also maybe good to keep track?? )

  } )
  # name the dfs in list
  names(all_event_lst)<-names(lst.data)
  # remove empty dfs frame list
  all_event_lst <- all_event_lst[lapply(all_event_lst,nrow)>0]

  # set key to be eid
  all_event_dt <- plyr::ldply(all_event_lst, data.frame) %>% data.table::as.data.table()
  all_event_dt$classification <- df.data.settings[match(all_event_dt$.id ,df.data.settings$datasource),]$classification
  if (nrow(all_event_dt) >0){
  data.table::setkey(all_event_dt,identifier)
  }
    return (all_event_dt)

}


# TODO composite phenotype from multiple sources get_all_events + TS
# simple count
# +time

