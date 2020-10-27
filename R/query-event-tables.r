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
get_all_events <- function (definition,lst.data=lst.data,df.data.settings,verbose=TRUE){

  # look up for all dataframes in list
  if(is.null(definition) | nrow(definition)==0 ){
    #message("no input given")
    return(NULL)
  }
  if(nrow(definition)>1){
    message("ERROR: Provide one definition")
    return(NULL)
  }

  if(verbose){
    message(paste("querying the following classifications: " ,paste(names(definition)[names(definition) %in% df.data.settings$classification],collapse=", ")))

    message("Filter records by thresholds specified in data setting")
  }
  all_event_lst<-lapply(names(lst.data), function(x) {
    classification=df.data.settings %>% dplyr::filter(datasource == x) %>% dplyr::pull(classification)
    datatype=df.data.settings %>% dplyr::filter(datasource == x) %>% dplyr::pull(datatype)
    codes <- to_datatype(strsplit(definition[,classification],split = ",")[[1]],datatype)
    dt<-lst.data[[x]][.(codes),nomatch=NULL] # nomatch is important, otherwise it will return row with NA if it didnt find the code (but which on the other hand may also maybe good to keep track?? )

    ############################################
    # min instance filter
    ############################################
    # # event=1 registry record /event=0 not real eventdate /event=2 self report eventdate
    # # for event =2 there is always an extra row event=0 (refer to conversion functions for self-reported field)
    # # some events do not have reported date hence only 1 row event=0
    # # keep events that are either 1)event ==1/2 [all events with event dates] 2) event=0 + id not in id_event2 (inds with actual event dates) [real event without eventdate reported]
    id_event2<-unique(dt[event==2,nomatch=NULL]$identifier)
    # this step only needed for the data type with event ==2
    if (length(id_event2)>0){
      dt<-dt[((!identifier %in% id_event2)|(event!=0)),nomatch=NULL]
    }
    ###########################################
    min.ins=df.data.settings %>% dplyr::filter(datasource == x) %>% dplyr::pull(minimum_instance)
    # lst.data[[x]]%>% dplyr::group_by(identifier)%>% dplyr::filter(n()>=min.ins)
    # this more memory / time efficient.....
    dt[,if(.N>=min.ins).SD,by=identifier,nomatch=NULL]

  })

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

  ############################################
  # min instance filter
  ############################################
  # if(verbose){
  #   message("Filter records by thresholds specified in data setting")
  # }
  # # event=1 registry record /event=0 not real eventdate /event=2 self report eventdate
  # # for event =2 there is always an extra row event=0 (refer to conversion functions for self-reported field)
  # # some events do not have reported date hence only 1 row event=0
  # # keep events that are either 1)event ==1/2 [all events with event dates] 2) event=0 + id not in id_event2 (inds with actual event dates) [real event without eventdate reported]
  # id_event2<-unique(all_event_dt[event==2,nomatch=NULL]$identifier)
  # all_event_dt<-all_event_dt[((!identifier %in% id_event2)|(event!=0)),nomatch=NULL]
  #
  # # first merge the datasource specific threshold values to df
  # all_event_dt$min.ins<- with(df.data.settings, minimum_instance[match(all_event_dt$.id,datasource)])
  #
  # # TODO how to make the message look better?
  # # #records by data type after filtering:c("tte.death.icd10.primary", "tte.death.icd10.secondary", "tte.hesin.icd10.primary", "tte.hesin.icd10.secondary", "tte.hesin.icd9.primary", "tte.hesin.icd9.secondary", "tte.hesin.oper4.primary", "tte.hesin.oper4.secondary", "tte.sr.20002")
  # # c(8900, 75114, 1234716, 7340, 1711, 337, 244658, 1604, 129186)
  # # too much information to read
  # # message(glue::glue("#records by data type before filtering: {glue::glue_collapse(all_event_dt%>% dplyr::group_by(.id)%>% count(),sep='\n')}"))
  # # print(all_event_dt%>% dplyr::group_by(.id)%>% count())
  #
  # # split df by .id, , group by identifier and filter records
  # all_event_dt<-plyr::ddply(.data=all_event_dt,.variables=".id",function(x) {
  #   # print(head(x$.id,1))
  #   # print(nrow(x))
  #   x%>% dplyr::group_by(identifier)%>% dplyr::filter(n()>=min.ins)
  #   # ins.min<-df.data.settings[df.data.settings$datasource==".id",]$minimum_instance
  #   # print(ins.min)
  #   # x%>%dplyr::group_by(identifier)%>% dplyr::filter(n()>ins.min)
  # })
  # all_event_dt<-as.data.table(all_event_dt)
  # message(glue::glue("#records by data type after filtering:{glue::glue_collapse(all_event_dt%>% dplyr::group_by(.id)%>% count(),sep='\n')}"))
  #############################################################################################################

}


