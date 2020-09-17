
to_datatype <- function(vct=c(),datatype){
  if(datatype=='numeric'){return(as.numeric(vct)) }
  if(datatype=='character'){return(as.character(vct)) }
  return(vct)
}


get_all_events <- function (definition,lst.data=lst.data,lst.data.settings){
  # definitions=dfDefinitions_processed_expanded[9,]
  # look up for all dataframes in list
  if(nrow(definition)>1){
    message("ERROR: Provide one definition")
    return(NULL)
  }
  message(paste("querying the following classifications: " ,paste(names(definition)[names(definition) %in% lst.data.settings$classification],collapse=", ")))
  
  all_event_lst<-lapply(names(lst.data), function(x) {
    classification=lst.data.settings %>% filter(datasource == x) %>% pull(classification)
    datatype=lst.data.settings %>% filter(datasource == x) %>% pull(datatype)
    codes <- to_datatype(strsplit(definition[,classification],split = ",")[[1]],datatype)
    lst.data[[x]][.(codes),nomatch=NULL] # nomatch is important, otherwise it will return row with NA if it didnt find the code (but which on the other hand may also maybe good to keep track?? )
    
  } )
  # name the dfs in list
  names(all_event_lst)<-names(lst.data)
  # remove empty dfs frame list 
  all_event_lst <- all_event_lst[lapply(all_event_lst,nrow)>0]
  
  # set key to be eid
  all_event_dt <- plyr::ldply(all_event_lst, data.frame) %>% as.data.table()
  all_event_dt$classification <- lst.data.settings[match(all_event_dt$.id ,lst.data.settings$datasource),]$classification
  setkey(all_event_dt,f.eid)    
  return (all_event_dt)
}
# 
# dfDefinitions_processed_expanded <- expand_dfDefinitions_processed(dfDefinitions_processed,lst.data.settings=lst.data.settings(),lst.counts = lst.counts)
# test<-get_all_events(dfDefinitions_processed_expanded[9,],lst) 
# key(test$tte.death.icd10.secondary)
# 
# 
# 
# # lst_datatables=lst
# colnames(dfDefinitions_processed)[c(1,11:28)]
# 


# TODO composite phenotype from multiple sources get_all_events + TS
# simple count
# +time

