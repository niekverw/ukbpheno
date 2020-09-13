#  columns not yet used in definitions
  # TS (touchscreen) take UKB columns
  # Minimum_Episode_duration
  # include_secondary - maybe not needed anymore


default_datatable_defCol_pair <- function() {
  data.frame(fread("datasource	classification	datatype	expand_codes	diagnosis	ignore.case
tte.sr.20002	n_20002	numeric	0	1	FALSE
tte.sr.20001	n_20001	numeric	0	1	FALSE
tte.sr.20004	n_20004	numeric	0	1	FALSE
sr.20003	n_20003	numeric	0	1	FALSE
tte.death.icd10.primary	ICD10	character	1	1	TRUE
tte.death.icd10.secondary	ICD10	character	1	2	TRUE
tte.hesin.oper3.primary	OPCS3	character	1	1	TRUE
tte.hesin.oper3.secondary	OPCS3	character	1	2	TRUE
tte.hesin.oper4.primary	OPCS4	character	1	1	TRUE
tte.hesin.oper4.secondary	OPCS4	character	1	2	TRUE
tte.hesin.icd10.primary	ICD10	character	1	1	TRUE
tte.hesin.icd10.secondary	ICD10	character	1	2	TRUE
tte.hesin.icd9.primary	ICD9	character	1	1	TRUE
tte.hesin.icd9.secondary	ICD9	character	1	2	TRUE
tte.gpclincal.read2	READ	character	0	2	FALSE
tte.gpclincal.read3	CTV3	character	0	2	FALSE
tte.gpscript.dmd.england	DMD	character	0	2	FALSE
tte.gpscript.bnf.england	BNF	character	0	2	FALSE
tte.gpscript.bnf.scotland	BNF	character	0	2	FALSE
tte.gpscript.read2.wales	READ	character	1	2	FALSE
ts	TS	character	0	1	TRUE"))
}

to_datatype <- function(vct=c(),datatype){
  if(datatype=='numeric'){return(as.numeric(vct)) }
  if(datatype=='character'){return(as.character(vct)) }
  return(vct)
}


process_ts_col <-function (definitions){
  # remove square bracket
  definitions$TS<-  gsub( " *\\[.*?\\] *", "", definitions$TS)
  # replace single equal sign = to == /larger/less symbol    ;\\b word boundary
  definitions$TS<-gsub("\\b[=]+\\b","==",definitions$TS,perl=TRUE)
  definitions$TS<-gsub("\\b[≥]\\b",">=",definitions$TS,perl=TRUE)
  definitions$TS<-gsub("\\b[≤]\\b","<=",definitions$TS,perl=TRUE)
  return(definitions)
}


get_all_events <- function (definition,lst.data=lst.data,datatable_defCol_pair=default_datatable_defCol_pair()){
  # definitions=dfDefinitions_processed_expanded[9,]
  # look up for all dataframes in list
  if(nrow(definition)>1){
    message("ERROR: Provide one definition")
    return(NULL)
  }
  message(paste("querying the following classifications: " ,paste(names(definition)[names(definition) %in% datatable_defCol_pair$classification],collapse=", ")))
  
  # TS needs cleaning
  definitions<-process_ts_col(definitions)
  
  all_event_lst<-lapply(names(lst.data), function(x) {
    classification=datatable_defCol_pair %>% filter(datasource == x) %>% pull(classification)
    datatype=datatable_defCol_pair %>% filter(datasource == x) %>% pull(datatype)
    codes <- to_datatype(strsplit(definition[,classification],split = ",")[[1]],datatype)
    lst.data[[x]][.(codes),nomatch=NULL] # nomatch is important, otherwise it will return row with NA if it didnt find the code (but which on the other hand may also maybe good to keep track?? )
    
  } )
  # name the dfs in list
  names(all_event_lst)<-names(lst.data)
  # remove empty dfs frame list 
  all_event_lst <- all_event_lst[lapply(all_event_lst,nrow)>0]
  
  # set key to be eid
  all_event_dt <- plyr::ldply(all_event_lst, data.frame) %>% as.data.table()
  all_event_dt$classification <- datatable_defCol_pair[match(all_event_dt$.id ,datatable_defCol_pair$datasource),]$classification
  setkey(all_event_dt,f.eid)    
  return (all_event_dt)
}
# 
# dfDefinitions_processed_expanded <- expand_dfDefinitions_processed(dfDefinitions_processed,datatable_defCol_pair=default_datatable_defCol_pair(),lst.counts = lst.counts)
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

