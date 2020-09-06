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
tte.gpscript.read2.wales	READ	character	1	2	FALSE"))
}


to_datatype <- function(vct=c(),datatype){
  if(datatype=='numeric'){return(as.numeric(vct)) }
  if(datatype=='character'){return(as.character(vct)) }
  return(vct)
}

get_all_events <- function (definitions,lst_dfs=lst,datatable_defCol_pair=default_datatable_defCol_pair()){
  # definitions=dfDefinitions_processed_expanded[9,]
  # look up for all dataframes in list
  all_event_lst<-lapply(names(lst_dfs), function(x) {
    classification=datatable_defCol_pair %>% filter(datasource == x) %>% pull(classification)
    datatype=datatable_defCol_pair %>% filter(datasource == x) %>% pull(datatype)
    codes <- to_datatype(strsplit(definitions[,classification],split = ",")[[1]],datatype)
    lst_dfs[[x]][.(codes),nomatch=NULL] # nomatch is important, otherwise it will return row with NA if it didnt find the code (but which on the other hand may also maybe good to keep track?? )
    
  } )
  # name the dfs in list
  names(all_event_lst)<-names(lst_dfs)
  # remove empty dfs frame list 
  all_event_lst <- all_event_lst[lapply(all_event_lst,nrow)>0]
  
  # set key to be eid
  all_event_dt <- plyr::ldply(all_event_lst, data.frame) %>% as.data.table()
  setkey(all_event_dt,f.eid)    
  return (all_event_dt)
}
# 
# dfDefinitions_processed_expanded <- expand_dfDefinitions_processed(dfDefinitions_processed,datatable_defCol_pair=default_datatable_defCol_pair(),lst.counts = lst.counts)
# test<-get_all_events(dfDefinitions_processed_expanded[9,],lst) #list of 11 dfs 
# 
# key(test$tte.death.icd10.secondary)
# 
# 
# 
# # lst_datatables=lst
# colnames(dfDefinitions_processed)[c(1,11:28)]
# 
# # TODO TS 


process_ts_cols<- function(trait,df_ukb,dfDefinitions_processed_=NULL){
  # default ,maybe parameter not needed at all? 
  if(is.null(dfDefinitions_processed_)) {
    dfDefinitions_processed_<- dfDefinitions_processed
  } 
  # look up the fields needed in ts
  ts_col<-dfDefinitions_processed[dfDefinitions_processed$TRAIT=="HxHrt",]$TS
  # parse fields
  ts_conditions<-unlist(strsplit(ts_col,","))
  # clean out extra characters that was not removed in process definition function
  # str_extract from stringr  , extract only the conditions
  ts_cols<-str_extract(ts_conditions,"\\d+[=|<|>][=]*\\d+")   # "20110=1"  "20107=1"  "20111<=1"
  # extract the fields for selectiion dfukb
  ts_fields <-str_extract(ts_cols,"\\d+") # "20110" "20107" "20111"
  #  immediately followed by non-digit \\D+  to avoid partial match of field number?
  ts_cols_ukb<-paste(unlist(ts_fields),"\\D+", sep="",collapse="|") # "20110\\D+|20107\\D+|20111\\D+"
  #  hard coded eid safe?
  ts_cols_ukb<-paste("eid",ts_cols_ukb,sep="|")
  df_ts <- df_ukb %>% select(matches(ts_cols_ukb))  
  
  # create a list to store the result
  ts_lst<-list()
  # for each field listed in ts
  for (col in ts_cols) {
    # parse the field and condition 
    field_prefix<-str_extract(col,"\\d+")
    cdn<-str_extract(col,"[=|<|>][=]*\\d+")
    # replace to logical equal
    cdn<-gsub("^={1}","==",cdn)
    # add dot to string for filtering
    cdn<-paste(".",cdn,sep="")
    # https://suzan.rbind.io/2018/02/dplyr-tutorial-3/#filtering-across-multiple-columns
    # field_prefix
    # names(df_ts)
    # add eid 
    cols_to_keep<-paste("eid",field_prefix,sep="|")
    
    # doesn't work anymore?eval(parse(text=cdn)))  why passing variable doesn't work ?
    
    test2<- df_ts %>% filter_at(vars(contains("20111")), any_vars(.<=1)) 
    #843 obs. of 133 variables
    test3<- df_ts %>% filter_at(vars(contains(field_prefix)), any_vars(.<=1)) 
    #843 obs. of 133 variables
    print(parse(text=cdn))
    #expression(.<=1)
    test4<- df_ts %>% filter_at(vars(contains(field_prefix)), any_vars(eval(parse(text=cdn))))
    #0 obs. of 133 variables
    test5<- df_ts %>% filter_at(vars(contains("20111")), any_vars(eval(parse(text=".<=1"))))
    #0 obs. of 133 variables
    
    #non dplyr version
    test6 <- apply(df_ts[,grep(field_prefix,names(df_ts)),with=FALSE], 1, function(x) paste(x[!is.na(x) & x != ""],collapse = ","))
    class(test6[1])
    length(test6[1])
    temp<-unlist(strsplit(test6[1],","))
    unlist(strsplit(StrRxCodes,","))
   temp<-as.numeric(temp)
    temp <=1
    
    # test2<- df_ts %>% {if("20111" %in% names(.)) filter(., a <= 1) else .}
    # %>% select(grep(cols_to_keep, names(df_ts)))
    
    
    # keep only columns with rows fulfilling condition and eid
    # TODO this doesn't work
    test<-test2 %>% select_if(~any(eval(parse(text=cdn)))|grepl("eid",names(.)) )
    
    # ts_lst[[paste("ts",field_prefix,sep=".")]] <-test
    

  }
  
  
  

  #  look up the codes from dfDefinitions_processed
  target_trait<-dfDefinitions_processed_[dfDefinitions_processed_$TRAIT==trait,]
  # codes are comma separated after processing, parse them into the search pattern
  #  [[]] for vector   "$" doesn't work for a variable (which needs to be evaluated first)
  codes_pat<-paste(unlist(strsplit(target_trait[[def_col_name]],",")), collapse="|")
  
  #dplyr NOT exact match i.e. I25 gives I251 ,I252,etc
  matched_rows <- filter(df, grepl(codes_pat, df$code))
  return (matched_rows)
  
}



# 
# ts_cols<-dfDefinitions_processed[dfDefinitions_processed$TRAIT=="Cad",]$TS
# ts_cols #"6150=1[3894]"
# ts_cols<-c(ts_cols,"615=1") #"6150=1[3894]" "615=1"  
# ts_cols_<-sub("=.*", "", ts_cols) # "6150"  "615"
# #  immediately followed by non-digit \\D+  to avoid partial match of field number?
# ts_cols<-paste(unlist(ts_cols_),"\\D+", sep="",collapse="|") #"6150\\D+|615\\D+"
# ts_cols
# df_ts <- dfukb %>% select(matches(ts_cols))  # 12 vars f.6150.x.x
# test <- dfukb %>% select(matches("6150|615")) #28 vars with also f.6153.x.x




#   # data.table version                  
# df[, grep("ABC", names(df)), with = FALSE]
#                     
# ts_conds<-dfDefinitions_processed[dfDefinitions_processed$TRAIT=="Cad",]$TS
# ts_conds<-sub("\\[.*", "", ts_conds)
# ts_conds<-c(ts_conds,"615=1")
# 
# value
# strsplit(ts_conds, "[=]")
# [[1]]
# 

# TODO composite phenotype from multiple sources get_all_events + TS
# simple count
# +time

