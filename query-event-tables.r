#  columns not yet used in definitions
  # TS (touchscreen) take UKB columns
  # Minimum_Episode_duration
  # include_secondary - maybe not needed anymore



default_datatable_defCol_pair <- function(){
  
  datatable_defCol_pair<- c(
    "tte.sr.20002" = "n_20002", # self reported non-cancer
    "tte.sr.20001" = "n_20001", # self reported cancer
    "tte.sr.20004" = "n_20004", # self reported operation
    "sr.20003" = "n_20003",     # self reported medication
    "tte.death.icd10.primary" = "ICD10",
    "tte.death.icd10.secondary" = "ICD10",
    "tte.hesin.oper3.primary" = "OPCS3",
    "tte.hesin.oper3.secondary" = "OPCS3",
    "tte.hesin.oper4.primary" = "OPCS4",
    "tte.hesin.oper4.secondary" = "OPCS4",
    "tte.hesin.icd10.primary" = "ICD10",
    "tte.hesin.icd10.secondary" = "ICD10",
    "tte.hesin.icd9.primary" = "ICD9",
    "tte.hesin.icd9.secondary" = "ICD9",
    "tte.gpclincal.read2" = "READ",
    "tte.gpclincal.read3" = "CTV3",
    "tte.gpscript.dmd.england" = "DMD",
    "tte.gpscript.bnf.england" = "BNF",
    "tte.gpscript.bnf.scotland" = "BNF", # 2 BNF may be needed in definition tables
    "tte.gpscript.read2.wales" = "READ"
    )
  
  return (datatable_defCol_pair)
  
  
}


# query single/specific table  
get_rows_with_matching_codes <- function(trait,df,def_col_name,dfDefinitions_processed_=NULL){
  # default ,maybe parameter not needed at all? 
  if(is.null(dfDefinitions_processed_)) {
    dfDefinitions_processed_<- dfDefinitions_processed
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
 



# query all tables
get_all_events <- function (trait,lst_dfs){
  # take the same data structure , a list of subsetted tables
  names_dfLst<-names(lst_dfs)
  
  all_event_lst <- list()
  col_name_map<-default_datatable_defCol_pair()
  

  # look up for all dataframes in list
  all_event_lst<-lapply(names(lst_dfs), function(x) {
    get_rows_with_matching_codes(trait=trait,lst_dfs[[x]],col_name_map[x])
    } )
  # name the dfs in list
  names(all_event_lst)<-names_dfLst
  # remove empty dfs frame list
  all_event_lst<-all_event_lst[map(all_event_lst, function(x) nrow(x)) > 0]
  # set key to be eid
  all_event_lst<-lapply(all_event_lst,function(x) {
    setkey(x,f.eid)    
  })
  
  
  return (all_event_lst)
}



test<-get_all_events("Cad",lst) #list of 11 dfs 

key(test$tte.death.icd10.secondary)



# lst_datatables=lst
colnames(dfDefinitions_processed)[c(1,11:28)]

# TODO TS 


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

