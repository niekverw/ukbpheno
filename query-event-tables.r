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
  
  return (all_event_lst)
}



test<-get_all_events("Cad",lst) #list of 11 dfs 

# lst_datatables=lst
colnames(dfDefinitions_processed)[c(1,11:28)]

# TS 
ts_cols<-dfDefinitions_processed[dfDefinitions_processed$TRAIT=="Cad",]$TS
ts_cols #"6150=1[3894]"
ts_cols<-c(ts_cols,"615=1") #"6150=1[3894]" "615=1"  
ts_cols_<-sub("=.*", "", ts_cols) # "6150"  "615"
#  immediately followed by non-digit \\D+  to avoid partial match of field number?
ts_cols<-paste(unlist(ts_cols_),"\\D+", sep="",collapse="|") #"6150\\D+|615\\D+"
ts_cols
df_ts <- dfukb %>% select(matches(ts_cols))  # 12 vars f.6150.x.x
test <- dfukb %>% select(matches("6150|615")) #28 vars with also f.6153.x.x
  # data.table version                  
df[, grep("ABC", names(df)), with = FALSE]
                    





