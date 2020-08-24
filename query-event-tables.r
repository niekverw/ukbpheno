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
    "tte.gpscript.bnf.scotland" = "BNF", # 2 BNF may be needed
    "tte.gpscript.read2.wales" = "READ"
    )
  
  return (datatable_defCol_pair)
  
  
}




# query single/specific table  
get_rows_with_matching_codes <- function(trait,datatable,dfDefinitions_processed_=NULL){
  # default ,maybe parameter not needed at all? 
  if(is.null(dfDefinitions_processed_)) {
    dfDefinitions_processed_<- dfDefinitions_processed
  } 
  
  
  # TODO here
  datatable_<-lst$tte.death.icd10.primary
  names(datatable_)
  def_col_name<-default_datatable_defCol_pair() 
  def_col_name<-default_datatable_defCol_pair()["tte.death.icd10.primary"]
  class(def_col_name)
  def_col_name
  dfDefinitions_processed_<-dfDefinitions_processed
  trait="Cad"
  #  look up the codes from dfDefinitions_processed
  target_trait<-dfDefinitions_processed_[dfDefinitions_processed_$TRAIT==trait,]
  # codes are comma separated after processing, parse them into the search pattern
  #  [[]] for vector   "$" doesn't work for a variable (which needs to be evaluated first)
  codes_pat<-paste(unlist(strsplit(target_trait[[def_col_name]],",")), collapse="|")
  
  #dplyr NOT exact match i.e. I25 gives I251 ,I252,etc
  matched_rows <- filter(datatable, grepl(codes_pat, datatable$code))
  return (matched_rows)
  
}
 

# query all tables
get_all_events <- function (trait,lst_datatables){
  # take the same data structure , a list of subsetted tables
  all_event_lst <- list()
  
  all_event_lst<-lapply(lst_datatables, function(x) {
    get_rows_with_matching_codes(trait=trait,x)
    } )
  return (all_event_lst)
}



# lst_datatables=lst
names(lst)
colnames(dfDefinitions_processed)[c(1,11:28)]










