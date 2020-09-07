


# TS per trait
process_ts_cols<- function(trait,df,dfDefinitions_processed=NULL){
  # default ,maybe parameter not needed at all? 
  if(is.null(dfDefinitions_processed)) {
    dfDefinitions_processed<- dfDefinitions_processed
  } 
  # look up the fields needed in ts
  trait<-"HxHrt"
  df<-dfukb

  # visit date code
  field_visit_date="53"
  # vector with name of the identifier col 
  identifierfield = names(df)[grepl("eid", names(df))]
  #  vector with names of all visit cols : "f.53.0.0" "f.53.1.0" "f.53.2.0" "f.53.3.0"
  visitdatefields = names(df)[grepl(paste0("[^0-9]",field_visit_date,"[^0-9]"), names(df))]
  visits = length(visitdatefields) #sum(grepl("53_", names(df)))
  
  # touchscreen col in processed definition table  # example "20110=1,20107=1[3894], 20111<=1" 
  ts_col<-dfDefinitions_processed[dfDefinitions_processed$TRAIT==trait,]$TS
  # parse touchscreen fields
  ts_conditions<-unlist(strsplit(ts_col,","))   #"20110=1"       "20107=1[3894]" "20111<=1" 
  
  df_out <-  matrix(ncol=5, nrow=0) # initiate output 
  
  # for each field listed in ts
  for (col in ts_conditions) {
    # parse the field and condition 
    col<-"20107=1[3894]"
    cdn<-str_extract(col,"[=|<|>][=]*\\d+")
    # replace one equal sign to logical equal
    cdn<-gsub("^={1}","==",cdn)
    
    field_ts_diagnosis<-str_extract(col,"\\d+")
    tsdiagnosisfields = names(df)[grepl( paste0("[^0-9]",field_ts_diagnosis,"[^0-9]"), names(df))]
    
    
    # optional col for age of diagnosis specified in bracket []   # https://stackoverflow.com/questions/52061753/r-capturing-string-inside-brackets
    # regmatches returns character(0) if there is no match
    age_diagnosis_col<-regmatches(col, regexpr("\\[\\K[^][]*", col,perl = TRUE))
    if(length(age_diagnosis_col)>0){tsdiagnosisdatefields = names(df)[grepl(age_diagnosis_col, names(df))]} else tsdiagnosisdatefields=NULL}
  
  
  # only need n_eid, visit dates and diag-codes + age-of-diag
  columns_to_keep = c(identifierfield,
                      visitdatefields,
                      tsdiagnosisfields, 
                      tsdiagnosisdatefields
  )
  # data.table - Setting with = FALSE disables the ability to refer to columns as if they are variables
  df_ <- df[,columns_to_keep,with=FALSE]
  df_$dummy <- NA
  
  for (v in 0:(visits-1)){ # for each visit, 
    message(paste0("querying visit ",v))
    diagfields = names(df_)[grepl(paste0("[^0-9]",field_ts_diagnosis,"[^0-9]",v),names(df_))]
    if(length(diagfields)==0){print(paste0("no data on visit ",v));next}
    # f.3894.v.0
    if(!is.null(tsdiagnosisdatefields)){diagdatefields = names(df_)[grepl(paste0("[^0-9]",age_diagnosis_col,"[^0-9]",v),names(df_))]}
    # f.53.v.0 
    visitdatefield = visitdatefields[v+1]
    for (i in 1:length(diagfields)) { # for each occurence of diagfield, find the corresponding age and convert it to  date - code and rbind() to df_out. 
      #agefield = paste0("age_",v)
      diagfield = diagfields[i]
      
      if(length(age_diagnosis_col)>0){ 
        # length(diagdatefields) will always be 1
        if( length(diagdatefields) == 1) { diagdatefield = diagdatefields } # in case of death (40001/2)
      } else {
        diagdatefield = "dummy" # e.g. incase of medication
      }
      
      # print(paste0((diagfield), " - ", diagdatefield))
      #  empty diagnosis column
      if(all(is.na(df_[,diagfield,with=FALSE]))){next}
      # diagfield example f.20107.0.3
      # for rows with non-empty current diagfield, select identifier,diagfield,diagdatefield,visitdatefield 
      df_sub <- df_[!is.na(get(diagfield) ),c(identifierfield,diagfield,diagdatefield,visitdatefield),with=FALSE]
      # find rows that fulfill the condition
      cdn_exp <-paste(diagfield,cdn,sep="") #"f.20107.0.3 ==1"
      df_sub<- df_sub %>% filter(eval((parse(text=cdn_exp))))
      # replace the diagfield content with the condition
      df_sub$diagfield<-cdn_exp
      # add visit instance
      df_sub$visit <- v
      names(df_sub) <- c(identifierfield,"code","eventdate","visitdate","visit")
      df_out <- rbind(df_out,as.matrix(df_sub[,c(identifierfield,"code","eventdate","visit","visitdate"),with=FALSE]))
    }
    
  }
  
  
  
  
  # 
  # cdn<-gsub("\\d+={1}","==",str_extract(col,"\\d+[=|<|>][=]*\\d+"))
  # cdn<-"20107==1"
  # 
  # field_prefix<-str_extract(col,"\\d+")
  # cdn<-str_extract(col,"[=|<|>][=]*\\d+")
  # # replace to logical equal
  # cdn<-gsub("^={1}","==",cdn)
  # # add dot to string for filtering
  # cdn<-paste(".",cdn,sep="")
  # # https://suzan.rbind.io/2018/02/dplyr-tutorial-3/#filtering-across-multiple-columns
  # # field_prefix
  # # names(df_ts)
  # # add eid 
  
  
  
  names(df_ts)[grepl(paste0("[^0-9]",field_prefix,"[^0-9]",v),names(df_ts))]
  
  
  
  
  
  message("convert to dataframe")
  # df_out contains all visits , each row in df_out is a event
  df_out <- data.table(df_out,stringsAsFactors=F)
  df_out <- df_out[, visitdate:=as.Date(visitdate)]
  
  
  #TODO add the  34 and 52  and change in nurseinterview as well
  # compute the event date from visitdate and age of
  df_out <- df_out[, eventdate:=as.numeric(eventdate)] ## as number.  age of diagnosis 
  df_out[eventdate <0,'eventdate']<-NA
  df_out$eventdate = df_out[,"visitdate"] - (df_out[,"eventdate"]*daysinyear)
  
}
  
#   
#   # doesn't work anymore?eval(parse(text=cdn)))  why passing variable doesn't work ?
#   test<- df_ts %>% filter(eval((parse(text=cdn))))
#   cdn
#   test2<- df_ts %>% filter(eval((parse(text=cdn)))) 
#   
#   #843 obs. of 133 variables
#   test3<- df_ts %>% filter_at(vars(contains(field_prefix)), any_vars(.<=1)) 
#   #843 obs. of 133 variables
#   print(parse(text=cdn))
#   #expression(.<=1)
#   test4<- df_ts %>% filter_at(vars(contains(field_prefix)), any_vars(eval(parse(text=cdn))))
#   #0 obs. of 133 variables
#   test5<- df_ts %>% filter_at(vars(contains("20111")), any_vars(eval(parse(text=".<=1"))))
#   #0 obs. of 133 variables
#   
#   #non dplyr version
#   test6 <- apply(df_ts[,grep(field_prefix,names(df_ts)),with=FALSE], 1, function(x) paste(x[!is.na(x) & x != ""],collapse = ","))
#   class(test6[1])
#   length(test6[1])
#   temp<-unlist(strsplit(test6[1],","))
#   unlist(strsplit(StrRxCodes,","))
#   temp<-as.numeric(temp)
#   temp <=1
#   
#   # test2<- df_ts %>% {if("20111" %in% names(.)) filter(., a <= 1) else .}
#   # %>% select(grep(cols_to_keep, names(df_ts)))
#   
#   
#   # keep only columns with rows fulfilling condition and eid
#   # TODO this doesn't work
#   test<-test2 %>% select_if(~any(eval(parse(text=cdn)))|grepl("eid",names(.)) )
#   
#   # ts_lst[[paste("ts",field_prefix,sep=".")]] <-test
#   
#   
# }
# 
# 
# 
# 
# #  look up the codes from dfDefinitions_processed
# target_trait<-dfDefinitions_processed_[dfDefinitions_processed_$TRAIT==trait,]
# # codes are comma separated after processing, parse them into the search pattern
# #  [[]] for vector   "$" doesn't work for a variable (which needs to be evaluated first)
# codes_pat<-paste(unlist(strsplit(target_trait[[def_col_name]],",")), collapse="|")
# 
# #dplyr NOT exact match i.e. I25 gives I251 ,I252,etc
# matched_rows <- filter(df, grepl(codes_pat, df$code))
# return (matched_rows)
# 
# }
