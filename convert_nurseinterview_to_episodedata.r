convert_nurseinterview_to_episodedata <- function(df,field_sr_diagnosis = "20002",field_sr_date = "20009",qc_treshold_year=10){
  # df
  # field_sr_diagnosis = "20002"
  # field_sr_date = "20009"
  # 
  daysinyear=365.25
  visits = sum(grepl("53_", names(df)))
  
  # only need n_eid, visit dates and diag-codes + age-of-diag
  df <- df[,c("n_eid",names(df)[grepl("53_", names(df))],
              names(df)[grepl(field_sr_diagnosis, names(df))], 
              names(df)[grepl(field_sr_date, names(df))]
  )]
  
  dfout <-  matrix(ncol=3, nrow=0) # initiate output 
  for (v in 0:(visits-1)){ # for each visit, 
    visitdatefield = paste0("ts_53_",v,"_0")
    df[,visitdatefield] <- as.Date (df[,visitdatefield],format="%d%b%Y")
    
    diagfields = names(df)[grepl(paste0("n_",field_sr_diagnosis,"_",v),names(df))]
    diagagefields = names(df)[grepl(paste0("n_",field_sr_date,"_",v),names(df))]

    for (i in 1:length(diagfields)){ # for each occurence of diagfield, find the corresponding age and convert it to  date - code and rbind() to dfout. 
      #agefield = paste0("age_",v)
      
      diagfield = diagfields[i]
      diagagefield = diagagefields[i]
      
      print(paste0((diagfield), "- ", diagagefield))
      
      dfsub <- df[!is.na(df[,diagfield] ) ,c("n_eid", visitdatefield,diagfield,diagagefield)]
      dfsub[,diagagefield][dfsub[,diagagefield] <0] <- NA      
      dfsub$eventdate = dfsub[,visitdatefield] - (dfsub[,diagagefield]*daysinyear)

      dfout <- rbind(dfout,as.matrix(dfsub[,c("n_eid","eventdate",diagfield)]))
    }
    
  }
  
  dfout <- as.data.frame(dfout,stringsAsFactors=F)
  names(dfout) <- c("n_eid","eventdate",paste0(field_sr_diagnosis))
  dfout$eventdate <- as.Date(dfout$eventdate,"%Y-%m-%d")
  dfout[,field_sr_diagnosis] <- as.integer(dfout[,field_sr_diagnosis])
  
  
  # deduplicate, min/max/mean/sd <- not very efficient? 
  library(dplyr)
  dfout_extrastats<- dfout %>% group_by(n_eid,!!as.name(field_sr_diagnosis)) %>%
    mutate(mindt = min(eventdate, na.rm = TRUE),maxdt = max(eventdate, na.rm = TRUE),meandt = mean(eventdate, na.rm = TRUE))
  
  dfout_extrastats$diffdt <- (dfout_extrastats$maxdt - dfout_extrastats$mindt)/daysinyear
  dfout_extrastats[dfout_extrastats$diffdt>qc_treshold_year ,"meandt"] <- NA
  
  dfout$eventdate <- dfout_extrastats$meandt
  dfout <- dfout[!dfout[,field_sr_diagnosis] %in% 99999,]
  
  dfout <- dfout[!duplicated(dfout),]
  
  head(dfout)
  return(dfout)
}
