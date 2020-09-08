get_stats_for_events <- function(all_event_dt){
  # show stats on codes
  stats.codes <- all_event_dt[, .(count=.N,sum.event = sum(event,na.rm = T),sum.epidur= sum(epidur,na.rm = T),median.epidur= median(epidur,na.rm = T),max.epidur= max(epidur,na.rm=T) ), keyby=list(f.eid,classification, code)]
  stats.codes <- stats.codes %>% dplyr::group_by(classification, code) %>% summarise(count=n() )
  stats.codes <- stats.codes %>% arrange(count)
  stats.codes$rank <- 1:nrow(stats.codes)
  p1 <- ggplot(stats.codes, aes(rank, count,label=code,color=classification)) + geom_point() + ylim(-((max(stats.codes$count))/3),NA)  + geom_text_repel(size =3,segment.size=0.5)
  p2 <- ggplot(stats.codes, aes(rank, count,label=code,color=classification)) + scale_y_continuous(trans='log2') + geom_point()   + geom_text_repel(size =3,segment.size=0.5)
  stats.codes.summary.table <- stats.codes
  stats.codes.summary.p <- ggarrange(p1,p2,nrow = 1, ncol = 2,common.legend = TRUE)
  
  # show co-occurences of classifications 
  stats.coocurrence <- table(all_event_dt[,c("f.eid" ,"classification")])#[1:10,]
  stats.coocurrence[stats.coocurrence>0] <-1
  mat <- crossprod(as.matrix(stats.coocurrence))     
  mat <- floor(t(mat * 100 / diag(mat)))                 # calculate the percentage
  diag(mat) <- NA
  stats.class.cooccur.table <- mat
  stats.class.cooccur.p <- pheatmap::pheatmap(mat,display_numbers=mat,cluster_cols = F,cluster_rows = F )
  
  # show co-occurences of codes
  stats.coocurrence <- table(all_event_dt[,c("f.eid" ,"code")])
  stats.coocurrence[stats.coocurrence>0] <-1
  mat <- crossprod(as.matrix(stats.coocurrence))     
  mat <- floor(t(mat * 100 / diag(mat)))                 # calculate the percentage
  diag(mat) <- NA
  stats.codes.cooccur.table <- mat
  stats.codes.cooccur.p <- pheatmap::pheatmap(mat,fontsize = 6,cluster_cols = F,cluster_rows = F)
  
  # filter on codes that with co-occurence of at least 10% to reduce sparseness and make clustering more informative. 
  stats.codes.cooccur.filtered.table <- stats.codes.cooccur.table[rowMaxs(stats.codes.cooccur.table,na.rm=T)>10,colMaxs(stats.codes.cooccur.table,na.rm=T)>10]
  stats.codes.cooccur.filtered.p <- pheatmap::pheatmap( stats.codes.cooccur.table.filtered  ,fontsize = 6)
  
  return(list(stats.codes.summary.table = stats.codes.summary.table,
         stats.codes.summary.p = stats.codes.summary.p,
         stats.class.cooccur.table = stats.class.cooccur.table,
         stats.class.cooccur.p = stats.class.cooccur.p,
         stats.codes.cooccur.table = stats.codes.cooccur.table,
         stats.codes.cooccur.p = stats.codes.cooccur.p,
         stats.codes.cooccur.filtered.p = stats.codes.cooccur.filtered.p))
}

### get prevalence (move this.. )
get_incidence_prevalence <- function(all_event_dt,reference_date,
                                     include_secondary_recurrence=TRUE,
                                     datatable_defCol_pair = default_datatable_defCol_pair(),
                                     return_dates=FALSE) {
  # event==0, event cannot be used forr age - of - diagnosiis or new events.
  # event==1, event can be used for any type of future events, as it is based on ICD10 type of data 
  # event==2, event can be used for any type of future events and age of diagnoses but only if there is no evidence in history, as it is self reported.
  
  # define window arround reference_dates? -- e.g. if you're diagnosed with XX, count follow-up events only after 20 days of no-events. 
  # define window to score the diagnosis on refereence date.
  # reference_dates # should be a vector of dates, with the identifier as name. 
  #reference_date = setNames(as.Date(as.character(dfukb$f.53.0.0),format="%Y-%m-%d"),dfukb$f.eid)
  
  ### Future: 1) count any new event==1 code (icd10 etc) 
  ###         2) count any primary event==1 if Hx==1 (recurrence) if 'include_secondary_recurrence'==TRUE
  ###         3) countany new event==2 (self report) age of diagnosis in future
  
  
  # 4139206
  df.ref <- data.frame(reference_date)
  df.ref$pid <- rownames(df.ref)
  
  df <- merge(all_event_dt,df.ref,by.x = 'f.eid' , by.y='pid',all.x = T) %>% arrange(eventdate) %>% as.data.table()
  df <- df %>% filter(!is.na(reference_date)) # comment out if missing f.eids is fixed. 
  
  df$days <- df$eventdate - df$reference_date
  setkey(df,days) # i don't know why, but setkey was alreaday on f.eid and cannot refresh..
  setkey(df,f.eid)
  
  if(include_secondary_recurrence){
    sources_recurrence_events <- datatable_defCol_pair %>% filter(diagnosis==1) %>% pull(datasource)
  } else {
    sources_recurrence_events <- datatable_defCol_pair$datasource
  } 
  ### History
  dfHx <- df[days<=0]
  Hx_days <- dfHx[event>0][, .(Hx_days=min(days,na.rm=T) ), keyby=list(f.eid)]
  dfHx[,Hx:=1]
  
  ### Future
  dfFu <- df[  ((days>(0) & event==1) & (!f.eid %in% dfHx$Hx)) |
               (days>(0) & event==1 & (f.eid %in% dfHx$Hx) & .id %in% sources_recurrence_events ) |
               (days>(0) & event==2  & (!f.eid %in% dfHx$Hx)) ]
  dfFu[,Fu:=1] #unique(dfFu$f.eid)
  Fu_days <- dfFu[,.(Fu_days= min(days,na.rm=T) ), keyby=list(f.eid)]
  
  ### age of diagnosis
  #system.time({ df %>% filter(event>0) %>% group_by(f.eid) %>% summarise(first_diagnosis_days=min(days)) }) # <- slow.. 
  #system.time({ df[df$event>0][,.(first_diagnosis_days=min(days,na.rm=T) ), by=f.eid] }) # <- fast..
  #df[df$event>0][df[, .I[which.max(days)], by=f.eid]$V1] # <- aanother way.. 
  first_diagnosis_days <- df[df$event>0][,.(first_diagnosis_days=min(days,na.rm=T) ), by=f.eid] # 
  
  ## some other stats, and to include event==0 individuals:
  stats <- all_event_dt[, .(count = .N,sum.epidur= sum(epidur,na.rm = T),median.epidur= median(epidur,na.rm = T),max.epidur= max(epidur,na.rm=T)), by = f.eid]
  stats[is.infinite(stats$max.epidur),]$max.epidur <-NA
  
  # test <- Reduce(function(...) merge(..., all = TRUE,by='f.eid'), list(df,
  #                                                                      Hx_days,
  #                                                                      Fu_days,
  #                                                                      unique(dfHx[,c("f.eid","Hx")]),
  #                                                                      unique(dfFu[,c("f.eid","Fu")]),
  #                                                                      first_diagnosis_days,
  #                                                                      stats)
  # )
  # View(test)
  
  all_event_dt.summary <- Reduce(function(...) merge(..., all = TRUE,by='f.eid'), list( 
  stats,
  Hx_days,
  Fu_days,
  unique(dfHx[,c("f.eid","Hx")]),
  unique(dfFu[,c("f.eid","Fu")]),
  first_diagnosis_days
  ))
    
  all_event_dt.summary[,Any:=1]
  
  
  return(all_event_dt.summary)
}
