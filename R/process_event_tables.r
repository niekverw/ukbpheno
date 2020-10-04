

#' Summary statistics of events for single phenotype
#'
#' Given a data.table with all events for a phenotype, get the summary of code occurence/occurence
#' @param all_event_dt data table containing all events
#' @return  a list of 3 summary tables and 4 plots
#' @keywords event stats
#' @export
#' @examples
#' all_event_dt <- get_all_events(dfDefinitions_processed_expanded[1,],lst.data,lst.data.settings)
#' get_stats_for_events(all_event_dt)
get_stats_for_events <- function(all_event_dt){
  # show stats on codes
  stats.codes <- all_event_dt[, .(count=.N,sum.event = sum(event,na.rm = T),sum.epidur= sum(epidur,na.rm = T),median.epidur= median(epidur,na.rm = T),max.epidur= max(epidur,na.rm=T) ), keyby=list(f.eid,classification, code)]
  stats.codes <- stats.codes %>% dplyr::group_by(classification, code) %>% summarise(count=n() )
  stats.codes <- stats.codes %>% arrange(count)
  stats.codes$rank <- 1:nrow(stats.codes)
  png("test.png")
  p1 <- ggplot(stats.codes, aes(rank, count,label=code,color=classification)) + geom_point() + ylim(-((max(stats.codes$count))/3),NA)  + geom_text_repel(size =3,segment.size=0.5)
  p1
  dev.off()
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

  # stats.class.cooccur.p <- pheatmap::pheatmap(mat,display_numbers=mat,cluster_cols = F,cluster_rows = F )
  # #########ggplot2 implementation of the heatmap##########################################################################
  longData<- reshape2::melt(mat)
  names(longData)<- c("Code_presence","Code_occur","%")
  stats.class.cooccur.p <- ggplot(longData, aes(x = Code_occur , y =Code_presence,fill=`%`)) +
    geom_tile()+
    scale_fill_distiller(palette = "RdYlBu",na.value = "grey85") +
    labs(x="Classfication of co-occurence", y="Classification present",title="Co-occurence of diagnosis by sources",caption = "[In the presence of row, column coexists with row by %]") +
    geom_text(aes(label = `%`))
 # ########################################################################################################################
  # show co-occurences of codes
  stats.coocurrence <- table(all_event_dt[,c("f.eid" ,"code")])
  stats.coocurrence[stats.coocurrence>0] <-1
  mat <- crossprod(as.matrix(stats.coocurrence))
  mat <- floor(t(mat * 100 / diag(mat)))                 # calculate the percentage
  diag(mat) <- NA
  stats.codes.cooccur.table <- mat


  # stats.codes.cooccur.p <- pheatmap::pheatmap(mat,fontsize = 6,cluster_cols = F,cluster_rows = F)
  # filter on codes that with co-occurence of at least 10% to reduce sparseness and make clustering more informative.
  stats.codes.cooccur.filtered.table <- stats.codes.cooccur.table[rowMaxs(stats.codes.cooccur.table,na.rm=T)>10,colMaxs(stats.codes.cooccur.table,na.rm=T)>10]
  # stats.codes.cooccur.filtered.p <- pheatmap::pheatmap( stats.codes.cooccur.filtered.table  ,fontsize = 6)
  # ########## ggplot implementation#################################################################################
  # Create ggplot version dendrogram from ggdendro
  code.dendro <- as.dendrogram(hclust(d = dist(x = stats.codes.cooccur.filtered.table)))
  ddata_x <- dendro_data(code.dendro)
  # to colour leaves by classifications
  lab_gp <- label(ddata_x)
  lab_gp$group <- stats.codes[match(lab_gp$label,stats.codes$code),]$classification
  stats.codes.cooccur.filtered.p.dendro <- ggplot(segment(ddata_x)) +
    geom_segment(aes(x=x, y=y+10, xend=xend, yend=yend+10)) +  geom_text(data=label(ddata_x),
                 aes(label=label, x=x, y=-5, colour=lab_gp$group),size =3,angle=45)   +theme(axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
           # legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())
  # ggplot version heatmap
  # code.order <- order.dendrogram(code.dendro)
  # TODO maker heatmap ordered like dendrogram?
  longData<- reshape2::melt(stats.codes.cooccur.filtered.table)
  names(longData)<- c("Code_presence","Code_occur","%")
  stats.codes.cooccur.filtered.p.heat <- ggplot(longData, aes(x = Code_occur , y =Code_presence,fill=`%`)) +
    geom_tile()+
    scale_fill_distiller(palette = "RdYlBu",na.value = "grey85") +
    labs(x="Code of co-occurence", y="Code present",title="Co-occurence of diagnosis code",caption = "[In the presence of row, column coexists with row by %]") + theme(axis.text.x = element_text(angle = 45))
###############################################################################################################


  return(list(stats.codes.summary.table = stats.codes.summary.table,
         stats.codes.summary.p = stats.codes.summary.p,
         stats.class.cooccur.table = stats.class.cooccur.table,
         stats.class.cooccur.p = stats.class.cooccur.p,
         stats.codes.cooccur.table = stats.codes.cooccur.table,
         # stats.codes.cooccur.p = stats.codes.cooccur.p,
         # stats.codes.cooccur.filtered.p = stats.codes.cooccur.filtered.p,
         stats.codes.cooccur.filtered.p.dendro = stats.codes.cooccur.filtered.p.dendro,
         stats.codes.cooccur.filtered.p.heat = stats.codes.cooccur.filtered.p.heat))
}



#' Get data for phenotype incidence and prevalence
#'
#' Given a data.table with all events for a phenotype and reference dates per individual, compute the time to event data for each individual. If no reference date is given then the date of first available event will be taken as reference date for each individual. ...
#' @param all_event_dt data table containing all events
#' @param lst.data.settings data frame containing data settings
#' @param reference_date reference dates for each individuals in the whole cohort as a named vector
#' @param include_secondary_recurrence
#' @param window_ref_days_include number of days around the reference date (visit) that should be used to indicates if individuals had the event on the reference date. Relevant if you want to know if participant took medication on the visit
#' @param window_fu_days_mask number of days that future events should not be counted
#' @return  a data table with individual with time to event information
#' @keywords time-to-event
#' @export
#' @examples
#' all_event_dt <- get_all_events(dfDefinitions_processed_expanded[1,],lst.data,lst.data.settings)
#' get_incidence_prevalence(all_event_dt,lst.data.settings, reference_date = setNames(as.Date(as.character(dfukb$f.53.0.0),format="%Y-%m-%d"),dfukb$f.eid))
### get prevalence (move this.. )
get_incidence_prevalence <- function(all_event_dt,
                                     lst.data.settings,
                                     reference_date=NULL,
                                     include_secondary_recurrence=FALSE,
                                     #return_dates=FALSE, # TODO: not working yet.
                                     window_ref_days_include=0,##  indicate number of days around the reference date (visit) that should be used to indicates if individuals had the event on the reference date. For example, relevant if you want to know if participant took medication on the visit
                                     window_fu_days_mask=0 ## indicates number of days that future events should not be counted; e.g. you could only count events after 10 days from the reference visit to avoid events related to the reference date. e.g. you could also use it to only count events after X years, in order to avoid assesment bias.
                                     ) {
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

  # If a reference date is given, it only returns individuals that match non-missing reference dates.
  #print(glue::glue("window_ref_days_include: {window_ref_days_include}"))
  #print(glue::glue("window_fu_days_mask: {window_fu_days_mask}"))

  # reference_date <- reference_date[!is.na(reference_date)]


  if(length(reference_date)==0){reference_date<-NULL}
  if(!is.null(reference_date)){
    df_referencedate <- data.table(reference_date)
    df_referencedate$f.eid <- names(reference_date)

    message(glue::glue("non missing reference_date: {length(reference_date)}"))
  } else{
    message("no reference_date given, taking the first available event as reference")
    df_referencedate <- all_event_dt[,.(reference_date= min(eventdate,na.rm = T)),by=f.eid]
  }
  if(include_secondary_recurrence){
    sources_recurrence_events <- lst.data.settings$datasource
  } else {
    sources_recurrence_events <- lst.data.settings %>% filter(diagnosis==1) %>% pull(datasource)
  }


  df <- merge(all_event_dt,df_referencedate,by = 'f.eid') %>% arrange(eventdate) %>% as.data.table()
  # df <- df %>% filter(!is.na(reference_date)) # comment out if missing f.eids is fixed.
  df$days <- df$eventdate - df$reference_date
  setkey(df,days) # i don't know why, but setkey was alreaday on f.eid and cannot refresh..
  setkey(df,f.eid)

  ### flag primary death records
  df$death.primary<-ifelse((df$.id %in% lst.data.settings[lst.data.settings$death,]$datasource)&(lst.data.settings[match(df$.id ,lst.data.settings$datasource),]$diagnosis==1),1,0)
  ### flag secondary death records
  df$death.secondary<-ifelse((df$.id %in% lst.data.settings[lst.data.settings$death,]$datasource)&(lst.data.settings[match(df$.id ,lst.data.settings$datasource),]$diagnosis==2),1,0)

  ### History
  dfHx <- df[days<=0]
  # if death record is the only record , they will appear here Hx=0
  Hx_days <- suppressWarnings(unique(dfHx[event>0][, .(Hx_days=min(days,na.rm=T) ,death.primary,death.secondary), keyby=list(f.eid)] ))
  dfHx[,Hx:=2]

  ### Future
  # window_fu_days_mask
  dfFu <- df[  ((days>(0+window_fu_days_mask) & event==1) & (!f.eid %in% dfHx$Hx)) |
               (days>(0+window_fu_days_mask) & event==1 & (f.eid %in% dfHx$Hx) & .id %in% sources_recurrence_events ) |
               (days>(0+window_fu_days_mask) & event==2  & (!f.eid %in% dfHx$Hx)) ]
  dfFu[,Fu:=2] #unique(dfFu$f.eid)
  # records if death succeed an diagnosis
  Fu_days <- suppressWarnings( unique(dfFu[,.(Fu_days= min(days,na.rm=T),death.primary,death.secondary ), keyby=list(f.eid)] ))
  ### age of diagnosis
  #system.time({ df %>% filter(event>0) %>% group_by(f.eid) %>% summarise(first_diagnosis_days=min(days)) }) # <- slow..
  #system.time({ df[df$event>0][,.(first_diagnosis_days=min(days,na.rm=T) ), by=f.eid] }) # <- fast..
  #df[df$event>0][df[, .I[which.max(days)], by=f.eid]$V1] # <- aanother way..
  #
  # first_diagnosis_days <- df[  ,.(first_diagnosis_days=min(days,na.rm=T)), by=f.eid]
  #
  #
  # first_diagnosis_days <- df[  ,.SD[which.min(days)], by=f.eid]
  #
  # df %>% filter(f.eid==1025336 )
  # first_diagnosis_days %>% filter(f.eid==1025336 )

  ### Data if participant had event/med  reference date (visit);+/- x day
  #TODO I think episodes for visitdate (event=0 rows) is messing up with this flag! i.e. almost all of them have this flag on
  dfRef <- df[df$eventdate>=(df$reference_date-window_ref_days_include) & df$eventdate<=(df$reference_date+window_ref_days_include),c("f.eid")]
  dfRef[,Ref:=2]


  ## some other stats, and to include event==0 individuals:
  # allow abesence of HES data , sr/death does not have epidur
  if(!all(unique(all_event_dt$event==0)) & !all(is.na(all_event_dt$epidur))){
    stats <- suppressWarnings( all_event_dt[, .(count = .N,sum.epidur= sum(epidur,na.rm = T),median.epidur= median(epidur,na.rm = T),max.epidur= max(epidur,na.rm=T)), by = f.eid] )
    stats[is.infinite(stats$max.epidur),]$max.epidur <-NA
  } else{
    stats <- all_event_dt[, .(count = .N,sum.epidur= NA,median.epidur= NA,max.epidur=NA), by = f.eid]
  }



  # test <- Reduce(function(...) merge(..., all = TRUE,by='f.eid'), list(df,
  #                                                                      Hx_days,
  #                                                                      Fu_days,
  #                                                                      unique(dfHx[,c("f.eid","Hx")]),
  #                                                                      unique(dfFu[,c("f.eid","Fu")]),
  #                                                                      unique(dfRef[,c("f.eid","Ref")]),
  #                                                                      first_diagnosis_days,
  #                                                                      df.referencedate
  #                                                                      stats)
  # )
  # View(test)
  all_event_dt.summary <- Reduce(function(...) merge(..., all = TRUE,by='f.eid'), list(
      stats,
      Hx_days,
      Fu_days,
      unique(dfHx[,c("f.eid","Hx")]),
      unique(dfFu[,c("f.eid","Fu")]),
      unique(dfRef[,c("f.eid","Ref")])
  ))

  # combine the values from 2 columns and keep 1
  all_event_dt.summary$death.primary.x<-fcoalesce(all_event_dt.summary$death.primary.x,all_event_dt.summary$death.primary.y)
  all_event_dt.summary$death.secondary.x<-fcoalesce(all_event_dt.summary$death.secondary.x,all_event_dt.summary$death.secondary.y)
  all_event_dt.summary<-all_event_dt.summary %>% select(f.eid, count ,sum.epidur, median.epidur, max.epidur,Hx_days,Fu_days, death.primary.x, death.secondary.x,Hx, Fu, Ref)
  names(all_event_dt.summary)<- c("f.eid", "count" ,"sum.epidur", "median.epidur", "max.epidur","Hx_days","Fu_days", "death.primary", "death.secondary","Hx", "Fu", "Ref")
  # NA are not true event (=1) records , which are not death records
  set(all_event_dt.summary,which(is.na(all_event_dt.summary$death.primary)),"death.primary",0)
  set(all_event_dt.summary,which(is.na(all_event_dt.summary$death.secondary)),"death.secondary",0)

  all_event_dt.summary$first_diagnosis_days <- pmin(all_event_dt.summary$Hx_days,all_event_dt.summary$Fu_days,na.rm = T)
  all_event_dt.summary[ Hx==2 & is.na(Hx_days),'first_diagnosis_days'] <- NA

  all_event_dt.summary <- merge(all_event_dt.summary,df_referencedate,by="f.eid")
  all_event_dt.summary[,Any:=2]
  all_event_dt.summary <- data.table(all_event_dt.summary)
  return(all_event_dt.summary)


}

#' Get case for a phenotype
#'
#' Given a phenotype, a list of episode data and reference dates per individual, identify cases and compute the time to event data for these individuals.  If no reference date is given then the date of first available event will be taken as reference date for each individual.
#' @param definitions phenotype/trait specified in definition table (a row in the table)
#' @param lst.data list of data table with all episode data
#' @param lst.data.settings data frame containing data settings
#' @param reference_date reference dates for each individuals in the whole cohort as a named vector

#' @return  a list of 2 data tables : all events for valid cases and an event summary containing time to event information for these individuals.
#' @keywords time-to-event
#' @export
#' @examples
#' get_cases(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="Nicm"), lst.data,lst.data.settings, reference_date=setNames(as.Date(as.character(dfukb$f.53.0.0),format="%Y-%m-%d"),dfukb$f.eid))
get_cases <- function(definitions,
                       lst.data,
                       lst.data.settings,
                       reference_date=NULL,
                       ...
                         ) {
  # define cases
  if(length(unique(definitions$TRAIT))>1){
    message("more than 1 TRAIT in definitions")
  }
  all_event_dt.Include_in_cases <- get_all_events(definitions %>% filter(Definitions =="Include_in_cases"),lst.data,lst.data.settings)   #MI
  all_event_dt.Include_in_cases.summary <- get_incidence_prevalence(all_event_dt = all_event_dt.Include_in_cases,lst.data.settings,
                                                                    reference_date = reference_date)
                                                                    #...)

  message(glue::glue("including {nrow(all_event_dt.Include_in_cases.summary)} cases"))
  all_event_dt.Exclude_from_cases <- get_all_events(definitions %>% filter(Definitions =="Exclude_from_cases"),lst.data,lst.data.settings)   #MI
  if(!is.null(all_event_dt.Exclude_from_cases)){
    exclude=all_event_dt.Include_in_cases.summary$f.eid %in% unique(all_event_dt.Exclude_from_cases$f.eid)
    message(glue::glue("excluding {sum(exclude)} cases"))
    set_to_na <- names(all_event_dt.Include_in_cases.summary)[!names(all_event_dt.Include_in_cases.summary) %in% c("f.eid","reference_date")]
    # the .summary needed to be flagged for the get_case_control() functoin
    all_event_dt.Include_in_cases.summary[exclude,(set_to_na):=-2]
    # remove these in all_event_dt
    all_event_dt.Include_in_cases<-all_event_dt.Include_in_cases[! (all_event_dt.Include_in_cases$f.eid %in% unique(all_event_dt.Exclude_from_cases$f.eid)),]
    message(glue::glue("{nrow(all_event_dt.Include_in_cases)} events fulfiling criteria in include_in_cases"))
  }
  return(list(all_event_dt.Include_in_cases=all_event_dt.Include_in_cases,
              all_event_dt.Include_in_cases.summary=all_event_dt.Include_in_cases.summary)
  )
}



######
#' Get case and controls for a phenotype
#'
#' Given a phenotype, a list of episode data and reference dates per individual, identify valid cases/controls and compute the time to event data for these individuals respectively.  If no reference date is given then the date of first available event will be taken as reference date for each individual.
#' @param definitions phenotype/trait specified in definition table (a row in the table)
#' @param lst.data list of data table with all episode data
#' @param lst.data.settings data frame containing data settings
#' @param reference_date reference dates for each individuals in the whole cohort as a named vector
#' @param lst.identifiers character vector listing the identifiers in the cohort. This is used to define controls if reference_date is not given

#' @return  a list of 2 data tables : all events for valid cases and an event summary containing time to event information for these individuals.
#' @keywords time-to-event
#' @export
#' @examples
#' get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="Nicm"), lst.data,lst.data.settings,  reference_date=setNames(as.Date(as.character(dfukb$f.53.0.0),format="%Y-%m-%d"),dfukb$f.eid))
get_cases_controls <- function (definitions,
                                 lst.data,
                                 lst.data.settings,
                                 reference_date=NULL,
                                 lst.identifiers=NULL # Used to define controls if reference_date is not given (NULL)
) {

  #reference_date = setNames(as.Date(as.character(dfukb$f.53.0.0),format="%Y-%m-%d"),dfukb$f.eid)

  reference_date <- reference_date[!is.na(reference_date)]
  reference_date <- reference_date[!is.na(names(reference_date))]
  if(is.null(reference_date)){
    message("reference_date=NULL, taking first occurence of case and all available identifiers (lst.data$all_identifiers)")
    if(length(lst.identifiers)<=1){
      message("ERROR: lst.identifiers is empty, please provide either reference_date or a list of lst.identifiers that should be used as population.")
      return(0)
    }
  }
  # define population
  all_event_dt.population <- get_all_events(definitions %>% filter(Definitions =="Study_population"),lst.data,lst.data.settings)   #MI

  if(!is.null(all_event_dt.population)) {
    # only consider event with real event date in study population, set those with visitdate to NA
    all_event_dt.population[all_event_dt.population$event==0,]$eventdate <-NA
    # get everyone , take the date of relevant events respectively
    all_event_dt.population.summary <- get_incidence_prevalence(all_event_dt = all_event_dt.population,lst.data.settings,reference_date = NULL,window_fu_days_mask = 15)
    reference_date = setNames(as.Date(as.character(all_event_dt.population.summary$reference_date),format="%Y-%m-%d"),all_event_dt.population.summary$f.eid)
    message(glue::glue("Population: {length(unique(all_event_dt.population.summary$f.eid))} individuals "))
  } else {
    message(glue::glue("Population: total non-missing reference date = {sum(!is.na(reference_date))}, total missing reference date= {sum(is.na(reference_date))} "))
  }

  # define cases
  cases <- get_cases(definitions,lst.data,lst.data.settings,reference_date,window_ref_days_include=0,window_fu_days_mask=0)
  all_event_dt.Include_in_cases.summary <- cases$all_event_dt.Include_in_cases.summary
  all_event_dt.Include_in_cases <- cases$all_event_dt.Include_in_cases
  # define exclude controls
  all_event_dt.Exclude_from_controls <- get_all_events(definitions %>% filter(Definitions =="Exclude_from_controls"),lst.data,lst.data.settings)   #MI
  ### define case & control
  if(is.null(reference_date)){
    reference_date = setNames(as.Date(rep(NA,length(lst.identifiers))),lst.identifiers)
  }

  df.casecontrol <- data.frame(reference_date=reference_date) %>% tibble::rownames_to_column('f.eid') %>% as.data.table()

  # exlude id in case_include & case_exclude from summary => potential control
  df.casecontrol <- df.casecontrol[!df.casecontrol$f.eid %in% all_event_dt.Include_in_cases.summary$f.eid,]
  df.casecontrol$reference_date <- as.Date(as.character(df.casecontrol$reference_date),format="%Y-%m-%d")

  # merge it with the cases
  df.casecontrol <- merge(df.casecontrol,all_event_dt.Include_in_cases.summary,by=c("f.eid","reference_date"),all=T)



  if(!is.null(all_event_dt.Exclude_from_controls)) {
    # !is.na(df.casecontrol$Any) for cases i.e. from all_event_dt.Include_in_cases.summary
    # here to exclude = 1)not a case  2) in exclude_from_control
    exclude=(is.na(df.casecontrol$Any) & (df.casecontrol$f.eid %in% all_event_dt.Exclude_from_controls$f.eid))
    message(glue::glue("excluding {sum(exclude)} controls"))
    # cols to be set to na
    set_to_na <- names(df.casecontrol)[!names(df.casecontrol) %in% c("f.eid","reference_date")]
    df.casecontrol[exclude,(set_to_na):=-1] #mark the non-control as -1
  }

  # 2 => case -2 => non-case  -1 => non-control , set control to 1
  df.casecontrol[is.na(df.casecontrol$Any),]$Any <- 1
  # as control has no event the event count is NA after merge, set it to 1
  df.casecontrol[is.na(df.casecontrol$count),]$count <- 1
  df.casecontrol[(is.na(df.casecontrol$Hx)&(df.casecontrol$Any == 1)),]$Hx <- 1
  df.casecontrol[(is.na(df.casecontrol$Fu)&(df.casecontrol$Any == 1)),]$Fu <- 1
  df.casecontrol[(is.na(df.casecontrol$Ref)&(df.casecontrol$Any == 1)),]$Ref <- 1
  # At this point NA is only possible in Hx/Fu/Ref for case_included i.e. Any==2
  # first consider Hx, cases either have event in past ==2 / future ==NA / unknown==NA [no ref date]
  # case that don't have event in the past but in future
  df.casecontrol[(is.na(df.casecontrol$Hx)&(df.casecontrol$Any == 2)&is.na(df.casecontrol$Hx_days)& (!is.na(df.casecontrol$Fu_days))),]$Hx <-  1
  # similar for Ref , cases either have event within window ==2 / no ==NA / unknown==NA [no ref date]
  # case that don't have event within window i.e. eventdate is known
  df.casecontrol[(is.na(df.casecontrol$Ref)&(df.casecontrol$Any == 2)&((!is.na(df.casecontrol$Hx_days)) | (!is.na(df.casecontrol$Fu_days)))),]$Ref <-  1
  # for Fu , cases either have event in past ==NA / future ==2/ unknown==NA [no ref date]
  # ******************but this is also affected by the sources_recurrence_events option in get_incidence_rate()?***************
  df.casecontrol[(is.na(df.casecontrol$Fu)&(df.casecontrol$Any == 2)& (!is.na(df.casecontrol$Hx_days))& (is.na(df.casecontrol$Fu_days))),]$Fu <-  1
  ############################################################################################################


  print(table(df.casecontrol$Any))

  return(list(df.casecontrol=df.casecontrol,
              all_event_dt.Include_in_cases=all_event_dt.Include_in_cases,
              all_event_dt.Include_in_cases.summary=all_event_dt.Include_in_cases.summary)
  )
}

#' Get survival data
#'
#' Given a phenotype and a list of episode data , compute the survival time after first diagnoisis in cases only
#' @param def phenotype/trait specified in definition table (a row in the table)
#' @param lst.data list of data table with all episode data
#' @param lst.data.settings data frame containing data settings
#' @param window_days_mask number of days that future events should not be counted since first diagnosis.Relevant when the death record is the only record i.e. t=0

#' @return  a data table summarizing time to event information for these individuals.
#' @keywords time-to-event
#' @export
#' @examples
#' get_survival_data(dfDefinitions_processed_expanded[14,],lst.data,lst.data.settings)
get_survival_data<-function(def,lst.data,
                             lst.data.settings,
                             include_secondary_recurrence=FALSE,
                             window_days_mask=0){
  # subset lst.data to get only death records
  # lst.data.death<-lst.data[grep("death", names(lst.data))]
  lst.data.death<-lst.data[lst.data.settings[match(names(lst.data),lst.data.settings$datasource),]$death]
  death_event_dt<-get_all_events(def,lst.data.death,lst.data.settings)

  # check consistency in the records w.r.t date
  discrepant_deaths<-death_event_dt%>%  group_by(f.eid) %>% summarize(max_date = max(eventdate, na.rm = TRUE),min_date = min(eventdate),same_date=(max_date==min_date))%>% filter(!same_date)

  if (nrow(discrepant_deaths)>0){
    message(glue::glue("Number of individuals have inconsistent death dates: {nrow(discrepant_deaths))}"))
    print(discrepant_deaths)

    # drop these individuals?
    lst.data.death$tte.death.icd10.primary<-lst.data.death$tte.death.icd10.primary[! lst.data.death$tte.death.icd10.primary$f.eid %in% discrepant_deaths$f.eid,]
    lst.data.death$tte.death.icd10.secondary<-lst.data.death$tte.death.icd10.secondary[! lst.data.death$tte.death.icd10.secondary$f.eid %in% discrepant_deaths$f.eid,]
  }

  # reference date is the first event date excluding the death records
  # SHOULD one take everything instead because that does reflect a survival event t=0 day....?
  # lst.data_nondeath<-lst.data[!lst.data.settings[match(names(lst.data),lst.data.settings$datasource),]$death]
  # all_evt_dt <- get_all_events(def,lst.data_nondeath,lst.data.settings)
  #
  all_evt_dt <- get_all_events(def,lst.data,lst.data.settings)

  df_referencedate <- all_evt_dt[,.(reference_date= min(eventdate,na.rm = T)),by=f.eid]

  death_event_dt.summary<-get_incidence_prevalence(death_event_dt,lst.data.settings,reference_date=setNames(df_referencedate$reference_date,df_referencedate$f.eid), include_secondary_recurrence,0, window_days_mask)
  # fu_days_mask masks FU event but not day==0 which is considered as Hx so apply the mask on the df

  # make the table cleaner
  death_event_dt.summary<-death_event_dt.summary %>% select(f.eid,first_diagnosis_days,reference_date,Any )
  names(death_event_dt.summary) <-  c("f.eid",	"days_after_diagnosis",	"reference_date",	"Any")

  death_event_dt.summary[death_event_dt.summary$days_after_diagnosis<window_days_mask ,"days_after_diagnosis"] <- NA

  return(death_event_dt.summary)

}
