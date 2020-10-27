

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
  p1 <-ggplot2::ggplot(stats.codes, ggplot2::aes(rank, count,label=code,color=classification)) + ggplot2::geom_point() + ggplot2::ylim(-((max(stats.codes$count))/3),NA)  + ggrepel::geom_text_repel(size =3,segment.size=0.5)

  p2 <- ggplot2::ggplot(stats.codes,ggplot2:: aes(rank, count,label=code,color=classification)) + ggplot2::scale_y_continuous(trans='log2') + ggplot2::geom_point()   + ggrepel::geom_text_repel(size =3,segment.size=0.5)
  stats.codes.summary.table <- stats.codes
  stats.codes.summary.p <- ggpubr::ggarrange(p1,p2,nrow = 1, ncol = 2,common.legend = TRUE)

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
  stats.class.cooccur.p <- ggplot2::ggplot(longData, ggplot2::aes(x = Code_occur , y =Code_presence,fill=`%`)) +
    ggplot2::geom_tile()+
    ggplot2::scale_fill_distiller(palette = "RdYlBu",na.value = "grey85") +
    ggplot2::labs(x="Classfication of co-occurence", y="Classification present",title="Co-occurence of diagnosis by sources",caption = "[In the presence of row, column coexists with row by %]") +
    ggplot2::geom_text(ggplot2::aes(label = `%`))
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
  stats.codes.cooccur.filtered.table <- stats.codes.cooccur.table[matrixStats::rowMaxs(stats.codes.cooccur.table,na.rm=T)>10,matrixStats::colMaxs(stats.codes.cooccur.table,na.rm=T)>10]
  # stats.codes.cooccur.filtered.p <- pheatmap::pheatmap( stats.codes.cooccur.filtered.table  ,fontsize = 6)
  # ########## ggplot implementation#################################################################################
  # Create ggplot version dendrogram from ggdendro

   code.dendro <- as.dendrogram(hclust(d = dist(x = stats.codes.cooccur.filtered.table)))
  ddata_x <-  ggdendro::dendro_data(code.dendro)
  # to colour leaves by classifications

  lab_gp <- ggdendro::label(ddata_x)
  lab_gp$group <- stats.codes[match(lab_gp$label,stats.codes$code),]$classification

  stats.codes.cooccur.filtered.p.dendro <- ggplot2::ggplot(ggdendro::segment(ddata_x)) +
    ggplot2::geom_segment(ggplot2::aes(x=x, y=y+10, xend=xend, yend=yend+10)) +  ggplot2::geom_text(data=ggdendro::label(ddata_x),
    ggplot2::aes(label=label, x=x, y=-5, colour=lab_gp$group),size =3,angle=45) +ggplot2::theme(axis.line=ggplot2::element_blank(),
         axis.text.x=ggplot2::element_blank(),
         axis.text.y=ggplot2::element_blank(),
         axis.ticks=ggplot2::element_blank(),
         axis.title.x=ggplot2::element_blank(),
         axis.title.y=ggplot2::element_blank(),
           # legend.position="none",
         panel.background=ggplot2::element_blank(),
         panel.border=ggplot2::element_blank(),
         panel.grid.major=ggplot2::element_blank(),
         panel.grid.minor=ggplot2::element_blank(),
         plot.background=ggplot2::element_blank())
  # ggplot version heatmap
  # code.order <- order.dendrogram(code.dendro)
  # TODO maker heatmap ordered like dendrogram?

  longData<- reshape2::melt(stats.codes.cooccur.filtered.table)
  names(longData)<- c("Code_presence","Code_occur","%")
  stats.codes.cooccur.filtered.p.heat <- ggplot2::ggplot(longData, ggplot2::aes(x = Code_occur , y =Code_presence,fill=`%`)) +
    ggplot2::geom_tile()+
    ggplot2::scale_fill_distiller(palette = "RdYlBu",na.value = "grey85") +
    ggplot2::labs(x="Code of co-occurence", y="Code present",title="Co-occurence of diagnosis code",caption = "[In the presence of row, column coexists with row by %]") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45))
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
    df_referencedate <-data.table::data.table(reference_date)
    df_referencedate$f.eid <- names(reference_date)

    message(glue::glue("non missing reference_date: {length(reference_date)}"))
  } else{
    message("no reference_date given, taking the first available event as reference")
    df_referencedate <- all_event_dt[,.(reference_date= min(eventdate,na.rm = T)),by=f.eid]
  }
  if(include_secondary_recurrence){
    sources_recurrence_events <- lst.data.settings$datasource
  } else {
    sources_recurrence_events <- lst.data.settings %>%  dplyr::filter(diagnosis==1) %>% dplyr::pull(datasource)
  }

  ###########################################################
  # in case of empty rows , the entire column will be cast to default of class NA i.e. logical
  # this could create error for downstream functions "Error in bmerge....Incompatible join types"
  col_num<-c("count","sum.epidur","median.epidur","max.epidur", "survival_days", "death.primary","death.any" , "Hx_days", "Fu_days" ,"Hx" , "Fu","Ref" , "first_diagnosis_days","Any")
  col_chr<-"f.eid"
  col_date<-"reference_date"
  # catch if there is no event, return an empty table which can be merged
  if (nrow(all_event_dt)==0){
    message("No event found.")
    results_cols<-c("f.eid","count","sum.epidur","median.epidur","max.epidur", "survival_days", "death.primary","death.any" , "Hx_days", "Fu_days" ,"Hx" , "Fu","Ref" , "first_diagnosis_days", "reference_date","Any")
    all_event_dt.summary <- setNames(data.table(matrix(nrow = 0, ncol = 16)),results_cols )
    all_event_dt.summary[, (col_num) := lapply(.SD, as.numeric), .SDcols = col_num]
    all_event_dt.summary[, (col_chr) := lapply(.SD, as.character), .SDcols = col_chr]
    all_event_dt.summary[, (col_date) := lapply(.SD, as.Date), .SDcols = col_date]
    return(all_event_dt.summary)
  }


  df <- merge(all_event_dt,df_referencedate,by = 'f.eid') %>% dplyr::arrange(eventdate) %>% data.table::as.data.table()

  # df <- df %>% filter(!is.na(reference_date)) # comment out if missing f.eids is fixed.
  df$days <- df$eventdate - df$reference_date

  data.table::setkey(df,days) # i don't know why, but setkey was alreaday on f.eid and cannot refresh..
  data.table::setkey(df,f.eid)


  #################################################################################
  # death
  dfDth<-df[df$.id %in% lst.data.settings[lst.data.settings$death,]$datasource,]

  if (nrow(dfDth) >0){
  # get death records
  # dfDth<-df[(df$.id %in% lst.data.settings[lst.data.settings$death,]$datasource),]
  # ### flag primary death records
  # dfDth$death.primary<-ifelse((lst.data.settings[match(dfDth$.id ,lst.data.settings$datasource),]$diagnosis==1),2,NA)
  # ### flag secondary death records
  # dfDth$death.secondary<-ifelse((lst.data.settings[match(dfDth$.id ,lst.data.settings$datasource),]$diagnosis==2),2,NA)
  ### this seems faster
  dfDth$death.primary <- NA
  dfDth$death.primary[lst.data.settings[match(dfDth$.id ,lst.data.settings$datasource),]$diagnosis==1] <- 2

  dfDth$death.secondary<- NA
  dfDth$death.secondary[lst.data.settings[match(dfDth$.id ,lst.data.settings$datasource),]$diagnosis==2] <- 2
  dfDth$death.secondary<-data.table::fcoalesce(dfDth$death.primary,dfDth$death.secondary)
  # dfDth
  dfDth<-dfDth[,c("f.eid","days","death.primary","death.secondary")]
  colnames(dfDth)<-c("f.eid","survival_days","death.primary","death.any")
  # in the case of duplicate death records,one per primary /secondary
  dfDth_extrastats <- suppressWarnings(dfDth[, .(mindy= min(survival_days,na.rm = T),maxdy= max(survival_days,na.rm = T),meandy= mean(survival_days,na.rm=T) ), keyby=list(f.eid)])
  id.diff.deathdt<-unique(dfDth_extrastats[dfDth_extrastats$mindy!=dfDth_extrastats$meandy|dfDth_extrastats$maxdy!=dfDth_extrastats$meandy,]$f.eid)

  if (length(id.diff.deathdt)>0){
  message(glue::glue("Inconsistent death records for these individuals: {glue::glue_collapse(id.diff.deathdt,sep = ',')}"))
  message("Mean survival day will be taken.")
  }
  # in case of multiple death records with different death dates , take mean
  dfDth<-stats::aggregate(x=dfDth[,!(names(dfDth) %in% c("f.eid")),with=FALSE], by=list(f.eid=dfDth$f.eid), mean, na.rm = TRUE)
  # Nan to NA
  dfDth[is.na(dfDth)]<-NA

  }else{
    # no records, take f.eid for the merge later
    dfDth<-unique(df[,"f.eid"])
    dfDth[,c("survival_days","death.primary","death.any")]<-as.numeric(NA)
  }
  ###############################################################################################################################

  ### History
  dfHx <- df[days<=0]
  Hx_days <- suppressWarnings(unique(dfHx[event>0][, .(Hx_days=min(days,na.rm=T)), keyby=list(f.eid)] ))
  dfHx[,Hx:=2]

  ### Future
  # window_fu_days_mask
  dfFu <- df[  ((days>(0+window_fu_days_mask) & event==1) & (!f.eid %in% dfHx$Hx)) |
               (days>(0+window_fu_days_mask) & event==1 & (f.eid %in% dfHx$Hx) & .id %in% sources_recurrence_events ) |
               (days>(0+window_fu_days_mask) & event==2  & (!f.eid %in% dfHx$Hx)) ]
  dfFu[,Fu:=2] #unique(dfFu$f.eid)
  Fu_days <- suppressWarnings( unique(dfFu[,.(Fu_days= min(days,na.rm=T)), keyby=list(f.eid)] ))
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

  # stats
  # Hx_days
  # Fu_days
  # dfDth
  # dfHx[,c("f.eid","Hx")]
  # dfHx
  # dfRef[,c("f.eid","Ref")]

  all_event_dt.summary <- Reduce(function(...) merge(..., all = TRUE,by='f.eid'), list(
      stats,
      unique(dfDth[ , c("f.eid","survival_days","death.primary","death.any")]),
      Hx_days,
      Fu_days,
      unique(dfHx[,c("f.eid","Hx")]),
      unique(dfFu[,c("f.eid","Fu")]),
      unique(dfRef[,c("f.eid","Ref")])
  ))


  all_event_dt.summary$first_diagnosis_days <- pmin(all_event_dt.summary$Hx_days,all_event_dt.summary$Fu_days,na.rm = T)
  all_event_dt.summary[ Hx==2 & is.na(Hx_days),'first_diagnosis_days'] <- NA

  all_event_dt.summary <- merge(all_event_dt.summary,df_referencedate,by="f.eid")
  all_event_dt.summary[,Any:=2]
  all_event_dt.summary <- data.table::data.table(all_event_dt.summary)

  ####################################type check######################################################
  all_event_dt.summary[, (col_num) := lapply(.SD, as.numeric), .SDcols = col_num]
  all_event_dt.summary[, (col_chr) := lapply(.SD, as.character), .SDcols = col_chr]
  all_event_dt.summary[, (col_date) := lapply(.SD, as.Date), .SDcols = col_date]
  #####################################################################################################
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
  if(length(definitions$TRAIT)==0){
    message("No TRAIT in definitions.Stop.")
    return(0)
  }


  all_event_dt.Include_in_cases <- get_all_events(definitions %>% dplyr::filter(Definitions =="Include_in_cases"),lst.data,lst.data.settings)   #MI
  all_event_dt.Include_in_cases.summary <- get_incidence_prevalence(all_event_dt = all_event_dt.Include_in_cases,lst.data.settings,
                                                                    reference_date = reference_date)
                                                                    #...)

  message(glue::glue("including {nrow(all_event_dt.Include_in_cases.summary)} cases"))
  all_event_dt.Exclude_from_cases <- get_all_events(definitions %>% dplyr::filter(Definitions =="Exclude_from_cases"),lst.data,lst.data.settings)   #MI
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

  if(length(definitions$TRAIT)==0){
    message("No TRAIT in definitions.Stop.")
    return(0)
  }

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
  all_event_dt.Exclude_from_controls <- get_all_events(definitions %>% dplyr::filter(Definitions =="Exclude_from_controls"),lst.data,lst.data.settings)
  ### define case & control
  if(is.null(reference_date)){
    reference_date = setNames(as.Date(rep(NA,length(lst.identifiers))),lst.identifiers)
  }

  # NOTE that if the rownames are not unique it will be discarded at data.frame() i.e. f.eid replaced by the running row number
  df.casecontrol <- data.frame(reference_date=reference_date) %>% tibble::rownames_to_column('f.eid') %>% data.table::as.data.table()
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
  discrepant_deaths<-death_event_dt%>%  dplyr::group_by(f.eid) %>% dplyr::summarize(max_date = max(eventdate, na.rm = TRUE),min_date = min(eventdate),same_date=(max_date==min_date))%>% dplyr::filter(!same_date)

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
  death_event_dt.summary<-death_event_dt.summary %>% dplyr::select(f.eid,first_diagnosis_days,reference_date,Any )
  names(death_event_dt.summary) <-  c("f.eid",	"days_after_diagnosis",	"reference_date",	"Any")

  death_event_dt.summary[death_event_dt.summary$days_after_diagnosis<window_days_mask ,"days_after_diagnosis"] <- NA

  return(death_event_dt.summary)

}
