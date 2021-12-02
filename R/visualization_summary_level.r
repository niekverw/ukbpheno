

#' return random colors of specified length
#'
#' Given a phenotype, a list of episode data and reference dates per individual, identify cases and control status as of *reference date" by source.
#' @param vec_leng phenotype/trait specified in definition table (a row in the table)
#' @param alpha list of data table with all episode data
#' @param cust_seed data frame containing data settings
#' @return  vector of rgb colour codes
#' @export
rand_col_vec<- function(vec_leng=5,alpha=0.75,cust_seed=NULL){
  if (! is.null(cust_seed)){
    set.seed(cust_seed)
  }
  # return(rgb(rgb_int[1],rgb_int[2],rgb_int[3],alpha))
  col_vec<-c()
  for (i in 1:vec_leng){

    rgb_int<-sample( 0:255 , 3 , replace=T)/255
    col_vec<-c(col_vec,do.call(rgb, as.list(c(rgb_int,alpha))))
    if (! is.null(cust_seed)){
      cust_seed<-cust_seed+2222*i
      set.seed(cust_seed)
    }
  }
  return(col_vec)

}


#############################################################
# data source trajectory over time
#############################################################
# timeline by eventdate
plot_disease_timeline_by_source <- function(definition,
                                            lst.data,
                                            df.data.settings,
                                            vct.identifiers
) {
  if(nrow(definition)==0){
    message("No definition is provided.Stop.")
    return(0)
  }
  if (!is.null(vct.identifiers)){
    lst.case_control <- get_cases_controls(definition, lst.data,df.data.settings, vct.identifiers=vct.identifiers,verbose=FALSE)
  }else{
    message("Both df.reference.dates and vct.identifiers are NULL, Please provide one of them.Exit.")
    return()
  }
  # for annotation , rounding to month
  dt_v0_min<-as.Date("2006-03-01")
  dt_v0_max<-as.Date("2010-10-01")
  dt_v1_min<-as.Date("2009-12-01")
  dt_v1_max<-as.Date("2013-06-01" )
  dt_v2_min<-as.Date("2014-05-01")
  dt_v2_max<-as.Date("2020-03-01")
  dt_v3_min<-as.Date("2019-06-01")
  dt_v3_max<-as.Date("2021-06-01")
  # get diease events
  lst.case_control <- get_cases_controls(definitions=definition, lst.data=lst.data,df.data.settings=df.data.settings,  vct.identifiers=vct.identifiers)
  all_event_dt<-lst.case_control$all_event_dt.Include_in_cases
  # # individuals without eventdates
  n_no_dt<-length(unique(all_event_dt[event==0]$identifier))


  # remove non-event  NOTE actual events without event date will be removed as well
  all_event_dt<-all_event_dt[event!=0]
  # get first date for each individual by source
  mindate<-all_event_dt[, .(mindt= min(eventdate,na.rm = T)), by = c('identifier','.id')]
  # to stitch the lines to final point
  max_dt<-max(mindate$mindt)
  # get a cumulative count with rank
  mindate[, cumcnt:= data.table::frank(mindt,ties.method='first'),by=.id]
  # create rows for stitching the line to final time
  dt_data_end<-mindate[, .(identifier=9999999,mindt=max_dt,cumcnt= max(cumcnt,na.rm=T)),by=.id
  ][, tail(.SD, 1) , by = .id  ]
  # rbind to actual dt
  mindate<-rbind(mindate,dt_data_end)
  # create a line of Any
  any_mindate<-all_event_dt[, .(mindt= min(eventdate,na.rm = T)), by = c('identifier')]
  any_mindate[, `:=`(.id="Any",cumcnt= data.table::frank(mindt,ties.method='first'))]
  # add Any
  mindate<-rbind(mindate,any_mindate)

  dt_data_end<-rbind(dt_data_end,any_mindate[any_mindate$cumcnt==max(any_mindate$cumcnt),])
  color_vec<-rand_col_vec(length(unique(mindate$.id)))
  # color_vec<-ggpubr::get_palette("rickandmorty",length(unique(mindate$.id)))
  plt_time<- ggplot2::ggplot( mindate,aes(x=mindt, y=cumcnt,colour=.id,group=.id)) +
    ggplot2::scale_color_manual(values = color_vec,name="source") +
    ggplot2::geom_line(alpha=1,size=1)+ggplot2::geom_point(size=0.8,alpha=0.1,shape=20) +ggplot2::xlab("Time")+
    ggplot2::ylab("Count")+
    # annotation window v0
    ggplot2::annotate("rect", xmin=dt_v0_min, xmax=dt_v0_max, ymin=0, ymax=Inf,alpha = .2,fill="#6395b6")+
    ggplot2::annotate(geom="text", x=as.Date(mean.Date(c(dt_v0_min,dt_v0_max))), y=max(dt_data_end$cumcnt)*1.12,label="Baseline visits",color="#213745") +
    # # annotation window v1
    ggplot2::annotate("rect", xmin=dt_v1_min, xmax=dt_v1_max, ymin=0, ymax=Inf,alpha = .2,fill="#f4b699")+
    ggplot2::annotate(geom="text", x=as.Date(mean.Date(c(dt_v1_min,dt_v1_max))), y=max(dt_data_end$cumcnt)*1.07,label="1st return",color="#7d300d") +
    # # annotation window v0
    ggplot2::annotate("rect", xmin=dt_v2_min, xmax=dt_v2_max, ymin=0, ymax=Inf,alpha = .2,fill="#d2baf8")+
    ggplot2::annotate(geom="text", x=as.Date(mean.Date(c(dt_v2_min,dt_v2_max))), y=max(dt_data_end$cumcnt)*1.1,label="Imaging visits",color="#2d0b66") +
    # # annotation window v0
    ggplot2::annotate("rect", xmin=dt_v3_min, xmax=dt_v3_max, ymin=0, ymax=Inf,alpha = .2,fill="#ebf179")+
    ggplot2::annotate(geom="text", x=dt_v3_min, y=max(dt_data_end$cumcnt)*1.05,label="1st return imaging",color="#777d0d") +
    # # TODO warn missing
    # ppl_miss<-paste0(n_no_dt, " cases without date of diagnosis not plotted")
    # ggplot2::annotate(geom="text", x=dt_v3_min, y=max(dt_data_end$cumcnt)*1.05,label="1st return imaging",color="#777d0d") +
    # # log2 scale not nice
    # scale_y_continuous(trans='log10')+
    ggpubr::theme_classic2(base_size = 12)+
    ggrepel::geom_text_repel( aes(label = cumcnt), data = dt_data_end,  fontface ="plain", nudge_y=max(dt_data_end$cumcnt)/50 , size = 4)

  # +ggplot2::geom_rect(aes(xmin=dt_v0_max, xmax=dt_v0_max, ymin=0, ymax=Inf),  alpha = .2,fill="#BDC7FA")

  return(plt_time)
}





#' Get case-control status by data sources
#'
#' Given a phenotype, a list of episode data and reference dates per individual, identify cases and control status as of *reference date" by source.
#' @param definitions phenotype/trait specified in definition table (a row in the table)
#' @param lst.data list of data table with all episode data
#' @param lst.data.settings data frame containing data settings
#' @param reference_date reference dates for each individuals in the whole cohort as a named vector

#' @return  a data table : case control status according to different sources and one column based on any of the sources.
#' @keywords time-to-event
#' @export
#' @examples
#' get_case_count_by_source(cancer_source,definition=dfDefinitions_processed_expanded %>% filter(TRAIT=="Nicm"), lst.data,lst.data.settings,lst.identifiers=dfukb$f.eid))
get_case_count_by_source <- function(definition,
                      lst.data,
                      df.data.settings,
                      df.reference.dates=NULL,
                      vct.identifiers=NULL,
                      standardize=TRUE
                      ) {
  if(nrow(definition)==0){
    message("No definition is provided.Stop.")
    return(0)
  }
  if (!is.null(df.reference.dates)){
    lst.case_control <- get_cases_controls(definition, lst.data,dfData.settings, df_reference_date=df.reference.dates,verbose=FALSE)
  }else if (!is.null(vct.identifiers)){
    lst.case_control <- get_cases_controls(definition, lst.data,dfData.settings, vct.identifiers=vct.identifiers,verbose=FALSE)
  }else{
      message("Both df.reference.dates and vct.identifiers are NULL, Please provide one of them.Exit.")
    return()
  }
  ###########################################################################
  # NOTE: this function depends on the get_cases_controls! modify accordingly
  ##########################################################################
  # ##########################################################################
  # exclude future events
  # otherwise the comparison between sources is not fair
  # ##########################################################################
  cases<- lst.case_control$df.casecontrol[lst.case_control$df.casecontrol$Any==2,]
  # only those ppl with only future events are excluded, because these ppl would be considered as control at the reference date
  ppl_future_only<-cases[cases$Fu==2 & cases$Any==2 & cases$Hx!=2,]$identifier
  if(nrow(definition)>1){
    # take first row for the message
    definition<-definition[1]
  }
  definition<-
  message(glue::glue("{definition$DESCRIPTION}: {length(ppl_future_only)} indivduals have events after reference dates and are not considered"))
  # discard future events
  lst.case_control$all_event_dt.Include_in_cases<-lst.case_control$all_event_dt.Include_in_cases[ (! lst.case_control$all_event_dt.Include_in_cases$identifier %in% ppl_future_only)]
  # get all available data sources from data
  all_sources<-unique(lst.case_control$all_event_dt.Include_in_cases$.id)
  # create new columns for each source
  for (source_name in all_sources ) {
  varname<-paste('Hx',source_name,sep='_')
  # create new columns for each source
  cases[, (varname)]<-as.numeric(NA)

  # lookup source from df with all episodes i.e. all_event_dt.Include_in_cases
  # if rows in include_in_case which originated from the target source, look up the eid and set to 2
  cases[[varname]][cases$identifier %in% lst.case_control$all_event_dt.Include_in_cases[lst.case_control$all_event_dt.Include_in_cases$.id %in% source_name,]$identifier] <-2
  }

  for (j in  paste0("Hx_",all_sources)){
    # for all new Hx_ columns , set everything to 0 for ease of counting
    # if Hx <0 , (excluded cases)
    set( cases,which( cases$Hx<=0),j,0)
    #  if the cell value was NA ->  control
    set( cases,which(is.na(cases[[j]])),j,0)
  }
  cases$Any[!cases$identifier %in% lst.case_control$all_event_dt.Include_in_cases$identifier] <-0

  cols<-c("Any",paste0("Hx_",all_sources))
  # alternatve to the block below
  # df_prop<-cases[ ,cols,with=FALSE]%>% dplyr::summarise(across(where(is.numeric),sum))
  # df_prop<-df_prop/2
  # df_prop<-data.table(t(df_prop),keep.rownames = TRUE)
  # colnames(df_prop)<-c("data_source","proportion")

  df_prop<-data.table(unlist(lapply(cases[ ,(cols)], function(y){
  v1 <- nrow(cases[cases[[y]]==2,])
  })))
  if (standardize){
      df_prop$V1<-df_prop$V1/max(df_prop$V1)
  colnames(df_prop)<-"proportion"
  }else{
      colnames(df_prop)<-"n"
  }
  df_prop$source<-c("Any",all_sources)
  # sort
  df_prop<-df_prop[order(df_prop$source),]
  return(df_prop)
}




#############################################################
# data source trajectory over varying reference dates
#############################################################



plot_source_proportion_over_time<-function(definition,lst.dfReferenceDate,lst.data,
                                           df.data.settings,standardize=TRUE){
  time.lab<-names(lst.dfReferenceDate)
  # lst.dfprop<-list()
  lst.dfprop<-lapply(lst.dfReferenceDate, function(x){
    df_prof<-get_case_count_by_source(definition=definition,lst.data,df.data.settings,df.reference.dates=x,standardize =standardize)
    df_prof
  })
  # print(lst.dfprop)
  # print(time.lab)
  # rename for merging
    lst.dfprop<-lapply(time.lab,function(x){
        # print(names(lst.dfprop[[x]]))
        # print(x)
        if(standardize){
        names(lst.dfprop[[x]])[names(lst.dfprop[[x]])=='proportion']<-x
      }else{
        names(lst.dfprop[[x]])[names(lst.dfprop[[x]])=='n']<-x
      }
    lst.dfprop[[x]]
  })

 # merge the data tables
  prop_df_time <- Reduce(
    function(x, y, ...) merge(x, y, by="source",all = TRUE, ...),
    lst.dfprop
  )
  color_vec<-rand_col_vec(nrow(prop_df_time))
  # long table for plotting
  prop_df_time<-reshape2::melt(prop_df_time)


  plt_source_over_time<-ggplot2::ggplot(prop_df_time, ggplot2::aes(x = variable, y = value, colour = source,group=source)) +
    ggplot2::scale_color_manual(values = color_vec) +
    ggplot2::geom_line()  + ggplot2::scale_y_continuous(breaks = pretty(1:max(prop_df_time$value),n=5))+
    ggplot2::geom_point() +ggplot2::xlab("Time point")+ggpubr::theme_pubclean(base_size = 12)

  if(standardize){
    plt_source_over_time<-plt_source_over_time + ggplot2::labs(y = "Proportion")
  }else{
    plt_source_over_time<-plt_source_over_time + ggplot2::labs(y = "Count")
  }

   return(plt_source_over_time)
}

# TODO dotplot !
