# library(ggplot2)
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
      cust_seed<-cust_seed+13*i
      set.seed(cust_seed)
    }
  }
  return(col_vec)

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
                      df_reference_dates=NULL,
                      vct_identifiers=NULL,
                      standardize=TRUE
                      ) {
  if(nrow(definition)==0){
    message("No definition is provided.Stop.")
    return(0)
  }
  if (!is.null(df_reference_dates)){
    lst.case_control <- get_cases_controls(definition, lst.data,dfData.settings, df_reference_date=df_reference_dates,verbose=FALSE)
  }else if (!is.null(vct_identifiers)){
    lst.case_control <- get_cases_controls(definition, lst.data,dfData.settings, vct.identifiers=vct_identifiers,verbose=FALSE)
  }else{
      message("Both df_reference_dates and vct_identifiers are NULL, Please provide one of them.Exit.")
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
#
# prop_af_v0<-get_case_count_by_source(definition=dfDefinitions_processed_expanded %>% filter(TRAIT==trait),lst.data,dfData.settings,df_reference_dates=df_reference_dt_v0,standardize = TRUE)
# prop_af_v2<-get_case_count_by_source(definition=dfDefinitions_processed_expanded %>% filter(TRAIT==trait),lst.data,dfData.settings,df_reference_dates=df_reference_dt_v2,standardize = TRUE)
#
#
#
# #############################################################
# # data source trajectory over time
# #############################################################
# # TODO:
# plot_source_proportion_over_time<-function(...){
#   prop_af_time <- Reduce(
#     function(x, y, ...) merge(x, y, by="source",all = TRUE, ...),
#     list(prop_af_v0,prop_af_v2,prop_af_today)
#   )
#
#
#   colnames(prop_af_time)<-c("source","baseline","v2","25-11-2021")
#   prop_af_time<-reshape2::melt(prop_af_time)
#
#   ggplot(prop_af_time, aes(x = variable, y = value, colour = source,group=source)) +
#     geom_line() +
#     geom_point()
#
#
# }
#
#
#
#
# #############################################################
# # radar plot  by data source
# #############################################################
# # TODO 1) a function to make the table 2) radar function taking 1) putput as input
# test_af<-get_case_count_by_source(definition=dfDefinitions_processed_expanded %>% filter(TRAIT==trait),lst.data,dfData.settings,vct_identifiers = dfukb$identifier)
#
#
# df_reference_dt_today<-df_reference_dt_v0
# df_reference_dt_today$f.53.0.0<-as.Date(as.character("2021-11-25"),format="%Y-%m-%d")
# df_reference_dt_today$f.53.0.0<-as.Date(df_reference_dt_today$f.53.0.0,format="%Y-%m-%d")
# as.Date(as.character(dfukb$f.53.0.0),format="%Y-%m-%d")
# rm(test_cad)
#
#
# prop_af_today<-get_case_count_by_source(definition=dfDefinitions_processed_expanded %>% filter(TRAIT==trait),lst.data,dfData.settings,df_reference_dates=df_reference_dt_today,standardize = TRUE)
# prop_cad_today<-get_case_count_by_source(definition=dfDefinitions_processed_expanded %>% filter(TRAIT=="Cad"),lst.data,dfData.settings,df_reference_dates=df_reference_dt_today,standardize = TRUE)
# prop_hf_today<-get_case_count_by_source(definition=dfDefinitions_processed_expanded %>% filter(TRAIT=="Hf"),lst.data,dfData.settings,df_reference_dates=df_reference_dt_today,standardize = TRUE)
# prop_hcm_today<-get_case_count_by_source(definition=dfDefinitions_processed_expanded %>% filter(TRAIT=="Hcm"),lst.data,dfData.settings,df_reference_dates=df_reference_dt_today,standardize = TRUE)
# prop_dcm_today<-get_case_count_by_source(definition=dfDefinitions_processed_expanded %>% filter(TRAIT=="dcm"),lst.data,dfData.settings,df_reference_dates=df_reference_dt_today,standardize = TRUE)
#
#
#
# make_radar_plot(df_prev,sort(cardiovascular_traits[[1]]),sources_vec =c("All","SR","GP","HESIN","Death") )
#
#   lst.proportion<-list()
#   for (trait in vct_diseases){
#     dt_disease_prop<-get_case_count_by_source(definition=dfDefinitions_processed_expanded %>% filter(TRAIT==trait))
#     lst.codemap[[cls]]<-fread(fmap)
#
#   }
#
#
#
# make_radar_plot<-function(df_prevalence,vct_diseases,sources_vec=c("All","SR","Cancer","HESIN","Death")){
#
#
#   prev_subset<-df_prevalence[,(names(df_prevalence) %in% vct_diseases)]
#   # prev_subset<-df_prevalence[,vct_diseases]
#
#   plot_data <- rbind(rep(max(prev_subset),ncol(prev_subset)) , rep(min(prev_subset),ncol(prev_subset)) , prev_subset[sources_vec,])
#
#   plot_data<-plot_data[,sort(colnames(prev_subset))]
#
#   disease_names<-dfDefinitions%>% filter(TRAIT %in% names(prev_subset))%>%select(TRAIT,DESCRIPTION) %>% arrange(TRAIT)
#
#   colnames(plot_data)<-disease_names$DESCRIPTION
#   colors_border=rand_col_vec(ncol(prev_subset),0.7,cust_seed=34)
#   colors_in=rand_col_vec(ncol(prev_subset),0.2,cust_seed=34)
#   radarchart( plot_data  , axistype=1 ,
#               #custom polygon
#               pcol=colors_border , pfcol=colors_in , plwd=2 , plty=1,
#               #custom the grid
#               cglcol="grey", cglty=1, axislabcol="grey", seg=4,caxislabels=round(seq(min(prev_subset),max(prev_subset),length.out=5),1), cglwd=1,
#               #custom labels
#               vlcex=1)
#   # Add a legend
#   legend(x=1.5, y=0.8, legend = rownames(plot_data[-c(1,2),]), bty = "n", pch=20 , col=colors_border , text.col = "black", cex=1.2, pt.cex=3)
# }
#
#
#
#
#
# # head(case_status)
# #
# # count(case_status,Any.sr)
# # count(case_status,Any.hes)
#
#
# #TODO
# # a function to show (number of case) unique to each source
#
# #add radarplots

