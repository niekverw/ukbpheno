

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
      cust_seed<-cust_seed+7*i
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
#' gc()
#' #
#' # prop_af_v0<-get_case_count_by_source(definition=dfDefinitions_processed_expanded %>% filter(TRAIT==trait),lst.data,dfData.settings,df.reference.dates=df_reference_dt_v0,standardize = TRUE)
#' # prop_af_v2<-get_case_count_by_source(definition=dfDefinitions_processed_expanded %>% filter(TRAIT==trait),lst.data,dfData.settings,df.reference.dates=df_reference_dt_v2,standardize = TRUE)
#' #
#' #
#' #
#############################################################
# data source trajectory over time
#############################################################

# lst_df_ref_dt<-list(Baseline=df_reference_dt_v0,Visit2=df_reference_dt_v2,Nov21=df_reference_dt_today)
plt_test<-plot_source_proportion_over_time(definition=dfDefinitions_processed_expanded %>% filter(TRAIT==trait),lst_df_ref_dt,lst.data,dfData.settings,FALSE)

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
    ggplot2::geom_point() +ggplot2::xlab("Time point")+ggpubr::theme_classic2(base_size = 12)

  if(standardize){
    plt_source_over_time<-plt_source_over_time + ggplot2::labs(y = "Proportion")
  }else{
    plt_source_over_time<-plt_source_over_time + ggplot2::labs(y = "Count")
  }

   return(plt_source_over_time)
}

#'
#' View(test)
#'
#'
#' #############################################################
#' # radar plot  --> dot plot by data source
#' #############################################################
#' # TODO 1) a function to make the table 2) radar function taking 1) putput as input
#' test_af<-get_case_count_by_source(definition=dfDefinitions_processed_expanded %>% filter(TRAIT==trait),lst.data,dfData.settings,vct.identifiers = dfukb$identifier)
#'
#'
#' df_reference_dt_today<-df_reference_dt_v0
#' df_reference_dt_today$f.53.0.0<-as.Date(as.character("2021-11-25"),format="%Y-%m-%d")
#' df_reference_dt_today$f.53.0.0<-as.Date(df_reference_dt_today$f.53.0.0,format="%Y-%m-%d")
#' as.Date(as.character(dfukb$f.53.0.0),format="%Y-%m-%d")
#' rm(test_cad)
#'
#'
#' dotplot_diseases_by_source<-function(lst.dfprop,disease_lab=NULL){
#'
#'
#' }
#'
#'
#' prop_af_today<-get_case_count_by_source(definition=dfDefinitions_processed_expanded %>% filter(TRAIT==trait),lst.data,dfData.settings,df.reference.dates=df_reference_dt_today,standardize = TRUE)
#' prop_cad_today<-get_case_count_by_source(definition=dfDefinitions_processed_expanded %>% filter(TRAIT=="Cad"),lst.data,dfData.settings,df.reference.dates=df_reference_dt_today,standardize = TRUE)
#' prop_hf_today<-get_case_count_by_source(definition=dfDefinitions_processed_expanded %>% filter(TRAIT=="Hf"),lst.data,dfData.settings,df.reference.dates=df_reference_dt_today,standardize = TRUE)
#' prop_hcm_today<-get_case_count_by_source(definition=dfDefinitions_processed_expanded %>% filter(TRAIT=="Hcm"),lst.data,dfData.settings,df.reference.dates=df_reference_dt_today,standardize = TRUE)
#' prop_dcm_today<-get_case_count_by_source(definition=dfDefinitions_processed_expanded %>% filter(TRAIT=="Dcm"),lst.data,dfData.settings,df.reference.dates=df_reference_dt_today,standardize = TRUE)
#' prop_DysLip_today<-get_case_count_by_source(definition=dfDefinitions_processed_expanded %>% filter(TRAIT=="HyperLip"),lst.data,dfData.settings,df.reference.dates=df_reference_dt_today,standardize = TRUE)
#' prop_Dm_today<-get_case_count_by_source(definition=dfDefinitions_processed_expanded %>% filter(TRAIT=="Dm"),lst.data,dfData.settings,df.reference.dates=df_reference_dt_today,standardize = TRUE)
#' disease_vec<-c('AF','CAD','DCM','HCM','HF','Hyperlipidemia','DM')
#'
#' lst_prop<-list(AF=prop_af_today,CAD=prop_cad_today,DCM=prop_dcm_today,HCM=prop_hcm_today,HF=prop_hf_today,Hyperlipidemia=prop_DysLip_today,DM=prop_Dm_today)
#'   # for (i in 1:length(lst_prop)){
#'   #   names(lst_prop[[i]])[names(lst_prop[[i]])=='proportion']<-disease_vec[[i]]
#'   #   print(lst_prop[[i]])
#'   #   }
#'
#'
#' lst_prop<-lapply(disease_vec,function(x){
#'   # names(x)[names(x)=='proportion']<-names(lst_prop)[names(lst_prop)==x]
#'   # print(colnames(lst_prop[[x]])[colnames(lst_prop[[x]])=='proportion'])
#'   # print(x)
#'   names(lst_prop[[x]])[names(lst_prop[[x]])=='proportion']<-x
#'   lst_prop[[x]]
#'   # names(lst.case_control$df.casecontrol)[names(lst.case_control$df.casecontrol) == paste(trait,visit,"f.eid",  sep = "_")]<-"f_eid"
#'   #
#'
#'   })
#'
#' lst_prop
#'
#'
#' gc()
#' class(names(lst_prop))
#'
#'
#' View(lst_prop$AF)
#' cardio_prop<- Reduce(
#'   function(x, y, ...) merge(x, y, by="source",all = TRUE, ...),
#'   lst_prop)
#'
#' df3<-reshape2::melt(cardio_prop)
#' df3[is.na(df3)]<-0
#' ggpubr::ggdotchart(
#'   df3, x = "variable", y = "value",
#'   group = "source", color = "source", palette = 'rickandmorty', ggtheme = ggpubr::theme_pubclean(base_size = 12),dot.size = 3.5,shape=18 ,                              # Large dot size
#'    size = 1,
#'   add = "segment", position = ggplot2::position_dodge(0.3),
#'   sorting = "descending"
#' )
#' # "rickandmorty"
#' ggpubr::show_point_shapes()
#'
#' df_plot<-cardio_prop[,2:6]
#' df_plot[is.na(df_plot)] <- 0
#' class(df_plot)
#' # colnames(plot_data)<-disease_names$DESCRIPTION
#' colors_border=rand_col_vec(nrow(df_plot),0.7,cust_seed=33)
#' colors_in=rand_col_vec(nrow(df_plot),0.2,cust_seed=34)
#' fmsb::radarchart( df_plot  , axistype=1 ,
#'             #custom polygon
#'             pcol=colors_border , pfcol=colors_in , plwd=2 , plty=1,
#'             #custom the grid
#'             cglcol="grey", cglty=1, axislabcol="grey", seg=4,caxislabels=round(seq(min(df_plot),max(df_plot),length.out=5),1), cglwd=1,
#'             #custom labels
#'             vlcex=1)
#' # Add a legend
#' legend(x=1.5, y=0.8, legend = cardio_prop[,1], bty = "n", pch=20 , col=colors_border , text.col = "black", cex=1.2, pt.cex=3)
#'
#' data <- as.data.frame(matrix( sample( 0:20 , 15 , replace=F) , ncol=5))
#' colnames(data) <- c("math" , "english" , "biology" , "music" , "R-coding" )
#' rownames(data) <- paste("mister" , letters[1:3] , sep="-")
#'
#' rownames(data[-c(1,2),])
#'
#' # plot_data <- rbind(rep(max(df_plot),ncol(df_plot),) , rep(min(df_plot),ncol(df_plot)) , df_plot)
#'
#' # make_radar_plot(df_prev,sort(cardiovascular_traits[[1]]),sources_vec =c("All","SR","GP","HESIN","Death") )
#' #
#' #   lst.proportion<-list()
#' #   for (trait in vct_diseases){
#' #     dt_disease_prop<-get_case_count_by_source(definition=dfDefinitions_processed_expanded %>% filter(TRAIT==trait))
#' #     lst.codemap[[cls]]<-fread(fmap)
#' #
#' #   }
#'
#'
#'
#' make_radar_plot<-function(df_prevalence,vct_diseases,sources_vec=c("All","SR","Cancer","HESIN","Death")){
#'
#'
#'   prev_subset<-df_prevalence[,(names(df_prevalence) %in% vct_diseases)]
#'   # prev_subset<-df_prevalence[,vct_diseases]
#'
#'   plot_data <- rbind(rep(max(prev_subset),ncol(prev_subset)) , rep(min(prev_subset),ncol(prev_subset)) , prev_subset[sources_vec,])
#'
#'   plot_data<-plot_data[,sort(colnames(prev_subset))]
#'
#'   disease_names<-dfDefinitions%>% filter(TRAIT %in% names(prev_subset))%>%select(TRAIT,DESCRIPTION) %>% arrange(TRAIT)
#'
#'   colnames(plot_data)<-disease_names$DESCRIPTION
#'   colors_border=rand_col_vec(ncol(prev_subset),0.7,cust_seed=34)
#'   colors_in=rand_col_vec(ncol(prev_subset),0.2,cust_seed=34)
#'   radarchart( plot_data  , axistype=1 ,
#'               #custom polygon
#'               pcol=colors_border , pfcol=colors_in , plwd=2 , plty=1,
#'               #custom the grid
#'               cglcol="grey", cglty=1, axislabcol="grey", seg=4,caxislabels=round(seq(min(prev_subset),max(prev_subset),length.out=5),1), cglwd=1,
#'               #custom labels
#'               vlcex=1)
#'   # Add a legend
#'   legend(x=1.5, y=0.8, legend = rownames(plot_data[-c(1,2),]), bty = "n", pch=20 , col=colors_border , text.col = "black", cex=1.2, pt.cex=3)
#' }
#'
#'
#'
#'
#'
#' # head(case_status)
#' #
#' # count(case_status,Any.sr)
#' # count(case_status,Any.hes)
#'
#'
#' #TODO
#' # a function to show (number of case) unique to each source
#'
#' #add radarplots
#'
