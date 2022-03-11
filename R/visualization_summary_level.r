

#' return random colors of specified length
#'
#' Given a phenotype, a list of episode data and reference dates per individual, identify cases and control status as of *reference date" by source.
#' @param vec_leng phenotype/trait specified in definition table (a row in the table)
#' @param alpha list of data table with all episode data
#' @param cust_seed data frame containing data settings
#' @return  vector of rgb colour codes
#' @export
rand_col_vec <- function(vec_leng = 5,
                         alpha = 0.75,
                         cust_seed = NULL) {
  if (!is.null(cust_seed)) {
    set.seed(cust_seed)
  }
  # return(rgb(rgb_int[1],rgb_int[2],rgb_int[3],alpha))
  col_vec <- c()
  for (i in 1:vec_leng) {
    rgb_int <- sample(0:255 , 3 , replace = T) / 255
    col_vec <- c(col_vec, do.call(rgb, as.list(c(rgb_int, alpha))))
    if (!is.null(cust_seed)) {
      cust_seed <- cust_seed + 2222 * i
      set.seed(cust_seed)
    }
  }
  return(col_vec)

}




#' Summary statistics of events for single phenotype
#'
#' Given a data.table with all events for a phenotype, get the summary of code occurence/occurence
#' @param all_event_dt data table containing all events
#' @return  a list of 3 summary tables and 4 plots
#' @keywords event stats
#' @export
#' @examples
#' all_event_dt <- get_all_events(dfDefinitions_processed_expanded[1,],lst.data,df.data.settings)
#' get_stats_for_events(all_event_dt)
get_stats_for_events <- function(all_event_dt,color_seed=117) {
  # show stats on codes
  stats.codes <-
    all_event_dt[, .(
      count = .N,
      sum.event = sum(event, na.rm = T),
      sum.epidur = sum(epidur, na.rm = T),
      median.epidur = median(epidur, na.rm = T),
      max.epidur = max(epidur, na.rm = T)
    ), keyby = list(identifier, classification, code)]
  stats.codes <-
    stats.codes %>% dplyr::group_by(classification, code) %>% summarise(count =
                                                                          n())
  stats.codes <- stats.codes %>% arrange(count)
  stats.codes$rank <- 1:nrow(stats.codes)

  color_vec <- rand_col_vec(length(unique(stats.codes$classification)),alpha = 1,cust_seed = color_seed)
  p1 <-
    ggplot2::ggplot(stats.codes,
                    ggplot2::aes(rank, count, label = code, color = classification)) + ggplot2::geom_point() + ggplot2::ylim(-((max(
                      stats.codes$count
                    )) / 3), NA)  + ggrepel::geom_text_repel(size =4.5, segment.size = 0.5) + ggpubr::theme_classic2(base_size = 16) +  ggplot2::scale_color_manual(values = color_vec)
  p2 <-
    ggplot2::ggplot(stats.codes,
                    ggplot2::aes(rank, count, label = code, color = classification)) + ggplot2::scale_y_continuous(trans =
                                                                                                                     'log2') + ggplot2::geom_point()   + ggpubr::theme_classic2(base_size = 16) + ggrepel::geom_text_repel(size = 4.5, segment.size =
                                                                                                                                                                                    0.5) +ggplot2::scale_color_manual(values = color_vec)
  stats.codes.summary.table <- stats.codes
  stats.codes.summary.p <-
    ggpubr::ggarrange(p1,
                      p2,
                      nrow = 1,
                      ncol = 2,
                      common.legend = TRUE)

  # show co-occurences of classifications
  stats.coocurrence <-
    table(all_event_dt[, c("identifier" , "classification")])#[1:10,]
  stats.coocurrence[stats.coocurrence > 0] <- 1

  mat <- crossprod(as.matrix(stats.coocurrence))
  mat <-
    floor(t(mat * 100 / diag(mat)))                 # calculate the percentage
  diag(mat) <- NA
  stats.class.cooccur.table <- mat

  # stats.class.cooccur.p <- pheatmap::pheatmap(mat,display_numbers=mat,cluster_cols = F,cluster_rows = F )
  # #########ggplot2 implementation of the heatmap##########################################################################
  longData <- reshape2::melt(mat)
  names(longData) <- c("Code_presence", "Code_occur", "%")
  stats.class.cooccur.p <-
    ggplot2::ggplot(longData,
                    ggplot2::aes(x = Code_occur , y = Code_presence, fill = `%`)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_distiller(palette = "RdYlBu", na.value = "grey85") +
    ggplot2::labs(
      x = "Classfication of co-occurence",
      y = "Classification present",
      title = "Co-occurence of diagnosis by sources",
      caption = "[% of column covered by row]"
    ) +
    ggplot2::geom_text(ggplot2::aes(label = `%`))
  # ########################################################################################################################
  # show co-occurences of codes
  stats.coocurrence <- table(all_event_dt[, c("identifier" , "code")])
  stats.coocurrence[stats.coocurrence > 0] <- 1
  mat <- crossprod(as.matrix(stats.coocurrence))
  mat <-
    floor(t(mat * 100 / diag(mat)))                 # calculate the percentage
  diag(mat) <- NA
  stats.codes.cooccur.table <- mat


  # stats.codes.cooccur.p <- pheatmap::pheatmap(mat,fontsize = 6,cluster_cols = F,cluster_rows = F)
  # filter on codes that with co-occurence of at least 10% to reduce sparseness and make clustering more informative.
  stats.codes.cooccur.filtered.table <-
    stats.codes.cooccur.table[matrixStats::rowMaxs(stats.codes.cooccur.table, na.rm =
                                                     T) > 10, matrixStats::colMaxs(stats.codes.cooccur.table, na.rm = T) > 10]
  # stats.codes.cooccur.filtered.p <- pheatmap::pheatmap( stats.codes.cooccur.filtered.table  ,fontsize = 6)
  # ########## ggplot implementation#################################################################################
  # Create ggplot version dendrogram from ggdendro

  code.dendro <-
    as.dendrogram(hclust(d = dist(x = stats.codes.cooccur.filtered.table)))
  ddata_x <-  ggdendro::dendro_data(code.dendro)
  # to colour leaves by classifications

  lab_gp <- ggdendro::label(ddata_x)
  lab_gp$group <-
    stats.codes[match(lab_gp$label, stats.codes$code), ]$classification

  stats.codes.cooccur.filtered.p.dendro <-
    ggplot2::ggplot(ggdendro::segment(ddata_x)) +
    ggplot2::geom_segment(ggplot2::aes(
      x = x,
      y = y + 10,
      xend = xend,
      yend = yend + 10
    )) +  ggplot2::geom_text(
      data = ggdendro::label(ddata_x),
      ggplot2::aes(
        label = label,
        x = x,
        y = -5,
        colour = lab_gp$group
      ),
      size = 3,
      angle = 45
    ) + ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.text.x =
        ggplot2::element_blank(),
      axis.text.y =
        ggplot2::element_blank(),
      axis.ticks =
        ggplot2::element_blank(),
      axis.title.x =
        ggplot2::element_blank(),
      axis.title.y =
        ggplot2::element_blank(),
      # legend.position="none",
      panel.background =
        ggplot2::element_blank(),
      panel.border =
        ggplot2::element_blank(),
      panel.grid.major =
        ggplot2::element_blank(),
      panel.grid.minor =
        ggplot2::element_blank(),
      plot.background =
        ggplot2::element_blank()
    )
  # ggplot version heatmap
  # code.order <- order.dendrogram(code.dendro)
  # TODO maker heatmap ordered like dendrogram?

  longData <- reshape2::melt(stats.codes.cooccur.filtered.table)
  names(longData) <- c("Code_presence", "Code_occur", "%")
  stats.codes.cooccur.filtered.p.heat <-
    ggplot2::ggplot(longData,
                    ggplot2::aes(x = Code_occur , y = Code_presence, fill = `%`)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_distiller(palette = "RdYlBu", na.value = "grey85") +
    ggpubr::theme_classic2(base_size = 16) +
    ggplot2::labs(
      x = "Code of co-occurence",
      y = "Code present",
      title = "Co-occurence of diagnosis code",
      caption = "[% of column covered by row]"
    ) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45))
  ###############################################################################################################


  return(
    list(
      stats.codes.summary.table = stats.codes.summary.table,
      stats.codes.summary.p = stats.codes.summary.p,
      stats.class.cooccur.table = stats.class.cooccur.table,
      stats.class.cooccur.p = stats.class.cooccur.p,
      stats.codes.cooccur.table = stats.codes.cooccur.table,
      # stats.codes.cooccur.p = stats.codes.cooccur.p,
      # stats.codes.cooccur.filtered.p = stats.codes.cooccur.filtered.p,
      stats.codes.cooccur.filtered.p.dendro = stats.codes.cooccur.filtered.p.dendro,
      stats.codes.cooccur.filtered.p.heat = stats.codes.cooccur.filtered.p.heat
    )
  )
}


#' Get case-control status (Any) by data sources
#'
#' Given a phenotype, a list of episode data and reference dates per individual, identify cases and control status as of *reference date" by source.
#' @param definition phenotype/trait specified in definition table (a row in the table)
#' @param lst.data list of data table with all episode data
#' @param df.data.settings data frame containing data settings
#' @param vct.identifiers character vector listing the identifiers in the cohort. This is used to define controls if reference_date is not given
#' @keywords event stats phenotype
#' @export
#' @examples
get_case_count_by_source <- function(definition,
                                            lst.data,
                                            df.data.settings,
                                            vct.identifiers,remove_event0=TRUE) {
  if (nrow(definition) == 0) {
    message("No definition is provided.Stop.")
    return(0)
  }
  if (!is.null(vct.identifiers)) {
    lst.case_control <- get_cases_controls( definition,lst.data,df.data.settings, vct.identifiers = vct.identifiers,verbose = FALSE )
  } else{
    message(
      "Parameter vct.identifiers is not provided, Please provide one of them.Exit."
    )
    return()
  }
  all_event_dt <- lst.case_control$all_event_dt.Include_in_cases
  #  individuals with eventdates
  id_dt<-all_event_dt[event != 0]$identifier
  # # individuals without eventdates
  n_no_dt <- length(unique(dplyr::setdiff(all_event_dt[event == 0]$identifier,id_dt)))
  if(remove_event0){
      # remove non-event  NOTE actual events without event date will be removed as well
  all_event_dt <- all_event_dt[event != 0]
  }

  # get first date for each individual by source
  mindate <-all_event_dt[, .(mindt = min(eventdate, na.rm = T)), by = c('identifier', '.id')]
  # to stitch the lines to final point
  max_dt <- max(mindate$mindt)
  # get a cumulative count with rank
  mindate[, cumcnt := data.table::frank(mindt, ties.method = 'first'), by =.id]
  # create rows for stitching the line to final time
  dt_data_end <-mindate[, .(identifier = 9999999, mindt = max_dt, cumcnt = max(cumcnt, na.rm = T)), by = .id][, tail(.SD, 1) , by = .id]
  # rbind to actual dt
  mindate <- rbind(mindate, dt_data_end)
  # create a line of Any
  any_mindate <-all_event_dt[, .(mindt = min(eventdate, na.rm = T)), by = c('identifier')]
  any_mindate[, `:=`(.id = "Any", cumcnt = data.table::frank(mindt, ties.method = 'first'))]
  # add Any
  mindate <- rbind(mindate, any_mindate)
  dt_data_end <-rbind(dt_data_end, any_mindate[any_mindate$cumcnt == max(any_mindate$cumcnt), ])

  return(list(mindate=mindate,dt_data_end=dt_data_end,n_no_dt=n_no_dt))
}


#' # data source trajectory over time
#'
#' Given a phenotype definition and list of data-tables with events stratified by sources, plot the case growth over time by sources
#' @param definitions phenotype/trait specified in definition table (a row in the table)
#' @param lst.data list of data table with all episode data
#' @param df.data.settings data frame containing data settings
#' @param vct.identifiers character vector listing the identifiers in the cohort. This is used to define controls if reference_date is not given
#' @return  ggplot object
#' @keywords event stats timeline phenotype
#' @export
#' @examples
plot_disease_timeline_by_source <- function(definition,
                                            lst.data,
                                            df.data.settings,
                                            vct.identifiers) {
  if (nrow(definition) == 0) {
    message("No definition is provided.Stop.")
    return(0)
  }
  if (!is.null(vct.identifiers)) {
    # get diease events

    lst.case_control <-
      get_cases_controls(
        definitions = definition,
        lst.data = lst.data,
        df.data.settings = df.data.settings,
        vct.identifiers = vct.identifiers,verbose=FALSE
      )
  } else{
    message(
      "Parameter vct.identifiers is not provided, Please provide one of them.Exit."
    )
    return()
  }
  # for annotation , rounding to month
  dt_v0_min <- as.Date("2006-03-01")
  dt_v0_max <- as.Date("2010-10-01")
  dt_v1_min <- as.Date("2009-12-01")
  dt_v1_max <- as.Date("2013-06-01")
  dt_v2_min <- as.Date("2014-05-01")
  dt_v2_max <- as.Date("2020-03-01")
  dt_v3_min <- as.Date("2019-06-01")
  dt_v3_max <- as.Date("2021-06-01")

  case_by_source<-get_case_count_by_source( definition = definition,
                                            lst.data = lst.data,
                                            df.data.settings = df.data.settings,
                                            vct.identifiers = vct.identifiers,remove_event0 = TRUE)
  mindate<-case_by_source$mindate
  dt_data_end<-case_by_source$dt_data_end
  n_no_dt<-case_by_source$n_no_dt


  color_vec <- rand_col_vec(length(unique(mindate$.id)),alpha=1,cust_seed = 477821469)
  # color_vec<-ggpubr::get_palette("rickandmorty",length(unique(mindate$.id)))
  plt_time <-
    ggplot2::ggplot(mindate, ggplot2::aes(
      x = mindt,
      y = cumcnt,
      colour = .id,
      group = .id
    )) +
    ggplot2::scale_color_manual(values = color_vec, name = "Source") +
    ggplot2::geom_line(alpha = 1, size = 1) + ggplot2::geom_point(size =
                                                                    0.8,
                                                                  alpha = 0.1,
                                                                  shape = 20) + ggplot2::xlab("Time of diagnosis") +
    ggplot2::ylab("Count") +
    # annotation window v0
    ggplot2::annotate(
      "rect",
      xmin = dt_v0_min,
      xmax = dt_v0_max,
      ymin = 0,
      ymax = Inf,
      alpha = .2,
      fill = "#6395b6"
    ) +
    ggplot2::annotate(
      geom = "text",
      x = as.Date(mean.Date(c(
        dt_v0_min, dt_v0_max
      ))),
      y = max(dt_data_end$cumcnt) * 1.15,
      label = "Baseline visits",
      color = "#213745"
    ) +
    # # annotation window v1
    ggplot2::annotate(
      "rect",
      xmin = dt_v1_min,
      xmax = dt_v1_max,
      ymin = 0,
      ymax = Inf,
      alpha = .2,
      fill = "#f4b699"
    ) +
    ggplot2::annotate(
      geom = "text",
      x = as.Date(mean.Date(c(
        dt_v1_min, dt_v1_max
      ))),
      y = max(dt_data_end$cumcnt) * 1.09,
      label = "1st repeat",
      color = "#7d300d"
    ) +
    # # annotation window v0
    ggplot2::annotate(
      "rect",
      xmin = dt_v2_min,
      xmax = dt_v2_max,
      ymin = 0,
      ymax = Inf,
      alpha = .2,
      fill = "#d2baf8"
    ) +
    ggplot2::annotate(
      geom = "text",
      x = as.Date(mean.Date(c(
        dt_v2_min, dt_v2_max
      ))),
      y = max(dt_data_end$cumcnt) * 1.12,
      label = "Imaging visit",
      color = "#2d0b66"
    ) +
    # # annotation window v0
    ggplot2::annotate(
      "rect",
      xmin = dt_v3_min,
      xmax = dt_v3_max,
      ymin = 0,
      ymax = Inf,
      alpha = .2,
      fill = "#ebf179"
    ) +
    ggplot2::annotate(
      geom = "text",
      x = dt_v3_min-25,
      y = max(dt_data_end$cumcnt) * 1.06,
      label = "1st repeat imaging",
      color = "#777d0d"
    ) +
    # # TODO warn missing
    # ppl_miss<-paste0(n_no_dt, " cases without date of diagnosis not plotted")
    ggplot2::labs(caption=paste0(n_no_dt, " cases without date of diagnosis not plotted"),size=16) +
    # # log2 scale not nice
    # scale_y_continuous(trans='log10')+
    ggpubr::theme_classic2(base_size = 18) +
    ggrepel::geom_text_repel(ggplot2::aes(label = cumcnt),
      data = dt_data_end,
      fontface = "plain",
      nudge_y = max(dt_data_end$cumcnt) / 50 ,
      size = 5
    )+ggplot2::theme(legend.text=ggplot2::element_text(size=ggplot2::rel(1)))

  # +ggplot2::geom_rect(aes(xmin=dt_v0_max, xmax=dt_v0_max, ymin=0, ymax=Inf),  alpha = .2,fill="#BDC7FA")

  return(plt_time)
}



#' Dotplot by data sources
#'
#' Given phenotypes, make dotplot by data source.
#' @param definitions phenotype/trait specified in definition table (a row in the table)
#' @param lst.data list of data table with all episode data
#' @param lst.data.settings data frame containing data settings
#' @param vct.identifiers character vector listing the identifiers in the cohort. This is used to define controls if reference_date is not given
#' @return  gg object
#' @keywords plot
#' @export
dotplot_diseases_by_source<-function(definitions,
                                     lst.data,
                                     df.data.settings,
                                     vct.identifiers = NULL,standardize=TRUE){

  if (nrow(definitions) == 0) {
    message("No definition is provided.Stop.")
    return(0)
  }
  vct_pheno_lab<-definitions[definitions$Definitions=="Include_in_cases",]$DESCRIPTION
  vct_pheno<-unique(definitions[definitions$Definitions=="Include_in_cases",]$TRAIT)

  vct.identifiers<-lst.harmonized.data$vct.identifiers
  lst.data=lst.harmonized.data$lst.data
  df.data.settings=dfData.settings
  if (is.null(vct.identifiers)) {
    message("Both df.reference.dates and vct.identifiers are NULL, Please provide one of them.Exit.")
    # TODO do i want to add the different reference dates options?
    return()
  }


    # use case 1: case distribution at present time
    lst.case_count<-lapply(vct_pheno,function(x){
    case_by_source<-get_case_count_by_source( definition = definitions[definitions$TRAIT==x,],
                                              lst.data = lst.data,
                                              df.data.settings = df.data.settings,
                                              vct.identifiers = vct.identifiers,remove_event0 = FALSE)
    # mindate<-case_by_source$mindate
    dt_data_end<-case_by_source$dt_data_end
    # n_no_dt<-case_by_source$n_no_dt
    dt_data_end
   })
    names(lst.case_count)<-vct_pheno_lab


   lst.case_count<-lapply(vct_pheno_lab,function(x){

    lst.case_count[[x]]<-lst.case_count[[x]][,c(".id","cumcnt")]
    if (standardize){
      lst.case_count[[x]]$cumcnt<-lst.case_count[[x]]$cumcnt/max(lst.case_count[[x]]$cumcnt)
      }
    names(lst.case_count[[x]])[names(lst.case_count[[x]])=='cumcnt']<-x
    lst.case_count[[x]]
    })

  pheno_cnt_table<- Reduce(
      function(x, y, ...) merge(x, y, by=".id",all = TRUE, ...),
      lst.case_count)

    pheno_cnt_long_table<-reshape2::melt(pheno_cnt_table)
    pheno_cnt_long_table[is.na(pheno_cnt_long_table)]<-0

     color_vec <- rand_col_vec(length(unique(pheno_cnt_long_table$.id)))

    dotplot<-ggpubr::ggdotchart(
      pheno_cnt_long_table, x = "variable", y = "value",
      group = ".id", color = ".id", palette =color_vec, ggtheme = ggpubr::theme_pubclean(base_size = 12),dot.size = 3.5,shape=18 ,                              # Large dot size
       size = 1,
      add = "segment", position = ggplot2::position_dodge(0.3),
      sorting = "descending"
    )+
    # "rickandmorty"
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,vjust=1))
    if(standardize){
      dotplot<-dotplot+ggplot2::labs(color="data source",x="Phenotype",y="proportion")
      }else{
        dotplot<-dotplot+ggplot2::labs(color="data source",x="Phenotype",y="count")
      }
    return(dotplot)

         }






#' Get case-control status (Hx) by data sources used for upset plot
#'
#' Given a phenotype, a list of episode data and reference dates per individual, identify cases and control status as of *reference date" by source.
#' @param definitions phenotype/trait specified in definition table (a row in the table)
#' @param lst.data list of data table with all episode data
#' @param lst.data.settings data frame containing data settings
#' @param reference_date reference dates for each individuals in the whole cohort as a named vector
#' @param vct.identifiers character vector listing the identifiers in the cohort. This is used to define controls if reference_date is not given
#' @return  a data table : case control status according to different sources and one column based on any of the sources.
#' @keywords time-to-event
#' @export
#' @examples
#' get_case_status_by_source(cancer_source,definition=dfDefinitions_processed_expanded %>% filter(TRAIT=="Nicm"), lst.data,lst.data.settings,lst.identifiers=dfukb$f.eid))
get_case_status_by_source <- function(definition,
                                     lst.data,
                                     df.data.settings,
                                     df.reference.dates = NULL,keep_any=FALSE,keep_id=FALSE) {
  if (nrow(definition) == 0) {
    message("No definition is provided.Stop.")
    return(0)
  }
  if (!is.null(df.reference.dates)) {
    lst.case_control <-
      get_cases_controls(
        definition,
        lst.data,
        df.data.settings = df.data.settings,
        df_reference_date = df.reference.dates,
        verbose = FALSE
      )
  } else{
    message(
      "Both df.reference.dates and vct.identifiers are NULL, Please provide one of them.Exit."
    )
    return()
  }
  ###########################################################################
  # NOTE: this function depends on the get_cases_controls! modify accordingly
  ##########################################################################
  # ##########################################################################
  # exclude future events
  # otherwise the comparison between sources is not fair
  # ##########################################################################

  cases <-
    lst.case_control$df.casecontrol[lst.case_control$df.casecontrol$Any == 2, ]
  # only those ppl with only future events are excluded, because these ppl would be considered as control at the reference date
  # ppl_future_only <-
    # cases[cases$Fu != 2 & cases$Any == 2 & cases$Hx == 2, ]$identifier
  if (nrow(definition) > 1) {
    # take first row for the message
    definition <- definition[1,]
  }

  # discard future events
    names(df.reference.dates)<-c("identifier","reference_date")

    lst.case_control$all_event_dt.Include_in_cases<-merge(lst.case_control$all_event_dt.Include_in_cases,df.reference.dates,by="identifier",all.x = TRUE,all.y=FALSE)
  lst.case_control$all_event_dt.Include_in_cases <-
    lst.case_control$all_event_dt.Include_in_cases[(
      lst.case_control$all_event_dt.Include_in_cases$eventdate <=lst.case_control$all_event_dt.Include_in_cases$reference_date
    )]
  message(glue::glue("{definition$DESCRIPTION}: {length(unique(lst.case_control$all_event_dt.Include_in_cases$identifier))} indivduals have events recorded by the corresponding reference dates"))



  # get all available data sources from data
  # all_sources <-
    # unique(lst.case_control$all_event_dt.Include_in_cases$.id)
  # all_sources <-
    # unique(lst.case_control$all_event_dt.Include_in_cases$classification)
  `Self-report (nurse interview)`<-c("tte.sr.20002","tte.sr.20001" ,"tte.sr.20004", "sr.20003")
  `Cancer registry`<-c("tte.cancer.icd10", "tte.cancer.icd9")
  `Death registry`<-c("tte.death.icd10.primary" ,"tte.death.icd10.secondary")
  `Hospital inpatient records` <-c("tte.hesin.oper3.primary" ,"tte.hesin.oper3.secondary","tte.hesin.oper4.primary","tte.hesin.oper4.secondary", "tte.hesin.icd10.primary","tte.hesin.icd10.secondary","tte.hesin.icd9.primary","tte.hesin.icd9.secondary")
  `Primary care` <-c("tte.gpclinical.read2","tte.gpclinical.read3","tte.gpscript.dmd.england","tte.gpscript.bnf.england","tte.gpscript.bnf.scotland","tte.gpscript.read2.wales")
  `Self-report (touchscreen)` <-c("ts")
  all_sources<-list(`Self-report (nurse interview)` = `Self-report (nurse interview)`,`Cancer registry`=`Cancer registry`,`Death registry`=`Death registry`,`Hospital inpatient records`=`Hospital inpatient records`, `Primary care`=`Primary care`,`Self-report (touchscreen)`=`Self-report (touchscreen)`)
  # create new columns for each source
  for (source_name in names(all_sources)) {
    # varname <- paste(source_name, 'Hx', sep = ' ')
    varname <-source_name

    # create new columns for each source
    cases[, (varname)] <- as.numeric(NA)
    # lookup source from df with all episodes i.e. all_event_dt.Include_in_cases
    # if rows in include_in_case which originated from the target source, look up the eid and set to 2
    cases[[varname]][cases$identifier %in% lst.case_control$all_event_dt.Include_in_cases[lst.case_control$all_event_dt.Include_in_cases$.id %in% all_sources[[source_name]], ]$identifier] <- 2
    # cases[[varname]][cases$identifier %in% lst.case_control$all_event_dt.Include_in_cases[lst.case_control$all_event_dt.Include_in_cases$classification %in% source_name, ]$identifier] <- 2

  }

  # for (j in  paste( names(all_sources),"Hx",sep=" ")) {
  for (j in  names(all_sources)) {

    # for all new Hx_ columns , set everything to 0 for ease of counting
    # if Hx <0 , (excluded cases)
    set(cases, which(cases$Hx <= 0), j, 0)
    #  if the cell value was NA ->  control
    set(cases, which(is.na(cases[[j]])), j, 0)
  }
  cases$Any[!cases$identifier %in% lst.case_control$all_event_dt.Include_in_cases$identifier] <-0
  # remove these future cases
  cases<-cases[cases$Any!=0]

  if(keep_any){
    # cols <- c("Any", paste( names(all_sources),"Hx",sep=" "))
    cols <- c("Any", names(all_sources))

  }else{
    # cols <- c(paste( names(all_sources),"Hx",sep=" "))
    cols <- c("Any", names(all_sources))

  }
  if(keep_id){
    cols <- c("identifier",cols)
  }

  # to 1/0 for upset plot
  cases<-cases[ ,cols,with=FALSE]
  #  to apply on (.SDcols) all columns except identifier, replace values in these columns non 2 to 0 and 2 to 1
  cases[,(names(cases)[!names(cases) %in%("identifier")]) := replace(.SD, .SD != 2, 0), .SDcols =names(cases)[!names(cases) %in%("identifier")]]
  cases[,(names(cases)[!names(cases) %in%("identifier")]) := replace(.SD, .SD == 2, 1), .SDcols =names(cases)[!names(cases) %in%("identifier")]]

  # cases[cases[!=2]<-0
  # cases[cases==2]<-1
  return(cases)
  # # alternatve to the block below
  # # df_prop<-cases[ ,cols,with=FALSE]%>% dplyr::summarise(across(where(is.numeric),sum))
  # # df_prop<-df_prop/2
  # # df_prop<-data.table(t(df_prop),keep.rownames = TRUE)
  # # colnames(df_prop)<-c("data_source","proportion")
  #
  # df_prop <- data.table(unlist(lapply(cases[, (cols)], function(y) {
  #   v1 <- nrow(cases[cases[[y]] == 2, ])
  # })))
  # if (standardize) {
  #   df_prop$V1 <- df_prop$V1 / max(df_prop$V1)
  #   colnames(df_prop) <- "proportion"
  # } else{
  #   colnames(df_prop) <- "n"
  # }
  # df_prop$source <- c("Any", all_sources)
  # # sort
  # df_prop <- df_prop[order(df_prop$source), ]
  # return(df_prop)
}


#' Make upset plot
#'
#' Given a phenotype and reference date , make an upset plot of cases (Hx) by the reference dates by source.
#' @param definitions phenotype/trait specified in definition table (a row in the table)
#' @param lst.data list of data table with all episode data
#' @param lst.data.settings data frame containing data settings
#' @param reference_date reference dates for each individuals in the whole cohort as a named vector
#' @return  upset plot
#' @keywords upset
#' @export
#' @examples
#' make_upsetplot(definition=dfDefinitions_processed_expanded%>%filter(TRAIT==trait),lst.harmonized.data$lst.data,df.data.settings,df.reference.dates = df_reference_dt_v0)
make_upsetplot <- function(definition,lst.data,
                                      df.data.settings,
                                      df.reference.dates) {
  # get the count
  dt_status_by_source<-get_case_status_by_source(definition,lst.data,df.data.settings,df.reference.dates)

  # make plot
  # dark green, mud yellow,brwon, blue green,light purple,bright light blue green
  cus_col<-c("#4d8076","#926c00","#4c1130","#087593","#9b89b3","#00c9a7")
  upset_plot<-ComplexUpset::upset(
    dt_status_by_source, c("Self-report (nurse interview)","Cancer registry","Death registry","Hospital inpatient records","Primary care","Self-report (touchscreen)"),min_size=0, base_annotations=list(
      'Intersection size'=ComplexUpset::intersection_size(
        text_colors=c(on_background='#000066', on_bar='#FF8C00'),text=list(size=5.5),fill="#333333"
      )
      # + ggplot2::annotate(
      #   geom='text', x=Inf, y=Inf,
      #   label=paste('Total:', nrow(case_status_by_classification)),
      #   vjust=1, hjust=1
      # )
      +  ggplot2::ylab('Intersection size')
      + ggplot2::theme(text=ggplot2::element_text(size=18),axis.text = ggplot2::element_text(size=18))
    ),
    width_ratio=0.25
    ,queries=list(
      ComplexUpset::upset_query(set="Self-report (nurse interview)", fill=cus_col[1]),
      ComplexUpset::upset_query(set="Cancer registry", fill=cus_col[2]),
      ComplexUpset::upset_query(set="Death registry", fill=cus_col[3]),
      ComplexUpset::upset_query(set="Hospital inpatient records", fill=cus_col[4]),
      ComplexUpset::upset_query(set="Primary care",fill=cus_col[5]),
      ComplexUpset::upset_query(set="Self-report (touchscreen)",fill=cus_col[6])
    ),
    encode_sets=FALSE,  # for annotate() to select the set by name disable encoding
    set_sizes=(
      ComplexUpset::upset_set_size()
      + ggplot2::geom_text(ggplot2::aes(label=..count..), hjust="inward", stat='count',size=5.5)
      # you can also add annotations on top of bars:
      # + annotate(geom='text', label='@', x='Drama', y=850, color='white', size=3)
      + ggplot2::expand_limits(y=nrow(dt_status_by_source))
      + ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90),text = ggplot2::element_text(size=16))
    ),themes=ComplexUpset::upset_default_themes(text=ggplot2::element_text(size=18))
  )
  return(upset_plot)
}





