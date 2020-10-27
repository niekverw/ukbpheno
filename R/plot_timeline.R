# library(scales)
# library(lubridate)



#' Get all episodes for a phenotype
#'
#' Given a phenotype and a list of episode data , extract events for this phenotype from all data sources
#'
#' @param df.data.settings data frame containing data settings
#' @param ind_all_event_dt a data table containing all events for the target individual
#' @param lst.data list of data tables with all episode data collapsed to separate datatables
#' @param ind_identifier identifier of the target individual for the timeline to be visualized
#' @return  a data table with all events
#' @keywords time-to-event
#' @export
#' @examples
#' lst.data.identifier<-lapply(lst.data,function(x) {x[, ("identifier") := lapply(.SD, as.numeric), .SDcols = "identifier"] })  set eid to numeric
#' lst.data.identifier<-lapply(lst.data,function(x) {setkey(x,identifier) }) # double check that everything has the same setkey.
#' plot_individual_timeline(df.data.settings,NULL,lst.data.identifier,ind_identifier="XXXXX")
plot_individual_timeline <- function(df.data.settings,ind_all_event_dt=NULL,lst.data=NULL,ind_identifier=XXXXX,plot_medication=FALSE,keep_classifications=NULL,font_size=5) {
  # credit:  https://benalexkeen.com/creating-a-timeline-graphic-using-r-and-ggplot2/
  # input: a collapsed datatable with ind_all_event_dt for one participant, e.g. for disease codes.
  # alternative input: lst.data and identifier, so that it generates the ind_all_event_dt based on all available data.
  # lst.data, can be keyed on identifier which is MUCH faster, we can create a functiion that checks the key and returns a new keyed object as global var if it iis not right.  lst.data.identifier<-lapply(lst.data,function(x) {setkey(x,identifier) }) # double check that everything has the same setkey.


  # ####################################################################################
  # Why to character? it throws an error when fetching the events for the individual
  # Error in bmerge .....Incompatible join types: x.identifier (double) and i.V1 (character)
  # identifier <- as.character(identifier)
  ind_identifier<-as.numeric(ind_identifier)
  # ####################################################################################
  if(is.null(ind_all_event_dt) & is.null(lst.data)){
    message("ind_all_event_dt and lst.data are null, please provide one")
    return(NULL)
  }
  if(is.null(ind_all_event_dt)){
    #  need to catch non-exisiting identifier first if this line is uncommented
    # lst.data<-lapply(lst.data,function(x) {data.table::setkey(x,identifier) })
    print("extract records for this individual")
      all_event_lst <- lapply(names(lst.data), function(x) {
          if(data.table::key(lst.data[[x]]) =="identifier"){
        lst.data[[x]] [ .(ind_identifier),nomatch=NULL] # nomatch is important
      } else{
        lst.data[[x]] [ lst.data[[x]]$identifier %in% ind_identifier]
      }
    } )
    names(all_event_lst)<-names(lst.data)
    # remove empty dfs frame list
    all_event_lst <- all_event_lst[lapply(all_event_lst,nrow)>0]
    # set key to be eid
    ind_all_event_dt <- plyr::ldply(all_event_lst, data.frame) %>% as.data.table()
    ind_all_event_dt$classification <- df.data.settings[match(ind_all_event_dt$.id ,df.data.settings$datasource),]$classification
    # if there is any event not belonging to the individual or if the whole table is empty which return *logical(0)* instead of FALSE
    if (any(!ind_all_event_dt$identifier %in% ind_identifier) | nrow(ind_all_event_dt) ==0 ){
      message(glue::glue("no data on {ind_identifier}"))
      return(0)
    }
    data.table::setkey(ind_all_event_dt,identifier)
    ######
  }

  ### get iit ini the right formatt.
  df <- ind_all_event_dt %>% dplyr::filter(identifier %in% ind_identifier) %>% as.data.frame()

  # Rationale for this: medication codes easily overwhelm the plot
  if (plot_medication==FALSE){
    df<-df[! df$classification %in% c("f.20003","READ2_drugs","DMD","BNF"),]
  }

  if(!is.null(keep_classifications)){
     df <- df[df$classification %in% keep_classifications, ]
  }
  df <- data.frame(month=month(df$eventdate),
            year=year(df$eventdate),
            code= df$code,
            event=df$event,
            classification=df$classification )

  df <- rbind(df %>% dplyr::filter (event==1) %>% dplyr::group_by(month,year,code,classification) %>% dplyr::mutate(dup=length(code)),
             df %>%  dplyr::filter (event==0 | event ==2)%>% dplyr::arrange(-event) %>% dplyr::distinct(code,classification, .keep_all = TRUE) %>% dplyr::group_by(month,year,code,classification) %>% dplyr::mutate(dup=length(code))
             )
  # change from factor to character
  df$code<- as.character(df$code)

  df[df$dup>1,]$code = paste0(df[df$dup>1,]$code,"(x",df[df$dup>1,]$dup,")")

  df <- unique(df)

  #############
  df$date <- with(df, lubridate::ymd(sprintf('%04d%02d%02d', year, month, 1)))
  df <- df[with(df, order(date)), ]
  # head(df)

  classification_levels <- unique(df$classification)
  # blue green yellow red   max allow 14 classification
  classification_colors <- c("#0070C0", "#00B050", "#FFB200", "#C00000","#7393B3","purple","pink","black","peru","darkblue","cyan4","seagreen","slateblue1","orangered1")
  df$classification <- factor(df$classification, levels=classification_levels, ordered=TRUE)

  positions <- c(0.5, -0.5, 1.0, -1.0, 1.5, -1.5)
  directions <- c(1, -1)

  line_pos <- data.frame(
    "date"=unique(df$date),
    "position"=rep(positions, length.out=length(unique(df$date))),
    "direction"=rep(directions, length.out=length(unique(df$date)))
  )

  df <- merge(x=df, y=line_pos, by="date", all = TRUE)
  df <- df[with(df, order(date, classification)), ]

  # head(df)

  text_offset <- 0.05

  df$month_count <- ave(df$date==df$date, df$date, FUN=cumsum)
  df$text_position <- (df$month_count * text_offset * df$direction) + df$position
  head(df)

  month_buffer <- 2

  month_date_range <- seq(min(df$date) - months(month_buffer), max(df$date) + months(month_buffer), by='month')
  month_format <- format(month_date_range, '%b')
  month_df <- data.frame(month_date_range, month_format)

  year_date_range <- seq(min(df$date) - months(month_buffer), max(df$date) + months(month_buffer), by='year')
  year_date_range <- as.Date(
    intersect(
      lubridate::ceiling_date(year_date_range, unit="year"),
      lubridate::floor_date(year_date_range, unit="year")
    ),  origin = "1970-01-01"
  )
  year_format <- format(year_date_range, '%Y')
  year_df <- data.frame(year_date_range, year_format)

  #### PLOT ####

  timeline_plot<-ggplot2::ggplot(df,ggplot2::aes(x=date,y=0, col=classification, label=code,size=font_size))
  timeline_plot<-timeline_plot+ggplot2::labs(col="Classifications",size=font_size)

  timeline_plot<-timeline_plot+ggplot2::scale_color_manual(values=classification_colors[1:length(classification_levels)], labels=classification_levels, drop = FALSE)
  timeline_plot<-timeline_plot+ggplot2::theme_classic()+ggplot2::theme(legend.title=ggplot2::element_text(size=font_size*3),
                                                             legend.text=ggplot2::element_text(size=font_size*3))

    # Plot horizontal black line for timeline
  timeline_plot<-timeline_plot+ggplot2::geom_hline(yintercept=0,
                                          color = "black", size=0.4)

  # Plot vertical segment lines for codes
  timeline_plot<-timeline_plot+ggplot2::geom_segment(data=df[df$month_count == 1,], ggplot2::aes(y=position,yend=0,xend=date), color='black', size=0.2)

  # Plot scatter points at zero and date
  timeline_plot<-timeline_plot+ggplot2::geom_point(ggplot2::aes(y=0), size=3)

  # Don't show axes, appropriately position legend
  timeline_plot<-timeline_plot+ggplot2::theme(axis.line.y=ggplot2::element_blank(),
                                     axis.text.y=ggplot2::element_blank(),
                                     axis.title.x=ggplot2::element_blank(),
                                     axis.title.y=ggplot2::element_blank(),
                                     axis.ticks.y=ggplot2::element_blank(),
                                     axis.text.x=ggplot2::element_blank(),
                                     axis.ticks.x =ggplot2::element_blank(),
                                     axis.line.x =ggplot2::element_blank(),
                                     legend.position = "bottom"
  )

  # Show text for each month
  #timeline_plot<-timeline_plot+geom_text(data=month_df, aes(x=month_date_range,y=-0.1,label=month_format),size=2.5,vjust=0.5, color='black', angle=90)
  # Show year text
  timeline_plot<-timeline_plot+ggplot2::geom_text(data=year_df, ggplot2::aes(x=year_date_range,y=-0.2,label=year_format, fontface="bold"),size=font_size, color='black',angle=90)
  # Show text for each code
  timeline_plot<-timeline_plot+ggplot2::geom_text(ggplot2::aes(y=text_position,label=code),size=font_size)
  print(timeline_plot)
  return(timeline_plot)
}
