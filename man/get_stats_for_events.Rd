% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization_summary_level.r
\name{get_stats_for_events}
\alias{get_stats_for_events}
\title{Summary statistics of events for single phenotype}
\usage{
get_stats_for_events(all_event_dt, color_seed = 117)
}
\arguments{
\item{all_event_dt}{data table containing all events}
}
\value{
a list of 3 summary tables and 4 plots
}
\description{
Given a data.table with all events for a phenotype, get the summary of code occurence/occurence
}
\examples{
all_event_dt <- get_all_events(dfDefinitions_processed_expanded[1,],lst.data,df.data.settings)
get_stats_for_events(all_event_dt)
}
\keyword{event}
\keyword{stats}
