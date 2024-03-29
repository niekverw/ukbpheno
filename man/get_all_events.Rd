% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query-event-tables.r
\name{get_all_events}
\alias{get_all_events}
\title{Get all episodes for a phenotype}
\usage{
get_all_events(
  definition,
  lst.data = lst.data,
  df.data.settings,
  verbose = TRUE
)
}
\arguments{
\item{definition}{phenotype/trait specified in definition table (a row in the table)}

\item{lst.data}{list of data table with all episode data collapsed to 1 datatable}

\item{df.data.settings}{data frame containing data settings}
}
\value{
a data table with all events
}
\description{
Given a phenotype and a list of episode data , extract events for this phenotype from all data sources
}
\examples{
get_all_events(dfDefinitions_processed_expanded[14,],lst.data,df.data.settings)
}
\keyword{time-to-event}
