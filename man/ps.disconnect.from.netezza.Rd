% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/psDataConnector.R
\name{ps.disconnect.from.netezza}
\alias{ps.disconnect.from.netezza}
\title{Disconnect from Netezza}
\usage{
ps.disconnect.from.netezza(conn)
}
\arguments{
\item{conn}{Database connection object created by ps.connect.to.netezza}
}
\description{
Close a connection to Netezza, freeing resources
}
\details{
This function closes an existing connection to the netezza appliance to free resources
}
\examples{
\dontrun{
ps.disconnect.from.netezza(nzconn)
}
}
