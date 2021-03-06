% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/psScore.R
\name{ps.score}
\alias{ps.score}
\title{Centaur Propensity Score Calculation}
\usage{
ps.score(data, covariates, ps.method = "glm", max.covariates = 200,
  max.twang = 30000, min.twang = 25, control.ratio = 5,
  random.seed = 43762116, odds.ratio = FALSE, lr.summary.file = NULL)
}
\arguments{
\item{data}{Data Frame - containing the dataset to be balanced. Must include treatment indicator as {0, 1} factor.}

\item{covariates}{Vector - containing the list of variable names to be included as confounding variables}

\item{ps.method}{String - name of the method to use for calculation of propensity scores. Options are
\code{glm} (default) and \code{twang}}

\item{max.covariates}{The maximum number of covariates that can be used to calculate propensity scores (default 200)}

\item{max.twang}{The maximum number of samples that can be used with the twang ps.method (default 30k)}

\item{min.twang}{The minimum number of covariates required to use the twang ps.method (default 25)}

\item{control.ratio}{The desired ratio of control to treatment subjects. A large imbalance between
control and treatment subjects can cause problems with algorithm convergence. If the control
population exceeds control.ratio times as many subjects as the treatment group, the control population
will be randomly sampled to the desired size.}

\item{random.seed}{Sets the random number generator seed, which determines how the control population is
sampled. Override the default (43762116) to generate a different sample of control subjects.}

\item{odds.ratio}{A logical argument indicating if the odds ratio and 95% confidence interval should be calculated.
This is an intensive calculation that can take a while.}

\item{lr.summary.file}{The file name where the logistic regression summary should be written to. If directory is
provided as part of the file name, it has to already exist. Default value is NULL. If not given, no file
will be written.}
}
\value{
Data Frame - containing the original dataset, trimmed of any incomplete cases with additional
         variable added for \code{ps_values} (the calculated propensity scores)
}
\description{
Performs propensity score Calculation.
}
\details{
This function calculates the propensity score based on the specified options.
In order for this function to work correctly, the \code{data} argument must be a data frame
containing the collection of confounding variables and a treatment indicator factor variable.
Additionally, the \code{covariates} argument must contain a vector of the variable names to
be included as covariates or confounding variables. See list of arguments for additional
options.
}
\examples{
\dontrun{
ps.score(myData, myCovariates)
ps.score(myData, myCovariates, ps.method = "twang")
}
}
