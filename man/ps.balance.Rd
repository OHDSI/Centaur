% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/psBalance.R
\name{ps.balance}
\alias{ps.balance}
\title{Propensity Score Population Balancing}
\usage{
ps.balance(data, covariates, estimand = "ATT", match.subjects = TRUE,
  match.exact = NULL, match.ratio = 1, caliper.sigma = 0,
  use.logit = FALSE, truncate.quantile = 0.95, truncate.method = "cap",
  max.matching = 50000)
}
\arguments{
\item{data}{Data Frame - containing the dataset with previously calculated PS. The data
frame must contain a treatment indicator variable called 'treat' and a propensity score value
called 'ps_values'.}

\item{covariates}{Vector, containing the variable names to be included as potential confounding variables}

\item{estimand}{String, specifying the desired estimand. Options are "ATT" (default) for the Average Treatment
Effect in the Treated, or "ATE" for the Average Treatment Effect.}

\item{match.subjects}{Boolean, indicating if matching should be used (default TRUE). This is only applicable
when the estimand is "ATT" as "ATE" can only be estimated via IPTW. If the estimand is "ATT",
SMRW will always be used to generate weights, but if match.subjects is set to TRUE, matches will
also be generated.}

\item{match.exact}{Vector, containing the list of covariate names to perform exact matching on.}

\item{match.ratio}{Number, indicating the match ratio of control:treat}

\item{caliper.sigma}{Number, indicating the width of the caliper to use in matching. A value of 0 indicates
that calipers will not be used. A non-zero value will turn on calipers}

\item{use.logit}{Boolean, indicating if the propensity score should be converted to logit before matching.
This is generally recommended when using calipers and leads to better balance in most cases.}

\item{truncate.quantile}{Number, indicating the upper quantile at which to apply weight trimming.}

\item{truncate.method}{String, indicating the approach to use to trim the dataset. Large weights can adversely
affect the ultimate balance of the population. Two approaches appear in the literature, capping weights
at a quantile value or dropping subjects from the dataset. The default for this parameter is "cap", which
will downwardly adjust any weights larger than the specified quantile to the value at that quantile.
Alternatively, set this parameter to "drop" to remove the subjects from the dataset completely. User will
be notified of how many subjects are lost in this step.}

\item{max.matching}{The maximum number of samples that can be used in matching (default 50k)}
}
\value{
psBalanceData - Object containing parameters used in balancing, along with the resulting data frame.
         The dataframe has additional variable(s) added for the weights and matches. When
         matching is used, the is_matched is [0,1] indicating if the subject was matched or
         not.
}
\description{
Calculates the weights and/or perform matching of subject to balance the population
}
\details{
This function performs propensity score based population balancing. The details of how
the population is balanced depend on the parameters specified by the user, including
the requested estimand.
}
\examples{
\dontrun{
ps.balance(myData, covariates)
ps.balance(myData, covariates, match.subjects = TRUE, match.exact = c("GENDER"))
}
}
