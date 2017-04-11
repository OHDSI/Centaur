
#' Trim Population
#'
#' Trims the population based on the propensity scores
#'
#' This function trims the supplied population based on propensity scores. Trimming can be
#' based on the overlapping regions, or a specified interquantile range.
#'
#' @param data Data Frame - containing the dataset with previously calculated PS. The data
#' frame must contain a treatment indicator variable called 'treat' and propensity score values called 'ps_values'.
#' @param trim.method String, specifying the method to use to trim the dataset. Available options
#'                      are "overlap" (default) to trim the non-overlapping tails of the control and
#'                      treatment distributions, or "quantile" to trim to an interquantil range.
#' @param trim.quantile Number, specifying the interquantile range to be trimmed. For example,
#'                      0.95 indicates that the dataset should be trimmed to the 0.025 / 0.975 interquantile range.
#' @param quantile.group Specifies the group on which to determine the quantile values. Options are "all" (default)
#'                      to determine the quantile values based on the propensity scores for the entire population, or
#'                      "treat" to determine the quantile values based on the propensity scores for only the treatment
#'                      group
#' @return Data Frame - trimmed dataset
#'
#' @examples
#' ps.trim(myData)
#' ps.trim(myData, trim.method = "quantile", trime.quantile = 0.90)
#' @export
ps.trim <- function(data, trim.method = "overlap", trim.quantile = 0.95, quantile.group = "all") {

  if (trim.method == "overlap") {
    min.trt <- min(data[data$treat == 1,]$ps_values)
    max.ctl <- max(data[data$treat == 0,]$ps_values)
    to.be.removed <- nrow(data[(data$treat == 0 & data$ps_values < min.trt) | (data$treat == 1 & data$ps_values > max.ctl),])

    print(paste(to.be.removed, " subjects to be removed by trimming tails of non-overlapping regions"))
    data <- data[(data$treat == 0 & data$ps_values >= min.trt) | (data$treat == 1 & data$ps_values <= max.ctl),]
  }
  else if (trim.method == "quantile") {
    if (trim.quantile < 0.5 | trim.quantile > 1) stop("Trim quantile must be [0.5:1] representing the interquantile range to be retained")

    trim.lower.quantile <- (1 - trim.quantile) / 2
    trim.upper.quantile <- 1 - trim.lower.quantile

    if (quantile.group == "all") {
      ps <- data$ps_values
    }
    else if (quantile.group == "treat") {
      ps <- data[data$treat == 1,]$ps_values
    }
    else {
      stop("Unrecognized value for propensity score quantile group. Options are \"all\" or \"treat\"")
    }

    q.upper <- quantile(ps, trim.upper.quantile)
    q.lower <- quantile(ps, trim.lower.quantile)
    to.be.removed <- nrow(data[data$ps_values < q.lower | data$ps_values > q.upper,])

    print(paste(to.be.removed, " subjects to be removed by trimming to ", trim.lower.quantile, "/", trim.upper.quantile, " quantile"))
    data <- data[data$ps_values > q.lower & data$ps_values < q.upper,]
  }
  else {
    stop(paste("Unrecognized trimming method specified:", trim.method))
  }

  return(data)
}
