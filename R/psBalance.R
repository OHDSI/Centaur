
#' Propensity Score Population Balancing
#'
#' Calculates the weights and/or perform matching of subject to balance the population
#'
#' This function performs propensity score based population balancing. The details of how
#' the population is balanced depend on the parameters specified by the user, including
#' the requested estimand.
#'
#' @param data Data Frame - containing the dataset with previously calculated PS. The data
#' frame must contain a treatment indicator variable called 'treat' and a propensity score value
#' called 'ps_values'.
#' @param covariates Vector, containing the variable names to be included as potential confounding variables
#' @param estimand String, specifying the desired estimand. Options are "ATT" (default) for the Average Treatment
#'          Effect in the Treated, or "ATE" for the Average Treatment Effect.
#' @param match.subjects Boolean, indicating if matching should be used (default TRUE). This is only applicable
#'          when the estimand is "ATT" as "ATE" can only be estimated via IPTW. If the estimand is "ATT",
#'          SMRW will always be used to generate weights, but if match.subjects is set to TRUE, matches will
#'          also be generated.
#' @param match.exact Vector, containing the list of covariate names to perform exact matching on.
#' @param match.ratio Number, indicating the match ratio of control:treat
#' @param caliper.sigma Number, indicating the width of the caliper to use in matching. A value of 0 indicates
#'          that calipers will not be used. A non-zero value will turn on calipers
#' @param use.logit Boolean, indicating if the propensity score should be converted to logit before matching.
#'          This is generally recommended when using calipers and leads to better balance in most cases.
#' @param truncate.quantile Number, indicating the upper quantile at which to apply weight trimming.
#' @param truncate.method String, indicating the approach to use to trim the dataset. Large weights can adversely
#'          affect the ultimate balance of the population. Two approaches appear in the literature, capping weights
#'          at a quantile value or dropping subjects from the dataset. The default for this parameter is "cap", which
#'          will downwardly adjust any weights larger than the specified quantile to the value at that quantile.
#'          Alternatively, set this parameter to "drop" to remove the subjects from the dataset completely. User will
#'          be notified of how many subjects are lost in this step.
#' @param max.matching The maximum number of samples that can be used in matching (default 50k)
#' @return psBalanceData - Object containing parameters used in balancing, along with the resulting data frame.
#'          The dataframe has additional variable(s) added for the weights and matches. When
#'          matching is used, the is_matched is [0,1] indicating if the subject was matched or
#'          not.
#'
#' @examples
#' \dontrun{
#' ps.balance(myData, covariates)
#' ps.balance(myData, covariates, match.subjects = TRUE, match.exact = c("GENDER"))
#' }
#' @export
ps.balance <- function(data, covariates, estimand = "ATT", match.subjects = TRUE, match.exact = NULL, match.ratio = 1,
                       caliper.sigma = 0, use.logit = FALSE, truncate.quantile = 0.95, truncate.method = "cap",
                       max.matching = 50000) {

  # Set up the data structure we'll use to package things up
  object <- structure(list(), class="psBalanceData")

  object[["covariates"]] <- covariates
  object[["estimand"]] <- estimand
  object[["exact_match"]] <- match.exact
  object[["match_ratio"]] <- match.ratio
  object[["caliper_sigma"]] <- caliper.sigma
  object[["use_logit"]] <- use.logit
  object[["truncate_quantile"]] <- truncate.quantile
  object[["truncate_method"]] <- truncate.method

  if (truncate.quantile > 1 | truncate.quantile < 0.5) stop("Weight truncation quantile must be [0.5, 1]")

  formula <- as.formula(paste("treat ~ ", paste(covariates, collapse = " + ")))
  object[["formula"]] <- formula

  if (estimand == "ATT") {

    # For ATT, evaluate weights via SMRW
    data$weights <- ifelse(data$treat == 1, 1, data$ps_values/(1-data$ps_values))

    # Truncate the dataset based on the calculated weights
    q <- quantile(data[data$treat == 0,]$weights, truncate.quantile)
    if (truncate.method == "cap") {
      data[data$treat == 0 & data$weights > q,]$weights = q
    }
    else if (truncate.method == "drop") {
      to.be.dropped <- nrow(data[data$treat == 0 & data$weights > q,])
      if (to.be.dropped > 0) {
        warning(paste(to.be.dropped, " subjects to be dropped by weight truncation"))
        data[data$treat == 0 & data$weights > q,]$weights = 0
      }
    }
    else {
      stop(paste("Unrecognized weight truncation method: ", truncate.method))
    }

    # Match subjects
    if (match.subjects) {

      if (nrow(data) > max.matching) stop(paste("Matching is only available for datasets with fewer than ", max.matching, " subjects"))

      ## Trim up the dataset based on the formula. Matchit doesn't like null values, even if we aren't using them
      d <- data[,names(data) %in% c(covariates, match.exact, "treat", "ps_values")]

      if(use.logit)
        psValues <- gtools::logit(data$ps_values)
      else
        psValues <- data$ps_values

      m <- MatchIt::matchit(formula, data = d, distance = psValues,
                            method = "nearest", m.order = "random", caliper = caliper.sigma,
                            ratio = match.ratio, exact = match.exact)
      data$is_matched <- m$weights

      object[["is_matched"]] <- m$weights
      object[["match_matrix"]] <- m$match.matrix
    }

    object[["data"]] <- data
    object[["weights"]] <- data$weights
  }
  else if (estimand == "ATE") {

    # For ATE, only IPTW weighting is available as an option
    data$weights <- ifelse(data$treat == 1, 1/data$ps_values, 1/(1-data$ps_values))

    # Truncate the dataset based on the calculated weights
    q <- quantile(data$weights, truncate.quantile)
    if (truncate.method == "cap") {
      data[data$weights > q,]$weights = q
    }
    else if (truncate.method == "drop") {
      to.be.dropped <- nrow(data[data$weights > q,])
      if (to.be.dropped > 0) {
        warning(paste(to.be.dropped, " subjects to be dropped by weight truncation"))
        data[data$weights > q,]$weights = 0
      }
    }
    else {
      stop(paste("Unrecognized weight truncation method: ", truncate.method))
    }

    object[["data"]] <- data
    object[["weights"]] <- data$weights
  }
  else {
    stop(paste("Unrecognized estimand: ", estimand))
  }

  return(object)
}
