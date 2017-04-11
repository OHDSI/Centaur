
#' Compute Strata
#'
#' Stratifies the popultion.
#'
#' This function stratifies a popultion based on an estimated propensity score. The stratum
#' of each subject is added to the dataset as a categorical value.
#'
#' @param data Data Frame - containing the dataset with previously calculated PS. The data frame
#' must contain propensity scores in a variable called 'ps_values'.
#' @param n Number of strata to divide the dataset into (default 5)
#' @return Vector - containing the categorical strata number for each subject
#'
#' @examples
#' ps.stratify(myData)
#' ps.stratify(myData, n = 10)
#' @export
ps.stratify <- function(data, n = 5) {

  return (dplyr::ntile(data$ps_values, n))
}


#' Evaluate a regression model
#'
#' Evaluate a regression outcome model for the stratified dataset.
#'
#' This function evaluates an outcome model for each strata and then pools the models
#' to evaluate the overall treatment effect
#'
#' @param data Data Frame - containing the stratified dataset. The data
#' frame must contain a treatment indicator variable called 'treat' and a propensity score value
#' called 'ps_values'.
#' @param strata Vector - containing the stratum that each subject is assigned to
#' @param outcome The outcome variable name
#' @param estimand String, the desired estimand ("ATT" or "ATE"), default (ATT).
#' @param covariates Vector, containing the names of the additional covariates to include in the outcome analysis
#' @param family - Model family, passed through to glm
#' @return An object of type "strata.model". This object contains the model information for each strata in object$models,
#'          as well as the pooled treatment effect information in object$treatment.effect.
#'
#' @examples
#' ps.strata.regression(myData, strata, "outcome")
#' @export
ps.strata.regression <- function(data, strata, outcome, estimand = "ATT", covariates = character(), family = binomial()) {

  # Test if the outcome variable is dichotomous...
  is.dichotomous <- all(data[,outcome] %in% 0:1)

  # We need to fit a regression model for each strata
  models <- data.frame(n = unique(strata));
  n <- nrow(models)

  # Create place holders for the computed information
  models$regression <- vector("list", n)
  models$treatment.effect <- vector("list", n)

  t.effect <- numeric(n)
  t.effect.var <- numeric(nrow(models))

  for (i in 1:n) {
    ds <- data[strata == models$n[i],]
    m <- ps.regression(ds, outcome, covariates = covariates, family = family)
    models$regression[[i]] <- broom::tidy(m)
    models$treatment.effect[[i]] <- m$treatment.effect
    t.effect[[i]] <- m$treatment.effect[[1]]
    if (is.dichotomous) {
      t.effect.var[[i]] <- exp(models$regression[[i]][["std.error"]][[2]]) * exp(models$regression[[i]][["std.error"]][[2]]) * nrow(ds)
    } else {
      t.effect.var[[i]] <- models$regression[[i]][["std.error"]][[2]] * models$regression[[i]][["std.error"]][[2]] * nrow(ds)
    }
  }

  # Calculate the weights for either ATT or ATE pooling
  if (estimand == "ATT") {
    nt <- nrow(data[data$treat == 1,])
    models$weights <- apply(models, 1, function(d) {
      return (nrow(data[strata == d$n & data$treat == 1,]) / nt)
    })
  }
  else if (estimand == "ATE") {
    nt <- nrow(data)
    models$weights <- apply(models, 1, function(d) {
      return (nrow(data[strata == d$n,]) / nt)
    })
  }
  else {
    stop("Unrecognized estimand. Allowed values are 'ATE' and 'ATT'")
  }

  # Calculate the pooled treatment effect
  pooled.treatment.effect <- weighted.mean(t.effect, models$weights)

  # Calculate the pooled variance
  pooled.te.var <- sum(models$weights * models$weights * t.effect.var)
  pooled.te.se <- sqrt(pooled.te.var / nrow(data))

  # Sort the strata to make them easier to interpret
  models <- models[order(models$n),]

  # Stack everything together into a return object
  o <- structure(list(), class="strata.model")
  o$models <- models
  o$treatment.effect <- data.frame(
                                    Treat.Effect = pooled.treatment.effect,
                                    SE.Treat.Effect = pooled.te.se,
                                    CI.Lower = pooled.treatment.effect - 1.96 * pooled.te.se,
                                    CI.Upper = pooled.treatment.effect + 1.96 * pooled.te.se
                                  )

  return(o)
}



#' Evaluate a survival model
#'
#' Evaluate a survival outcome model for the stratified  dataset.
#'
#' This function evaluates an outcome model for each stratum and then pools the models
#' to evaluate the overall treatment effect
#'
#' @param data Data Frame - containing the stratified dataset. The data
#' frame must contain a treatment indicator variable called 'treat' and a propensity score value
#' called 'ps_values'.
#' @param strata Vector - containing the stratum that each subject is assigned to
#' @param outcome The outcome variable name. The variable must contain a precomputed survival object
#' @param estimand String, the desired estimand ("ATT" or "ATE"), default (ATT).
#' @param covariates Vector, containing the names of the additional covariates to include in the outcome analysis
#' @return An object of type "strata.model". This object contains the model information for each strata in object$models,
#'          as well as the pooled treatment effect information in object$treatment.effect.
#'
#' @examples
#' ps.strata.survival(myData, myData$strata, "MISURVIVAL")
#' @export
ps.strata.survival <- function(data, strata, outcome, estimand = "ATT", covariates = character()) {

  # We need to fit a regression model for each strata
  models <- data.frame(n = unique(strata));
  n <- nrow(models)

  # Create place holders for the computed information
  models$km <- vector("list", n)
  models$cph <- vector("list", n)

  t.effect <- numeric(n)
  t.effect.var <- numeric(nrow(models))

  for (i in 1:n) {
    # Get the subset of data in this strata
    ds <- data[strata == models$n[i],]
    # Use the Centaur survival model
    m <- ps.survival(ds, outcome, covariates = covariates)

    # Organize and store the model information...
    models$km[[i]] <- broom::tidy(m$km)
    models$cph[[i]] <- broom::tidy(m$cph)

    # Grab some information to use in computing the pooled treatment effect (log-hazard ratios)
    t.effect[[i]] <- models$cph[[i]][1,"estimate"]
    t.effect.var[[i]] <- models$cph[[i]][1,"std.error"] * models$cph[[i]][1,"std.error"] * nrow(ds)
  }

  # Calculate the weights for either ATT or ATE pooling
  if (estimand == "ATT") {
    nt <- nrow(data[data$treat == 1,])
    models$weights <- apply(models, 1, function(d) {
      return (nrow(data[strata == d$n & data$treat == 1,]) / nt)
    })
  }
  else if (estimand == "ATE") {
    nt <- nrow(data)
    models$weights <- apply(models, 1, function(d) {
      return (nrow(data[strata == d$n,]) / nt)
    })
  }
  else {
    stop("Unrecognized estimand. Allowed values are 'ATE' and 'ATT'")
  }

  # Calculate the pooled treatment effect
  pooled.treatment.effect <- weighted.mean(t.effect, models$weights)

  # Calculate the pooled variance
  pooled.te.var <- sum(models$weights * models$weights * t.effect.var)
  pooled.te.se <- sqrt(pooled.te.var / nrow(data))

  # Sort the strata to make them easier to interpret
  models <- models[order(models$n),]

  # Stack everything together into a return object
  o <- structure(list(), class="strata.model")
  o$models <- models
  o$treatment.effect <- data.frame(
    Treat.Effect = pooled.treatment.effect,
    SE.Treat.Effect = pooled.te.se,
    CI.Lower = pooled.treatment.effect - 1.96 * pooled.te.se,
    CI.Upper = pooled.treatment.effect + 1.96 * pooled.te.se
  )
  o$treatment.effect <- exp(o$treatment.effect)

  return(o)
}



#' Strata Balance
#'
#' Visualize the covariate balance within each stratum
#'
#' This function creates a bar or box-and-whisker plot visualizing the summary statistics
#' for each covariate in each of the strata
#'
#' @param data Data Frame, containing the dataset. The data
#' frame must contain a treatment indicator variable called 'treat'.
#' @param strata Vector, containing the strata for the dataset
#' @param covariates Vector, containing the covariate variable names to include
#' @return NULL
#'
#' @examples
#' ps.strata.balance(myData, myData$strata, covariates)
#' @export
ps.strata.balance <- function(data, strata, covariates) {

  # Attach the strata to the data frame
  data$ps.strata <- strata

  # sort of strata - gets things neatly sorted for the graphs below
  data <- data[order(data$ps.strata),]
  names <- unique(data$ps.strata)
  names <- sort(names)

  # compute the range of ps values for each strata
  stratum_min = vector()
  stratum_max = vector()
  for(name in names) {
    stratum_min <- c(stratum_min, min(data$ps_values[data$ps.strata  == name]))
    stratum_max <- c(stratum_max, max(data$ps_values[data$ps.strata  == name]))
  }


  ud <- names * 3
  at <- c(ud-0.5, ud+0.5)
  at <- sort(at)

  for (i in 1:length(covariates)) {

    # Test if the covariate variable is dichotomous...
    if (all(data[,covariates[i]] %in% 0:1)) {

      # For dichotomous, we'll do a simple bar chart
      agg <- aggregate(data[,covariates[i]], by=list(data$ps.strata, data$treat), FUN=mean)
      m <- matrix(data = agg$x, nrow = 2, ncol = length(names))
      barplot(m, horiz = TRUE, beside = TRUE, main=paste("Probability of ", covariates[i], " Per Stratum"), names.arg=paste(format(round(stratum_min, 2),nsmall = 2), " - ", format(round(stratum_max, 2), nsmall = 2)),  ylab = "Propensity", col=topo.colors(2, .2))
      legend('topright', c("Control", "Treatment"), fill=topo.colors(2, 0.3))
    }
    else {

      # For continuous, we'll do a simple box-and-whisker chart
      boxplot(as.formula(paste(covariates[i], "~ treat * ps.strata")), data = data, horizontal=TRUE, at=at, yaxt="n", main=paste("Distribution of ", covariates[i], " Per Stratum"), ylab = "Propensity",  col=topo.colors(2, .3))
      axis(2, at=ud, labels=paste(format(round(stratum_min, 2), nsmall = 2), " - ", format(round(stratum_max, 2), nsmall = 2)))
      legend('topright', c("Control", "Treatment"), fill=topo.colors(2, 0.3))
    }
  }
}
