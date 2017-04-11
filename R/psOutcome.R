# Outcome analysis functions for propensity score balanced populations.

#' Outcome Regression Model
#'
#' Fits a regression model to the specified outcome variable
#'
#' This function uses the glm package to fit a regression model to the specified outcome variable. By default, the
#' regression model will fit outcome ~ treat + ps_values. Inclusion of the PS values in the outcome model is
#' recommended by literature producing a "doubly robust" analysis. In addition, any unbalanced covariates be
#' included in the vector of covariate names parameter. These will also be included in the regression model. Weights
#' for the regression model can be specified, or default to use the calculated weights in the data frame.
#'
#' @param data Data frame, containing the dataset to be analyzed. The data
#' frame must contain a treatment indicator variable called 'treat' and propensity score values called 'ps_values'.
#' @param outcome String, containing the outcome variable name to be analyzed
#' @param covariates Vector, containing the set of covariate variable names to include in the regression
#' @param family Model family, passed through to glm. Defaults to binomial() - see ?glm for additional documentation
#' @param w Vector, containing the subject weights. Defaults to equally weighted.
#'          If analysis of matched data is desired, set this value to myData$is_matched.
#' @return Object, containing fitted model values. In addition to standard glm/lm output, the treatment effect
#'          is appended to the model object as model$treatment.effect. For dichotomous outcome variables, this is
#'          the odds ratio with confidence interval.
#'
#' @examples
#' ps.regression(myData, "outcome")
#' ps.regression(myData, "outcome", covariates, w = myData$is_matched)
#' @export
ps.regression <- function(data, outcome, covariates = character(), family = binomial(), w = NULL) {

  # Construct the model formula
  covariates <- c("treat", "ps_values", covariates)
  formula <- as.formula(paste(outcome, paste(covariates, collapse = " + "), sep = " ~ "))

  if (is.null(w)) w = rep.int(1, nrow(data))

  # Test if the outcome variable is dichotomous...
  is.dichotomous <- all(data[,outcome] %in% 0:1)

  # Regress
  if (is.dichotomous) {
    model <- glm(formula, data, family = family, weights = w, na.action = na.exclude)
  }
  else {
    model <- lm(formula, data, weights = w, na.action = na.exclude)
  }

  # Evaluate the treatment effect
  model$treatment.effect <- cbind(coef(model), confint(model))
  if (is.dichotomous) {
    model$treatment.effect <- exp(model$treatment.effect)
  }
  colnames(model$treatment.effect) <- c("Treat. Eff.", "2.5 %", "97.5 %")
  model$treatment.effect <- model$treatment.effect["treat",]

  return(model)
}


#' Outcome Survival Model
#'
#' Performs survival analysis of an outcome variable
#'
#' This function uses the survival package to do basic survival analysis for an outcome variable. The outcome
#' variable should be previously computed as a survival object (see ps.compute.survival). The analysis will
#' automatically perform survival::survfit and survival::coxph.
#'
#' @param data Data frame, containing the dataset to be analyzed. The data
#' frame must contain a treatment indicator variable called 'treat' and propensity score values called 'ps_values'.
#' @param outcome String, containing the outcome variable name to be analyzed. The variable must be a Surv object
#' @param covariates Vector, containing the set of covariate variable names to include in the model
#' @param w Vector, containing the subject weights. Defaults to equal weighting.
#'          If analysis of matched data is desired, set this value to myData$is_matched.
#' @return Object, containing fitted model values. In addition to standard glm/lm output, the treatment effect
#'          is appended to the model object as model$treatment.effect. For dichotomous outcome variables, this is
#'          the odds ratio with confidence interval.
#'
#' @examples
#' ps.survival(myData, "outcome")
#' @export
ps.survival <- function(data, outcome, covariates = character(), w = NULL) {

  if (!survival::is.Surv(data[,outcome])) stop("Outcome variable for survival analysis must be Survival Object. Use ps.compute.survival()")

  # Construct the model formula
  covariates <- c("treat", "ps_values", covariates)
  formula <- as.formula(paste(outcome, paste(covariates, collapse = " + "), sep = " ~ "))

  if (is.null(w)) w = rep.int(1, nrow(data))

  # survival analysis does not support weights = 0, filter the dataset down (particularly for matched datasets where weights are [0,1])
  data <- data[w > 0,]
  w <- w[w > 0]

  km.analysis <- survival::survfit(as.formula(paste(outcome, "treat", sep = " ~ ")), data = data, weights = w, na.action = na.exclude)
  cph.analysis <- survival::coxph(formula, data = data, weights = w, na.action = na.exclude)

  model <- list(km = km.analysis, cph = cph.analysis)
  return(model)
}


#' Compute Survival Objects
#'
#' Compute survival information using the survival package to be used in a response model
#'
#' This function leverages the survival package Surv function to compute survival information to be
#' used as an outcome variable in outcome analysis. See survival::Surv documentation for more information.
#'
#' @param data Data frame, containing the dataset to be analyzed
#' @param outcome.dates String, containing the variable name with outcome dates
#' @param index.dates String, containing the variable name for the index dates
#' @param obs.end.dates String, containing the variable name for observation end dates (default = "OBSEND")
#' @param death.dates String, containing the variable name for death dates (default = "DEATHDATE")
#' @return Vector, containing computed survival objects based on the outcome dates specified
#'
#' @examples
#' ps.compute.survival(myData, "outcome")
#' @export
ps.compute.survival <- function(data, outcome.dates, index.dates,
                                obs.end.dates = "OBSEND",
                                death.dates = "DEATHDATE") {

  #-----------------------------------------------------------
  # Do some basic QA to confirm that things are defined properly

  # Ensure that the input is a data frame
  if (class(data) != class(data.frame())) stop("The supplied data must be of the class 'data.frame'")

  # Ensure that outcomeDates is not null, and it is in myData
  if(is.null(outcome.dates)) stop("outcome.dates must be supplied")
  if(!(outcome.dates %in% colnames(data))) stop("The data frame must have a column named matching outcome.dates")
  if(!(index.dates %in% colnames(data))) stop("The data frame must have a column named matching index.dates")

  # Calculate the follow up time for each subject - including variety of censoring events
  data$ps.survival.time <- NA

  # First find the time from index to observation end...
  if (obs.end.dates %in% colnames(data)) {
    data$ps.survival.time <- difftime(data[,obs.end.dates], data[,index.dates], units = "days")
  }

  # Censor for death...
  if (death.dates %in% colnames(data)) {
    idx <- !is.na(data[,death.dates])
    data$ps.survival.time[idx] <- difftime(data[idx,death.dates], data[idx,index.dates], units = "days")
  }

  # Now find outcome event dates
  idx <- !is.na(data[,outcome.dates])
  data$ps.survival.time[idx] <- difftime(data[idx,outcome.dates], data[idx,index.dates], units = "days")

  # Truncate negative event times to 0
  data$ps.survival.time <- ifelse(data$ps.survival.time < 0, 0, data$ps.survival.time)

  return(survival::Surv(as.double(data$ps.survival.time), ifelse(is.na(data[,outcome.dates]), 0, 1)))
}


#' Kaplan-Meier Plot
#'
#' Creates a basic Kaplan-Meier plot for the survival analysis
#'
#' This function creates a basic Kaplan-Meier plot for the supplied survival analysis. Basic assumptions are
#' included about extent of axes, labels, etc.
#'
#' @param km.analysis A survival analysis object. See ps.survival()
#' @param title Chart title (default "Kaplan-Meier Curves")
#' @return NULL
#'
#' @examples
#' ps.km.plot(model$km)
#' @export
ps.km.plot <- function(km.analysis, title = "Kaplan-Meier Curves") {

  plot(km.analysis, cex = 0, yscale = 100, ymin = 0.9, conf.int = FALSE, col = 1:2)
  title(title, xlab = "Days", ylab = "Survival Rate (%)")
  legend(x = 150, y = .95, legend = c("Control", "Treatment"), col = 1:2, lty = 1)
}
