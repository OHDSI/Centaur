
#
# Various AUC functions for PS diagnostics
#


#' Accuracy
#'
#' Calculates the accuracy curve of calculated propensity scores.
#'
#' This function uses the 'AUC' R package to calculate the accuracy of the
#' propensity scores. For this to work properly, the argument \code{...} should
#' be a data frame containing the dataset with previously calculated propensity
#' scores in a variable called ps_values.
#'
#' @param data Data Frame - containing the dataset with previously calculated PS. The data
#' frame must contain a treatment indicator variable called 'treat' and a propensity score value
#' called 'ps_values'.
#' @return Calculated accuracy of propensity scores
#'
#' @examples
#' ps.accuracy(myData)
#' @export
ps.accuracy <- function(data) {

  return (AUC::accuracy(data$ps_values, factor(data$treat)))
}



#' Sensitivity
#'
#' Calculates the sensitivity curve of calculated propensity scores.
#'
#' This function uses the 'AUC' R package to calculate the sensitivity of the
#' propensity scores. For this to work properly, the argument \code{...} should
#' be a data frame containing the dataset with previously calculated propensity
#' scores in a variable called ps_values.
#'
#' @param data Data Frame - containing the dataset with previously calculated PS. The data
#' frame must contain a treatment indicator variable called 'treat' and a propensity score value
#' called 'ps_values'.
#' @return Calculated sensitivity of propensity scores
#'
#' @examples
#' ps.sensitivity(myData)
#' @export
ps.sensitivity <- function(data) {

  return (AUC::sensitivity(data$ps_values, factor(data$treat)))
}



#' Specificity
#'
#' Calculates the specificity curve of calculated propensity scores.
#'
#' This function uses the 'AUC' R package to calculate the specificity of the
#' propensity scores. For this to work properly, the argument \code{...} should
#' be a data frame containing the dataset with previously calculated propensity
#' scores in a variable called ps_values.
#'
#' @param data Data Frame - containing the dataset with previously calculated PS. The data
#' frame must contain a treatment indicator variable called 'treat' and a propensity score value
#' called 'ps_values'.
#' @return Calculated specificity of propensity scores
#'
#' @examples
#' ps.specificity(myData)
#' @export
ps.specificity <- function(data) {

  return (AUC::specificity(data$ps_values, factor(data$treat)))
}



#' ROC
#'
#' Generates a ROC curve for a set of calculated propensity scores.
#'
#' This function uses the 'AUC' R package to generate a ROC curve of the
#' propensity scores. For this to work properly, the argument \code{...} should
#' be a data frame containing the dataset with previously calculated propensity
#' scores in a variable called ps_values.
#'
#' @param data Data Frame - containing the dataset with previously calculated PS. The data
#' frame must contain a treatment indicator variable called 'treat' and a propensity score value
#' called 'ps_values'.
#' @return NULL
#'
#' @examples
#' ps.roc(myData)
#' @export
ps.roc <- function(data) {

  plot(AUC::roc(data$ps_values, factor(data$treat)))
  title(main="Propensity Score ROC")
}

#' AUC
#'
#' Calculates the area under the curve for a set of calculated propensity scores.
#'
#' This function uses the 'AUC' R package to generate a AUC value of the
#' propensity scores. For this to work properly, the argument \code{...} should
#' be a data frame containing the dataset with previously calculated propensity
#' scores in a variable called ps_values.
#'
#' @param data Data Frame - containing the dataset with previously calculated PS. The data
#' frame must contain a treatment indicator variable called 'treat' and a propensity score value
#' called 'ps_values'.
#' @return The area under the curve
#'
#' @examples
#' ps.auc(myData)
#' @export
ps.auc <- function(data) {

  return (AUC::auc(AUC::roc(data$ps_values, factor(data$treat))))
}
