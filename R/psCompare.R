
#' Propensity Score Comparison
#'
#' Generates a graph used to compare the distribution of propensity score values in
#' the cotnrol and treatment groups
#'
#' This function uses sm.density.compare to plot an overlapping density distribution plot
#' of the propensity scores for the two groups allowing one to visually assess the area
#' of common support
#'
#' @param data Data Frame - containing the dataset with previously calculated PS. The data
#' frame must contain a treatment indicator variable called 'treat' and a propensity score value
#' called 'ps_values'.
#' @return NULL
#'
#' @examples
#' \dontrun{
#' ps.compare(myData)
#' }
#' @export
ps.compare <- function(data) {
  sm::sm.density.compare(data$ps_values, data$treat, xlab="Propensity Score")
  title(main="Propensity Score Distribution")
  legend("topleft", c("control", "treat"), fill=c(2,3))
}
