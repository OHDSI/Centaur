#
# Population balance assessment utility functions
#


#' Simple population summary of matched samples
#'
#' Summarize the matched samples.
#'
#' This function does a simple summary counting of the matched vs unmatched samples.
#' Print the resulting table to see counts of the original population, vs the
#' populations that result from previously completed matching. (Not really useful
#' for weighting as all subjects will have a non-zero weight)
#'
#' @param data Data Frame - containing the dataset with previously calculated is_matched. The data
#' frame must contain a treatment indicator variable called 'treat' and a variable indicating whether
#' the subject was matched called 'is_matched'.
#' @return Matrix - containing the summary counts of treat and control in original and matched
#'                populations
#'
#' @examples
#' ps.matched.summary(myData)
#' @export
ps.matched.summary <- function(data) {

  sizes <- matrix(ncol = 2, nrow = 2)
  colnames(sizes) <- c('Control', 'Treated')
  rownames(sizes) <- c('All Data', 'Matched')
  sizes[1,1] <- sum(data$treat == 0)
  sizes[1,2] <- sum(data$treat == 1)
  sizes[2,1] <- sum(data$treat == 0 & data$is_matched != 0)
  sizes[2,2] <- sum(data$treat == 1 & data$is_matched != 0)

  return (sizes)
}


#' Matched Sample Scatterplot
#'
#' Simple visual representation of matched samples
#'
#' This function creates a simple scatter plot of matched and unmatched samples. Works best for
#' smaller sample sizes.
#'
#' @param data Data Frame - containing the dataset with previously calculated matches.  The data
#' frame must contain a treatment indicator variable called 'treat' and a variable indicating whether
#' the subject was matched called 'is_matched'.
#' @return NULL
#'
#' @examples
#' ps.matched.scatterplot(myData)
#' @export
ps.matched.scatterplot <- function(data) {

  # Drop each sample into categories for matched / unmatched and treated / control
  cats <- ifelse(data$is_matched == 0, ifelse(data$treat == 1, 3, 0), ifelse(data$treat == 1, 2, 1))

  # Draw a simple scatter plot of the data
  plot(jitter(cats, 0.5) ~ data$ps_values, type = "p", xlab = "Propensity Score",
       ylab = "", main = "Distribution of Propensity Scores",
       yaxt="n", ylim=c(-0.3, 3.7))

  x <- min(data$ps_values) + ((max(data$ps_values)-min(data$ps_values)) / 2)

  # Add labels to each category of samples
  text(x, y = c(0.5, 1.5, 2.5, 3.5),
       labels=c("Unmatched Control Units",
                "Matched Control Units",
                "Matched Treatment Units",
                "Unmatched Treatment Units"))
}


#' Matched Sample Violinplot
#'
#' Simple visual representation of matched samples
#'
#' This function creates a violin plot of matched and unmatched samples. This plot shows the
#' kernel density for the two populations - with an inner box and whisker plot of summary statistics.
#' This plot tends to work better for large populations than the matched scatter plot.
#'
#' @param data Data Frame - containing the dataset with previously calculated matches. The data
#' frame must contain a treatment indicator variable called 'treat' and a variable indicating whether
#' the subject was matched called 'is_matched'.
#' @return NULL
#'
#' @examples
#' ps.matched.violinplot(myData)
#' @export
ps.matched.violinplot <- function(data) {

  # Drop each sample into categories for matched / unmatched and treated / control
  cats <- ifelse(data$is_matched == 0, ifelse(data$treat == 1, 3, 0), ifelse(data$treat == 1, 2, 1))

  x0 <- data$ps_values[cats == 0]
  x1 <- data$ps_values[cats == 1]
  x2 <- data$ps_values[cats == 2]
  x3 <- data$ps_values[cats == 3]

  if (length(x0) > 0 && length(x3) > 0) {
    vioplot::vioplot(x0, x1, x2, x3,
            col="gold",
            names=(c("Unmatched Control", "Matched Control", "Matched Treatment", "Unmatched Treatment")))
  }
  else if (length(x0) > 0) {
    vioplot::vioplot(x0, x1, x2,
            col="gold",
            names=(c("Unmatched Control", "Matched Control", "Matched Treatment")))
  }
  else if (length(x3) > 0) {
    vioplot::vioplot(x1, x2, x3,
            col="gold",
            names=(c("Matched Control", "Matched Treatment", "Unmatched Treatment")))
  }
  else {
    vioplot::vioplot(x1, x2,
            col="gold",
            names=(c("Matched Control", "Matched Treatment")))
  }

  title("Distribution of Propensity Scores")
}


#' Matched Sample Box and Whisker
#'
#' Simple visual representation of matched samples
#'
#' This function creates a box and whisker plot of matched and unmatched samples.
#'
#' @param data Data Frame - containing the dataset with previously calculated matches. The data
#' frame must contain a treatment indicator variable called 'treat' and a variable indicating whether
#' the subject was matched called 'is_matched'.
#' @return NULL
#'
#' @examples
#' ps.matched.boxplot(myData)
#' @export
ps.matched.boxplot <- function(data) {

  # Drop each sample into categories for matched / unmatched and treated / control
  cats <- ifelse(data$is_matched == 0,
                 ifelse(data$treat == 1, "Unmatched Treatment", "Unmatched Control"),
                 ifelse(data$treat == 1, "Matched Treatment", "Matched Control"))

  cats <- factor(cats, levels = c("Unmatched Control", "Matched Control", "Matched Treatment", "Unmatched Treatment"))

  boxplot(data$ps_values~cats,data=data, main="Distribution of Propensity Scores",
                   xlab="Population", ylab="Propensity Score")

}


#' Covariate Summary Matrix
#'
#' Calculate the summary statistics for covariates
#'
#' This function calculates a variety of summary statistics and organizes them into a matrix.
#' Included summary statistics are the mean and sd in the treatment and control populations,
#' the difference in means between the treat and control populations and the standardized
#' difference of the means between the populations.
#'
#' @param data Data Frame - containing the dataset with previously calculated weights or matches. The data
#' frame must contain a treatment indicator variable called 'treat'.
#' @param covariates Vector - containing the set of covariates for which to evaluate statistics.
#' @param weights Vector, containing the weights to use in evaluating statistics for population. Defaults to
#'        equal weighting.
#' @param outputs String - indicates which covariate type(s) to output. Accepted inputs are "all" (Default), "dichotomous", and "continuous"
#' @return Matrix - containing summary statistics for covariates
#'
#' @examples
#' ps.covariate.statistics(myData, covariates)
#' ps.covariate.statistics(myData, covariates, weights = myData$weights)
#' ps.covariate.statistics(myData, covariates, weights = myData$is_matched)
#' ps.covariate.statistics(myData, covariates, outputs = "dichotomous")
#' @export
ps.covariate.statistics <- function(data, covariates, weights = NULL, outputs = c("all", "dichotomous", "continuous")) {

  outputs <- match.arg(outputs)

  # Identify which features to keep
  # TODO: This can be a lengthy check - might be wise to extract and add this as a parameter to this function
  d.variables <- names(data)[names(data) %in% covariates & apply(data, 2, function(x) { all(x %in% 0:1) })]
  n <- length(covariates)

  if (outputs == "dichotomous") {
    covariates <- covariates[covariates %in% d.variables]
    n <- length(covariates)
    columnNames <- c("Ctrl N", "Treat N", "N Diff", "Std Diff", "Ctrl %", "Treat %")
  } else if (outputs == "continuous") {
    covariates <- covariates[!(covariates %in% d.variables)]
    n <- length(covariates)
    columnNames <- c("Ctrl Mean", "Treat Mean", "Mean Diff", "Std Diff", "Ctrl StdDev", "Treat StdDev")
  } else {
    columnNames <- c("Ctrl N/Mean", "Treat N/Mean", "N/Mean Diff", "Std Diff", "Ctrl %/StdDev", "Treat %/StdDev")
  }

  # Set up the return data structure
  summary <- matrix(ncol = 6, nrow = n)
  colnames(summary) <- columnNames
  rownames(summary) <- covariates

  if (is.null(weights)) weights <- rep.int(1, nrow(data))
  data$w <- weights

  # Pre-allocate our variables
  controlMeans <- numeric(length = n)
  treatMeans <- numeric(length = n)
  controlStdDevs <- numeric(length = n)
  treatStdDevs <- numeric(length = n)
  stdDiffs <- vector(length = n)
  meanDiffs <- numeric(length = n)

  for(i in 1:n)
  {
    var <- covariates[i]

    # First check if this is a dichotomous variable
    if (covariates[i] %in% d.variables)
    {
      controlMeansDichot <- SDMTools::wt.mean(data[data$treat == 0, var], data[data$treat == 0, "w"])
      treatMeansDichot <- SDMTools::wt.mean(data[data$treat == 1, var], data[data$treat == 1, "w"])

      # Calculate standardized difference for dichotomous variable
      stdDiffs[i] <- (treatMeansDichot - controlMeansDichot) / sqrt((treatMeansDichot * (1 - treatMeansDichot) + controlMeansDichot * (1 - controlMeansDichot)) / 2) * 100

      # Count the number of occurances for dichotomous variables (rather than mean)
      controlMeans[i] <- sum(data[(data$treat == 0 & data$w != 0), var])
      treatMeans[i] <- sum(data[(data$treat == 1 & data$w != 0), var])

      # Calculate the weighted prevalence for dichotomous variables (rather than std.dev)
      controlStdDevs[i] <- 100*controlMeansDichot
      treatStdDevs[i] <- 100*treatMeansDichot

    } else {

      # Calculate mean for continuous variables
      controlMeans[i] <- SDMTools::wt.mean(data[data$treat == 0, var], data[data$treat == 0, "w"])
      treatMeans[i] <- SDMTools::wt.mean(data[data$treat == 1, var], data[data$treat == 1, "w"])
      meanDiffs[i] <- treatMeans[i] - controlMeans[i]

      # Calculate the standard deviations for continuous variables
      controlStdDevs[i] <- SDMTools::wt.sd(data[data$treat == 0, var], data[data$treat == 0, "w"])
      treatStdDevs[i] <- SDMTools::wt.sd(data[data$treat == 1, var], data[data$treat == 1, "w"])

      # Calculate standardized difference for continuous variables
      stdDiffs[i] <- (treatMeans[i] - controlMeans[i]) / sqrt((controlStdDevs[i]^2 + treatStdDevs[i]^2) / 2) * 100

    }
  }

  summary[,1] <- controlMeans
  summary[,2] <- treatMeans
  summary[,3] <- summary[,2] - summary[,1]
  summary[,4] <- stdDiffs
  summary[,5] <- controlStdDevs
  summary[,6] <- treatStdDevs

  return(summary)
}


#' Percent Reduction in Standardized Differences of Means
#'
#' Calculate the percent reduction in the standardized differences of the means for
#' covariates
#'
#' This function uses the calculated covariate summary information for a balanced and
#' total / unbalanced population and computes the percent reduction in the standardized
#' differences of the means
#'
#' @param balanced Matrix, containing the summary statistics for covariates in the balanced population
#' @param unbalanced Matrix, containing the summary statistics for covariates in the unbalanced
#'          / total population
#' @param abs Boolean, indicating if the absolute percent reduction should be calculated. Default (FALSE)
#' @return Matrix, containing the percent reduction in standardized differences of means for all covariates
#'
#' @examples
#' ps.percent.reduction(summary.balanced, summary.unbalanced)
#' @export
ps.percent.reduction <- function(balanced, unbalanced, abs = FALSE) {

  # Compare standard differences
  stdDiffChange <- matrix(ncol = 1, nrow = nrow(balanced))
  colnames(stdDiffChange) <- c("% Red. Std. Diff.")
  rownames(stdDiffChange) <- rownames(balanced)

  if (abs) {
    stdDiffChange[,1] <- (abs(unbalanced[,"Std Diff"]) - abs(balanced[,"Std Diff"]))/abs(unbalanced[,"Std Diff"]) * 100
  }
  else {
    stdDiffChange[,1] <- (unbalanced[,"Std Diff"] - balanced[,"Std Diff"])/unbalanced[,"Std Diff"] * 100
  }

  return(stdDiffChange)
}


#' Propensity Score Balance
#'
#' Summarizes the balance of the propensity scores in treatment and control populations
#'
#' This function calculates key summary statistics based on the propensity scores in the treatment
#' and control populations. For a balanced population, one should observe similar mean values of the
#' propensity scores in the two populations. Similarly to covariate diagnostics, one should observe a
#' decrease in the standardized difference of the mean propensity score after balancing. Finally, one
#' should see the ratio of variances in the propensity score move closer to 1 after balancing.
#'
#' @param data Data frame, containing the balanced population, including weights and/or matches. The data
#' frame must contain a treatment indicator variable called 'treat'.
#' @param weights Vector, containing the weights to use in assessing the balanced population (defaults to unweighted).
#' @return Matrix, containing the summary statistics for the calculated propensity score
#'
#' @examples
#' ps.score.summary(myData)
#' ps.score.summary(myData, weights = myData$is_matched)
#' @export
ps.score.summary <- function(data, weights = NULL) {


  summary <- matrix(ncol = 6, nrow = 2)
  colnames(summary) <- c("Ctrl PS Mean", "Ctrl PS Var", "Treat PS Mean", "Treat PS Var", "Std Diff", "Ratio of Var")
  rownames(summary) <- c("Original Population", "Balanced Population")

  if (is.null(weights)) weights = rep.int(1, nrow(data))
  data$w <- weights

  summary[1,1] <- mean(data$ps_values[data$treat == 0])
  summary[1,2] <- var(data$ps_values[data$treat == 0])
  summary[1,3] <- mean(data$ps_values[data$treat == 1])
  summary[1,4] <- var(data$ps_values[data$treat == 1])
  summary[2,1] <- SDMTools::wt.mean(data$ps_values[data$treat == 0 & data$w != 0], wt = data$w[data$treat == 0 & data$w != 0])
  summary[2,2] <- SDMTools::wt.var(data$ps_values[data$treat == 0 & data$w != 0], wt = data$w[data$treat == 0 & data$w != 0])
  summary[2,3] <- SDMTools::wt.mean(data$ps_values[data$treat == 1 & data$w != 0], wt = data$w[data$treat == 1 & data$w != 0])
  summary[2,4] <- SDMTools::wt.var(data$ps_values[data$treat == 1 & data$w != 0], wt = data$w[data$treat == 1 & data$w != 0])

  summary[,5] <- (summary[,3] - summary[,1]) / sqrt((summary[,4]+summary[,2]) / 2)
  summary[,6] <- summary[,4] / summary[,2]

  return(summary)
}


#' QQ Plots
#'
#' Create QQ plots for each of the continuous covariates
#'
#' This function identifies each of the continuous covariates and creates a QQ plot, which allows one
#' to visually inspect the similarity of distributions of covariates in the control and treatment populations
#'
#' @param data Data frame, containing the balanced population data. The data
#' frame must contain a treatment indicator variable called 'treat'.
#' @param covariates Vector, containing the list of covariate names
#' @param weights Vector, containing the weights to use in balanced population calculations. Defaults to unweighted.
#' @return NULL
#'
#' @examples
#' ps.covariate.qq(myData, covariates)
#' ps.covariate.qq(myData, covariates, myData$is_matched)
#' @export
ps.covariate.qq <- function(data, covariates, weights = NULL) {

  if (is.null(weights)) weights = rep.int(1, nrow(data))
  data$w <- weights

  # Remove any omitted samples
  bd <- data[data$w != 0,]

  d.variables <- names(data)[names(data) %in% covariates & apply(data, 2, function(x) { all(x %in% 0:1) })]
  c.variables <- names(data)[names(data) %in% covariates & !(names(data) %in% d.variables)]

  for (i in 1:length(c.variables)) {

    # Data must be sorted by the covariate of interest
    bd <- bd[order(bd[[c.variables[i]]]),]

    x <- findInterval(as.numeric(unlist(bd[c.variables[i]])), Hmisc::wtd.quantile(bd[bd$treat == 0, c.variables[i]], bd[bd$treat == 0, "w"], seq(0, 1, length = 201)), all.inside = T) / 2
    y <- findInterval(as.numeric(unlist(bd[c.variables[i]])), Hmisc::wtd.quantile(bd[bd$treat == 1, c.variables[i]], bd[bd$treat == 1, "w"], seq(0, 1, length = 201)), all.inside = T) / 2

    qqplot(x, y, plot.it = TRUE,
           main=paste(c.variables[i], " - QQ Plot"),
           xlab="Control",
           ylab="Treatment")
    abline(0, 1)
  }
}


#' Covariate Standard Differences of Means Plot
#'
#' Visualization of the covariate standardized differences of means
#'
#' This function creates a plot of each covariate's standardized differences of means before and after balancing
#'
#' @param balanced Matrix, containing the summary statistics for covariates in the balanced population
#' @param unbalanced Matrix, containing the summary statistics for covariates in the unbalanced
#'          / total population
#' @param abs Boolean, indicating if the absolute percent reduction should be calculated. Default (FALSE)
#' @return NULL
#'
#' @examples
#' ps.covariate.biasplot(summary.balanced, summary.unbalanced)
#' ps.covariate.biasplot(summary.balanced, summary.unbalanced, abs = TRUE)
#' @export
ps.covariate.biasplot <- function(balanced, unbalanced, abs = FALSE) {

  pre <- unbalanced[,"Std Diff"]
  post <- balanced[,"Std Diff"]
  vars <- rownames(balanced)

  # Remove infinite values that may occur
  pre[!is.finite(pre)] <- NA
  post[!is.finite(post)] <- NA

  x.label = "Standardized Difference of Means"

  if (abs) {
    pre <- abs(pre)
    post <- abs(post)
    x.label = "Absolute Standardized Difference of Means"
  }

  opar <- par()
  par(mar=c(5, 8, 4, 2) + 0.1)

  low <- min(min(pre, na.rm = T), min(post, na.rm = T))
  high <- max(max(pre, na.rm = T), max(post, na.rm = T))
  xlim = c(min(low, -30), max(high, 30))
  if (abs) {
    xlim = c(0, max(high, 30))
  }

  yticks <- seq(1, length(vars))
  plot(pre, yticks,
       main="Covariate Std. Diff. Reduction",
       xlab=x.label,
       ylab="",
       yaxt="n",
       col="red",
       xlim=xlim)
  points(post, yticks, pch=4)
  axis(2, at=yticks, labels=vars, las=2)

  abline(v=10, lty=2, col="gray")
  if (!abs) {
    abline(v=0)
    abline(v=-10, lty=2, col="gray")
  }

  par(opar)
}


#' Covariate Residuals Plot
#'
#' A graphical balance diagnostic examining the ratio of residuals orthogonal to the propensity scores
#'
#' This function creates a plot of the residuals orthogonal to the propensity scores. Each covariate is
#' regressed on the propensity scores for the control and treatment populations. The variance of the residuals
#' of this regression are then compared. For a balanced population, the ratio of variances should be close to
#' 1. An imbalance in the variances indicates that there is additional variance potentially attributed to the
#' covariate in one of the populations that has not been captured in the propensity scores.
#'
#' @param data Data Frame, containing the dataset. The data
#' frame must contain a treatment indicator variable called 'treat'.
#' @param covariates Vector, containing the list of covariates to include in the plot
#' @param weights Vector, containing the weights to use in assessment of balanced population. Defaults to unweighted
#' @param dichotomous Boolean value indicating if dichotomous covariates should be included in the calculation or not
#' @return NULL
#'
#' @examples
#' ps.covariate.residualplot(myData, covariates)
#' @export
ps.covariate.residualplot <- function(data, covariates, weights = NULL, dichotomous = FALSE) {

  if (is.null(weights)) weights = rep.int(1, nrow(data))
  data$w <- weights

  # Check dichotomous input
  if (!is.logical(dichotomous)) {
    cat("\nValue entered for dichotomous is not a valid logical value. Will default to FALSE\n\n")
    dichotomous <- FALSE
  }

  # Identify which features are dichotomous
  if (dichotomous == FALSE) {
    d.variables <- names(data)[names(data) %in% covariates & apply(data, 2, function(x) { all(x %in% 0:1) })]
    covariates <- subset(covariates, !(covariates %in% d.variables))
  }

  # Analysis is only on balanced data
  bd <- data[data$w != 0,]
  ratios <- vector(length = length(covariates))

  for (i in 1:length(covariates)) {

    model.t <- lm(as.formula(paste(covariates[i], "ps_values", sep="~")), bd[bd$treat == 1,], weights=bd[bd$treat == 1,]$w)
    model.c <- lm(as.formula(paste(covariates[i], "ps_values", sep="~")), bd[bd$treat == 0,], weights=bd[bd$treat == 0,]$w)

    resid.t <- residuals(model.t)
    resid.c <- residuals(model.c)

    var.t <- var(resid.t)
    var.c <- var(resid.c)

    ratios[i] <- ifelse(var.c != 0, var.t / var.c, NA)
  }

  opar <- par()
  if (dichotomous == TRUE) {
    par(mar=c(8, 8, 4, 2) + 0.1)
  } else {
    par(mar=c(5, 8, 4, 2) + 0.1)
  }

  xlim= c(0, max(2.2, max(ratios, na.rm=TRUE)))

  yticks <- seq(1, length(covariates))
  plot(ratios, yticks,
       main="Covariate Ratio of Variance Orthogonal to Propensity Scores",
       xlab="Ratio of Variance",
       ylab="",
       yaxt="n",
       col="red",
       xlim=xlim)

  if (dichotomous == TRUE) {
    mtext(text="Warning: low frequency events will not be meaningful.", side = 1, cex = 1, line = 5, col = "red")
  }

  axis(2, at=yticks, labels=covariates, las=2)

  abline(v=0.5, lty=2, col="gray")
  abline(v=1)
  abline(v=2, lty=2, col="gray")

  suppressWarnings(par(opar))
}
