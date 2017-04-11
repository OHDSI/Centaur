#' Centaur Propensity Score Calculation
#'
#' Performs propensity score Calculation.
#'
#' This function calculates the propensity score based on the specified options.
#' In order for this function to work correctly, the \code{data} argument must be a data frame
#' containing the collection of confounding variables and a treatment indicator factor variable.
#' Additionally, the \code{covariates} argument must contain a vector of the variable names to
#' be included as covariates or confounding variables. See list of arguments for additional
#' options.
#'
#' @param data Data Frame - containing the dataset to be balanced. Must include treatment indicator as {0, 1} factor.
#' @param covariates Vector - containing the list of variable names to be included as confounding variables
#' @param ps.method String - name of the method to use for calculation of propensity scores. Options are
#'        \code{glm} (default) and \code{twang}
#' @param max.covariates The maximum number of covariates that can be used to calculate propensity scores (default 200)
#' @param max.twang The maximum number of samples that can be used with the twang ps.method (default 30k)
#' @param min.twang The minimum number of covariates required to use the twang ps.method (default 25)
#' @param control.ratio The desired ratio of control to treatment subjects. A large imbalance between
#'        control and treatment subjects can cause problems with algorithm convergence. If the control
#'        population exceeds control.ratio times as many subjects as the treatment group, the control population
#'        will be randomly sampled to the desired size.
#' @param random.seed Sets the random number generator seed, which determines how the control population is
#'        sampled. Override the default (43762116) to generate a different sample of control subjects.
#' @param odds.ratio A logical argument indicating if the odds ratio and 95% confidence interval should be calculated.
#'        This is an intensive calculation that can take a while.
#' @param lr.summary.file The file name where the logistic regression summary should be written to. If directory is
#'        provided as part of the file name, it has to already exist. Default value is NULL. If not given, no file
#'        will be written.
#' @return Data Frame - containing the original dataset, trimmed of any incomplete cases with additional
#'          variable added for \code{ps_values} (the calculated propensity scores)
#'
#' @examples
#' ps.score(myData, myCovariates)
#' ps.score(myData, myCovariates, ps.method = "twang")
#' @export
#'
ps.score <- function(data,
                     covariates,
                     ps.method = "glm",
                     max.covariates = 200,
                     max.twang = 30000,
                     min.twang = 25,
                     control.ratio = 5,
                     random.seed = 43762116,
                     odds.ratio = FALSE,
                     lr.summary.file = NULL) {

  #-----------------------------------------------------------
  # Do some basic QA to confirm that things are defined properly

  # Ensure that the input is a data frame
  if (class(data) != class(data.frame())) {
    stop("The supplied data must be of the class 'data.frame'")
  }

  # Ensure that the data frame has a column called "treat"
  if (!("treat" %in% colnames(data))) {
    stop("The data frame must have a column named 'treat'")
  }

  # Check odds.ratio input
  if (!is.logical(odds.ratio)) {
    cat("\nValue entered for odds.ratio is not a valid logical value. Will default to FALSE\n\n")
    odds.ratio <- FALSE
  }

  if (length(covariates) > max.covariates) {
    stop(paste("Too many covariates. ps.balance can handle up to ", max.covariates, " covariates"))
  }

  # Check the initial balance of the dataset
  treat.count <- nrow(data[data$treat == 1,])
  control.count <- nrow(data[data$treat == 0,])
  control.count.desired <- floor(treat.count * control.ratio)
  if (control.count > control.count.desired) {
    warning(paste("Number of control subjects exceeds recommended ratio of control to treat.\n", control.count.desired, "control subjects will be randomly sampled from the population."))
    # Take all the treatment subjects, then sample the desired number of control subjects, and combine back into a single data frame
    set.seed(random.seed)
    data <- rbind(data[data$treat == 1,], data[sample(which(data$treat == 0), control.count.desired),])
  }

  #-----------------------------------------------------------
  # Digest the data and build the general model formula
  #
  #d.variables <- names(data)[names(data) %in% covariates & apply(data, 2, function(x) { all(x %in% 0:1) })]
  #c.variables <- names(data)[names(data) %in% covariates & !(names(data) %in% d.variables)]

  formula <- as.formula(paste("treat ~ ", paste(covariates, collapse = " + ")))

  # Clean the data - remove entries where any data fields are NA
  to.be.removed <- nrow(data) - sum(complete.cases(data[,names(data) %in% covariates]))
  if (to.be.removed > 0) {
    warning(paste("Dataset contains incomplete cases. ", to.be.removed, " subjects will be removed."), immediate. = TRUE)
    data <- data[rowSums(is.na(data[,names(data) %in% covariates]))<1,]
  }


  # Calculate the Propensity Score
  if (ps.method == "glm") {

    # Prepare summary matrix
    numCovariates <- length(covariates)
    logisticSummary <- matrix(ncol = 7, nrow = (numCovariates+1))
    formattedSummary <- matrix(ncol = 7, nrow = (numCovariates+1))
    colnames(logisticSummary) <- c("Coeff. ", "Odds Ratio", "Std. Err.", "z", "P>|z|", "2.5%", "97.5%")
    rownames(logisticSummary) <- union(c("(Intercept)"), covariates)
    colnames(formattedSummary) <- colnames(logisticSummary)
    rownames(formattedSummary) <- rownames(logisticSummary)

    # Calculate the PS via logistic regression
    logisticModel = glm(formula, data = data, family = binomial(), control = list(maxit = 100));
    if (odds.ratio) {
      oddsRatio <- exp(cbind(OR = coef(logisticModel), confint(logisticModel)))
    }

    data$ps_values <-  predict(logisticModel, type = "response")

    # Covariates handled by the logistic regression
    logisticCovariates <- rownames(summary(logisticModel)$coefficients)

    # Populate logisticSummary matrix according to the desired layout
    # This matrix is kept around in case we need to output it at some point in the future.
    for (i in 1:(numCovariates+1)) {

      covariateName <- rownames(logisticSummary)[i]

      # Only care about covariates handled by the logistic regression
      if (covariateName %in% logisticCovariates) {
        logisticSummary[i,1] <- summary(logisticModel)$coefficients[covariateName, 1]
        logisticSummary[i,3] <- summary(logisticModel)$coefficients[covariateName, 2]
        logisticSummary[i,4] <- summary(logisticModel)$coefficients[covariateName, 3]
        logisticSummary[i,5] <- summary(logisticModel)$coefficients[covariateName, 4]
        if (odds.ratio) {
          logisticSummary[i,2] <- oddsRatio[covariateName,1]
          logisticSummary[i,6] <- oddsRatio[covariateName,2]
          logisticSummary[i,7] <- oddsRatio[covariateName,3]
        }
      } else {
        logisticSummary[i,1] <- NA
      }
    }

    # This matrix is constructed purely for reporting aesthetics
    for (i in 1:(numCovariates+1)) {
      for (j in 1:7)  {

        # If we have an NA, skip it
        if (is.na(logisticSummary[i,j])){
          next
        }

        # Conditional formatting of the columns
        if (j == 5) {
          formattedSummary[i,j] <- format(logisticSummary[i,j], scientific=TRUE, digits = 4)
        } else {
          if (abs(logisticSummary[i,j]) < 0.00001) {
            formattedSummary[i,j] <- format(logisticSummary[i,j], scientific=TRUE, digits = 3)
          } else {
            formattedSummary[i,j] <- format(logisticSummary[i,j], digits = 4, scientific=FALSE)
          }
        }
      }
    }

    # Remove the unpopulated columns
    if (!odds.ratio) {
      formattedSummary<-formattedSummary[,-c(2,6,7)] #delete odds ratio columns
    }

    print(formattedSummary, quote = FALSE, print.gap = 3)

    # Output to file
    if (!is.null(lr.summary.file)){

      #Get the directory name, check if it exists
      output.directory <- dirname(lr.summary.file)
      if (dir.exists(output.directory)){

        # Safely write the file
        tryCatch(
          {
          output.file <- file(lr.summary.file, open = "wt")
          sink(output.file)
          print(formattedSummary, quote = FALSE, print.gap = 3)
          },
          error = function(e) {
            cat("\n\nError writing output file... skipping file output\n")
          },
          finally = {
            # Close any sinks left open
            while(sink.number() > 0) {
              sink()
            }
            close(output.file)
          })
      } else {
        cat("\n\nOutput directory does not exist... skipping file output\n")
      }
    }

  }
  else if (ps.method == "twang") {

    # Validate that Twang is appropriate given the dataset size
    #
    if (length(covariates) < min.twang) {
      stop("Too few covariates. Twang generalized boosted method performs poorly with few covariates. Please use glm instead")
    }
    if (nrow(data) > max.twang) {
      stop(paste("Too many cases. Twang is unavailable for greater than ", max.twang, " cases. Please use glm instead"))
    }

    # Load the twang library and calculate the PS.
    data$ps_values <-  twang::ps(formula,
                                  data = data,
                                  n.trees = 5000,
                                  interaction.depth = 2,
                                  shrinkage = 0.01,
                                  bag.fraction = 1,
                                  perm.test.iters=0,
                                  print.level = 2,
                                  iterlim = 1000,
                                  verbose = FALSE,
                                  estimand = "ATT",
                                  stop.method = "ks.mean")$ps[,1]
  }
  else {
    stop(paste("Unrecognized PS calculation method: ", ps.method))
  }

  return(data)
}
