
#' Cohort Method Data
#'
#' Reformat a data frame into a CohortMethodData object for use with OHDSI package
#'
#' This function takes a data frame containing an observational study dataset and reformats
#' into an object of type CohortMethodData that can be used with the OHDSI CohortMethod library.
#'
#' @param data Data Frame - containing the dataset to reformat.
#' @param outcomes - outcome vector, or list of outcome vectors to include in data object. For TTE outcomes,
#' include event times (NA for no event), for dichtomous data include 1/0 indicator
#' @param covariates Vector of covariate names to include in the data object
#' @param treat.name Name of the treatment indicator [0,1] variable (default "treat")
#' @param person.id.name Name of the subject id variable (default "PERSON_ID")
#' @return CohortMethodData - reformatted data object
#'
#' @examples
#' ps.getCohortMethodData(myData)
#' @export
ps.getCohortMethodData <- function(data, outcomes, covariates, treat.name = "treat", person.id.name = "PERSON_ID") {

  cId = 2

  # outcomes can be a list or a list of lists...
  if (length(outcomes[[1]]) == 1) outcomes <- list(outcomes)

  # Add row id's used for indexing
  data$rowId <- seq(1, nrow(data))

  # Create the cohorts object - this defines the treatment and control group subjects
  cohorts <- data.frame(
                          rowId = data$rowId,
                          treatment = data[,treat.name],
                          personId = data[,person.id.name]
                          # TODO (timeToObsPeriodEnd, timeToCohortEnd)
                       )

  # Create the outcomes object - this defines the outcome variables (outcomes)
  outcomeData <- data.frame(rowId = numeric(), outcomeId = numeric(), timeToEvent = numeric())

  # Loop over each outcome variable
  for (i in 1:length(outcomes)) {
    mask <- !is.na(outcomes[[i]]) & outcomes[[i]] != 0
    outcomeData <- rbind(outcomeData,
                         data.frame(rowId = data[mask,]$rowId, outcomeId = rep(i + cId, sum(mask)), timeToEvent = outcomes[[i]][mask], y = outcomes[[i]][mask])
                    )
  }
  outcomeData[order(outcomeData$rowId),]
  cId = cId + length(outcomes)

  # Create the covariates object
  cov <- data.frame(rowId = numeric(), covariateId = numeric(), covariateValue = numeric())
  covariateRef <- data.frame(covariateId = numeric(), covariateName = character(), analysisId = numeric(), conceptId = numeric())

  # Loop over each of the covariates specified
  for (i in 1:length(covariates)) {
    # Determine the type of data - dichotomous indicator variables are stored in a tall-table format
    if (all(unique(data[,covariates[[i]]] %in% c(0, 1)))) {
      # For binary, indicator variables we take only the rows that are flagged as TRUE
      cRows <- data[data[,covariates[[i]]] == 1,]$rowId
      covariateValue <- rep(1, length(cRows))
    }
    else {
      cRows <- data$rowId
      covariateValue <- data[,covariates[[i]]]
    }

    cov <- rbind(cov, data.frame(rowId = cRows, covariateId = rep(i + cId, length(cRows)), covariateValue = covariateValue))
    covariateRef <- rbind(covariateRef, data.frame(covariateId = i + cId, covariateName = covariates[i], analysisId = 1, conceptId = i + cId))
  }
  cov[order(cov$covariateId, cov$rowId),]

  # Create an empty exclude data object
  exclude = data.frame(rowId = numeric(), outcomeId = numeric())

  # Assemble a metadata object
  metaData <- structure(list())
  metaData$targetId <- 1
  metaData$comparatorId <- 2
  metaData$outcomeIds <- unique(outcomeData$outcomeId)

  # Evaluate the counts
  n <- nrow(data)
  nt <- length(cohorts$treatment)
  metaData$counts <- data.frame(treatment = c(0, 1), notExcludedCount = c(n-nt, nt))

  # Assemble everything into the cohort method data object
  object <- structure(list(), class="cohortMethodData")

  object[["outcomes"]] <- ff::as.ffdf(outcomeData)
  object[["cohorts"]] <- ff::as.ffdf(cohorts)
  object[["covariates"]] <- ff::as.ffdf(cov)
  object[["exclude"]] <- ff::as.ffdf(data.frame(rowId = c(NA), outcomeId = c(NA)))
  object[["covariateRef"]] <- ff::as.ffdf(covariateRef)
  object[["metaData"]] <- metaData

  return (object)
}

