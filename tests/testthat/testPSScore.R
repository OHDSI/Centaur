library(Centaur)
context("Propensity Score Computation")

test_that("ps.score adds ps_values to the data", {
  set.seed(20)
n = 2000
mydata = {}
mydata$AGE = c(rnorm(n, 65, 30), rnorm(n, 25, 30))
mydata$WEIGHT = c(rnorm(n, 200, 30), rnorm(n, 120, 30))
mydata$BP = c(rnorm(n, 140, 5), rnorm(n, 120, 5))
mydata$treat = c(matrix(1, n, 1), matrix(0,n,1))
mydata = as.data.frame(mydata);

covar = names(mydata);
covar = covar[1:3]
mydata$treat = c(matrix(1, n, 1), matrix(0,n,1))
res = ps.score(mydata, covar, ps.method = "glm")

## Verify that GLM was used
expect_true("ps_values" %in% colnames(res))})



context("PS comparison")
mydata <- data.frame(ps_values =rnorm(100),
                     treat = rep(1:2, rep(50,2)))
comp <- ps.compare(mydata)
expect_is(comp, "list")
