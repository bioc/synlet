context("Test output of bscore and t-test")

test_that("b-score is correctly calculated", {
  bscore.res <- sapply("P001", bScore, exampleDat, 
    control = "control", treatment = "treatment", simplify = FALSE)
  bscore.res.1 <- bscore.res$P001
  expect_equal(as.numeric(bscore.res.1[1, ]), c(-0.3888287, -0.6491217, 
    0.1497299, 1.5912198, 1.0046037, -1.5310097), tolerance = 1e-3)
  bscore.ttest  <- sapply(names(bscore.res), tTest, bscore.res, 
      numTreat = 3, numCont = 3, simplify = FALSE, USE.NAMES = TRUE)
  expect_equal(as.numeric(bscore.ttest$P001[1, ]),
      c(0.5706962, -0.6510114, 0.8836586), tolerance = 1e-3)

})




