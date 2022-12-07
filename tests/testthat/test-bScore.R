context("Test output of bscore and t-test")

test_that("b-score is correctly calculated", {
  res <- sapply("P001", bScore, example_dt, control = "control",
                       treatment = "treatment", simplify = FALSE)
  res_1 <- res$P001
  expect_equal(as.numeric(res_1[1, ]), c(-0.379059117,  1.4344181, -0.8484694,
                                         -0.1088711, -1.14224946, -0.0854374), tolerance = 1e-3)
  ttest_res  <- tTest(res_1, 3, 3)
  expect_equal(as.numeric(ttest_res[1, ]), c(0.55670002,  0.51448255, 0.8761181), tolerance = 1e-3)
})
