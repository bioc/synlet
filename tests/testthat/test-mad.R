context("Test median +- k*MAD")

test_that("Correct calculation of median +- k*MAD", {
  res <- madSelect("P001", example_dt, control = "control", treatment = "treatment")
  expect_equal(as.numeric(res[1, 2:4]), c(0.9427558, 1.0201578, 1.0821018), tolerance = 1e-3)
  expect_equal(rownames(res)[1], "A2M si1")
})
