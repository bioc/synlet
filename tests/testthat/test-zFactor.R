context("Test calculation of z factor")

test_that("Correct calculation of z factor", {
  res <- zFactor(exampleDat, negativeCon = "scrambled control si1",
    positiveCon = "PLK1 si1")
  expect_equal(as.numeric(res[1, ]), c(-10.63422, -41.36317), tolerance = 1e-3)
})








