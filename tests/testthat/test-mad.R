context("Test median +- k*MAD")

test_that("Correct calculation of median +- k*MAD", {
  madSelection <- sapply("P001", madSelect, exampleDat,
    control = "control", treatment = "treatment", simplify = FALSE)
  expect_equal(as.numeric(madSelection$P001[1, 2:4]),
      c(0.8966, 0.8981399, 1.001757), tolerance = 1e-3)
  expect_equal(as.character(madSelection$P001[1, 5]), "No")
  expect_equal(rownames(madSelection$P001)[1], "AAK1 si3")
})






