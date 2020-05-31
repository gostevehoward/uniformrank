require(testthat)
require(uniformrank)

test_that("scores and test statistics are calculated correctly",{
  result <- get_scores_and_T_by_score(c(-6, -1, 8, 4), score_fns("wsrt"))
  expect_equal(result$scores, (1:4) / 5)
  expect_equal(result$T, c(4/5, 4/5, 6/5, 6/5))
})

test_that("scores and test statistics work with tied data",{
  result <- get_scores_and_T_by_score(c(5, 1, -5, 5, 10, 5, -20), score_fns("wsrt"))
  expect_equal(result$scores, c(1, 3.5, 3.5, 3.5, 3.5, 6, 7) / 8)
  expect_equal(result$T, c(0, 6, 6, 6, 6, 16.5, 17.5) / 8)
})

