test_that("blblm works", {
  x <- blblm(mpg ~ wt * disp * drat, data= mtcars, m = 3, B = 100)
  expect_s3_class(x, "blblm")
  colm <- coef(x)
  expect_equal(length(x), 2)
})
