test_that("blbglm works", {
  glm <- blbglm(Species ~ Sepal.Length * Sepal.Width, data = iris[1:100,],
                   m = 3, B = 100, family = binomial)
  expect_s3_class(glm, "blbglm")
  glm2 <- coef.blblm(glm)
  expect_equal(length(glm2), 4)
})