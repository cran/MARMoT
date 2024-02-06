test_that("ASB works", {
  expect_length(ASB(data = MARMoT_data,
                   confounders = c("race", "age"), treatment = "hospital"), 2)
})
