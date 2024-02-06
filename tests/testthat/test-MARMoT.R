test_that("multiplication works", {
  expect_length(MARMoT(data = MARMoT_data, confounders = c("race", "age"),
                       treatment = "hospital", n.cores = 1), 4)
})
