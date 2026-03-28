test_that("kbsd runs and returns correct structure", {

  data <- data.frame(
    x = rnorm(50),
    y = rnorm(50)
  )

  int_data_list <- list(
    data,
    transform(data, y = 1)
  )

  disthalf_vec <- c(x = 1, y = 1)

  result <- kbsd(data, int_data_list, disthalf_vec)

  expect_true(is.data.frame(result))
  expect_true(all(c("intervention", "diagnostic", "observation") %in% names(result)))
  expect_true(nrow(result) > 0)

})
