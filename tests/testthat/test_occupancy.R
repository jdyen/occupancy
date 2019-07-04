context("occupancy")

expect_ok <- function(expr)
  expect_error(expr, NA)

test_that("occupancy models can be fitted in all variations", {

  # test full models
  expect_ok(occupancy(response ~ occ_predictor1 + occ_predictor2 + (1 | occ_random1) + (1 | occ_random2),
                      ~ detect_predictor1 + detect_predictor2 + (1 | detect_random1),
                      site_id = "site", survey_id = "survey",
                      data = occupancy_data, jags_settings = list(n_iter = 50, n_burnin = 10)))
  
  # test without detection random effects
  expect_ok(occupancy(response ~ occ_predictor1 + occ_predictor2 + (1 | occ_random1) + (1 | occ_random2),
                      ~ detect_predictor1 + detect_predictor2,
                      site_id = "site", survey_id = "survey",
                      data = occupancy_data, jags_settings = list(n_iter = 50, n_burnin = 10)))

  # test without occupancy random effects
  expect_ok(occupancy(response ~ occ_predictor1 + occ_predictor2,
                      ~ detect_predictor1 + detect_predictor2 + (1 | detect_random1),
                      site_id = "site", survey_id = "survey",
                      data = occupancy_data, jags_settings = list(n_iter = 50, n_burnin = 10)))

  # test without any random effects
  expect_ok(occupancy(response ~ occ_predictor1 + occ_predictor2,
                      ~ detect_predictor1 + detect_predictor2,
                      site_id = "site", survey_id = "survey",
                      data = occupancy_data, jags_settings = list(n_iter = 50, n_burnin = 10)))
  
  # test without detection fixed effects
  expect_ok(occupancy(response ~ occ_predictor1 + occ_predictor2 + (1 | occ_random1) + (1 | occ_random2),
                      ~ (1 | detect_random1),
                      site_id = "site", survey_id = "survey",
                      data = occupancy_data, jags_settings = list(n_iter = 50, n_burnin = 10)))
  
  # test without occupancy fixed effects
  expect_ok(occupancy(response ~ (1 | occ_random1) + (1 | occ_random2),
                      ~ detect_predictor1 + detect_predictor2 + (1 | detect_random1),
                      site_id = "site", survey_id = "survey",
                      data = occupancy_data, jags_settings = list(n_iter = 50, n_burnin = 10)))
  
  # test without any fixed effects
  expect_ok(occupancy(response ~ (1 | occ_random1) + (1 | occ_random2),
                      ~ (1 | detect_random1),
                      site_id = "site", survey_id = "survey",
                      data = occupancy_data, jags_settings = list(n_iter = 50, n_burnin = 10)))
  
  # test intercept only models
  expect_ok(occupancy(response ~ 1,
                      ~ 1,
                      site_id = "site", survey_id = "survey",
                      data = occupancy_data, jags_settings = list(n_iter = 50, n_burnin = 10)))
  
})

test_that("methods work as expected", {
  
  # test full models
  mod_tmp <- occupancy(response ~ occ_predictor1 + occ_predictor2 + (1 | occ_random1) + (1 | occ_random2),
                       ~ detect_predictor1 + detect_predictor2 + (1 | detect_random1),
                       site_id = "site", survey_id = "survey",
                       data = occupancy_data, jags_settings = list(n_iter = 50, n_burnin = 10))
  
  # list of methods at this stage, add tests of dims and output types
  expect_ok(coef(mod_tmp))
  expect_ok(summary(mod_tmp))
  expect_ok(fitted(mod_tmp))
  expect_ok(r2_calc(mod_tmp))
  expect_ok(calculate_metrics(mod_tmp))
  expect_ok(plot(mod_tmp))
  expect_ok(plot_pr_occ(mod_tmp))
  expect_ok(plot_pr_detect(mod_tmp))
  
})

# errors on data in wrong format
test_that("methods work as expected", {
  
  # test full models
  mod_tmp <- occupancy(response ~ occ_predictor1 + occ_predictor2 + (1 | occ_random1) + (1 | occ_random2),
                       ~ detect_predictor1 + detect_predictor2 + (1 | detect_random1),
                       site_id = "site", survey_id = "survey",
                       data = occupancy_data, jags_settings = list(n_iter = 50, n_burnin = 10))
  
  # list of methods at this stage, add tests of dims and output types
  expect_ok(coef(mod_tmp))
  expect_ok(summary(mod_tmp))
  expect_ok(fitted(mod_tmp))
  expect_ok(r2_calc(mod_tmp))
  expect_ok(calculate_metrics(mod_tmp))
  expect_ok(plot(mod_tmp))
  expect_ok(plot_pr_occ(mod_tmp))
  expect_ok(plot_pr_detect(mod_tmp))
  
})

# predictions
test_that("models can predict to fitted and new data", {
  
  # test full models
  mod_tmp <- occupancy(response ~ occ_predictor1 + occ_predictor2 + (1 | occ_random1) + (1 | occ_random2),
                       ~ detect_predictor1 + detect_predictor2 + (1 | detect_random1),
                       site_id = "site", survey_id = "survey",
                       data = occupancy_data, jags_settings = list(n_iter = 50, n_burnin = 10))
  
  # predictions to in-sample data
  expect_ok(predict(mod_tmp))
  
  # predictions to new data
  idx <- sample(seq_len(100), size = 200, replace = TRUE)
  expect_ok(predict(mod_tmp, newdata = occupancy_data[idx, ]))
  
  # errors if variable missing
  idx <- sample(seq_len(100), size = 200, replace = TRUE)
  idy <- sample(seq_len(ncol(occupancy_data)), size = 5, replace = FALSE)
  expect_error(predict(mod_tmp, newdata = occupancy_data[idx, idy]))
  
})

# spatial predictions
test_that("models can predict to new spatial data", {
  
  # test full models
  mod_tmp <- occupancy(response ~ occ_predictor1 + occ_predictor2 + (1 | occ_random1) + (1 | occ_random2),
                       ~ detect_predictor1 + detect_predictor2 + (1 | detect_random1),
                       site_id = "site", survey_id = "survey",
                       data = occupancy_data, jags_settings = list(n_iter = 50, n_burnin = 10))
  
  # predictions to new data
  expect_ok(spatial_predict(mod_tmp, newdata = raster_tmp))
  
})

# validate
test_that("methods work as expected", {
  
  # test full models
  mod_tmp <- occupancy(response ~ occ_predictor1 + occ_predictor2 + (1 | occ_random1) + (1 | occ_random2),
                       ~ detect_predictor1 + detect_predictor2 + (1 | detect_random1),
                       site_id = "site", survey_id = "survey",
                       data = occupancy_data, jags_settings = list(n_iter = 50, n_burnin = 10))
  
  # list of methods at this stage, add tests of dims and output types
  expect_ok(coef(mod_tmp))
  expect_ok(summary(mod_tmp))
  expect_ok(fitted(mod_tmp))
  expect_ok(r2_calc(mod_tmp))
  expect_ok(calculate_metrics(mod_tmp))
  expect_ok(plot(mod_tmp))
  expect_ok(plot_pr_occ(mod_tmp))
  expect_ok(plot_pr_detect(mod_tmp))
  
})