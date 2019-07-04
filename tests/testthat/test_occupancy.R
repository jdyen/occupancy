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

# validate

# plot, coef

# errors on data in wrong format

# predictions

# spatial predictions
