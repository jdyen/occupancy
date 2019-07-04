### occupancy: an R package for fitting occupancy-detection models with JAGS.

`occupancy` lets you write your model with a standard R mixed-model formula, and provides several methods for summarising, visualising, validating, and predicting from fitted models.

You can install the current version of the package (0.0.1) from GitHub:

``` r
devtools::install_github("jdyen/occupancy")
```


### System requirements
The `occupancy' package uses the JAGS software to fit all models. JAGS will need to be downloaded and installed separately. You can find details about JAGS installation [here](http://mcmc-jags.sourceforge.net/).


### Example analysis
The package includes some simulated data, and you can fit an occupancy-detection model to these data with the following code:
```
# load the occupancy package
library(occupancy)

# fit a model
mod <- occupancy(response ~ occ_predictor1 + occ_predictor2 +
                   (1 | occ_random1) + (1 | occ_random2),
                 ~ detect_predictor1 + detect_predictor2 + 
                   (1 | detect_random1),
                 site_id = "site",
                 survey_id = "survey",
                 data = occupancy_data,
                 jags_settings = list(n_iter = 1000, n_burnin = 500))

# plot the model coefficients
plot(mod)

# extract model coefficients
coef(mod)

# calculate model fit statistics
calculate_metrics(mod)

# cross validate a fitted model
validate(mod, n_cv = 5)
```

[![build status](https://travis-ci.org/jdyen/occupancy.svg?branch=master)](https://travis-ci.org/jdyen/occupancy) [![codecov.io](https://codecov.io/github/jdyen/occupancy/coverage.svg?branch=master)](https://codecov.io/github/jdyen/occupancy?branch=master) [![license](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)