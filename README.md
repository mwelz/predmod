# predmod
Predictive modeling of heterogeneous treatment effects. Implements risk modeling, effect modeling, ratio ratio estimation,  and enables estimation via double/debiased machine learning and causal random forests. Risk modeling and effetc modeling are described in detail in the [PATH statement](https://www.acpjournals.org/doi/full/10.7326/M18-3668?rfr_dat=cr_pub++0pubmed&url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org) of Kent et al. (2020, Annals of Internal Medicine). This package is developed by [Max Welz](https://mwelz.github.io/), [Kevin ten Haaf](https://employees.publichealthrotterdam.com/profile/kevin-ten-haaf/), and [Andreas Alfons](https://personal.eur.nl/alfons/).

**This package is in development and we cannot yet guarantee stability or correctness.** In case of questions, please get in touch with Max Welz (`welz <at> ese <dot> eur <dot> nl`).

## Installation
To install `predmod`, use the function `install_github()` of the `devtools` package.

```R
install.packages("devtools")
devtools::install_github("mwelz/predmod")
```

## Example

```R
# load package
library("predmod", quietly = TRUE)
set.seed(1)

# generate data for 200 individuals for whom we observe 5 variables
n <- 200
p <- 5

# random treatment assignment
w <- rbinom(n, 1, 0.5)

# covariates
X <- matrix(runif(n * p), n)
colnames(X) <- paste0("var",1:p)

# mortality probability (heterogeneous; driven by X1 and X2)
pi <- plogis(X[,1] + X[,2] + w * X[,2])

# observed outcome 
status <- rbinom(n, 1, pi)

# run effect model
EM <- effect_model(X = X, status = status, w = w)

# make calibration plot
calibration_plot(EM)

# run risk model
RM <- risk_model(X = X, status = status, w = w)
calibration_plot(RM)

```
