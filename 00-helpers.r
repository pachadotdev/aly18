try(dir.create("results"), silent = TRUE)

if (!require("data.table")) install.packages("data.table", repos = "https://cran.r-project.org")
if (!require("capybara")) install.packages("capybara", repos = "https://cran.r-project.org")

library(data.table)
library(capybara)

dt <- fread("ge_ppml_data.csv", integer64 = "numeric")
sigma <- 7
ref_country <- "0-DEU"

# adds exp(fixed effects) to dt, normalized so the importer fixed effect of the reference country is zero
add_fes <- function(dt, fit, exp_name, imp_name) {
  fe_imp <- fit$fixed_effects$importer
  fe_exp <- fit$fixed_effects$exporter + fe_imp[ref_country]
  fe_imp <- fe_imp - fe_imp[ref_country]

  dt[data.table(exporter = names(fe_exp), fe = exp(fe_exp)), on = "exporter", (exp_name) := i.fe]
  dt[data.table(importer = names(fe_imp), fe = exp(fe_imp)), on = "importer", (imp_name) := i.fe]
}

# OMR and IMR are power transforms of the price indexes, the paper reports the change in the price indexes
price_change <- function(new, old) ((new / old)^(1 / (1 - sigma)) - 1) * 100

# tight tolerances, otherwise the fixed effects are only accurate to ~1e-4 and the error is amplified in the
# multilateral resistances
ctrl <- list(return_fe = TRUE, dev_tol = 1e-12, center_tol = 1e-12)
