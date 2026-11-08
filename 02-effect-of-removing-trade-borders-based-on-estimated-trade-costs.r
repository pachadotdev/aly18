library(data.table)
library(capybara)
library(aly2018)

sigma <- 7
ref_country <- "0-DEU"

# tight tolerances, otherwise the fixed effects are only accurate to ~1e-4 and the error is amplified in the
# multilateral resistances
ctrl <- list(return_fe = TRUE, dev_tol = 1e-12, center_tol = 1e-12)

# step 1: baseline ----

fit_bln <- fepoisson(
  trade ~ ln_dist + cntg + brdr | exporter + importer,
  data = aly2018,
  control = ctrl
)

summary(fit_bln)

# uncomment to see the robust std error

# fit_bln <- fepoisson(
#   trade ~ ln_dist + cntg + brdr | exporter + importer | pair,
#   data = aly2018,
#   control = ctrl
# )

# summary(fit_bln)

# trade costs, in the counterfactual there are no borders and the internal trade costs remain the same
b <- coef(fit_bln)

aly2018 <- copy(aly2018)

aly2018[, trade_bln := predict(fit_bln)]
aly2018[, t_ij_bln := exp(b["ln_dist"] * ln_dist + b["cntg"] * cntg + b["brdr"] * brdr)]
aly2018[, t_ij_ctrf := fifelse(exporter == importer, t_ij_bln, exp(b["ln_dist"] * ln_dist + b["cntg"] * cntg))]
aly2018[, t_ij_ctrf_1 := log(t_ij_ctrf)]

# baseline indexes, equations (7) and (8)
add_fes(aly2018, fit_bln, "exp_fes_bln", "imp_fes_bln")

aly2018[, omr_bln := output * expndr_deu / exp_fes_bln]
aly2018[, imr_bln := expndr / (imp_fes_bln * expndr_deu)]

aly2018[, rGDP_bln := sum((output / imr_bln^(1 / (1 - sigma)))[exporter == importer]), by = exporter]
aly2018[, tot_exp_bln := sum(trade_bln[exporter != importer]), by = exporter]

# step 2: conditional ----

# equation (11), the trade costs are fixed with an offset and only the fixed effects are estimated
fit_cdl <- fepoisson(
  trade ~ 0 | exporter + importer,
  data = aly2018,
  offset = ~t_ij_ctrf_1,
  control = ctrl
)

aly2018[, trade_cdl := predict(fit_cdl)]
add_fes(aly2018, fit_cdl, "exp_fes_cdl", "imp_fes_cdl")

aly2018[, omr_cdl := output * expndr_deu / exp_fes_cdl]
aly2018[, imr_cdl := expndr / (imp_fes_cdl * expndr_deu)]

aly2018[, rGDP_cdl := sum((output / imr_cdl^(1 / (1 - sigma)))[exporter == importer]), by = exporter]
aly2018[, tot_exp_cdl := sum(trade_cdl[exporter != importer]), by = exporter]

# step 3: full endowment ----

# output and expenditure change with the factory-gate prices, we iterate until the prices converge
# expenditure is a constant share of output (phi)
# the columns ending in _old store the values from the previous iteration, initially those of the conditional step

aly2018[, phi := (expndr / output)[exporter == importer], by = importer]

aly2018[, trade_pred := trade_cdl]
aly2018[, c("fes_old", "expndr_deu_old", "omr_old", "imr_old") := .(exp_fes_cdl, expndr_deu, omr_cdl, imr_cdl)]
aly2018[, fes_dom_old := exp_fes_cdl[exporter == importer], by = importer]

# factory-gate prices of each exporter (p_exp) and of each importer as a producer (p_imp), and changes in the
# multilateral resistances
aly2018[, p_exp := (exp_fes_cdl / exp_fes_bln)^(1 / (1 - sigma))]
aly2018[, p_imp := (fes_dom_old / exp_fes_bln[exporter == importer])^(1 / (1 - sigma)), by = importer]
aly2018[, c("omr_ch", "imr_ch") := 1]

# converged when the standard deviation and the maximum of the change in p_exp are below 0.001
p_exp_old <- 0
sd_diff <- 1
max_diff <- 1
iter <- 1

while (sd_diff > 1e-3 | max_diff > 1e-3) {
  iter <- iter + 1

  # equation (14), rescale the trade flows and estimate equation (11)
  aly2018[, trade_iter := trade_pred * p_exp * p_imp / (omr_ch * imr_ch)]

  fit_iter <- fepoisson(
    trade_iter ~ 0 | exporter + importer,
    data = aly2018,
    offset = ~t_ij_ctrf_1,
    control = ctrl
  )

  aly2018[, trade_pred := predict(fit_iter)]
  add_fes(aly2018, fit_iter, "exp_fes", "imp_fes")
  aly2018[, fes_dom := exp_fes[exporter == importer], by = importer]

  # update output and expenditure
  aly2018[, output_iter := sum(trade_pred), by = exporter]
  aly2018[, expndr_iter := (phi * output_iter)[exporter == importer], by = importer]
  aly2018[, expndr_deu_iter := sum(trade_pred[importer == ref_country])]

  # convergence criteria with the change in p_exp from the previous iteration
  diffs <- aly2018$p_exp - p_exp_old
  sd_diff <- sd(diffs)
  max_diff <- abs(max(diffs))
  p_exp_old <- aly2018$p_exp

  # update prices
  aly2018[, p_exp := ((exp_fes / fes_old) / (expndr_deu_iter / expndr_deu_old))^(1 / (1 - sigma))]
  aly2018[, p_imp := ((fes_dom / fes_dom_old) / (expndr_deu_iter / expndr_deu_old))^(1 / (1 - sigma))]

  # update the multilateral resistances, equations (7) and (8)
  aly2018[, omr_iter := output_iter * expndr_deu_iter / exp_fes]
  aly2018[, imr_iter := expndr_iter / (imp_fes * expndr_deu_iter)]

  aly2018[, c("omr_ch", "imr_ch") := .(omr_iter / omr_old, imr_iter / imr_old)]

  aly2018[, c("fes_old", "fes_dom_old", "expndr_deu_old", "omr_old", "imr_old") := .(
    exp_fes, fes_dom, expndr_deu_iter, omr_iter, imr_iter
  )]

  cat("Iteration", iter, "- sd:", sd_diff, "- max:", max_diff, "\n")
}

# full endowment indexes, as in the original Stata code the expenditure of the reference country
# is phi * output from the last iteration, and not the sum of the predicted flows used in the loop
aly2018[, expndr_deu_full := expndr_iter[importer == ref_country][1]]

aly2018[, p_full := ((exp_fes / exp_fes_bln) / (expndr_deu_full / expndr_deu))^(1 / (1 - sigma))]
aly2018[, output_full := p_full * output]
aly2018[, expndr_full := (phi * output_full)[exporter == importer], by = importer]

aly2018[, omr_full := output_full * expndr_deu_full / exp_fes]
aly2018[, imr_full := expndr_iter / (imp_fes * expndr_deu_full)]

aly2018[, rGDP_full := sum((output_full / imr_full^(1 / (1 - sigma)))[exporter == importer]), by = exporter]

aly2018[, trade_full := output_full * expndr_full * t_ij_ctrf / (imr_full * omr_full)]
aly2018[, tot_exp_full := sum(trade_full[exporter != importer]), by = exporter]

# step 4: indexes ----

# the exporter indexes are constant for each exporter and the importer indexes for each importer, so we keep
# one observation of each
indexes_exp <- aly2018[, .(
  tot_exp_bln = tot_exp_bln[1], tot_exp_cdl = tot_exp_cdl[1], tot_exp_full = tot_exp_full[1],
  rGDP_bln = rGDP_bln[1], rGDP_cdl = rGDP_cdl[1], rGDP_full = rGDP_full[1],
  omr_bln = omr_bln[1], omr_cdl = omr_cdl[1], omr_full = omr_full[1],
  p_full = p_full[1]
), by = .(country = exporter)]

indexes_imp <- aly2018[, .(
  imr_bln = imr_bln[1], imr_cdl = imr_cdl[1], imr_full = imr_full[1]
), by = .(country = importer)]

indexes <- merge(indexes_exp, indexes_imp, by = "country")
indexes[, country := fifelse(country == ref_country, "DEU", country)]
setorder(indexes, country)

# changes in percentage as reported in the paper
indexes_final <- indexes[, .(
  country,
  exports_cdl = (tot_exp_cdl / tot_exp_bln - 1) * 100,
  exports_full = (tot_exp_full / tot_exp_bln - 1) * 100,
  rGDP_cdl = (rGDP_cdl / rGDP_bln - 1) * 100,
  rGDP_full = (rGDP_full / rGDP_bln - 1) * 100,
  imr_cdl = price_change(imr_cdl, imr_bln),
  imr_full = price_change(imr_full, imr_bln),
  omr_cdl = price_change(omr_cdl, omr_bln),
  omr_full = price_change(omr_full, omr_bln),
  prices_full = (p_full - 1) * 100
)]

cols <- setdiff(names(indexes_final), "country")
indexes_final[, (cols) := lapply(.SD, round, 2), .SDcols = cols]

indexes_final

fwrite(indexes_final, "results/table1-effect-of-removing-trade-borders-based-on-estimated-trade-costs.csv")

# step 5: validation with the non-linear solver ----

# solves the systems of the original MATLAB code (ge_ppml_indexes.m) directly, without iterating on regressions
indexes_nls <- nls_geppml(aly2018, sigma = sigma, ref_country = ref_country)

indexes_nls <- indexes_nls[, lapply(.SD, round, 2), by = country]

cols <- setdiff(names(indexes_final), "country")

m <- merge(indexes_final, indexes_nls, by = "country", suffixes = c("_it", "_nls"))

diff_dt <- m[, c(list(country = country),
                 setNames(lapply(cols, function(v) abs(get(paste0(v, "_it")) - get(paste0(v, "_nls")))), cols))]

# max gap per index
diff_dt[, lapply(.SD, max), .SDcols = cols]
